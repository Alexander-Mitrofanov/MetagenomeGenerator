#!/usr/bin/env python3
from __future__ import annotations
"""BLASTN-based filtering of endogenous viral elements (EVEs) in non-viral sequences.

Aligns non-viral genomes (bacteria, archaea, plasmid) against a viral database,
parses hits, and produces per-sequence intervals to exclude when chunking so that
EVE regions are not used as negative training data.

For proper prophage/EVE detection, the viral DB should be comprehensive (e.g. all
RefSeq viral from a snapshot). Use build-viral-db to create it, then pass
--viral-db or --viral-reference-fasta to blastn-filter.
"""

import hashlib
import json
import logging
import subprocess
import sys
import time
from pathlib import Path

from .genome_layout import get_viral_fasta_paths, get_nonviral_fasta_paths

logger = logging.getLogger(__name__)


def _progress(msg: str) -> None:
    """Line to stdout with timestamp; always flush (important when stdout is piped to tee)."""
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}", flush=True)


# BLAST+ outfmt 6: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
BLAST_COLS = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore".split()
BLAST_TASK_CHOICES = {"blastn", "megablast", "dc-megablast"}


def _merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Merge overlapping or adjacent intervals. Input and output 0-based half-open [start, end)."""
    if not intervals:
        return []
    sorted_i = sorted(intervals)
    merged = [list(sorted_i[0])]
    for s, e in sorted_i[1:]:
        if s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return [(a, b) for a, b in merged]


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _blast_db_files(db_prefix: Path) -> list[Path]:
    files = sorted(p for p in db_prefix.parent.glob(f"{db_prefix.name}.*") if p.is_file())
    if not files:
        raise FileNotFoundError(f"No BLAST DB files found for prefix: {db_prefix}")
    return files


def viral_db_fingerprint(db_prefix: Path, *, log_each_file: bool = False) -> str:
    """Stable aggregate fingerprint over all BLAST DB files for db_prefix."""
    files = _blast_db_files(db_prefix)
    h = hashlib.sha256()
    n = len(files)
    for i, fp in enumerate(files, start=1):
        if log_each_file:
            sz = fp.stat().st_size
            _progress(f"EVE cache: hashing viral DB part {i}/{n}: {fp.name} ({sz:,} bytes)…")
        h.update(f"{fp.name}:{_sha256_file(fp)}\n".encode("utf-8"))
    return h.hexdigest()


def write_viral_db_manifest(
    db_prefix: Path,
    manifest_path: Path,
    *,
    snapshot_timestamp: str | None = None,
    snapshot_date: str | None = None,
    source_snapshot: str | None = None,
    viral_accession_count: int | None = None,
) -> Path:
    """Write manifest JSON for a viral BLAST DB with per-file checksums + aggregate fingerprint."""
    files = _blast_db_files(db_prefix)
    file_checksums = {fp.name: _sha256_file(fp) for fp in files}
    fp_hash = viral_db_fingerprint(db_prefix)
    payload = {
        "created_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "db_prefix": str(db_prefix),
        "db_name": db_prefix.name,
        "db_directory": str(db_prefix.parent),
        "snapshot_timestamp": snapshot_timestamp or "",
        "snapshot_date": snapshot_date or "",
        "source_snapshot": source_snapshot or "",
        "viral_accession_count": int(viral_accession_count) if viral_accession_count is not None else None,
        "aggregate_sha256": fp_hash,
        "file_sha256": file_checksums,
    }
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    with manifest_path.open("w") as f:
        json.dump(payload, f, indent=2)
    return manifest_path


def verify_viral_db(
    db_prefix: Path,
    *,
    manifest_path: Path | None = None,
    required_aggregate_sha256: str | None = None,
) -> None:
    """Validate DB exists and optionally matches manifest and/or explicit aggregate fingerprint."""
    files = _blast_db_files(db_prefix)
    _ = files  # non-empty assertion via helper
    if manifest_path is not None:
        if not manifest_path.exists():
            raise FileNotFoundError(f"Viral DB manifest not found: {manifest_path}")
        with manifest_path.open() as f:
            manifest = json.load(f)
        file_sha = manifest.get("file_sha256", {})
        for fp in files:
            expected = file_sha.get(fp.name)
            if not expected:
                raise ValueError(f"Manifest missing checksum for DB file: {fp.name}")
            observed = _sha256_file(fp)
            if observed != expected:
                raise ValueError(
                    f"Viral DB checksum mismatch for {fp.name}: expected {expected}, got {observed}"
                )
        expected_agg = manifest.get("aggregate_sha256")
        observed_agg = viral_db_fingerprint(db_prefix)
        if expected_agg and observed_agg != expected_agg:
            raise ValueError(
                f"Viral DB aggregate checksum mismatch: expected {expected_agg}, got {observed_agg}"
            )
    if required_aggregate_sha256:
        observed_agg = viral_db_fingerprint(db_prefix)
        if observed_agg != required_aggregate_sha256:
            raise ValueError(
                f"Viral DB fingerprint mismatch: expected {required_aggregate_sha256}, got {observed_agg}"
            )


def build_viral_db(viral_fasta: Path, db_dir: Path) -> Path:
    """Run makeblastdb on viral FASTA. Returns path to DB prefix (no extension)."""
    db_dir.mkdir(parents=True, exist_ok=True)
    db_prefix = db_dir / "viral_db"
    # Do not capture stdout/stderr: long-running BLAST tools can fill pipe buffers and deadlock.
    subprocess.run(
        ["makeblastdb", "-in", str(viral_fasta), "-out", str(db_prefix), "-dbtype", "nucl"],
        check=True,
    )
    return db_prefix


def run_blastn(
    query_fasta: Path,
    db_prefix: Path,
    out_tsv: Path,
    *,
    evalue: float = 1e-5,
    perc_identity: float = 70.0,
    max_target_seqs: int = 5,
    num_threads: int = 4,
    task: str = "dc-megablast",
) -> None:
    """Run blastn (query = non-viral, db = viral). Writes tabular output to out_tsv."""
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    if task not in BLAST_TASK_CHOICES:
        raise ValueError(f"Invalid BLAST task: {task}. Choose from {sorted(BLAST_TASK_CHOICES)}")
    fmt = "6 " + " ".join(BLAST_COLS)
    cmd = [
        "blastn",
        "-query", str(query_fasta),
        "-db", str(db_prefix),
        "-out", str(out_tsv),
        "-outfmt", fmt,
        "-evalue", str(evalue),
        "-perc_identity", str(perc_identity),
        "-max_target_seqs", str(max_target_seqs),
        "-num_threads", str(num_threads),
        "-task", task,
    ]
    # Do not capture stdout/stderr: blastn may write substantial progress to stderr; pipes can fill
    # and block the child while the parent waits — looks like an infinite hang.
    sys.stdout.flush()
    sys.stderr.flush()
    subprocess.run(cmd, check=True)
    sys.stdout.flush()
    sys.stderr.flush()


def _file_signature(path: Path) -> dict[str, int | str]:
    st = path.stat()
    return {
        "path": str(path.resolve()),
        "size": int(st.st_size),
        "mtime_ns": int(st.st_mtime_ns),
    }


def _load_json_if_exists(path: Path) -> dict | None:
    if not path.exists():
        return None
    with path.open() as f:
        return json.load(f)


EVE_QUERY_STORE_SCHEMA = 2


def _query_file_signature(path: Path) -> str:
    """Stable signature for a query FASTA (resolved path + size + mtime). Symlinks to the same file share a signature."""
    st = path.stat()
    raw = f"{path.resolve()}|{st.st_size}|{st.st_mtime_ns}".encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def _eve_blast_params_key(
    *,
    evalue: float,
    perc_identity: float,
    max_target_seqs: int,
    task: str,
) -> dict[str, float | int | str]:
    return {
        "evalue": float(evalue),
        "perc_identity": float(perc_identity),
        "max_target_seqs": int(max_target_seqs),
        "task": str(task),
    }


def _eve_params_bucket_id(params: dict[str, float | int | str]) -> str:
    blob = json.dumps(params, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(blob).hexdigest()[:28]


def _try_load_eve_query_entry(
    entry_path: Path,
    *,
    query_stem: str,
    query_sig: str,
    db_fingerprint: str,
    params: dict[str, float | int | str],
) -> dict[tuple[str, str], list[tuple[int, int]]] | None:
    data = _load_json_if_exists(entry_path)
    if not data or data.get("schema") != EVE_QUERY_STORE_SCHEMA:
        return None
    if data.get("query_stem") != query_stem or data.get("query_sig") != query_sig:
        return None
    if data.get("db_fingerprint") != db_fingerprint or data.get("params") != params:
        return None
    out: dict[tuple[str, str], list[tuple[int, int]]] = {}
    for k, v in data.get("intervals", {}).items():
        parts = k.split("\t", 1)
        if len(parts) != 2:
            continue
        out[(parts[0], parts[1])] = [tuple(x) for x in v]  # type: ignore[misc]
    return out


def _save_eve_query_entry(
    entry_path: Path,
    *,
    query_stem: str,
    query_sig: str,
    db_fingerprint: str,
    params: dict[str, float | int | str],
    intervals: dict[tuple[str, str], list[tuple[int, int]]],
) -> None:
    entry_path.parent.mkdir(parents=True, exist_ok=True)
    serial = {f"{a}\t{b}": [list(t) for t in ivs] for (a, b), ivs in intervals.items()}
    payload = {
        "schema": EVE_QUERY_STORE_SCHEMA,
        "created_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "query_stem": query_stem,
        "query_sig": query_sig,
        "db_fingerprint": db_fingerprint,
        "params": params,
        "intervals": serial,
    }
    with entry_path.open("w") as f:
        json.dump(payload, f, indent=0)


def parse_blastn_tabular(tsv_path: Path) -> dict[str, list[tuple[int, int]]]:
    """Parse BLAST outfmt 6 TSV. Returns dict: qseqid -> list of (start, end) 0-based half-open intervals."""
    by_qseqid: dict[str, list[tuple[int, int]]] = {}
    with tsv_path.open() as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 8:
                continue
            qseqid = parts[0]
            qstart, qend = int(parts[6]), int(parts[7])
            start = qstart - 1
            end = qend
            if start > end:
                start, end = end - 1, start
            by_qseqid.setdefault(qseqid, []).append((start, end))
    for qseqid in by_qseqid:
        by_qseqid[qseqid] = _merge_intervals(by_qseqid[qseqid])
    return by_qseqid


def _concat_viral_fasta(genome_dir: Path, out_fasta: Path) -> Path:
    """Concatenate viral FASTAs (virus/ or viral_*.fasta) in genome_dir into a single FASTA for makeblastdb."""
    from Bio import SeqIO
    viral_files = get_viral_fasta_paths(genome_dir)
    if not viral_files:
        raise FileNotFoundError(f"No viral FASTAs found in {genome_dir} (virus/ or viral_*.fasta)")
    out_fasta.parent.mkdir(parents=True, exist_ok=True)
    with out_fasta.open("w") as out_handle:
        n = 0
        for fp in viral_files:
            for rec in SeqIO.parse(fp, "fasta"):
                rec.id = f"{fp.stem}_{n}_{rec.id}"
                rec.description = ""
                SeqIO.write(rec, out_handle, "fasta")
                n += 1
    return out_fasta


def run_blastn_nonviral(
    nonviral_fasta_paths: list[Path],
    viral_fasta: Path,
    out_dir: Path,
    *,
    evalue: float = 1e-5,
    perc_identity: float = 70.0,
    max_target_seqs: int = 5,
    num_threads: int = 4,
    task: str = "dc-megablast",
) -> dict[tuple[str, str], list[tuple[int, int]]]:
    """Build viral DB, run BLASTN for each non-viral FASTA, return EVE intervals per (file_stem, qseqid)."""
    out_dir.mkdir(parents=True, exist_ok=True)
    blast_dir = out_dir / "blastn"
    blast_dir.mkdir(parents=True, exist_ok=True)

    db_prefix = build_viral_db(viral_fasta, blast_dir)
    logger.info("BLASTN filter: viral DB built at %s", db_prefix)
    _progress(f"Viral BLAST DB: {db_prefix}")

    eve_intervals: dict[tuple[str, str], list[tuple[int, int]]] = {}
    nq = len(nonviral_fasta_paths)
    for qi, qf in enumerate(nonviral_fasta_paths, start=1):
        out_tsv = blast_dir / f"{qf.stem}.blastn.tsv"
        logger.info("BLASTN filter: query=%s -> %s (evalue=%s, perc_identity=%s)",
                    qf.name, out_tsv.name, evalue, perc_identity)
        _progress(f"BLASTN ({qi}/{nq}) START {qf.name} -> {out_tsv.name} (task={task}, threads={num_threads})")
        t0 = time.perf_counter()
        run_blastn(
            qf,
            db_prefix,
            out_tsv,
            evalue=evalue,
            perc_identity=perc_identity,
            max_target_seqs=max_target_seqs,
            num_threads=num_threads,
            task=task,
        )
        _progress(f"BLASTN ({qi}/{nq}) DONE  {qf.name} in {time.perf_counter() - t0:.1f}s")
        by_id = parse_blastn_tabular(out_tsv)
        for qseqid, intervals in by_id.items():
            key = (qf.stem, qseqid)
            eve_intervals[key] = intervals
        if by_id:
            logger.info("BLASTN filter: %s -> %d sequences with EVE hits (excluded from chunking)",
                        qf.stem, len(by_id))

    eve_json = out_dir / "eve_intervals.json"
    serializable = {f"{k[0]}\t{k[1]}": v for k, v in eve_intervals.items()}
    with eve_json.open("w") as f:
        json.dump(serializable, f, indent=0)
    logger.info("BLASTN filter: eve_intervals.json written; %d sequences with EVE hits total", len(eve_intervals))
    _progress(f"EVE intervals: {eve_json} ({len(eve_intervals)} sequences with hits)")
    return eve_intervals


def run_blastn_from_dirs(
    genome_dir: Path,
    out_dir: Path,
    *,
    evalue: float = 1e-5,
    perc_identity: float = 70.0,
    max_target_seqs: int = 5,
    num_threads: int = 4,
    task: str = "dc-megablast",
    viral_reference_fasta: Path | None = None,
    viral_db_prefix: Path | None = None,
    viral_db_manifest: Path | None = None,
    required_viral_db_sha256: str | None = None,
    reuse_cache: bool = True,
    eve_query_store: Path | None = None,
) -> dict[tuple[str, str], list[tuple[int, int]]]:
    """Run BLASTN: non-viral (bacteria/, archaea/, plasmid/) vs viral DB.

    Viral DB source (first that is set wins):
    - viral_db_prefix: use this existing BLAST DB (e.g. from build-viral-db).
    - viral_reference_fasta: build DB from this FASTA (e.g. comprehensive RefSeq viral).
    - else: concatenate virus/ from genome_dir (current run's viral set only; limited for EVE detection).

    Per-query EVE cache: ``eve_query_store`` (default ``out_dir / "eve_query_store"``) holds one JSON
    per non-viral FASTA (keyed by resolved path + size + mtime, viral DB fingerprint, and BLAST params).
    New genomes are BLASTed and merged; existing entries are reused when ``reuse_cache`` is True.

    For proper prophage/EVE detection, use a comprehensive viral reference (build-viral-db from snapshot
    or a full RefSeq viral download) and pass --viral-reference-fasta or --viral-db.
    """
    nonviral = get_nonviral_fasta_paths(genome_dir)
    if not nonviral:
        raise FileNotFoundError(f"No non-viral genome FASTAs in {genome_dir} (bacteria/, archaea/, plasmid/ or flat)")

    out_dir.mkdir(parents=True, exist_ok=True)
    blast_dir = out_dir / "blastn"
    blast_dir.mkdir(parents=True, exist_ok=True)
    eve_json = out_dir / "eve_intervals.json"
    store_root = Path(eve_query_store) if eve_query_store is not None else (out_dir / "eve_query_store")

    if viral_db_prefix is not None:
        db_prefix = Path(viral_db_prefix)
        nhr = db_prefix.with_suffix(".nhr") if db_prefix.suffix else db_prefix.parent / (db_prefix.name + ".nhr")
        if not nhr.exists():
            raise FileNotFoundError(f"BLAST DB not found at {viral_db_prefix} (expected {nhr})")
        verify_viral_db(
            db_prefix,
            manifest_path=viral_db_manifest,
            required_aggregate_sha256=required_viral_db_sha256,
        )
        logger.info("BLASTN filter: using existing viral DB at %s", db_prefix)
        _progress(f"Using existing viral BLAST DB: {db_prefix}")
    elif viral_reference_fasta is not None:
        vref = Path(viral_reference_fasta)
        if not vref.exists():
            raise FileNotFoundError(f"Viral reference FASTA not found: {vref}")
        db_prefix = build_viral_db(vref, blast_dir)
        logger.info("BLASTN filter: viral DB built from reference at %s", db_prefix)
        _progress(f"Viral BLAST DB built from reference: {vref}")
    else:
        viral_concat = out_dir / "viral_concat.fasta"
        _concat_viral_fasta(genome_dir, viral_concat)
        db_prefix = build_viral_db(viral_concat, blast_dir)
        logger.info("BLASTN filter: viral DB built from genome_dir virus/ at %s", db_prefix)
        _progress(
            f"Viral BLAST DB: {db_prefix} (from genome_dir virus/; for better EVE detection use "
            "--viral-reference-fasta or build-viral-db)"
        )

    logger.info("BLASTN filter: computing viral DB fingerprint for per-query EVE store")
    _progress("Computing viral DB fingerprint for EVE query store (not running blastn yet)…")
    db_fingerprint = viral_db_fingerprint(db_prefix, log_each_file=True)
    _progress("Viral DB fingerprint complete.")
    params = _eve_blast_params_key(
        evalue=evalue,
        perc_identity=perc_identity,
        max_target_seqs=max_target_seqs,
        task=task,
    )
    params_bucket = _eve_params_bucket_id(params)
    store_bucket = store_root / db_fingerprint[:28] / params_bucket
    store_bucket.mkdir(parents=True, exist_ok=True)
    _progress(f"EVE query store: {store_root} (bucket …/{db_fingerprint[:12]}…/{params_bucket[:12]}…)")

    eve_intervals: dict[tuple[str, str], list[tuple[int, int]]] = {}
    nq = len(nonviral)
    n_cache_hit = 0
    for qi, qf in enumerate(nonviral, start=1):
        out_tsv = blast_dir / f"{qf.stem}.blastn.tsv"
        qsig = _query_file_signature(qf)
        entry_path = store_bucket / f"{qsig}.json"
        cached: dict[tuple[str, str], list[tuple[int, int]]] | None = None
        if reuse_cache:
            cached = _try_load_eve_query_entry(
                entry_path,
                query_stem=qf.stem,
                query_sig=qsig,
                db_fingerprint=db_fingerprint,
                params=params,
            )
        if cached is not None:
            n_cache_hit += 1
            eve_intervals.update(cached)
            _progress(f"BLASTN ({qi}/{nq}) CACHE HIT {qf.name} ({len(cached)} keyed interval rows)")
            continue

        logger.info("BLASTN filter: query=%s -> %s (evalue=%s, perc_identity=%s)",
                    qf.name, out_tsv.name, evalue, perc_identity)
        _progress(f"BLASTN ({qi}/{nq}) START {qf.name} -> {out_tsv.name} (task={task}, threads={num_threads})")
        t0 = time.perf_counter()
        run_blastn(
            qf,
            db_prefix,
            out_tsv,
            evalue=evalue,
            perc_identity=perc_identity,
            max_target_seqs=max_target_seqs,
            num_threads=num_threads,
            task=task,
        )
        _progress(f"BLASTN ({qi}/{nq}) DONE  {qf.name} in {time.perf_counter() - t0:.1f}s")
        by_id = parse_blastn_tabular(out_tsv)
        partial: dict[tuple[str, str], list[tuple[int, int]]] = {}
        for qseqid, intervals in by_id.items():
            key = (qf.stem, qseqid)
            eve_intervals[key] = intervals
            partial[key] = intervals
        if by_id:
            logger.info("BLASTN filter: %s -> %d sequences with EVE hits (excluded from chunking)",
                        qf.stem, len(by_id))
        _save_eve_query_entry(
            entry_path,
            query_stem=qf.stem,
            query_sig=qsig,
            db_fingerprint=db_fingerprint,
            params=params,
            intervals=partial,
        )

    serializable = {f"{k[0]}\t{k[1]}": v for k, v in eve_intervals.items()}
    with eve_json.open("w") as f:
        json.dump(serializable, f, indent=0)
    logger.info(
        "BLASTN filter: eve_intervals.json written; %d sequences with hits; %d/%d queries from store",
        len(eve_intervals), n_cache_hit, nq,
    )
    _progress(
        f"EVE intervals: {eve_json} ({len(eve_intervals)} sequences with hits); "
        f"store hits {n_cache_hit}/{nq}"
    )
    return eve_intervals


def load_eve_intervals(eve_json: Path) -> dict[tuple[str, str], list[tuple[int, int]]]:
    """Load EVE intervals from run_blastn_nonviral output (eve_intervals.json)."""
    with eve_json.open() as f:
        data = json.load(f)
    return {tuple(k.split("\t")): [tuple(iv) for iv in v] for k, v in data.items()}


def export_eve_regions_fasta(
    nonviral_fasta_paths: list[Path],
    eve_intervals: dict[tuple[str, str], list[tuple[int, int]]],
    out_fasta: Path,
    *,
    min_interval_length: int = 1,
) -> int:
    """Export EVE/provirus intervals as FASTA.

    Writes one FASTA record per interval with header:
      <file_stem>|<qseqid>|<start1>-<end1>

    where coordinates are 1-based inclusive on the original query sequence.
    """
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    if min_interval_length < 1:
        raise ValueError("min_interval_length must be >= 1")

    out_fasta.parent.mkdir(parents=True, exist_ok=True)

    n_written = 0
    with out_fasta.open("w") as out_handle:
        for fp in nonviral_fasta_paths:
            by_qseqid: dict[str, list[tuple[int, int]]] = {}
            for (stem, qseqid), intervals in eve_intervals.items():
                if stem != fp.stem:
                    continue
                if intervals:
                    by_qseqid[qseqid] = intervals
            if not by_qseqid:
                continue

            for rec in SeqIO.parse(fp, "fasta"):
                qseqid = rec.id
                intervals = by_qseqid.get(qseqid)
                if not intervals:
                    continue
                seq_str = str(rec.seq)
                for (start0, end0) in intervals:
                    if end0 <= start0:
                        continue
                    if (end0 - start0) < min_interval_length:
                        continue
                    sub = seq_str[start0:end0]
                    start1 = start0 + 1
                    end1 = end0
                    rid = f"{fp.stem}|{qseqid}|{start1}-{end1}"
                    out_rec = SeqRecord(Seq(sub), id=rid, description="")
                    SeqIO.write(out_rec, out_handle, "fasta")
                    n_written += 1

    logger.info("BLASTN filter: exported %d EVE intervals to %s", n_written, out_fasta)
    print(f"Exported EVE intervals FASTA: {out_fasta} ({n_written} intervals)")
    return n_written


def chunk_overlaps_eve(chunk_start: int, chunk_end: int, eve_intervals: list[tuple[int, int]]) -> bool:
    """True if [chunk_start, chunk_end) overlaps any EVE interval (0-based half-open)."""
    for s, e in eve_intervals:
        if chunk_start < e and chunk_end > s:
            return True
    return False


def _snapshot_timestamp_to_date(ts: str) -> str:
    """Extract YYYY-MM-DD from snapshot timestamp (ISO8601). Fallback to today if unparseable."""
    if not ts:
        return time.strftime("%Y-%m-%d", time.gmtime())
    part = ts.split("T")[0] if "T" in ts else ts[:10]
    if len(part) >= 10 and part[4] == "-" and part[7] == "-":
        return part
    return time.strftime("%Y-%m-%d", time.gmtime())


def run_build_viral_db(accessions_file: Path, output_dir: Path) -> Path:
    """Download all viral genomes from snapshot, concatenate to one FASTA, build BLAST DB.

    Output is written to a dated subfolder (snapshot date) under output_dir, e.g.
    output_dir/viral_db_YYYY-MM-DD/, so the viral DB is stamped with the same date
    as the snapshot. Use the resulting DB with blastn-filter --viral-db for proper
    prophage/EVE detection. Returns path to BLAST DB prefix.
    """
    from .download_genomes import (
        ACCESSIONS_KEY_ARCHAEA,
        ACCESSIONS_KEY_BACTERIAL,
        ACCESSIONS_KEY_PLASMID,
        ACCESSIONS_KEY_TIMESTAMP,
        ACCESSIONS_KEY_VIRAL,
        get_accession_lists_from_data,
        load_accessions,
        download_genomes,
    )

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    data = load_accessions(accessions_file)
    _bacterial, viral_ids, _archaea, _plasmid = get_accession_lists_from_data(data)
    if not viral_ids:
        raise ValueError(f"No viral accessions in {accessions_file}")

    ts = data.get(ACCESSIONS_KEY_TIMESTAMP, time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))
    date_str = _snapshot_timestamp_to_date(ts)
    work_dir = output_dir / f"viral_db_{date_str}"
    work_dir.mkdir(parents=True, exist_ok=True)

    viral_only = {
        ACCESSIONS_KEY_TIMESTAMP: ts,
        ACCESSIONS_KEY_BACTERIAL: [],
        ACCESSIONS_KEY_VIRAL: viral_ids,
        ACCESSIONS_KEY_ARCHAEA: [],
        ACCESSIONS_KEY_PLASMID: [],
    }
    temp_json = work_dir / "_viral_only.json"
    with temp_json.open("w") as f:
        json.dump(viral_only, f, indent=0)

    logger.info("Building viral reference DB (snapshot date %s): downloading %d viral genomes", date_str, len(viral_ids))
    print(f"Viral DB folder (snapshot date {date_str}): {work_dir}")
    print(f"Downloading {len(viral_ids):,} viral genomes from snapshot...")
    download_genomes(0, 0, work_dir, accessions_file=temp_json)
    temp_json.unlink(missing_ok=True)

    viral_ref_fasta = work_dir / "viral_ref.fasta"
    _concat_viral_fasta(work_dir, viral_ref_fasta)
    db_dir = work_dir / "blastn_db"
    db_prefix = build_viral_db(viral_ref_fasta, db_dir)
    manifest_path = work_dir / "viral_db_manifest.json"
    write_viral_db_manifest(
        db_prefix,
        manifest_path,
        snapshot_timestamp=ts,
        snapshot_date=date_str,
        source_snapshot=str(accessions_file),
        viral_accession_count=len(viral_ids),
    )
    logger.info("Viral BLAST DB built at %s (snapshot date %s)", db_prefix, date_str)
    print(f"Viral reference BLAST DB: {db_prefix}")
    print(f"Viral DB manifest: {manifest_path}")
    print(f"Viral DB aggregate SHA256: {viral_db_fingerprint(db_prefix)}")
    print(f"Use with: metagenome-generator blastn-filter --genome-dir <nonviral_dir> --out-dir <out> --viral-db {db_prefix}")
    return db_prefix
