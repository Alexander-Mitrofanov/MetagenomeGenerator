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

import json
import logging
import subprocess
import time
from pathlib import Path

from .genome_layout import get_viral_fasta_paths, get_nonviral_fasta_paths

logger = logging.getLogger(__name__)

# BLAST+ outfmt 6: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
BLAST_COLS = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore".split()


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


def build_viral_db(viral_fasta: Path, db_dir: Path) -> Path:
    """Run makeblastdb on viral FASTA. Returns path to DB prefix (no extension)."""
    db_dir.mkdir(parents=True, exist_ok=True)
    db_prefix = db_dir / "viral_db"
    subprocess.run(
        ["makeblastdb", "-in", str(viral_fasta), "-out", str(db_prefix), "-dbtype", "nucl"],
        check=True,
        capture_output=True,
        text=True,
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
) -> None:
    """Run blastn (query = non-viral, db = viral). Writes tabular output to out_tsv."""
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
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
    ]
    subprocess.run(cmd, check=True, capture_output=True, text=True)


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
) -> dict[tuple[str, str], list[tuple[int, int]]]:
    """Build viral DB, run BLASTN for each non-viral FASTA, return EVE intervals per (file_stem, qseqid)."""
    out_dir.mkdir(parents=True, exist_ok=True)
    blast_dir = out_dir / "blastn"
    blast_dir.mkdir(parents=True, exist_ok=True)

    db_prefix = build_viral_db(viral_fasta, blast_dir)
    logger.info("BLASTN filter: viral DB built at %s", db_prefix)
    print(f"Viral BLAST DB: {db_prefix}")

    eve_intervals: dict[tuple[str, str], list[tuple[int, int]]] = {}
    for qf in nonviral_fasta_paths:
        out_tsv = blast_dir / f"{qf.stem}.blastn.tsv"
        logger.info("BLASTN filter: query=%s -> %s (evalue=%s, perc_identity=%s)",
                    qf.name, out_tsv.name, evalue, perc_identity)
        print(f"BLASTN: {qf.name} -> {out_tsv.name}")
        run_blastn(
            qf,
            db_prefix,
            out_tsv,
            evalue=evalue,
            perc_identity=perc_identity,
            max_target_seqs=max_target_seqs,
        )
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
    print(f"EVE intervals: {eve_json} ({len(eve_intervals)} sequences with hits)")
    return eve_intervals


def run_blastn_from_dirs(
    genome_dir: Path,
    out_dir: Path,
    *,
    evalue: float = 1e-5,
    perc_identity: float = 70.0,
    max_target_seqs: int = 5,
    viral_reference_fasta: Path | None = None,
    viral_db_prefix: Path | None = None,
) -> dict[tuple[str, str], list[tuple[int, int]]]:
    """Run BLASTN: non-viral (bacteria/, archaea/, plasmid/) vs viral DB.

    Viral DB source (first that is set wins):
    - viral_db_prefix: use this existing BLAST DB (e.g. from build-viral-db).
    - viral_reference_fasta: build DB from this FASTA (e.g. comprehensive RefSeq viral).
    - else: concatenate virus/ from genome_dir (current run's viral set only; limited for EVE detection).

    For proper prophage/EVE detection, use a comprehensive viral reference (build-viral-db from snapshot
    or a full RefSeq viral download) and pass --viral-reference-fasta or --viral-db.
    """
    nonviral = get_nonviral_fasta_paths(genome_dir)
    if not nonviral:
        raise FileNotFoundError(f"No non-viral genome FASTAs in {genome_dir} (bacteria/, archaea/, plasmid/ or flat)")

    out_dir.mkdir(parents=True, exist_ok=True)
    blast_dir = out_dir / "blastn"
    blast_dir.mkdir(parents=True, exist_ok=True)

    if viral_db_prefix is not None:
        db_prefix = Path(viral_db_prefix)
        nhr = db_prefix.with_suffix(".nhr") if db_prefix.suffix else db_prefix.parent / (db_prefix.name + ".nhr")
        if not nhr.exists():
            raise FileNotFoundError(f"BLAST DB not found at {viral_db_prefix} (expected {nhr})")
        logger.info("BLASTN filter: using existing viral DB at %s", db_prefix)
        print(f"Using existing viral BLAST DB: {db_prefix}")
    elif viral_reference_fasta is not None:
        vref = Path(viral_reference_fasta)
        if not vref.exists():
            raise FileNotFoundError(f"Viral reference FASTA not found: {vref}")
        db_prefix = build_viral_db(vref, blast_dir)
        logger.info("BLASTN filter: viral DB built from reference at %s", db_prefix)
        print(f"Viral BLAST DB built from reference: {vref}")
    else:
        viral_concat = out_dir / "viral_concat.fasta"
        _concat_viral_fasta(genome_dir, viral_concat)
        db_prefix = build_viral_db(viral_concat, blast_dir)
        logger.info("BLASTN filter: viral DB built from genome_dir virus/ at %s", db_prefix)
        print(f"Viral BLAST DB: {db_prefix} (from genome_dir virus/; for better EVE detection use --viral-reference-fasta or build-viral-db)")

    eve_intervals: dict[tuple[str, str], list[tuple[int, int]]] = {}
    for qf in nonviral:
        out_tsv = blast_dir / f"{qf.stem}.blastn.tsv"
        logger.info("BLASTN filter: query=%s -> %s (evalue=%s, perc_identity=%s)",
                    qf.name, out_tsv.name, evalue, perc_identity)
        print(f"BLASTN: {qf.name} -> {out_tsv.name}")
        run_blastn(
            qf,
            db_prefix,
            out_tsv,
            evalue=evalue,
            perc_identity=perc_identity,
            max_target_seqs=max_target_seqs,
        )
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
    print(f"EVE intervals: {eve_json} ({len(eve_intervals)} sequences with hits)")
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
    logger.info("Viral BLAST DB built at %s (snapshot date %s)", db_prefix, date_str)
    print(f"Viral reference BLAST DB: {db_prefix}")
    print(f"Use with: metagenome-generator blastn-filter --genome-dir <nonviral_dir> --out-dir <out> --viral-db {db_prefix}")
    return db_prefix
