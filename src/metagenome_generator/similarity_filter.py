#!/usr/bin/env python3
from __future__ import annotations
"""Filter sequences by pairwise similarity using BLASTN.

Used to ensure no two sequences in the final set exceed a similarity threshold
(e.g. 90% identity over sufficient alignment length). Sequences that are too
similar to already-kept sequences are removed; the dataset is then enriched
by generating more candidates until the target count is reached or max rounds.
"""

import logging
import subprocess
import tempfile
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger(__name__)

# BLAST outfmt 6 columns we use: qseqid, sseqid, pident, length
BLAST_SIM_COLS = "qseqid sseqid pident length"


def _build_blast_db(fasta_path: Path, db_dir: Path, *, prefix: str | None = None) -> Path:
    """Run makeblastdb on fasta_path. Returns path to DB prefix (no extension).

    The output prefix defaults to ``<fasta stem>_db`` inside ``db_dir`` so that
    rebuilding the DB from the same FASTA intentionally overwrites it (the
    within-call ``filter_by_similarity`` flow relies on this), while DBs built
    from different FASTAs in the same ``db_dir`` (e.g. first-batch self-BLAST
    vs. the growing kept-set DB, or concurrent callers sharing ``work_dir``)
    get distinct prefixes and no longer clobber each other. Pass ``prefix`` to
    override (e.g. for tests or bespoke pipelines).
    """
    db_dir.mkdir(parents=True, exist_ok=True)
    base = prefix if prefix is not None else f"{fasta_path.stem}_db"
    db_prefix = db_dir / base
    subprocess.run(
        [
            "makeblastdb",
            "-in", str(fasta_path),
            "-out", str(db_prefix),
            "-dbtype", "nucl",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    return db_prefix


def _run_blastn_similarity(
    query_fasta: Path,
    db_prefix: Path,
    out_tsv: Path,
    *,
    perc_identity: float = 90.0,
    max_target_seqs: int = 5,
    num_threads: int = 4,
    use_megablast: bool = True,
) -> None:
    """Run blastn query vs db; write tabular output (qseqid, sseqid, pident, length).

    use_megablast=True (default): -task megablast for much faster high-identity (≥90%%) search.
    num_threads: parallel threads for BLAST (default 4).
    """
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    fmt = "6 " + BLAST_SIM_COLS
    cmd = [
        "blastn",
        "-query", str(query_fasta),
        "-db", str(db_prefix),
        "-out", str(out_tsv),
        "-outfmt", fmt,
        "-perc_identity", str(perc_identity),
        "-max_target_seqs", str(max_target_seqs),
        "-num_threads", str(num_threads),
    ]
    if use_megablast:
        cmd.extend(["-task", "megablast"])
    subprocess.run(cmd, check=True, capture_output=True, text=True)


def _parse_similar_hits(
    tsv_path: Path,
    query_lengths: dict[str, int],
    pident_threshold: float,
    min_coverage: float,
) -> set[str]:
    """Parse BLAST tabular output; return set of qseqid that have a hit with pident>=threshold and length>=min_coverage*query_len."""
    similar_ids: set[str] = set()
    if not tsv_path.exists():
        return similar_ids
    with tsv_path.open() as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 4:
                continue
            qseqid, _sseqid, pident_str, length_str = parts[0], parts[1], parts[2], parts[3]
            try:
                pident = float(pident_str)
                length = int(length_str)
            except (ValueError, TypeError):
                continue
            qlen = query_lengths.get(qseqid)
            if qlen is None or qlen == 0:
                continue
            if pident >= pident_threshold and length >= min_coverage * qlen:
                similar_ids.add(qseqid)
    return similar_ids


def _parse_similar_neighbors(
    tsv_path: Path,
    query_lengths: dict[str, int],
    pident_threshold: float,
    min_coverage: float,
) -> dict[str, set[str]]:
    """Parse all-vs-all BLAST output into a symmetric adjacency dict.

    Used for within-batch deduplication when seeding the kept set: any pair (q, s)
    with pident >= threshold and alignment length >= min_coverage * qlen is
    considered a near-duplicate. Self-hits (q == s) are ignored.
    """
    neighbors: dict[str, set[str]] = {qid: set() for qid in query_lengths}
    if not tsv_path.exists():
        return neighbors
    with tsv_path.open() as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 4:
                continue
            qseqid, sseqid, pident_str, length_str = parts[0], parts[1], parts[2], parts[3]
            if qseqid == sseqid:
                continue
            try:
                pident = float(pident_str)
                length = int(length_str)
            except (ValueError, TypeError):
                continue
            qlen = query_lengths.get(qseqid)
            if qlen is None or qlen == 0:
                continue
            if pident >= pident_threshold and length >= min_coverage * qlen:
                neighbors.setdefault(qseqid, set()).add(sseqid)
                neighbors.setdefault(sseqid, set()).add(qseqid)
    return neighbors


def filter_by_similarity(
    records: list[SeqRecord],
    target_count: int,
    *,
    similarity_threshold: float = 90.0,
    min_coverage: float = 0.8,
    oversample_factor: float = 2.0,
    batch_size: int = 500,
    work_dir: Path | None = None,
    max_refill_rounds: int | None = None,
    num_threads: int = 4,
    use_megablast: bool = True,
) -> tuple[list[SeqRecord], dict]:
    """Filter records so no two kept sequences are >= similarity_threshold similar (BLASTN pident over min_coverage of query).

    Single-pass over ``records``: the first batch is deduplicated against itself,
    then each subsequent batch is BLAST'd against the growing ``kept`` set and
    only non-similar sequences are retained. Returns
    ``(filtered list of up to target_count records, stats dict)`` with keys
    ``generated``, ``removed``, ``kept``, ``rounds``, ``warning``.

    ``max_refill_rounds`` is accepted for backward compatibility only — the
    previous outer loop was dead code (candidates were cleared after the first
    round), so it never triggered. Refill is the caller's responsibility:
    generate more candidates and call :func:`filter_candidates_against_kept`.
    ``oversample_factor`` is likewise a caller-side hint and is not used here.
    """
    _ = max_refill_rounds, oversample_factor  # intentionally unused (see docstring)

    stats = {"generated": len(records), "removed": 0, "kept": 0, "rounds": 1, "warning": None}
    if not records:
        return [], stats
    if target_count <= 0:
        return [], stats

    kept: list[SeqRecord] = []
    candidates = list(records)

    use_dir = work_dir or Path(tempfile.mkdtemp(prefix="simfilter_"))
    use_dir.mkdir(parents=True, exist_ok=True)
    kept_fasta = use_dir / "kept.fasta"
    db_dir = use_dir / "db"
    db_dir.mkdir(parents=True, exist_ok=True)

    db_prefix: str | None = None
    removed_count = 0
    i = 0
    while i < len(candidates) and len(kept) < target_count:
        batch = candidates[i : i + batch_size]
        i += batch_size
        if not batch:
            break
        query_lengths = {rec.id: len(rec.seq) for rec in batch}
        batch_fasta = use_dir / "batch.fasta"
        SeqIO.write(batch, batch_fasta, "fasta")

        if db_prefix is None:
            # First batch: dedup within the batch before seeding `kept`. Previously
            # the entire first batch was appended unchecked, which allowed pairs of
            # near-duplicates from the same batch to both remain in the output
            # and violate the "no two kept sequences >= threshold similar" contract.
            db_prefix = _build_blast_db(batch_fasta, db_dir)
            self_tsv = use_dir / "batch_self.tsv"
            # For pairwise dedup we need every significant neighbor, not just the top
            # 5 (default). Cap at 500 to bound memory for very large first batches.
            self_max_targets = max(5, min(len(batch), 500))
            _run_blastn_similarity(
                batch_fasta,
                db_prefix,
                self_tsv,
                perc_identity=similarity_threshold,
                max_target_seqs=self_max_targets,
                num_threads=num_threads,
                use_megablast=use_megablast,
            )
            neighbors = _parse_similar_neighbors(
                self_tsv, query_lengths, similarity_threshold, min_coverage
            )
            kept_ids: set[str] = set()
            for rec in batch:
                if any(nid in kept_ids for nid in neighbors.get(rec.id, ())):
                    removed_count += 1
                    continue
                if len(kept) >= target_count:
                    break
                kept.append(rec)
                kept_ids.add(rec.id)
            # Rebuild DB from the deduplicated kept set for subsequent batches.
            if kept:
                SeqIO.write(kept, kept_fasta, "fasta")
                db_prefix = _build_blast_db(kept_fasta, db_dir)
            continue

        out_tsv = use_dir / "batch_vs_kept.tsv"
        _run_blastn_similarity(
            batch_fasta,
            db_prefix,
            out_tsv,
            perc_identity=similarity_threshold,
            max_target_seqs=5,
            num_threads=num_threads,
            use_megablast=use_megablast,
        )
        similar_ids = _parse_similar_hits(
            out_tsv, query_lengths, similarity_threshold, min_coverage
        )
        for rec in batch:
            if rec.id in similar_ids:
                removed_count += 1
                continue
            if len(kept) >= target_count:
                break
            kept.append(rec)
        if kept:
            SeqIO.write(kept, kept_fasta, "fasta")
            db_prefix = _build_blast_db(kept_fasta, db_dir)

    stats["removed"] = removed_count
    stats["kept"] = len(kept)
    if len(kept) < target_count:
        logger.warning(
            "Similarity filter: kept %d < target %d after single pass; caller may refill via filter_candidates_against_kept",
            len(kept), target_count,
        )
        stats["warning"] = (
            f"Dataset could not be fully created: kept {len(kept)} < target {target_count} "
            f"(not enough unique sequences after filtering at {similarity_threshold}% identity)."
        )

    result = kept[:target_count] if len(kept) >= target_count else kept
    stats["kept"] = len(result)
    return result, stats


def filter_candidates_against_kept(
    candidates: list[SeqRecord],
    kept: list[SeqRecord],
    *,
    similarity_threshold: float = 90.0,
    min_coverage: float = 0.8,
    batch_size: int = 2000,
    work_dir: Path | None = None,
    num_threads: int = 4,
    use_megablast: bool = True,
) -> list[SeqRecord]:
    """Return those candidates that are not >= similarity_threshold similar to any sequence in kept.
    Used for train-test split: remove from test any sequence similar to train.

    Optimizations: large batch_size (fewer BLAST runs), megablast task (faster for high identity), num_threads."""
    if not kept:
        return list(candidates)
    if not candidates:
        return []
    use_dir = work_dir or Path(tempfile.mkdtemp(prefix="simfilter_"))
    use_dir.mkdir(parents=True, exist_ok=True)
    kept_fasta = use_dir / "kept.fasta"
    SeqIO.write(kept, kept_fasta, "fasta")
    db_dir = use_dir / "db"
    db_prefix = _build_blast_db(kept_fasta, db_dir)
    passing: list[SeqRecord] = []
    for i in range(0, len(candidates), batch_size):
        batch = candidates[i : i + batch_size]
        query_lengths = {rec.id: len(rec.seq) for rec in batch}
        batch_fasta = use_dir / "batch.fasta"
        SeqIO.write(batch, batch_fasta, "fasta")
        out_tsv = use_dir / "batch_vs_kept.tsv"
        _run_blastn_similarity(
            batch_fasta, db_prefix, out_tsv,
            perc_identity=similarity_threshold, max_target_seqs=5,
            num_threads=num_threads, use_megablast=use_megablast,
        )
        similar_ids = _parse_similar_hits(
            out_tsv, query_lengths, similarity_threshold, min_coverage
        )
        for rec in batch:
            if rec.id not in similar_ids:
                passing.append(rec)
    return passing


_FASTQ_SUFFIXES = {".fastq", ".fq"}


def _format_from_suffix(path: Path) -> str:
    """Return 'fastq' if the file extension is .fastq/.fq, else 'fasta'.

    Matches the detection used elsewhere in the CLI (see
    split-metagenome-train-test) so FASTQ inputs are not silently mis-parsed.
    """
    return "fastq" if path.suffix.lower() in _FASTQ_SUFFIXES else "fasta"


def filter_test_against_train(
    train_fasta: Path,
    test_fasta: Path,
    output_fasta: Path,
    *,
    similarity_threshold: float = 90.0,
    min_coverage: float = 0.8,
    work_dir: Path | None = None,
    batch_size: int = 2000,
    num_threads: int = 4,
) -> tuple[int, int]:
    """Load train and test metagenomes; remove from test any read >= similarity_threshold similar to train; write filtered test.

    Input formats are inferred from file extensions: ``.fastq``/``.fq`` are read
    as FASTQ (Phred qualities preserved), anything else as FASTA. The output is
    written in the test input's format so per-base qualities survive the round
    trip. Previously both inputs were always parsed as FASTA, which silently
    produced empty results when FASTQ files were passed.

    Returns (n_removed, n_kept). Use for temporal split: after building train and test metagenomes separately,
    run this to drop test reads that are highly similar to train (e.g. different strains of same species).
    Requires BLAST+.
    """
    train_fmt = _format_from_suffix(train_fasta)
    test_fmt = _format_from_suffix(test_fasta)
    out_fmt = test_fmt  # Preserve Phred qualities when test input is FASTQ.

    train_records = list(SeqIO.parse(train_fasta, train_fmt))
    test_records = list(SeqIO.parse(test_fasta, test_fmt))
    output_fasta.parent.mkdir(parents=True, exist_ok=True)

    if not train_records:
        SeqIO.write(test_records, output_fasta, out_fmt)
        return 0, len(test_records)
    if not test_records:
        # Empty output; still write a well-formed empty file.
        output_fasta.write_text("")
        return 0, 0
    filtered = filter_candidates_against_kept(
        test_records,
        train_records,
        similarity_threshold=similarity_threshold,
        min_coverage=min_coverage,
        batch_size=batch_size,
        work_dir=work_dir,
        num_threads=num_threads,
        use_megablast=True,
    )
    n_removed = len(test_records) - len(filtered)
    SeqIO.write(filtered, output_fasta, out_fmt)
    return n_removed, len(filtered)
