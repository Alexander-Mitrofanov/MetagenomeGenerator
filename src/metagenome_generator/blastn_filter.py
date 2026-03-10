#!/usr/bin/env python3
from __future__ import annotations
"""BLASTN-based filtering of endogenous viral elements (EVEs) in non-viral sequences.

Aligns non-viral genomes (bacteria, archaea, plasmid) against a viral database,
parses hits, and produces per-sequence intervals to exclude when chunking so that
EVE regions are not used as negative training data.
"""

import json
import logging
import subprocess
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
) -> dict[tuple[str, str], list[tuple[int, int]]]:
    """Run BLASTN: non-viral (bacteria/, archaea/, plasmid/ or flat) vs viral (virus/ or viral_*.fasta).
    Builds viral DB from concatenated viral FASTAs. Returns EVE intervals and writes eve_intervals.json."""
    viral_concat = out_dir / "viral_concat.fasta"
    _concat_viral_fasta(genome_dir, viral_concat)
    nonviral = get_nonviral_fasta_paths(genome_dir)
    if not nonviral:
        raise FileNotFoundError(f"No non-viral genome FASTAs in {genome_dir} (bacteria/, archaea/, plasmid/ or flat)")
    return run_blastn_nonviral(
        nonviral,
        viral_concat,
        out_dir,
        evalue=evalue,
        perc_identity=perc_identity,
        max_target_seqs=max_target_seqs,
    )


def load_eve_intervals(eve_json: Path) -> dict[tuple[str, str], list[tuple[int, int]]]:
    """Load EVE intervals from run_blastn_nonviral output (eve_intervals.json)."""
    with eve_json.open() as f:
        data = json.load(f)
    return {tuple(k.split("\t")): [tuple(iv) for iv in v] for k, v in data.items()}


def chunk_overlaps_eve(chunk_start: int, chunk_end: int, eve_intervals: list[tuple[int, int]]) -> bool:
    """True if [chunk_start, chunk_end) overlaps any EVE interval (0-based half-open)."""
    for s, e in eve_intervals:
        if chunk_start < e and chunk_end > s:
            return True
    return False
