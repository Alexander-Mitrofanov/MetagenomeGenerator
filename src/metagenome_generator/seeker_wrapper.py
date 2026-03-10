#!/usr/bin/env python3
from __future__ import annotations
"""Run Seeker (predict-metagenome) on a metagenome FASTA and write predictions + phage subset.

Seeker requires sequences of at least 200 nt. This module filters short reads if needed,
calls the CLI (via conda or direct), parses output, and writes TSV + phage FASTA.
"""

import subprocess
from pathlib import Path

from Bio import SeqIO

# Seeker's model asserts each sequence has at least this many bases.
SEEKER_MIN_LENGTH = 200


def run_seeker(
    input_path: Path,
    output_dir: Path,
    *,
    threshold: float = 0.5,
    conda_env: str | None = "seeker",
    predictions_tsv: str | None = None,
    phage_fasta: str | None = None,
    min_length: int = SEEKER_MIN_LENGTH,
) -> tuple[Path, Path]:
    """Run Seeker on input FASTA; write predictions TSV and phage FASTA. Return (predictions_path, phage_path)."""
    output_dir.mkdir(parents=True, exist_ok=True)
    if not input_path.exists():
        raise FileNotFoundError(f"Input does not exist: {input_path}")

    predictions_name = (
        predictions_tsv
        if predictions_tsv is not None
        else f"{input_path.stem}.seeker_predictions.tsv"
    )
    phage_fasta_name = (
        phage_fasta
        if phage_fasta is not None
        else f"{input_path.stem}.seeker_phage.fasta"
    )
    predictions_path = output_dir / predictions_name
    phage_fasta_path = output_dir / phage_fasta_name

    # Build input for Seeker: filter to reads >= min_length. Only write filtered file if some are dropped.
    filtered_input = input_path
    total = 0
    kept = 0
    for rec in SeqIO.parse(str(input_path), "fasta"):
        total += 1
        if len(rec.seq) >= min_length:
            kept += 1

    if kept < total:
        tmp_filtered = (
            output_dir / f"{input_path.stem}.seeker_filtered_min{min_length}.fasta"
        )
        with tmp_filtered.open("w") as out_handle:
            for rec in SeqIO.parse(str(input_path), "fasta"):
                if len(rec.seq) >= min_length:
                    SeqIO.write(rec, out_handle, "fasta")
        filtered_input = tmp_filtered
        print(
            f"Seeker input filter: kept {kept}/{total} reads "
            f"(min_length={min_length}) -> {filtered_input}"
        )
    elif total > 0:
        print(f"Seeker input: all {total} reads >= {min_length} nt, using input as-is.")

    cmd = (
        ["conda", "run", "-n", conda_env, "predict-metagenome", str(filtered_input)]
        if conda_env
        else ["predict-metagenome", str(filtered_input)]
    )
    proc = subprocess.run(cmd, text=True, capture_output=True)
    if proc.returncode != 0:
        details = (proc.stderr or proc.stdout or "").strip()
        raise RuntimeError(details or "Seeker failed")

    rows: list[tuple[str, str, float]] = []
    for line in proc.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) != 3:
            continue
        name, pred, score_str = parts
        if name == "name" and pred == "prediction":
            continue
        try:
            score = float(score_str)
        except ValueError:
            continue
        rows.append((name, pred, score))

    if not rows:
        raise RuntimeError("Seeker produced no parseable predictions")

    with predictions_path.open("w") as f:
        f.write("name\tprediction\tscore\n")
        for name, pred, score in rows:
            f.write(f"{name}\t{pred}\t{score}\n")

    phage_ids = {name for (name, _pred, score) in rows if score >= threshold}
    with phage_fasta_path.open("w") as out_handle:
        for rec in SeqIO.parse(str(filtered_input), "fasta"):
            if rec.id in phage_ids:
                SeqIO.write(rec, out_handle, "fasta")

    print(f"Seeker predictions TSV: {predictions_path}")
    print(f"Seeker phage FASTA (score >= {threshold}): {phage_fasta_path}")
    return predictions_path, phage_fasta_path
