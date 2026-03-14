#!/usr/bin/env python3
"""Structured benchmark recipe: fixed N per category, optional replicates.

Samples genomes from an accession snapshot (no NCBI search), runs download + chunk
for each replicate with a distinct seed so benchmarks are comparable and reproducible.
"""

from __future__ import annotations

import random
from pathlib import Path
from typing import Callable

from .download_genomes import (
    get_accession_lists_from_data,
    load_accessions,
    save_accessions,
)
from .chunk_genomes import build_metagenome
from .genome_layout import validate_genome_dir


def sample_accessions_from_snapshot(
    accessions_file: Path,
    n_bacterial: int,
    n_viral: int,
    n_archaea: int = 0,
    n_plasmid: int = 0,
    seed: int = 42,
) -> tuple[list[str], list[str], list[str], list[str]]:
    """Sample without replacement from each category in the snapshot.

    Returns (bacterial_ids, viral_ids, archaea_ids, plasmid_ids). If a category
    has fewer than the requested N, all are returned (no error).
    """
    data = load_accessions(accessions_file)
    bacterial, viral, archaea, plasmid = get_accession_lists_from_data(data)

    rng = random.Random(seed)

    def take(n: int, lst: list[str]) -> list[str]:
        if not lst or n <= 0:
            return []
        if n >= len(lst):
            return list(lst)
        return list(rng.sample(lst, n))

    return (
        take(n_bacterial, bacterial),
        take(n_viral, viral),
        take(n_archaea, archaea),
        take(n_plasmid, plasmid),
    )


def run_benchmark_recipe(
    accessions_file: Path,
    output_dir: Path,
    per_category: int,
    replicates: int = 1,
    *,
    seed: int = 42,
    n_archaea: int = 0,
    n_plasmid: int = 0,
    sequence_length: int = 250,
    reads_per_organism: int = 1000,
    output_fasta_name: str = "metagenome.fasta",
    progress_callback: Callable[[int, int, str], None] | None = None,
    output_fastq: bool = False,
) -> list[Path]:
    """Run a structured benchmark: per_category genomes per group, R replicates.

    For each replicate i (1..replicates):
      - Sample per_category bacterial, per_category viral (and n_archaea, n_plasmid)
        from the snapshot with seed = seed + i.
      - Write sampled accessions to output_dir/replicate_XXX/accessions_sampled.json.
      - Download those genomes into output_dir/replicate_XXX/downloaded/.
      - Chunk into output_dir/replicate_XXX/<output_fasta_name> with seed = seed + i.

    Returns the list of metagenome FASTA paths (one per replicate).
    progress_callback(replicate_index, total_replicates, message) is called for progress.
    """
    if per_category < 1:
        raise ValueError("per_category must be >= 1")
    if replicates < 1:
        raise ValueError("replicates must be >= 1")

    from .download_genomes import download_genomes

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    result_paths: list[Path] = []

    for i in range(replicates):
        rep_num = i + 1
        rep_name = f"replicate_{rep_num:03d}"
        rep_dir = output_dir / rep_name
        rep_dir.mkdir(parents=True, exist_ok=True)
        downloaded_dir = rep_dir / "downloaded"
        out_fasta = rep_dir / output_fasta_name

        if progress_callback:
            progress_callback(rep_num, replicates, f"Sampling replicate {rep_num}/{replicates}")

        rep_seed = seed + i
        bacterial, viral, archaea, plasmid = sample_accessions_from_snapshot(
            accessions_file,
            n_bacterial=per_category,
            n_viral=per_category,
            n_archaea=n_archaea,
            n_plasmid=n_plasmid,
            seed=rep_seed,
        )
        sampled_path = rep_dir / "accessions_sampled.json"
        save_accessions(sampled_path, bacterial, viral, archaea, plasmid)

        if progress_callback:
            progress_callback(rep_num, replicates, f"Downloading replicate {rep_num}/{replicates}")
        download_genomes(
            len(bacterial),
            len(viral),
            downloaded_dir,
            num_archaea=len(archaea),
            num_plasmid=len(plasmid),
            accessions_file=sampled_path,
        )

        if progress_callback:
            progress_callback(rep_num, replicates, f"Chunking replicate {rep_num}/{replicates}")
        ok, err = validate_genome_dir(downloaded_dir)
        if not ok:
            raise RuntimeError(f"Replicate {rep_num}: {err}")
        build_metagenome(
            downloaded_dir,
            out_fasta,
            sequence_length,
            reads_per_organism,
            seed=rep_seed,
            output_fastq=output_fastq,
        )
        written = out_fasta.with_suffix(".fastq") if output_fastq else out_fasta
        result_paths.append(written)

    return result_paths
