#!/usr/bin/env python3
"""Structured benchmark recipe: fixed N per category, optional replicates.

Samples genomes from an accession snapshot (no NCBI search), runs download + chunk
for each replicate with a distinct seed so benchmarks are comparable and reproducible.
"""

from __future__ import annotations

import json
import random
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Callable

from Bio import SeqIO

from .download_genomes import (
    ACCESSIONS_KEY_ARCHAEA,
    ACCESSIONS_KEY_BACTERIAL,
    ACCESSIONS_KEY_PLASMID,
    ACCESSIONS_KEY_VIRAL,
    get_accession_lists_from_data,
    load_accessions,
    save_accessions,
)
from .chunk_genomes import build_metagenome, split_train_test_and_write
from .genome_layout import (
    ARCHAEA_DIR,
    BACTERIA_DIR,
    PLASMID_DIR,
    VIRUS_DIR,
    validate_genome_dir,
)


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


def _sample_ids_with_optional_overlap(
    rng: random.Random,
    ids: list[str],
    n: int,
    *,
    exclude: set[str] | None = None,
) -> list[str]:
    """Sample up to n ids, preferably excluding `exclude` when possible.

    If not enough ids are available without overlap, the remainder is sampled
    from the full pool (overlap becomes unavoidable).
    """
    if n <= 0:
        return []
    exclude = exclude or set()
    available = [x for x in ids if x not in exclude]
    if len(available) >= n:
        return list(rng.sample(available, n))
    chosen = list(available)
    remaining = n - len(chosen)
    if remaining <= 0:
        return chosen
    pool = [x for x in ids if x not in chosen]
    if not pool:
        return chosen
    chosen.extend(list(rng.sample(pool, min(remaining, len(pool)))))
    return chosen


def _write_concat_fasta_and_lengths(accessions: list[str], cache_cat_dir: Path, out_fasta: Path) -> dict[str, int]:
    """Concatenate cached per-accession FASTA files into a multi-FASTA and return qseqid->length."""
    out_fasta.parent.mkdir(parents=True, exist_ok=True)
    qlens: dict[str, int] = {}
    # Write atomically-ish: if the file exists, overwrite it.
    with out_fasta.open("w") as w:
        for acc in accessions:
            fp = cache_cat_dir / f"{acc}.fasta"
            if not fp.exists():
                raise FileNotFoundError(f"Cached genome FASTA not found: {fp}")
            rec = SeqIO.read(fp, "fasta")
            qid = rec.id.split()[0]
            qlens[qid] = len(rec.seq)
            # Append raw FASTA for speed (stream, not read_text).
            with fp.open("r") as src:
                shutil.copyfileobj(src, w)
    return qlens


def _build_blast_db_from_concat_fasta(fasta_path: Path, db_prefix: Path) -> None:
    """Create a nucleotide BLAST DB at `db_prefix` using makeblastdb."""
    db_prefix.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [
            "makeblastdb",
            "-in",
            str(fasta_path),
            "-out",
            str(db_prefix),
            "-dbtype",
            "nucl",
        ],
        check=True,
        capture_output=True,
        text=True,
    )


def _run_blastn_megablast(query_fasta: Path, db_prefix: Path, out_tsv: Path, *, perc_identity: float, num_threads: int) -> None:
    """Run blastn megablast query vs db; write outfmt6 with qseqid/sseqid/pident/length."""
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [
            "blastn",
            "-query",
            str(query_fasta),
            "-db",
            str(db_prefix),
            "-out",
            str(out_tsv),
            "-outfmt",
            "6 qseqid sseqid pident length",
            "-task",
            "megablast",
            "-perc_identity",
            str(perc_identity),
            "-max_target_seqs",
            "5",
            "-num_threads",
            str(num_threads),
        ],
        check=True,
        capture_output=True,
        text=True,
    )


def _count_similar_queries_in_blast_tsv(
    tsv_path: Path,
    qlens: dict[str, int],
    *,
    perc_identity: float,
    min_coverage: float,
) -> int:
    """Count how many query records have at least one acceptable BLAST hit."""
    if not tsv_path.exists():
        return 0
    similar: set[str] = set()
    with tsv_path.open() as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            qseqid, _sseqid, pident_str, length_str = parts[0], parts[1], parts[2], parts[3]
            try:
                pident = float(pident_str)
                alen = int(length_str)
            except ValueError:
                continue
            qlen = qlens.get(qseqid)
            if not qlen:
                continue
            if pident >= perc_identity and alen >= min_coverage * qlen:
                similar.add(qseqid)
    return len(similar)


def _ensure_cache_has_accessions(
    cache_root: Path,
    *,
    bacterial: list[str],
    viral: list[str],
    archaea: list[str],
    plasmid: list[str],
    sample_seed: int,
) -> None:
    """Download only missing genomes into the on-disk cache."""
    cache_root.mkdir(parents=True, exist_ok=True)
    cache_dirs = {
        "bacteria": cache_root / BACTERIA_DIR,
        "virus": cache_root / VIRUS_DIR,
        "archaea": cache_root / ARCHAEA_DIR,
        "plasmid": cache_root / PLASMID_DIR,
    }
    # Ensure base dirs exist so existence checks work consistently.
    for d in cache_dirs.values():
        d.mkdir(parents=True, exist_ok=True)

    def _missing(accessions: list[str], cat_dir: Path) -> list[str]:
        return [acc for acc in accessions if not (cat_dir / f"{acc}.fasta").exists()]

    missing_b = _missing(bacterial, cache_dirs["bacteria"])
    missing_v = _missing(viral, cache_dirs["virus"])
    missing_a = _missing(archaea, cache_dirs["archaea"])
    missing_p = _missing(plasmid, cache_dirs["plasmid"])

    if not (missing_b or missing_v or missing_a or missing_p):
        return

    # Download only missing accessions, preserving exact lists (no sampling).
    tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False)
    tmp_path = Path(tmp.name)
    try:
        json.dump(
            {
                ACCESSIONS_KEY_BACTERIAL: missing_b,
                ACCESSIONS_KEY_VIRAL: missing_v,
                ACCESSIONS_KEY_ARCHAEA: missing_a,
                ACCESSIONS_KEY_PLASMID: missing_p,
            },
            tmp,
            indent=2,
        )
    finally:
        tmp.close()

    from .download_genomes import download_genomes

    download_genomes(
        0,
        0,
        cache_root,
        accessions_file=tmp_path,
        sample_seed=sample_seed,
    )

    tmp_path.unlink(missing_ok=True)


def _copy_cached_genomes_to_downloaded(
    cache_root: Path,
    downloaded_root: Path,
    *,
    bacterial: list[str],
    viral: list[str],
    archaea: list[str],
    plasmid: list[str],
) -> None:
    """Copy cached per-accession FASTA files into replicate downloaded/."""
    downloaded_root.mkdir(parents=True, exist_ok=True)
    for subdir, ids in [
        (BACTERIA_DIR, bacterial),
        (VIRUS_DIR, viral),
        (ARCHAEA_DIR, archaea),
        (PLASMID_DIR, plasmid),
    ]:
        out_dir = downloaded_root / subdir
        out_dir.mkdir(parents=True, exist_ok=True)
        src_dir = cache_root / subdir
        for acc in ids:
            src = src_dir / f"{acc}.fasta"
            if not src.exists():
                raise FileNotFoundError(f"Cache missing for accession {acc}: {src}")
            shutil.copy2(src, out_dir / f"{acc}.fasta")


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
    train_test_split: float = 80.0,
    train_test_similarity_threshold: float = 90.0,
    similarity_min_coverage: float = 0.8,
    train_test_blast_threads: int = 4,
    train_test_blast_batch_size: int = 2000,
    diversity_max_attempts: int = 3,
    diversity_blast_perc_identity: float = 90.0,
    diversity_blast_min_coverage: float = 0.8,
    diversity_blast_threads: int = 4,
    progress_callback: Callable[[int, int, str], None] | None = None,
    output_fastq: bool = False,
) -> list[Path]:
    """Run a structured benchmark with maximally diverse replicates.

    For each replicate i (1..replicates):
    1) Greedily select a candidate genome accession set per category that minimizes
       genome-level BLAST similarity (megablast, perc_identity/min_coverage) against
       the genomes used in already-selected replicates.
    2) Download/cache the selected genomes.
    3) Generate reads from the selected genomes and split them into train/test
       using the existing BLAST-based train-vs-test similarity filtering.

    Returns the list of test FASTA/FASTQ paths (one per replicate).
    progress_callback(replicate_index, total_replicates, message) is called for progress.
    """
    if per_category < 1:
        raise ValueError("per_category must be >= 1")
    if replicates < 1:
        raise ValueError("replicates must be >= 1")
    if train_test_split <= 0 or train_test_split > 100:
        raise ValueError("--train-test-split must be in (0,100]")
    if diversity_max_attempts < 1:
        raise ValueError("--diversity-max-attempts must be >= 1")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    result_test_paths: list[Path] = []

    # Load snapshot once; we will sample repeatedly.
    snapshot = load_accessions(accessions_file)
    bacterial_all, viral_all, archaea_all, plasmid_all = get_accession_lists_from_data(snapshot)

    # On-disk cache to avoid repeated downloads across candidate attempts.
    cache_root = output_dir / ".genome_cache"
    cache_root.mkdir(parents=True, exist_ok=True)

    prev_selected_b = set()
    prev_selected_v = set()
    prev_selected_a = set()
    prev_selected_p = set()

    output_stem = Path(output_fasta_name).stem
    ext = "fastq" if output_fastq else "fasta"

    for i in range(replicates):
        rep_num = i + 1
        rep_name = f"replicate_{rep_num:03d}"
        rep_dir = output_dir / rep_name
        rep_dir.mkdir(parents=True, exist_ok=True)
        downloaded_dir = rep_dir / "downloaded"

        # If rerun, avoid mixing old files with new ones.
        shutil.rmtree(downloaded_dir, ignore_errors=True)

        work_dir = rep_dir / ".diversity_work"
        work_dir.mkdir(parents=True, exist_ok=True)

        if progress_callback:
            progress_callback(rep_num, replicates, f"Selecting diverse genome set ({rep_num}/{replicates})")

        rep_seed = seed + i
        best_candidate: tuple[list[str], list[str], list[str], list[str]] | None = None
        best_penalty: float = float("inf")

        # Build BLAST DBs from previously accepted replicate genomes (one per category).
        db_dir = work_dir / "blast_dbs"
        db_dir.mkdir(parents=True, exist_ok=True)

        def _build_prev_db(prev_ids: set[str], cat_subdir: Path, db_prefix: Path) -> Path | None:
            if not prev_ids:
                return None
            prev_fasta = db_prefix.with_suffix(".fasta")
            _write_concat_fasta_and_lengths(sorted(prev_ids), cat_subdir, prev_fasta)
            _build_blast_db_from_concat_fasta(prev_fasta, db_prefix)
            return db_prefix

        prev_db_b = _build_prev_db(prev_selected_b, cache_root / BACTERIA_DIR, db_dir / "prev_bacterial_db")
        prev_db_v = _build_prev_db(prev_selected_v, cache_root / VIRUS_DIR, db_dir / "prev_viral_db")
        prev_db_a = _build_prev_db(prev_selected_a, cache_root / ARCHAEA_DIR, db_dir / "prev_archaea_db")
        prev_db_p = _build_prev_db(prev_selected_p, cache_root / PLASMID_DIR, db_dir / "prev_plasmid_db")

        # Candidate attempts: sample different genome sets and pick the one with smallest penalty.
        for attempt in range(diversity_max_attempts):
            if progress_callback:
                progress_callback(rep_num, replicates, f"  Candidate {attempt + 1}/{diversity_max_attempts} scoring")

            cand_rng = random.Random(rep_seed + attempt * 10_000)
            exclude_b = prev_selected_b
            exclude_v = prev_selected_v
            exclude_a = prev_selected_a
            exclude_p = prev_selected_p

            cand_b = _sample_ids_with_optional_overlap(cand_rng, bacterial_all, per_category, exclude=exclude_b)
            cand_v = _sample_ids_with_optional_overlap(cand_rng, viral_all, per_category, exclude=exclude_v)
            cand_a = _sample_ids_with_optional_overlap(cand_rng, archaea_all, n_archaea, exclude=exclude_a)
            cand_p = _sample_ids_with_optional_overlap(cand_rng, plasmid_all, n_plasmid, exclude=exclude_p)

            # Ensure candidates exist in cache before scoring.
            _ensure_cache_has_accessions(
                cache_root,
                bacterial=cand_b,
                viral=cand_v,
                archaea=cand_a,
                plasmid=cand_p,
                sample_seed=seed,
            )

            penalty = 0

            # Score bacterial category.
            if prev_db_b is not None and cand_b:
                q_fasta = work_dir / f"cand_{attempt}_bacterial_query.fasta"
                qlens = _write_concat_fasta_and_lengths(cand_b, cache_root / BACTERIA_DIR, q_fasta)
                out_tsv = work_dir / f"cand_{attempt}_bacterial.tsv"
                _run_blastn_megablast(
                    q_fasta,
                    prev_db_b,
                    out_tsv,
                    perc_identity=diversity_blast_perc_identity,
                    num_threads=diversity_blast_threads,
                )
                penalty += _count_similar_queries_in_blast_tsv(
                    out_tsv,
                    qlens,
                    perc_identity=diversity_blast_perc_identity,
                    min_coverage=diversity_blast_min_coverage,
                )

            # Score viral category.
            if prev_db_v is not None and cand_v:
                q_fasta = work_dir / f"cand_{attempt}_viral_query.fasta"
                qlens = _write_concat_fasta_and_lengths(cand_v, cache_root / VIRUS_DIR, q_fasta)
                out_tsv = work_dir / f"cand_{attempt}_viral.tsv"
                _run_blastn_megablast(
                    q_fasta,
                    prev_db_v,
                    out_tsv,
                    perc_identity=diversity_blast_perc_identity,
                    num_threads=diversity_blast_threads,
                )
                penalty += _count_similar_queries_in_blast_tsv(
                    out_tsv,
                    qlens,
                    perc_identity=diversity_blast_perc_identity,
                    min_coverage=diversity_blast_min_coverage,
                )

            # Score archaea category.
            if prev_db_a is not None and cand_a:
                q_fasta = work_dir / f"cand_{attempt}_archaea_query.fasta"
                qlens = _write_concat_fasta_and_lengths(cand_a, cache_root / ARCHAEA_DIR, q_fasta)
                out_tsv = work_dir / f"cand_{attempt}_archaea.tsv"
                _run_blastn_megablast(
                    q_fasta,
                    prev_db_a,
                    out_tsv,
                    perc_identity=diversity_blast_perc_identity,
                    num_threads=diversity_blast_threads,
                )
                penalty += _count_similar_queries_in_blast_tsv(
                    out_tsv,
                    qlens,
                    perc_identity=diversity_blast_perc_identity,
                    min_coverage=diversity_blast_min_coverage,
                )

            # Score plasmid category.
            if prev_db_p is not None and cand_p:
                q_fasta = work_dir / f"cand_{attempt}_plasmid_query.fasta"
                qlens = _write_concat_fasta_and_lengths(cand_p, cache_root / PLASMID_DIR, q_fasta)
                out_tsv = work_dir / f"cand_{attempt}_plasmid.tsv"
                _run_blastn_megablast(
                    q_fasta,
                    prev_db_p,
                    out_tsv,
                    perc_identity=diversity_blast_perc_identity,
                    num_threads=diversity_blast_threads,
                )
                penalty += _count_similar_queries_in_blast_tsv(
                    out_tsv,
                    qlens,
                    perc_identity=diversity_blast_perc_identity,
                    min_coverage=diversity_blast_min_coverage,
                )

            if progress_callback:
                progress_callback(rep_num, replicates, f"  Candidate penalty={penalty}")

            if penalty < best_penalty:
                best_penalty = penalty
                best_candidate = (cand_b, cand_v, cand_a, cand_p)
                if best_penalty == 0:
                    break

        if best_candidate is None:
            raise RuntimeError("Failed to select a candidate genome set (unexpected)")

        bacterial_sel, viral_sel, archaea_sel, plasmid_sel = best_candidate

        sampled_path = rep_dir / "accessions_sampled.json"
        save_accessions(sampled_path, bacterial_sel, viral_sel, archaea_sel, plasmid_sel)

        if progress_callback:
            progress_callback(rep_num, replicates, f"Materializing genomes for replicate {rep_num}/{replicates}")

        _ensure_cache_has_accessions(
            cache_root,
            bacterial=bacterial_sel,
            viral=viral_sel,
            archaea=archaea_sel,
            plasmid=plasmid_sel,
            sample_seed=seed,
        )

        _copy_cached_genomes_to_downloaded(
            cache_root,
            downloaded_dir,
            bacterial=bacterial_sel,
            viral=viral_sel,
            archaea=archaea_sel,
            plasmid=plasmid_sel,
        )

        ok, err = validate_genome_dir(downloaded_dir)
        if not ok:
            raise RuntimeError(f"Replicate {rep_num}: {err}")

        if progress_callback:
            progress_callback(rep_num, replicates, f"Chunking and splitting reads for replicate {rep_num}/{replicates}")

        dummy_out = rep_dir / output_fasta_name
        _n, records = build_metagenome(
            downloaded_dir,
            dummy_out,
            sequence_length,
            reads_per_organism,
            seed=rep_seed,
            output_fastq=output_fastq,
            return_records=True,
        )

        split_train_test_and_write(
            records,
            train_test_split,
            rep_seed,
            rep_dir,
            output_stem=output_stem,
            similarity_threshold=train_test_similarity_threshold,
            similarity_min_coverage=similarity_min_coverage,
            work_dir=rep_dir / ".train_test_sim_work",
            blast_batch_size=train_test_blast_batch_size,
            blast_num_threads=train_test_blast_threads,
            write_fastq=output_fastq,
        )

        test_path = rep_dir / f"{output_stem}_test.{ext}"
        result_test_paths.append(test_path)

        # Update "previously accepted" genomes for diversity scoring in later replicates.
        prev_selected_b.update(bacterial_sel)
        prev_selected_v.update(viral_sel)
        prev_selected_a.update(archaea_sel)
        prev_selected_p.update(plasmid_sel)

    return result_test_paths
