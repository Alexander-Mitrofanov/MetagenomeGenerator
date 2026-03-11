#!/usr/bin/env python3
from __future__ import annotations
"""Chunk genome FASTA file(s) into fixed- or variable-length reads and write a metagenome FASTA.

Supports task-specific training data: fixed-length (e.g. 250 nt) or variable-length
contigs (e.g. 300–2000 bp uniform) to simulate metagenomic data.
Optional similarity filtering: oversample, then keep only sequences not ≥90% similar to already-kept; refill if needed.
Exposes `build_metagenome` for programmatic use and a CLI when run directly.
"""

import argparse
import logging
import random
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .genome_layout import (
    ARCHAEA_PREFIX,
    BACTERIA_PREFIX,
    PLASMID_PREFIX,
    VIRUS_PREFIX,
    iter_genome_fastas,
)

logger = logging.getLogger(__name__)

# Category name from file prefix stem (e.g. bacterial_1 -> bacterial)
_PREFIX_TO_CATEGORY: list[tuple[str, str]] = [
    ("bacterial", BACTERIA_PREFIX),
    ("viral", VIRUS_PREFIX),
    ("archaea", ARCHAEA_PREFIX),
    ("plasmid", PLASMID_PREFIX),
]


def _category_from_prefix(prefix: str) -> str:
    """Return category name for a genome file prefix (bacterial_1 -> bacterial)."""
    for category, p in _PREFIX_TO_CATEGORY:
        if prefix.startswith(p) or prefix == p.rstrip("_"):
            return category
    return "other"


def _compute_read_limits(
    prefixes_and_files: list[tuple[str, Path]],
    base_reads_per_file: int,
    *,
    abundance_profile: dict[str, float] | None = None,
    abundance_distribution: str | None = None,
    seed: int | None = None,
) -> list[int]:
    """Return one read limit per file. If no abundance, all limits = base_reads_per_file."""
    n = len(prefixes_and_files)
    if n == 0:
        return []
    if abundance_profile:
        limits = []
        for prefix, _ in prefixes_and_files:
            cat = _category_from_prefix(prefix)
            w = abundance_profile.get(cat, 1.0)
            limits.append(max(1, int(base_reads_per_file * w)))
        return limits
    if abundance_distribution == "exponential":
        rng = random.Random(seed)
        weights = [rng.expovariate(1.0) for _ in range(n)]
        total = sum(weights)
        if total <= 0:
            return [base_reads_per_file] * n
        scale = n / total
        limits = [max(1, int(base_reads_per_file * w * scale)) for w in weights]
        return limits
    return [base_reads_per_file] * n


BASES = "ACGT"


def _apply_mutations(
    seq_str: str,
    substitution_rate: float,
    indel_rate: float,
    rng: random.Random,
) -> str:
    """Apply substitutions and optionally indels to a sequence (ACGT only).

    substitution_rate: per-base probability of substituting to a random different base.
    indel_rate: per-base probability of an indel (50% insert random base, 50% delete).
    With indels, output length can change. Uses rng for reproducibility.
    """
    if substitution_rate <= 0 and indel_rate <= 0:
        return seq_str
    seq_str = seq_str.upper()
    out: list[str] = []
    for c in seq_str:
        if c not in BASES:
            out.append(c)
            continue
        # Indel (before the base)
        if indel_rate > 0 and rng.random() < indel_rate:
            if rng.random() < 0.5:
                out.append(rng.choice(BASES))
            else:
                continue  # delete this base
        # Substitution
        if substitution_rate > 0 and rng.random() < substitution_rate:
            others = [b for b in BASES if b != c]
            out.append(rng.choice(others) if others else c)
        else:
            out.append(c)
    return "".join(out)


def _apply_mutations_to_record(
    rec: SeqRecord,
    substitution_rate: float,
    indel_rate: float,
    rng: random.Random,
) -> SeqRecord:
    """Return a new SeqRecord with mutations applied to the sequence. Preserves id and description."""
    new_seq = _apply_mutations(str(rec.seq), substitution_rate, indel_rate, rng)
    return SeqRecord(Seq(new_seq), id=rec.id, description=rec.description)


def chunk_sequence(record, prefix: str, chunk_size: int, yield_coords: bool = False):
    """Split a sequence into non-overlapping fixed-size chunks.

    Only full-length chunks are emitted (i.e. trailing remainder is dropped).
    If yield_coords, yields (rec, start, end) for EVE overlap checks; else yields rec.
    """
    seq = record.seq
    idx = 0
    for i in range(0, len(seq) - chunk_size + 1, chunk_size):
        sub = seq[i : i + chunk_size]
        rec = SeqRecord(
            sub,
            id=f"{prefix}_chunk_{idx}",
            description=f"len={len(sub)} from {record.id}",
        )
        if yield_coords:
            yield rec, i, i + chunk_size
        else:
            yield rec
        idx += 1


def chunk_sequence_variable(
    record,
    prefix: str,
    min_len: int,
    max_len: int,
    reads_per_organism: int | None,
    rng: random.Random | None = None,
    yield_coords: bool = False,
):
    """Split a sequence into non-overlapping contigs of random length in [min_len, max_len] (uniform).

    Simulates variable-length metagenomic contigs (e.g. 300–2000 bp).
    If yield_coords, yields (rec, start, end); else yields rec.
    """
    rng = rng or random.Random()
    seq = record.seq
    idx = 0
    pos = 0
    while pos < len(seq):
        length = rng.randint(min_len, max_len)
        end = min(pos + length, len(seq))
        if end - pos < min_len:
            break
        sub = seq[pos:end]
        rec = SeqRecord(
            sub,
            id=f"{prefix}_contig_{idx}",
            description=f"len={len(sub)} from {record.id}",
        )
        if yield_coords:
            yield rec, pos, end
        else:
            yield rec
        idx += 1
        pos = end
        if reads_per_organism is not None and idx >= reads_per_organism:
            break


def get_file_stats(
    input_path: Path,
    sequence_length: int,
    *,
    min_length: int | None = None,
    max_length: int | None = None,
) -> list[tuple[str, Path, int, int]]:
    """Scan input (file or dir of *.fasta) and return (prefix, path, total_bases, max_reads) per file.

    For fixed-length: max_reads = total_bases // sequence_length.
    For variable-length (min_length and max_length set): uses mean length for estimate.
    """
    if input_path.is_file():
        files = [(input_path.stem, input_path)]
    else:
        files = iter_genome_fastas(input_path)

    effective_len = sequence_length
    if min_length is not None and max_length is not None:
        effective_len = (min_length + max_length) // 2

    result: list[tuple[str, Path, int, int]] = []
    for prefix, fp in files:
        total_bases = sum(len(rec.seq) for rec in SeqIO.parse(fp, "fasta"))
        max_reads = total_bases // effective_len if total_bases >= effective_len else 0
        result.append((prefix, fp, total_bases, max_reads))
    return result


def _is_allowed_sequence(rec: SeqRecord, allow_ambiguous: bool) -> bool:
    """True if sequence is acceptable under the ambiguity policy.

    When allow_ambiguous is False, only A/C/G/T (case-insensitive) are allowed;
    any other character (e.g. N, R, Y, etc.) causes the sequence to be rejected.
    """
    if allow_ambiguous:
        return True
    seq_str = str(rec.seq).upper()
    return all(ch in "ACGT" for ch in seq_str)


def _collect_chunks_for_file(
    fp: Path,
    sequence_length: int,
    reads_per_organism: int | None,
    prefix: str,
    *,
    min_length: int | None = None,
    max_length: int | None = None,
    seed: int | None = None,
    eve_intervals: dict[tuple[str, str], list[tuple[int, int]]] | None = None,
    allow_ambiguous: bool = True,
    substitution_rate: float = 0.0,
    indel_rate: float = 0.0,
    mutation_rng: random.Random | None = None,
) -> list[SeqRecord]:
    """Collect chunks for a single FASTA file (fixed- or variable-length).

    Skips chunks overlapping EVE if eve_intervals is given.
    If substitution_rate or indel_rate > 0, applies mutations to each chunk (use mutation_rng for reproducibility).
    """
    from .blastn_filter import chunk_overlaps_eve

    chunks: list[SeqRecord] = []
    rng = random.Random(seed) if seed is not None else None
    use_eve = eve_intervals is not None
    apply_mutations = (substitution_rate > 0 or indel_rate > 0) and mutation_rng is not None

    for record in SeqIO.parse(fp, "fasta"):
        key = (prefix, record.id)
        intervals = eve_intervals.get(key, []) if eve_intervals else []

        if min_length is not None and max_length is not None:
            for item in chunk_sequence_variable(
                record, prefix, min_length, max_length, reads_per_organism, rng=rng, yield_coords=use_eve
            ):
                if use_eve:
                    rec, start, end = item
                    if chunk_overlaps_eve(start, end, intervals):
                        continue
                else:
                    rec = item
                if not _is_allowed_sequence(rec, allow_ambiguous):
                    continue
                if apply_mutations:
                    rec = _apply_mutations_to_record(rec, substitution_rate, indel_rate, mutation_rng)
                chunks.append(rec)
        else:
            for i, item in enumerate(chunk_sequence(record, prefix, sequence_length, yield_coords=use_eve)):
                if reads_per_organism is not None and i >= reads_per_organism:
                    break
                if use_eve:
                    rec, start, end = item
                    if chunk_overlaps_eve(start, end, intervals):
                        continue
                else:
                    rec = item
                if not _is_allowed_sequence(rec, allow_ambiguous):
                    continue
                if apply_mutations:
                    rec = _apply_mutations_to_record(rec, substitution_rate, indel_rate, mutation_rng)
                chunks.append(rec)
        break
    return chunks


def _collect_chunks_from_multirecord_fasta(
    fp: Path,
    sequence_length: int,
    reads_per_organism: int | None,
    prefix_base: str,
    *,
    min_length: int | None = None,
    max_length: int | None = None,
    seed: int | None = None,
    allow_ambiguous: bool = True,
    substitution_rate: float = 0.0,
    indel_rate: float = 0.0,
) -> list[SeqRecord]:
    """Collect chunks from a multi-record FASTA (e.g. metavirome contigs). No EVE filtering.

    Each record is chunked with prefix prefix_base_0, prefix_base_1, ... Same length/balance
    and mutation logic as main input; use seed for reproducibility.
    """
    chunks: list[SeqRecord] = []
    use_mutations = substitution_rate > 0 or indel_rate > 0
    for rec_idx, record in enumerate(SeqIO.parse(fp, "fasta")):
        prefix = f"{prefix_base}_{rec_idx}"
        file_seed = (seed + 20000 + rec_idx) if seed is not None else None
        rng = random.Random(file_seed) if file_seed is not None else None
        file_mutation_rng = random.Random(seed + 30000 + rec_idx) if use_mutations else None
        apply_mutations = use_mutations and file_mutation_rng is not None
        if min_length is not None and max_length is not None:
            for item in chunk_sequence_variable(
                record, prefix, min_length, max_length, reads_per_organism, rng=rng, yield_coords=False
            ):
                rec = item
                if not _is_allowed_sequence(rec, allow_ambiguous):
                    continue
                if apply_mutations:
                    rec = _apply_mutations_to_record(rec, substitution_rate, indel_rate, file_mutation_rng)
                chunks.append(rec)
        else:
            for i, rec in enumerate(chunk_sequence(record, prefix, sequence_length, yield_coords=False)):
                if reads_per_organism is not None and i >= reads_per_organism:
                    break
                if not _is_allowed_sequence(rec, allow_ambiguous):
                    continue
                if apply_mutations:
                    rec = _apply_mutations_to_record(rec, substitution_rate, indel_rate, file_mutation_rng)
                chunks.append(rec)
    return chunks


def build_metagenome(
    input_path: Path,
    out_path: Path,
    sequence_length: int,
    reads_per_organism: int | None,
    *,
    min_length: int | None = None,
    max_length: int | None = None,
    seed: int | None = None,
    cap_total_reads: int | None = None,
    eve_intervals: dict[tuple[str, str], list[tuple[int, int]]] | None = None,
    allow_ambiguous: bool = True,
    substitution_rate: float = 0.0,
    indel_rate: float = 0.0,
    filter_similar: bool = False,
    similarity_threshold: float = 90.0,
    similarity_min_coverage: float = 0.8,
    oversample_factor: float = 2.0,
    similarity_work_dir: Path | None = None,
    max_refill_rounds: int = 3,
    return_records: bool = False,
    extra_viral_fasta: Path | None = None,
    abundance_profile: dict[str, float] | None = None,
    abundance_distribution: str | None = None,
) -> int | tuple[int, list[SeqRecord]]:
    """Build a metagenome FASTA from input_path. Fixed-length or variable-length (min_length–max_length) contigs.

    If cap_total_reads is set, randomly downsample to that many reads (for balancing negative to positive).
    If eve_intervals is set, chunks overlapping those intervals (EVE regions) are excluded.
    If substitution_rate or indel_rate > 0, applies mutations to each chunk (seed used for reproducibility; default 42).
    If filter_similar is True: generate more chunks than needed (oversample), filter out sequences that are
    >= similarity_threshold (default 90%%) similar to any already-kept sequence, then refill until target or max rounds.
    If return_records is True, do not write to out_path and return (count, list[SeqRecord]) for train-test split.
    If extra_viral_fasta is set, chunks from that FASTA (multi-record, e.g. metavirome contigs) are merged in
    with prefix extra_viral_0, extra_viral_1, ... using the same length and mutation settings.
    If abundance_profile is set (e.g. {"bacterial": 0.5, "viral": 2.0}), scale reads per file by category.
    If abundance_distribution == "exponential", assign per-genome weights from Exp(1) (use seed for reproducibility).
    """
    from .similarity_filter import filter_by_similarity, filter_candidates_against_kept

    use_mutations = substitution_rate > 0 or indel_rate > 0
    if extra_viral_fasta is not None and not extra_viral_fasta.is_file():
        raise FileNotFoundError(f"extra_viral_fasta not found or not a file: {extra_viral_fasta}")
    if use_mutations and seed is None:
        seed = 42
    mutation_rng = random.Random(seed) if use_mutations else None

    if input_path.is_file():
        prefixes_and_files = [(input_path.stem, input_path)]
    else:
        prefixes_and_files = iter_genome_fastas(input_path)
    num_files = len(prefixes_and_files)

    target_count: int | None = None
    if cap_total_reads is not None:
        target_count = cap_total_reads
    elif reads_per_organism is not None and num_files > 0:
        target_count = reads_per_organism * num_files

    # When filter_similar we may oversample: ask for more per file when not capped (e.g. not balanced)
    reads_per_organism_gen = reads_per_organism
    if filter_similar and target_count is not None and target_count > 0:
        oversample_per_file = max(1, int(target_count * oversample_factor) // num_files)
        reads_per_organism_gen = max(reads_per_organism or 0, oversample_per_file)
        if reads_per_organism_gen == 0:
            reads_per_organism_gen = reads_per_organism

    if reads_per_organism_gen is None:
        read_limits: list[int | None] = [None] * num_files
    else:
        read_limits = _compute_read_limits(
            prefixes_and_files,
            reads_per_organism_gen,
            abundance_profile=abundance_profile,
            abundance_distribution=abundance_distribution,
            seed=seed,
        )

    all_chunks: list[SeqRecord] = []
    for i, (prefix, fp) in enumerate(prefixes_and_files):
        file_seed = (seed + i) if seed is not None else None
        file_mutation_rng = random.Random(seed + 10000 + i) if mutation_rng is not None else None
        per_file_limit = read_limits[i]
        chunks = _collect_chunks_for_file(
            fp,
            sequence_length,
            per_file_limit,
            prefix,
            min_length=min_length,
            max_length=max_length,
            seed=file_seed,
            eve_intervals=eve_intervals,
            allow_ambiguous=allow_ambiguous,
            substitution_rate=substitution_rate,
            indel_rate=indel_rate,
            mutation_rng=file_mutation_rng,
        )
        all_chunks.extend(chunks)

    if extra_viral_fasta is not None:
        extra_seed = seed if seed is not None else 42
        extra_chunks = _collect_chunks_from_multirecord_fasta(
            extra_viral_fasta,
            sequence_length,
            reads_per_organism_gen,
            "extra_viral",
            min_length=min_length,
            max_length=max_length,
            seed=extra_seed,
            allow_ambiguous=allow_ambiguous,
            substitution_rate=substitution_rate,
            indel_rate=indel_rate,
        )
        all_chunks.extend(extra_chunks)
        logger.info("Extra viral FASTA: added %d chunks from %s", len(extra_chunks), extra_viral_fasta)

    if filter_similar and target_count is not None and target_count > 0:
        work_dir = similarity_work_dir or (out_path.parent / ".simfilter_work")
        kept, sim_stats = filter_by_similarity(
            all_chunks,
            target_count,
            similarity_threshold=similarity_threshold,
            min_coverage=similarity_min_coverage,
            oversample_factor=oversample_factor,
            work_dir=work_dir,
            max_refill_rounds=max_refill_rounds,
        )
        logger.info(
            "Similarity filter: generated=%d, removed=%d, kept=%d, rounds=%d",
            sim_stats.get("generated", 0), sim_stats.get("removed", 0),
            sim_stats.get("kept", 0), sim_stats.get("rounds", 1),
        )
        if sim_stats.get("warning"):
            logger.warning("Similarity filter: %s", sim_stats["warning"])
        # Refill: if we have fewer than target, generate more chunks (different seed) and filter against kept
        rng = random.Random(seed)
        refill_round = 0
        while len(kept) < target_count and refill_round < max_refill_rounds:
            refill_round += 1
            refill_seed = (seed + 1000 + refill_round) if seed is not None else None
            more_chunks: list[SeqRecord] = []
            for i, (prefix, fp) in enumerate(prefixes_and_files):
                file_seed = (refill_seed + i) if refill_seed is not None else None
                file_mutation_rng = random.Random(refill_seed + 10000 + i) if mutation_rng is not None else None
                chunks = _collect_chunks_for_file(
                    fp,
                    sequence_length,
                    read_limits[i],
                    prefix,
                    min_length=min_length,
                    max_length=max_length,
                    seed=file_seed,
                    eve_intervals=eve_intervals,
                    allow_ambiguous=allow_ambiguous,
                    substitution_rate=substitution_rate,
                    indel_rate=indel_rate,
                    mutation_rng=file_mutation_rng,
                )
                more_chunks.extend(chunks)
            if not more_chunks:
                break
            additions = filter_candidates_against_kept(
                more_chunks, kept,
                similarity_threshold=similarity_threshold,
                min_coverage=similarity_min_coverage,
                work_dir=work_dir,
            )
            kept.extend(additions)
            if len(kept) >= target_count:
                break
            if not additions:
                logger.warning(
                    "Similarity refill round %d: no new unique sequences; kept %d < target %d",
                    refill_round, len(kept), target_count,
                )
                break
        all_chunks = kept[:target_count] if len(kept) >= target_count else kept
        if len(all_chunks) < target_count:
            logger.warning(
                "Dataset could not be fully created: %d sequences after filtering (target %d). Consider relaxing similarity threshold or adding more source genomes.",
                len(all_chunks), target_count,
            )
    else:
        if cap_total_reads is not None and len(all_chunks) > cap_total_reads:
            rng = random.Random(seed)
            all_chunks = rng.sample(all_chunks, cap_total_reads)

    if return_records:
        return len(all_chunks), all_chunks
    count = SeqIO.write(all_chunks, out_path, "fasta")
    return int(count)


def split_train_test_and_write(
    records: list[SeqRecord],
    train_pct: float,
    seed: int | None,
    output_dir: Path,
    output_stem: str,
    *,
    similarity_threshold: float = 90.0,
    similarity_min_coverage: float = 0.8,
    work_dir: Path | None = None,
    blast_batch_size: int = 2000,
    blast_num_threads: int = 4,
) -> tuple[int, int]:
    """Split records into train (train_pct%%) and test; remove from test any sequence similar to train; write train and test FASTAs.

    Similarity check: BLAST (megablast task, faster for high identity) of test vs train DB in batches;
    test sequences with a hit >= similarity_threshold over min_coverage of length are dropped.
    Returns (n_train, n_test_after_filter). Output files: {output_stem}_train.fasta, {output_stem}_test.fasta in output_dir.
    """
    from .similarity_filter import filter_candidates_against_kept

    if not records:
        return 0, 0
    train_pct = max(0.0, min(100.0, train_pct))
    n_train_target = max(1, int(len(records) * train_pct / 100.0))
    n_test_target = len(records) - n_train_target
    if n_train_target <= 0 or n_test_target <= 0:
        logger.warning("Train-test split: not enough sequences (%d); writing all to train", len(records))
        train_path = output_dir / f"{output_stem}_train.fasta"
        SeqIO.write(records, train_path, "fasta")
        (output_dir / f"{output_stem}_test.fasta").write_text("")
        return len(records), 0

    rng = random.Random(seed)
    shuffled = list(records)
    rng.shuffle(shuffled)
    train_list = shuffled[:n_train_target]
    test_list = shuffled[n_train_target:]

    work_dir = work_dir or (output_dir / ".train_test_sim_work")
    work_dir.mkdir(parents=True, exist_ok=True)
    test_filtered = filter_candidates_against_kept(
        test_list,
        train_list,
        similarity_threshold=similarity_threshold,
        min_coverage=similarity_min_coverage,
        work_dir=work_dir,
        batch_size=blast_batch_size,
        num_threads=blast_num_threads,
        use_megablast=True,
    )
    n_removed = len(test_list) - len(test_filtered)
    if n_removed:
        logger.info(
            "Train-test split: removed %d sequences from test (similar to train at >=%.1f%%); test size %d -> %d",
            n_removed, similarity_threshold, len(test_list), len(test_filtered),
        )

    output_dir.mkdir(parents=True, exist_ok=True)
    train_path = output_dir / f"{output_stem}_train.fasta"
    test_path = output_dir / f"{output_stem}_test.fasta"
    SeqIO.write(train_list, train_path, "fasta")
    SeqIO.write(test_filtered, test_path, "fasta")
    return len(train_list), len(test_filtered)


def _cli(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="Chunk genome FASTAs into fixed-length reads and write a metagenome FASTA."
    )
    parser.add_argument(
        "--input",
        type=Path,
        required=True,
        help="Input FASTA file or directory of genome FASTAs (e.g. output from download_genomes.py)",
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Output filename for the metagenome FASTA (e.g. metagenome_250nt.fasta)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("output"),
        help="Directory to write the output file. Default: output",
    )
    parser.add_argument(
        "--sequence-length",
        type=int,
        default=250,
        help="Fixed read length in nucleotides (used when not using variable-length). Default: 250",
    )
    parser.add_argument(
        "--min-contig-length",
        type=int,
        default=None,
        help="Min contig length for variable-length mode (with --max-contig-length). e.g. 300",
    )
    parser.add_argument(
        "--max-contig-length",
        type=int,
        default=None,
        help="Max contig length for variable-length mode (with --min-contig-length). e.g. 2000",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed for variable-length chunking (reproducibility).",
    )
    parser.add_argument(
        "--reads-per-organism",
        type=int,
        default=None,
        help="Max reads to extract per organism (per input file). Default: all reads",
    )
    parser.add_argument(
        "--balanced",
        action="store_true",
        help="Use equal reads per file: set reads-per-organism to the minimum max reads across all files (so each file contributes the same number of reads).",
    )
    parser.add_argument(
        "--cap-total-reads",
        type=int,
        default=None,
        help="Downsample to at most this many reads (for balancing e.g. non-viral to match viral count). Use --seed for reproducibility.",
    )
    parser.add_argument(
        "--eve-intervals",
        type=Path,
        default=None,
        help="Path to eve_intervals.json from blastn-filter step; chunks overlapping EVE regions are excluded.",
    )
    parser.add_argument(
        "--filter-similar",
        action="store_true",
        help="Filter out chunks that are >=90%% similar to any already-kept chunk; oversample then refill to reach target.",
    )
    parser.add_argument(
        "--similarity-threshold",
        type=float,
        default=90.0,
        help="Max allowed similarity (BLASTN pident); chunks above this vs. kept set are removed. Default: 90",
    )
    parser.add_argument(
        "--similarity-min-coverage",
        type=float,
        default=0.8,
        help="Min fraction of query length that must align for a hit to count as similar. Default: 0.8",
    )
    parser.add_argument(
        "--oversample-factor",
        type=float,
        default=2.0,
        help="When --filter-similar: generate up to this many times target count before filtering. Default: 2.0",
    )
    parser.add_argument(
        "--forbid-ambiguous",
        action="store_true",
        help="Discard chunks containing ambiguous bases (non-ACGT, e.g. N). By default, such chunks are kept.",
    )
    parser.add_argument(
        "--train-test-split",
        type=float,
        default=None,
        metavar="PCT",
        help="Split output into train and test: PCT%% train, rest test (e.g. 80). Test sequences similar to train are removed. Writes {output_stem}_train.fasta and {output_stem}_test.fasta.",
    )
    parser.add_argument(
        "--train-test-similarity-threshold",
        type=float,
        default=90.0,
        help="Max BLASTN percent identity for train-test: remove from test if similar to train above this. Default: 90",
    )
    args = parser.parse_args(argv)

    if not args.input.exists():
        parser.error(f"--input does not exist: {args.input}")

    eve_intervals = None
    if args.eve_intervals is not None:
        if not args.eve_intervals.exists():
            parser.error(f"--eve-intervals file not found: {args.eve_intervals}")
        from .blastn_filter import load_eve_intervals
        eve_intervals = load_eve_intervals(args.eve_intervals)
        print(f"Loaded EVE intervals for {len(eve_intervals)} sequences (excluding overlapping chunks).")

    use_variable = args.min_contig_length is not None and args.max_contig_length is not None
    if use_variable and (args.min_contig_length < 1 or args.max_contig_length < args.min_contig_length):
        parser.error("--min-contig-length and --max-contig-length must be positive and min <= max")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    out_path = args.output_dir / args.output

    effective_length = args.sequence_length
    min_len, max_len = None, None
    if use_variable:
        min_len, max_len = args.min_contig_length, args.max_contig_length
        effective_length = (min_len + max_len) // 2
        print(f"Variable-length contigs: {min_len}–{max_len} bp (uniform)")

    reads_per_organism = args.reads_per_organism
    if args.balanced:
        stats = get_file_stats(
            args.input,
            effective_length,
            min_length=min_len,
            max_length=max_len,
        )
        if not stats:
            parser.error("No FASTA files found under --input")
        for prefix, _fp, total_bases, max_reads in stats:
            print(f"  {prefix}: {total_bases} bp -> ~{max_reads} reads")
        min_reads = min(max_reads for (_, _, _, max_reads) in stats)
        reads_per_organism = min_reads
        print(f"  Balanced: using {reads_per_organism} reads per file.")

    do_train_test_split = getattr(args, "train_test_split", None) is not None
    allow_ambiguous = not getattr(args, "forbid_ambiguous", False)
    result = build_metagenome(
        args.input,
        out_path,
        args.sequence_length,
        reads_per_organism,
        min_length=min_len,
        max_length=max_len,
        seed=args.seed,
        cap_total_reads=args.cap_total_reads,
        eve_intervals=eve_intervals,
        allow_ambiguous=allow_ambiguous,
        filter_similar=getattr(args, "filter_similar", False),
        similarity_threshold=getattr(args, "similarity_threshold", 90.0),
        similarity_min_coverage=getattr(args, "similarity_min_coverage", 0.8),
        oversample_factor=getattr(args, "oversample_factor", 2.0),
        return_records=do_train_test_split,
    )
    if do_train_test_split:
        _count, records = result
        output_stem = Path(args.output).stem
        n_train, n_test = split_train_test_and_write(
            records,
            args.train_test_split,
            args.seed,
            args.output_dir,
            output_stem,
            similarity_threshold=getattr(args, "train_test_similarity_threshold", 90.0),
            similarity_min_coverage=0.8,
        )
        train_path = args.output_dir / f"{output_stem}_train.fasta"
        test_path = args.output_dir / f"{output_stem}_test.fasta"
        print(f"Train-test split ({args.train_test_split}% train): wrote {n_train} to {train_path}, {n_test} to {test_path}")
    else:
        count = result
        if args.cap_total_reads is not None and count == args.cap_total_reads:
            print(f"Wrote {count} sequences (capped) to {out_path}")
        else:
            print(f"Wrote {count} sequences to {out_path}")


def main() -> None:
    _cli()


if __name__ == "__main__":
    main()
