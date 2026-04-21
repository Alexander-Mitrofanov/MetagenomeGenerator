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
import math
import random
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .genome_layout import (
    ARCHAEA_DIR,
    ARCHAEA_PREFIX,
    BACTERIA_DIR,
    BACTERIA_PREFIX,
    PLASMID_DIR,
    PLASMID_PREFIX,
    VIRUS_DIR,
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
    """Return category name for a genome file prefix (bacterial_1 -> bacterial; legacy)."""
    for category, p in _PREFIX_TO_CATEGORY:
        if prefix.startswith(p) or prefix == p.rstrip("_"):
            return category
    return "other"


def _category_from_path(path: Path) -> str:
    """Return category name from parent directory (bacteria/, virus/, archaea/, plasmid/)."""
    parent = path.parent.name.lower()
    if parent == BACTERIA_DIR:
        return "bacteria"
    if parent == VIRUS_DIR:
        return "virus"
    if parent == ARCHAEA_DIR:
        return "archaea"
    if parent == PLASMID_DIR:
        return "plasmid"
    return "other"


def _apply_viral_taxonomy_balance(
    prefixes_and_files: list[tuple[str, Path]],
    read_limits: list[int | None],
    viral_taxonomy: dict[str, str],
    *,
    prefix_to_max_reads: dict[str, int] | None = None,
    input_path: Path | None = None,
    sequence_length: int = 0,
    min_length: int | None = None,
    max_length: int | None = None,
) -> list[int | None]:
    """Overwrite read_limits for viral files so each taxonomy group contributes equally.

    viral_taxonomy: prefix -> group (e.g. viral_1 -> Herpesviridae). Viral files not in
    the dict get limit unchanged. Non-viral files are unchanged.
    When prefix_to_max_reads is provided, stats are not recomputed (avoids double get_file_stats).
    """
    if prefix_to_max_reads is None:
        if input_path is None:
            return read_limits
        stats = get_file_stats(
            input_path,
            sequence_length,
            min_length=min_length,
            max_length=max_length,
        )
        if not stats:
            return read_limits
        prefix_to_max_reads = {p: max_r for p, _fp, _tb, max_r in stats}
    # Group virus prefixes by taxonomy
    groups: dict[str, list[str]] = {}
    viral_prefixes: list[str] = []
    for prefix, fp in prefixes_and_files:
        if _category_from_path(fp) != "virus":
            continue
        viral_prefixes.append(prefix)
        group = viral_taxonomy.get(prefix, "unknown")
        groups.setdefault(group, []).append(prefix)
    if not groups:
        return read_limits
    # Warn when the taxonomy JSON fails to cover the viral prefixes — a common
    # silent failure mode is an accession-format mismatch (e.g. JSON keyed by
    # versioned IDs like "NC_001798.2" while file stems are unversioned, or
    # the reverse). Without this warning, every virus collapses into the
    # "unknown" bucket and balancing becomes a no-op.
    if viral_prefixes:
        unknown_viral = len(groups.get("unknown", []))
        unknown_frac = unknown_viral / len(viral_prefixes)
        if unknown_frac >= 0.5:
            missing_sample = [p for p in viral_prefixes if viral_taxonomy.get(p, "unknown") == "unknown"][:3]
            sample_keys = list(viral_taxonomy.keys())[:3]
            logger.warning(
                "Viral taxonomy balance: %d/%d viral prefixes (%.1f%%) are unmapped and fell into 'unknown'. "
                "Check that the keys in --viral-taxonomy match file stems. "
                "Unmapped example prefixes: %s. Example taxonomy keys: %s.",
                unknown_viral, len(viral_prefixes), 100.0 * unknown_frac,
                missing_sample, sample_keys,
            )
    # Target per group = min over groups of (sum of max_reads in group)
    total_per_group = {}
    for group, prefs in groups.items():
        total_per_group[group] = sum(prefix_to_max_reads.get(p, 0) for p in prefs)
    target_total = min(total_per_group.values()) if total_per_group else 0
    if target_total <= 0:
        return read_limits
    # Per-file limit in each group so group total = target_total
    result = list(read_limits)
    for i, (prefix, fp) in enumerate(prefixes_and_files):
        if _category_from_path(fp) != "virus" or i >= len(result):
            continue
        group = viral_taxonomy.get(prefix, "unknown")
        prefs = groups.get(group, [])
        if not prefs:
            continue
        max_r = prefix_to_max_reads.get(prefix, 0)
        per_file = max(1, target_total // len(prefs))
        result[i] = min(max_r, per_file) if max_r else result[i]
    return result


def _compute_coverage_read_limits(
    prefixes_and_files: list[tuple[str, Path]],
    coverage: float,
    effective_read_length: int,
    *,
    coverage_cv: float = 0.0,
    abundance_profile: dict[str, float] | None = None,
    seed: int | None = None,
) -> list[int]:
    """Return per-file read counts derived from a target coverage.

    For each input file, compute ``ceil(total_bp * per_file_coverage /
    read_length)``. When ``coverage_cv > 0``, per-file coverage is drawn from
    a log-normal with mean ≈ ``coverage`` and coefficient of variation
    ``coverage_cv`` so different genomes get different depths (models real
    metagenomic sequencing, where coverage is uneven across organisms).
    Returns at least 1 read per file (genomes shorter than the read length
    still contribute one attempted read which is skipped downstream if it
    doesn't fit; ensures no genome is silently dropped).
    """
    if effective_read_length <= 0:
        raise ValueError("effective_read_length must be > 0")
    rng = random.Random(seed)
    limits: list[int] = []
    # LogNormal parameters: sigma chosen so the coefficient of variation of
    # exp(N(mu, sigma)) equals ``coverage_cv``; mu set so E[X] == 1 (a pure
    # multiplier applied to ``coverage``). Formulas:
    #   CV = sqrt(exp(sigma^2) - 1)  =>  sigma^2 = log(CV^2 + 1)
    #   E[X] = exp(mu + sigma^2 / 2) == 1  =>  mu = -sigma^2 / 2
    if coverage_cv > 0:
        sigma2 = math.log(coverage_cv * coverage_cv + 1.0)
        sigma = math.sqrt(sigma2)
        mu = -sigma2 / 2.0
    else:
        sigma = 0.0
        mu = 0.0
    for _prefix, fp in prefixes_and_files:
        total_bp = sum(len(rec.seq) for rec in SeqIO.parse(fp, "fasta"))
        cov = coverage
        if sigma > 0:
            cov = coverage * math.exp(rng.gauss(mu, sigma))
        if abundance_profile is not None:
            cat = _category_from_path(fp)
            cov *= abundance_profile.get(cat, 1.0)
        n_reads = max(1, int(math.ceil(total_bp * cov / effective_read_length)))
        limits.append(n_reads)
    return limits


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
        for prefix, fp in prefixes_and_files:
            cat = _category_from_path(fp)
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

# Illumina-like position-dependent error: low at 5', higher at 3' (quality drop)
ILLUMINA_BASE_ERROR = 0.001   # ~0.1% at start of read
ILLUMINA_END_ERROR = 0.012    # ~1.2% toward end (typical for 250 bp)

# 3rd-generation sequencing error profiles (per-base rates). Values are
# representative of contemporary (~2023–2025) platforms and are intentionally
# conservative so synthetic datasets remain usable for classifier training.
# Nanopore R10: substitution-dominated overall error ~4–5%, indels mostly at
# homopolymers; PacBio HiFi (CCS) is very accurate (~0.3% total); PacBio CLR is
# high-error with indel parity.
NANOPORE_SUB_RATE = 0.025
NANOPORE_INS_RATE = 0.015
NANOPORE_DEL_RATE = 0.015
NANOPORE_HOMOPOLYMER_MULT = 2.5   # inflate error rate inside homopolymer runs (>=3 bases)
NANOPORE_HOMOPOLYMER_MIN_RUN = 3

PACBIO_HIFI_SUB_RATE = 0.002
PACBIO_HIFI_INS_RATE = 0.0005
PACBIO_HIFI_DEL_RATE = 0.0005

PACBIO_CLR_SUB_RATE = 0.013
PACBIO_CLR_INS_RATE = 0.055
PACBIO_CLR_DEL_RATE = 0.055

THIRD_GEN_ERROR_MODELS = ("nanopore", "pacbio-hifi", "pacbio-clr")
ALL_ERROR_MODELS = ("illumina",) + THIRD_GEN_ERROR_MODELS

# Contig-quality presets inspired by real metagenome benchmarks:
# low/medium/high bins map to different length ranges and weights.
CONTIG_QUALITY_PRESETS: dict[str, list[tuple[int, int, float]]] = {
    "realistic": [(300, 1200, 0.55), (1200, 5000, 0.35), (5000, 15000, 0.10)],
    "high-quality": [(1000, 4000, 0.20), (4000, 12000, 0.50), (12000, 30000, 0.30)],
    "low-quality": [(300, 900, 0.75), (900, 2500, 0.20), (2500, 6000, 0.05)],
}
CONTIG_QUALITY_PRESET_NAMES = tuple(sorted(CONTIG_QUALITY_PRESETS.keys()))


def normalize_train_split_percent(train_split: float) -> float:
    """Normalize train split to percentage points.

    Accepts either:
    - percentage in [0, 100], e.g. 80
    - fraction in (0, 1], e.g. 0.8 (interpreted as 80%)

    The exact value ``1.0`` is ambiguous (could mean "100%" or "1%"); we follow
    the fraction branch and return ``100.0``, but emit a warning so users who
    meant "1%" realize their mistake. Pass ``1`` (int) or ``1.01``/``0.99`` to
    avoid the warning.
    """
    if train_split == 1.0 and isinstance(train_split, float):
        logger.warning(
            "train_test_split=1.0 is ambiguous; interpreting as 100%% (full train). "
            "If you meant 1%%, pass 0.01. If you meant 100%%, pass 100 to silence this warning."
        )
    if 0.0 < train_split <= 1.0:
        return train_split * 100.0
    return train_split


def _contig_profile_components(name: str) -> list[tuple[int, int, float]]:
    if name not in CONTIG_QUALITY_PRESETS:
        raise ValueError(f"Unknown contig quality profile: {name}. Choose one of {CONTIG_QUALITY_PRESET_NAMES}")
    return CONTIG_QUALITY_PRESETS[name]


def contig_quality_profile_mean_length(name: str) -> int:
    comps = _contig_profile_components(name)
    weighted = sum(((lo + hi) / 2.0) * w for (lo, hi, w) in comps)
    return max(1, int(round(weighted)))


def _weighted_length_from_profile(rng: random.Random, components: list[tuple[int, int, float]]) -> int:
    x = rng.random()
    acc = 0.0
    for lo, hi, w in components:
        acc += w
        if x <= acc:
            return rng.randint(lo, hi)
    lo, hi, _w = components[-1]
    return rng.randint(lo, hi)


def _apply_illumina_like_errors(
    seq_str: str,
    rng: random.Random,
    base_error: float = ILLUMINA_BASE_ERROR,
    end_error: float = ILLUMINA_END_ERROR,
) -> str:
    """Apply position-dependent substitution errors (Illumina-like: higher error toward 3' end).

    No indels (Illumina indel rate is very low). Uses rng for reproducibility.
    """
    if base_error <= 0 and end_error <= 0:
        return seq_str
    seq_str = seq_str.upper()
    L = len(seq_str)
    if L == 0:
        return seq_str
    out: list[str] = []
    for i, c in enumerate(seq_str):
        if c not in BASES:
            out.append(c)
            continue
        # Linear increase from 5' to 3'
        t = i / max(1, L - 1)
        rate = base_error + (end_error - base_error) * t
        if rate > 0 and rng.random() < rate:
            others = [b for b in BASES if b != c]
            out.append(rng.choice(others) if others else c)
        else:
            out.append(c)
    return "".join(out)


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


def _apply_third_gen_errors(
    seq_str: str,
    rng: random.Random,
    *,
    sub_rate: float,
    ins_rate: float,
    del_rate: float,
    homopolymer_mult: float = 1.0,
    homopolymer_min_run: int = 3,
) -> str:
    """Apply long-read (nanopore / PacBio) error patterns to a sequence.

    - Substitutions, insertions, and deletions are applied per-base with
      independent rates.
    - When ``homopolymer_mult > 1`` the per-base rates are multiplied inside
      homopolymer runs of length >= ``homopolymer_min_run`` to emulate the
      well-known nanopore basecaller weakness there.
    - Output length can differ from input (indel-dependent).
    - Uses ``rng`` for reproducibility.
    """
    if sub_rate <= 0 and ins_rate <= 0 and del_rate <= 0:
        return seq_str
    seq_str = seq_str.upper()
    L = len(seq_str)
    if L == 0:
        return seq_str

    # Precompute homopolymer-run lengths so each base knows whether it sits
    # inside a qualifying run. Walking once in O(L) beats recomputing inside
    # the hot loop, which matters for long reads (~10 kb typical nanopore).
    run_lengths: list[int] = [0] * L
    if homopolymer_mult > 1.0 and homopolymer_min_run >= 2:
        i = 0
        while i < L:
            j = i
            while j + 1 < L and seq_str[j + 1] == seq_str[i]:
                j += 1
            run_len = j - i + 1
            for k in range(i, j + 1):
                run_lengths[k] = run_len
            i = j + 1

    out: list[str] = []
    for i, c in enumerate(seq_str):
        if c not in BASES:
            out.append(c)
            continue
        mult = (
            homopolymer_mult
            if (homopolymer_mult > 1.0 and run_lengths[i] >= homopolymer_min_run)
            else 1.0
        )
        eff_ins = min(1.0, ins_rate * mult)
        eff_del = min(1.0, del_rate * mult)
        eff_sub = min(1.0, sub_rate * mult)
        # Insertion (before the current base) — duplicates the current base with
        # probability 0.5 (typical long-read bias) else a random base.
        if eff_ins > 0 and rng.random() < eff_ins:
            out.append(c if rng.random() < 0.5 else rng.choice(BASES))
        # Deletion
        if eff_del > 0 and rng.random() < eff_del:
            continue
        # Substitution
        if eff_sub > 0 and rng.random() < eff_sub:
            others = [b for b in BASES if b != c]
            out.append(rng.choice(others) if others else c)
        else:
            out.append(c)
    return "".join(out)


def _introduce_chimeras(
    records: list[SeqRecord],
    chimera_rate: float,
    rng: random.Random,
) -> tuple[list[SeqRecord], int]:
    """Replace a fraction of records with two-parent chimeras.

    Returns (new_records, n_chimeras). Each chimera glues the 5' half of one
    parent to the 3' half of a different parent (random pick), tagged with
    ``chimera=parentA|parentB`` in the description for traceability. The
    overall list length is preserved so downstream read counts match user
    expectation (``chimera_rate`` = fraction of emitted reads that are
    chimeric). Needs at least 2 records and ``chimera_rate > 0`` to do
    anything.
    """
    if chimera_rate <= 0 or len(records) < 2:
        return records, 0
    chimera_rate = min(1.0, chimera_rate)
    n = len(records)
    n_chimeras = int(round(n * chimera_rate))
    if n_chimeras == 0:
        return records, 0
    indices = list(range(n))
    rng.shuffle(indices)
    to_replace = set(indices[:n_chimeras])
    out: list[SeqRecord] = list(records)
    for idx in to_replace:
        parent_a = records[idx]
        # Pick a different parent; fall back to any parent if degenerate
        candidates = [i for i in range(n) if i != idx]
        parent_b_idx = rng.choice(candidates)
        parent_b = records[parent_b_idx]
        seq_a = str(parent_a.seq)
        seq_b = str(parent_b.seq)
        split_a = len(seq_a) // 2
        split_b = len(seq_b) // 2
        new_seq = seq_a[:split_a] + seq_b[split_b:]
        chimera = SeqRecord(
            Seq(new_seq),
            id=f"chimera_{idx}",
            description=f"chimera parents={parent_a.id}|{parent_b.id}",
        )
        # Preserve Phred quality when upstream records have it, so FASTQ still
        # works. We concatenate the matching halves of parents' qualities.
        qa = parent_a.letter_annotations.get("phred_quality")
        qb = parent_b.letter_annotations.get("phred_quality")
        if qa is not None and qb is not None:
            chimera.letter_annotations["phred_quality"] = list(qa[:split_a]) + list(qb[split_b:])
        out[idx] = chimera
    return out, n_chimeras


def _introduce_pcr_duplicates(
    records: list[SeqRecord],
    duplicate_rate: float,
    rng: random.Random,
) -> tuple[list[SeqRecord], int]:
    """Append PCR-like duplicates to the record list.

    Each record is independently duplicated with probability
    ``duplicate_rate``; duplicates get ``_dup`` appended to their id and
    ``pcr_duplicate=true`` in the description. Returns (new_records,
    n_duplicates_added). The duplicate is bit-identical to the source; callers
    that want slight mutations on duplicates should apply an error model
    first (before this call).
    """
    if duplicate_rate <= 0 or not records:
        return records, 0
    duplicate_rate = min(1.0, duplicate_rate)
    out = list(records)
    dup_counter = 0
    for rec in records:
        if rng.random() < duplicate_rate:
            dup_counter += 1
            dup = SeqRecord(
                Seq(str(rec.seq)),
                id=f"{rec.id}_dup",
                description=(f"{rec.description} pcr_duplicate=true" if rec.description else "pcr_duplicate=true"),
            )
            q = rec.letter_annotations.get("phred_quality")
            if q is not None:
                dup.letter_annotations["phred_quality"] = list(q)
            out.append(dup)
    return out, dup_counter


def _apply_error_model_to_record(
    rec: SeqRecord,
    error_model: str,
    rng: random.Random,
) -> SeqRecord:
    """Apply platform-specific error model to a record. Preserves id and description.

    Supported models:

    - ``illumina`` — position-dependent substitution only (short-read).
    - ``nanopore`` — substitutions + indels with homopolymer inflation.
    - ``pacbio-hifi`` — low-error substitution-biased (CCS-like).
    - ``pacbio-clr`` — high-error, indel-heavy continuous long reads.
    """
    if error_model == "illumina":
        new_seq = _apply_illumina_like_errors(str(rec.seq), rng)
    elif error_model == "nanopore":
        new_seq = _apply_third_gen_errors(
            str(rec.seq),
            rng,
            sub_rate=NANOPORE_SUB_RATE,
            ins_rate=NANOPORE_INS_RATE,
            del_rate=NANOPORE_DEL_RATE,
            homopolymer_mult=NANOPORE_HOMOPOLYMER_MULT,
            homopolymer_min_run=NANOPORE_HOMOPOLYMER_MIN_RUN,
        )
    elif error_model == "pacbio-hifi":
        new_seq = _apply_third_gen_errors(
            str(rec.seq),
            rng,
            sub_rate=PACBIO_HIFI_SUB_RATE,
            ins_rate=PACBIO_HIFI_INS_RATE,
            del_rate=PACBIO_HIFI_DEL_RATE,
            homopolymer_mult=1.0,
        )
    elif error_model == "pacbio-clr":
        new_seq = _apply_third_gen_errors(
            str(rec.seq),
            rng,
            sub_rate=PACBIO_CLR_SUB_RATE,
            ins_rate=PACBIO_CLR_INS_RATE,
            del_rate=PACBIO_CLR_DEL_RATE,
            homopolymer_mult=1.0,
        )
    else:
        return rec
    return SeqRecord(Seq(new_seq), id=rec.id, description=rec.description)


def _illumina_phred_at_position(i: int, length: int) -> int:
    """Return Phred quality (0–41) for position i (0-based) in a read of length, matching our Illumina-like error profile."""
    if length <= 0:
        return 30
    t = i / max(1, length - 1)
    rate = ILLUMINA_BASE_ERROR + (ILLUMINA_END_ERROR - ILLUMINA_BASE_ERROR) * t
    rate = max(1e-10, min(1.0, rate))
    q = -10.0 * math.log10(rate)
    return int(min(41, max(0, round(q))))


def _phred_from_rate(rate: float) -> int:
    """Convert a per-base error rate to a flat Phred quality score (clamped 0..60)."""
    rate = max(1e-10, min(1.0, rate))
    q = -10.0 * math.log10(rate)
    return int(min(60, max(0, round(q))))


# Approximate total per-base error rate used to back out a flat Phred score for
# 3rd-generation platforms. Keep in sync with the ``*_SUB_RATE/INS_RATE/DEL_RATE``
# constants above; we sum rather than compound because real basecallers report
# error as an additive "total". Illumina keeps the position-dependent curve.
_ERROR_MODEL_FLAT_RATE = {
    "nanopore": NANOPORE_SUB_RATE + NANOPORE_INS_RATE + NANOPORE_DEL_RATE,
    "pacbio-hifi": PACBIO_HIFI_SUB_RATE + PACBIO_HIFI_INS_RATE + PACBIO_HIFI_DEL_RATE,
    "pacbio-clr": PACBIO_CLR_SUB_RATE + PACBIO_CLR_INS_RATE + PACBIO_CLR_DEL_RATE,
}


def add_illumina_qualities_to_record(rec: SeqRecord, error_model: str | None = None) -> SeqRecord:
    """Add ``letter_annotations['phred_quality']`` to a SeqRecord.

    For backward compatibility the default (``error_model=None`` or
    ``"illumina"``) keeps the position-dependent Illumina-like profile. For
    3rd-generation models (``nanopore``, ``pacbio-hifi``, ``pacbio-clr``) a
    flat Phred score derived from the model's total per-base error rate is
    written — this is a simplification, but matches the long-read convention
    of reporting mean read quality rather than position-dependent qualities.
    Returns the same record, modified in place.
    """
    L = len(rec.seq)
    if error_model in _ERROR_MODEL_FLAT_RATE:
        q = _phred_from_rate(_ERROR_MODEL_FLAT_RATE[error_model])
        rec.letter_annotations["phred_quality"] = [q] * L
    else:
        rec.letter_annotations["phred_quality"] = [_illumina_phred_at_position(i, L) for i in range(L)]
    return rec


def chunk_sequence(
    record,
    prefix: str,
    chunk_size: int,
    yield_coords: bool = False,
    *,
    random_chunk_start: bool = False,
    rng: random.Random | None = None,
    reads_per_organism: int | None = None,
):
    """Split a sequence into non-overlapping fixed-size simulated reads.

    Only full-length reads are emitted (trailing remainder is dropped).
    Record id: {prefix}_read_{idx}; description includes start/end (0-based) for traceability.
    If yield_coords, yields (rec, start, end) for EVE overlap checks; else yields rec.

    When random_chunk_start is True, chooses a random offset in [0, chunk_size) so
    reads start at e.g. 17..1017 instead of always 0..1000, while still keeping
    reads non-overlapping (step = chunk_size).

    If ``reads_per_organism`` is set, emission stops after that many reads so the
    cap matches :func:`chunk_sequence_variable` and
    :func:`chunk_sequence_quality_profile`. Callers may still apply their own
    cap on top (e.g. a per-file budget across multi-record FASTAs); this
    internal cap just prevents direct callers from accidentally consuming the
    whole genome when they only wanted a few reads.
    """
    seq = record.seq
    start_offset = 0
    if random_chunk_start:
        if rng is None:
            rng = random.Random()
        start_offset = rng.randint(0, max(0, chunk_size - 1)) if chunk_size > 0 else 0

    idx = 0
    start = start_offset
    while start + chunk_size <= len(seq):
        if reads_per_organism is not None and idx >= reads_per_organism:
            break
        start, end = start, start + chunk_size
        sub = seq[start:end]
        rec = SeqRecord(
            sub,
            id=f"{prefix}_read_{idx}",
            description=f"start={start} end={end}",
        )
        if yield_coords:
            yield rec, start, end
        else:
            yield rec
        idx += 1
        start += chunk_size


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
    Record id: {prefix}_contig_{idx}; description includes start/end (0-based) for traceability.
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
        start = pos
        sub = seq[start:end]
        rec = SeqRecord(
            sub,
            id=f"{prefix}_contig_{idx}",
            description=f"start={start} end={end}",
        )
        if yield_coords:
            yield rec, start, end
        else:
            yield rec
        idx += 1
        pos = end
        if reads_per_organism is not None and idx >= reads_per_organism:
            break


def chunk_sequence_paired(
    record,
    prefix: str,
    read_length: int,
    insert_size: int,
    insert_size_sd: float,
    reads_per_organism: int | None,
    rng: random.Random | None = None,
):
    """Yield (R1, R2) SeqRecord pairs simulating paired-end sequencing.

    The genome is walked in non-overlapping steps of ``insert_size`` (same
    non-overlapping contract as :func:`chunk_sequence`). For each step, a
    fragment length ``L`` is drawn from ``Normal(insert_size, insert_size_sd)``
    and clamped to ``[read_length, insert_size + 3*insert_size_sd]`` (or
    simply ``read_length`` if ``insert_size_sd <= 0``). R1 is the 5' end of
    the fragment, R2 is the reverse-complement of the 3' end.

    Record ids use the ``/1`` and ``/2`` Illumina convention. Descriptions
    carry ``start``/``end``/``insert`` for traceability.

    ``reads_per_organism`` caps the number of emitted **pairs**, matching the
    semantics of the other chunker helpers (1 "read" in the cap corresponds
    to 1 fragment / 1 pair in paired mode).
    """
    rng = rng or random.Random()
    seq = record.seq
    L = len(seq)
    if read_length <= 0 or insert_size <= 0:
        return
    if insert_size < read_length:
        # Each fragment must fit two non-overlapping reads; otherwise skip.
        return
    max_frag = insert_size + int(round(3 * insert_size_sd)) if insert_size_sd > 0 else insert_size
    step = insert_size
    idx = 0
    pos = 0
    while pos + max_frag <= L:
        if reads_per_organism is not None and idx >= reads_per_organism:
            break
        if insert_size_sd > 0:
            frag_len = int(round(rng.gauss(insert_size, insert_size_sd)))
            frag_len = max(read_length, min(frag_len, max_frag))
        else:
            frag_len = insert_size
        if pos + frag_len > L:
            pos += step
            continue
        r1_start = pos
        r1_end = pos + read_length
        r2_end = pos + frag_len
        r2_start = r2_end - read_length
        r1_seq = seq[r1_start:r1_end]
        # BioPython's reverse_complement handles any ambiguous bases safely.
        r2_seq = seq[r2_start:r2_end].reverse_complement()
        r1 = SeqRecord(
            r1_seq,
            id=f"{prefix}_pair_{idx}/1",
            description=f"start={r1_start} end={r1_end} insert={frag_len} mate=R1",
        )
        r2 = SeqRecord(
            r2_seq,
            id=f"{prefix}_pair_{idx}/2",
            description=f"start={r2_start} end={r2_end} insert={frag_len} mate=R2",
        )
        yield r1, r2
        idx += 1
        pos += step


def chunk_sequence_quality_profile(
    record,
    prefix: str,
    profile_name: str,
    reads_per_organism: int | None,
    rng: random.Random | None = None,
    yield_coords: bool = False,
):
    """Split sequence into non-overlapping contigs using a quality-profile length mixture."""
    rng = rng or random.Random()
    seq = record.seq
    idx = 0
    pos = 0
    components = _contig_profile_components(profile_name)
    min_len = min(lo for (lo, _hi, _w) in components)
    while pos < len(seq):
        length = _weighted_length_from_profile(rng, components)
        end = min(pos + length, len(seq))
        if end - pos < min_len:
            break
        start = pos
        sub = seq[start:end]
        rec = SeqRecord(
            sub,
            id=f"{prefix}_contig_{idx}",
            description=f"start={start} end={end}",
        )
        if yield_coords:
            yield rec, start, end
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
    contig_quality_profile: str | None = None,
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
    if contig_quality_profile is not None:
        effective_len = contig_quality_profile_mean_length(contig_quality_profile)
    elif min_length is not None and max_length is not None:
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


def _embed_tax_in_description(
    rec: SeqRecord,
    cat: str,
    prefix: str,
    viral_taxonomy: dict[str, str] | None,
) -> SeqRecord:
    """Append ``tax=<group>`` to the record description for gold-standard labels.

    Non-viral categories get ``tax=<category>``. Viral records use
    ``viral_taxonomy.get(prefix, "unknown")`` so the same JSON produced by the
    ``viral-taxonomy`` subcommand feeds both balancing and header embedding.

    The original description is preserved (``start=… end=… tax=…``) so
    downstream parsers that already look for ``start=`` / ``end=`` keep
    working.
    """
    if cat == "virus":
        group = (viral_taxonomy or {}).get(prefix, "unknown")
    else:
        group = cat
    existing = rec.description or ""
    # Avoid duplicating the id that BioPython prepends when it writes the file;
    # we only append once per record.
    rec.description = f"{existing} tax={group}".strip()
    return rec


def _collect_chunks_for_file(
    fp: Path,
    sequence_length: int,
    reads_per_organism: int | None,
    prefix: str,
    *,
    min_length: int | None = None,
    max_length: int | None = None,
    contig_quality_profile: str | None = None,
    seed: int | None = None,
    eve_intervals: dict[tuple[str, str], list[tuple[int, int]]] | None = None,
    allow_ambiguous: bool = True,
    substitution_rate: float = 0.0,
    indel_rate: float = 0.0,
    error_model: str | None = None,
    mutation_rng: random.Random | None = None,
    random_chunk_start: bool = False,
    embed_taxonomy: bool = False,
    viral_taxonomy: dict[str, str] | None = None,
) -> list[SeqRecord]:
    """Collect chunks for a single FASTA file (fixed- or variable-length).

    ``reads_per_organism`` is a **per-file** cap shared across every record in a
    multi-record FASTA (e.g. a multi-segment virus stored as one file). All
    three chunker helpers now honour the cap consistently:

    - ``chunk_sequence`` (fixed-length) reads the outer budget via
      ``reads_per_organism`` passed through on each record invocation; this
      file-level loop also shortcuts as soon as the budget is exhausted, so no
      record past the limit is opened.
    - ``chunk_sequence_variable`` and ``chunk_sequence_quality_profile``
      enforce the same cap internally.

    Chunk IDs are disambiguated across records using a ``_seg{N}`` infix for
    the 2nd record onward (single-record files keep the plain
    ``{prefix}_read_{idx}`` naming), so no two chunks from the same file share
    a header.

    Skips chunks overlapping EVE if ``eve_intervals`` is given. If
    ``error_model`` is set (e.g. ``'illumina'``) platform-specific errors are
    applied; else if ``substitution_rate`` or ``indel_rate`` > 0, uniform
    mutations are applied.
    """
    from .blastn_filter import chunk_overlaps_eve

    chunks: list[SeqRecord] = []
    cat = _category_from_path(fp)
    rng = random.Random(seed) if seed is not None else None
    use_eve = eve_intervals is not None
    apply_platform = error_model and mutation_rng is not None
    apply_mutations = not apply_platform and (substitution_rate > 0 or indel_rate > 0) and mutation_rng is not None

    # reads_per_organism is a per-file cap (see CLI help "per input file"). When a
    # FASTA contains multiple records (multi-segment viruses, multi-contig drafts,
    # user-supplied multi-record FASTAs), iterate every record but share the budget
    # across them instead of processing only the first. Record index is only
    # appended to the prefix from the second record onward so single-record
    # NCBI downloads keep their historical chunk IDs (``{accession}_read_{idx}``).
    for rec_idx, record in enumerate(SeqIO.parse(fp, "fasta")):
        if reads_per_organism is not None and len(chunks) >= reads_per_organism:
            break
        remaining = None if reads_per_organism is None else reads_per_organism - len(chunks)
        record_prefix = prefix if rec_idx == 0 else f"{prefix}_seg{rec_idx}"
        key = (prefix, record.id)
        intervals = eve_intervals.get(key, []) if eve_intervals else []

        if contig_quality_profile is not None:
            for item in chunk_sequence_quality_profile(
                record, record_prefix, contig_quality_profile, remaining, rng=rng, yield_coords=use_eve
            ):
                if use_eve:
                    rec, start, end = item
                    if chunk_overlaps_eve(start, end, intervals):
                        continue
                else:
                    rec = item
                # Ensure FASTA/FASTQ headers start with the category class label.
                rec.id = f"{cat}_{rec.id}"
                if not _is_allowed_sequence(rec, allow_ambiguous):
                    continue
                if apply_platform:
                    rec = _apply_error_model_to_record(rec, error_model, mutation_rng)
                elif apply_mutations:
                    rec = _apply_mutations_to_record(rec, substitution_rate, indel_rate, mutation_rng)
                if embed_taxonomy:
                    rec = _embed_tax_in_description(rec, cat, prefix, viral_taxonomy)
                chunks.append(rec)
                if reads_per_organism is not None and len(chunks) >= reads_per_organism:
                    break
        elif min_length is not None and max_length is not None:
            for item in chunk_sequence_variable(
                record, record_prefix, min_length, max_length, remaining, rng=rng, yield_coords=use_eve
            ):
                if use_eve:
                    rec, start, end = item
                    if chunk_overlaps_eve(start, end, intervals):
                        continue
                else:
                    rec = item
                # Ensure FASTA/FASTQ headers start with the category class label.
                rec.id = f"{cat}_{rec.id}"
                if not _is_allowed_sequence(rec, allow_ambiguous):
                    continue
                if apply_platform:
                    rec = _apply_error_model_to_record(rec, error_model, mutation_rng)
                elif apply_mutations:
                    rec = _apply_mutations_to_record(rec, substitution_rate, indel_rate, mutation_rng)
                if embed_taxonomy:
                    rec = _embed_tax_in_description(rec, cat, prefix, viral_taxonomy)
                chunks.append(rec)
                if reads_per_organism is not None and len(chunks) >= reads_per_organism:
                    break
        else:
            for i, item in enumerate(
                chunk_sequence(
                    record,
                    record_prefix,
                    sequence_length,
                    yield_coords=use_eve,
                    random_chunk_start=random_chunk_start,
                    rng=rng,
                )
            ):
                if remaining is not None and i >= remaining:
                    break
                if use_eve:
                    rec, start, end = item
                    if chunk_overlaps_eve(start, end, intervals):
                        continue
                else:
                    rec = item
                # Ensure FASTA/FASTQ headers start with the category class label.
                rec.id = f"{cat}_{rec.id}"
                if not _is_allowed_sequence(rec, allow_ambiguous):
                    continue
                if apply_platform:
                    rec = _apply_error_model_to_record(rec, error_model, mutation_rng)
                elif apply_mutations:
                    rec = _apply_mutations_to_record(rec, substitution_rate, indel_rate, mutation_rng)
                if embed_taxonomy:
                    rec = _embed_tax_in_description(rec, cat, prefix, viral_taxonomy)
                chunks.append(rec)
                if reads_per_organism is not None and len(chunks) >= reads_per_organism:
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
    contig_quality_profile: str | None = None,
    seed: int | None = None,
    allow_ambiguous: bool = True,
    substitution_rate: float = 0.0,
    indel_rate: float = 0.0,
    error_model: str | None = None,
    category_label: str = "virus",
    random_chunk_start: bool = False,
    embed_taxonomy: bool = False,
    viral_taxonomy: dict[str, str] | None = None,
) -> list[SeqRecord]:
    """Collect chunks from a multi-record FASTA (e.g. metavirome contigs). No EVE filtering.

    Each record is chunked with prefix prefix_base_0, prefix_base_1, ... Same length/balance
    and mutation/error-model logic as main input; use seed for reproducibility.
    """
    chunks: list[SeqRecord] = []
    use_platform = bool(error_model)
    use_mutations = not use_platform and (substitution_rate > 0 or indel_rate > 0)
    for rec_idx, record in enumerate(SeqIO.parse(fp, "fasta")):
        prefix = f"{prefix_base}_{rec_idx}"
        file_seed = (seed + 20000 + rec_idx) if seed is not None else None
        rng = random.Random(file_seed) if file_seed is not None else None
        file_mutation_rng = random.Random(seed + 30000 + rec_idx) if (use_platform or use_mutations) else None
        apply_platform_here = use_platform and file_mutation_rng is not None
        apply_mutations = use_mutations and file_mutation_rng is not None
        if contig_quality_profile is not None:
            for item in chunk_sequence_quality_profile(
                record, prefix, contig_quality_profile, reads_per_organism, rng=rng, yield_coords=False
            ):
                rec = item
                rec.id = f"{category_label}_{rec.id}"
                if not _is_allowed_sequence(rec, allow_ambiguous):
                    continue
                if apply_platform_here:
                    rec = _apply_error_model_to_record(rec, error_model, file_mutation_rng)
                elif apply_mutations:
                    rec = _apply_mutations_to_record(rec, substitution_rate, indel_rate, file_mutation_rng)
                if embed_taxonomy:
                    rec = _embed_tax_in_description(rec, category_label, prefix, viral_taxonomy)
                chunks.append(rec)
        elif min_length is not None and max_length is not None:
            for item in chunk_sequence_variable(
                record, prefix, min_length, max_length, reads_per_organism, rng=rng, yield_coords=False
            ):
                rec = item
                rec.id = f"{category_label}_{rec.id}"
                if not _is_allowed_sequence(rec, allow_ambiguous):
                    continue
                if apply_platform_here:
                    rec = _apply_error_model_to_record(rec, error_model, file_mutation_rng)
                elif apply_mutations:
                    rec = _apply_mutations_to_record(rec, substitution_rate, indel_rate, file_mutation_rng)
                if embed_taxonomy:
                    rec = _embed_tax_in_description(rec, category_label, prefix, viral_taxonomy)
                chunks.append(rec)
        else:
            for i, rec in enumerate(
                chunk_sequence(
                    record,
                    prefix,
                    sequence_length,
                    yield_coords=False,
                    random_chunk_start=random_chunk_start,
                    rng=rng,
                )
            ):
                if reads_per_organism is not None and i >= reads_per_organism:
                    break
                if not _is_allowed_sequence(rec, allow_ambiguous):
                    continue
                rec.id = f"{category_label}_{rec.id}"
                if apply_platform_here:
                    rec = _apply_error_model_to_record(rec, error_model, file_mutation_rng)
                elif apply_mutations:
                    rec = _apply_mutations_to_record(rec, substitution_rate, indel_rate, file_mutation_rng)
                if embed_taxonomy:
                    rec = _embed_tax_in_description(rec, category_label, prefix, viral_taxonomy)
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
    contig_quality_profile: str | None = None,
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
    viral_taxonomy_json: Path | None = None,
    balance_viral_by_taxonomy: bool = False,
    error_model: str | None = None,
    output_fastq: bool = False,
    write_abundance: bool = False,
    random_chunk_start: bool = False,
    embed_taxonomy: bool = False,
    chimera_rate: float = 0.0,
    pcr_duplicate_rate: float = 0.0,
    coverage: float | None = None,
    coverage_cv: float = 0.0,
) -> int | tuple[int, list[SeqRecord]]:
    """Build a metagenome FASTA from input_path. Fixed-length or variable-length contigs.

    If cap_total_reads is set, randomly downsample to that many reads (for balancing negative to positive).
    If eve_intervals is set, reads/contigs overlapping those intervals (EVE regions) are excluded.
    If error_model is set (e.g. 'illumina'), applies platform-specific position-dependent errors; else if substitution_rate or indel_rate > 0, applies uniform mutations (seed used for reproducibility; default 42).
    If output_fastq is True, writes FASTQ instead of FASTA: applies Illumina-like errors and adds position-dependent Phred quality scores per base (output path suffix becomes .fastq).
    If contig_quality_profile is set, contig lengths are sampled from preset low/medium/high quality strata.
    If random_chunk_start is True, fixed-length chunks (reads) are generated starting at a random offset (still non-overlapping) rather than always starting at position 0.
    If filter_similar is True: generate more reads than needed (oversample), filter out sequences that are
    >= similarity_threshold (default 90%%) similar to any already-kept sequence, then refill until target or max rounds.
    If return_records is True, do not write to out_path and return (count, list[SeqRecord]) for train-test split.
    If extra_viral_fasta is set, reads from that FASTA (multi-record, e.g. metavirome contigs) are merged in
    with prefix extra_viral_0, extra_viral_1, ... using the same length and mutation settings.
    If abundance_profile is set (e.g. {"bacteria": 0.5, "virus": 2.0}), scale reads per file by category.
    If abundance_distribution == "exponential", assign per-genome weights from Exp(1) (use seed for reproducibility).
    If viral_taxonomy_json and balance_viral_by_taxonomy, viral read limits are set so each taxonomy group (e.g. family) contributes equally.
    If write_abundance is True, write a tab-separated file next to the output (stem_abundance.txt) with columns: genome_id (prefix), read_count, proportion (ground-truth composition for benchmarking).
    """
    from .similarity_filter import filter_by_similarity, filter_candidates_against_kept

    if contig_quality_profile is not None and (min_length is not None or max_length is not None):
        raise ValueError("Use either contig_quality_profile OR min_length/max_length, not both.")
    if output_fastq and not error_model:
        error_model = "illumina"
    use_mutations = substitution_rate > 0 or indel_rate > 0 or bool(error_model)
    if extra_viral_fasta is not None and not extra_viral_fasta.is_file():
        raise FileNotFoundError(f"extra_viral_fasta not found or not a file: {extra_viral_fasta}")
    if (use_mutations or output_fastq or random_chunk_start) and seed is None:
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

    if coverage is not None and coverage > 0:
        # Coverage-driven mode: derive per-file read counts from genome length
        # and coverage. This overrides ``reads_per_organism`` (which is still
        # honoured as an upper cap per file after the coverage calculation so
        # users can safely combine ``--coverage`` with ``--reads-per-organism``
        # as a safety valve). When ``filter_similar`` is on we multiply by the
        # oversample factor here so the downstream filter has enough
        # candidates to choose from.
        eff_read_len = sequence_length
        if contig_quality_profile is not None:
            eff_read_len = contig_quality_profile_mean_length(contig_quality_profile)
        elif min_length is not None and max_length is not None:
            eff_read_len = (min_length + max_length) // 2
        eff_coverage = coverage * oversample_factor if filter_similar else coverage
        read_limits = _compute_coverage_read_limits(
            prefixes_and_files,
            eff_coverage,
            max(1, eff_read_len),
            coverage_cv=coverage_cv,
            abundance_profile=abundance_profile,
            seed=seed,
        )
        if reads_per_organism is not None:
            read_limits = [min(rl, reads_per_organism) for rl in read_limits]
    elif reads_per_organism_gen is None:
        read_limits: list[int | None] = [None] * num_files
    else:
        read_limits = _compute_read_limits(
            prefixes_and_files,
            reads_per_organism_gen,
            abundance_profile=abundance_profile,
            abundance_distribution=abundance_distribution,
            seed=seed,
        )

    # Load viral taxonomy once when needed for either balancing or
    # ``--embed-taxonomy`` header tagging, so we don't parse the JSON twice.
    viral_tax: dict[str, str] | None = None
    if viral_taxonomy_json and viral_taxonomy_json.is_file() and (balance_viral_by_taxonomy or embed_taxonomy):
        from .viral_taxonomy import load_viral_taxonomy
        viral_tax = load_viral_taxonomy(viral_taxonomy_json) or None

    if viral_tax and balance_viral_by_taxonomy:
        # Reuse single get_file_stats result for balance (avoids re-scanning all FASTAs)
        stats = get_file_stats(
            input_path,
            sequence_length,
            min_length=min_length,
            max_length=max_length,
            contig_quality_profile=contig_quality_profile,
        )
        prefix_to_max_reads = {p: max_r for p, _fp, _tb, max_r in stats} if stats else {}
        read_limits = _apply_viral_taxonomy_balance(
            prefixes_and_files,
            read_limits,
            viral_tax,
            prefix_to_max_reads=prefix_to_max_reads,
        )
        logger.info("Viral taxonomy balance: applied for %d viral prefixes", len(viral_tax))

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
            contig_quality_profile=contig_quality_profile,
            seed=file_seed,
            eve_intervals=eve_intervals,
            allow_ambiguous=allow_ambiguous,
            substitution_rate=substitution_rate,
            indel_rate=indel_rate,
            error_model=error_model,
            mutation_rng=file_mutation_rng,
            random_chunk_start=random_chunk_start,
            embed_taxonomy=embed_taxonomy,
            viral_taxonomy=viral_tax,
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
            contig_quality_profile=contig_quality_profile,
            seed=extra_seed,
            allow_ambiguous=allow_ambiguous,
            substitution_rate=substitution_rate,
            indel_rate=indel_rate,
            error_model=error_model,
            category_label="virus",
            random_chunk_start=random_chunk_start,
            embed_taxonomy=embed_taxonomy,
            viral_taxonomy=viral_tax,
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
                    contig_quality_profile=contig_quality_profile,
                    seed=file_seed,
                    eve_intervals=eve_intervals,
                    allow_ambiguous=allow_ambiguous,
                    substitution_rate=substitution_rate,
                    indel_rate=indel_rate,
                    error_model=error_model,
                    mutation_rng=file_mutation_rng,
                    random_chunk_start=random_chunk_start,
                    embed_taxonomy=embed_taxonomy,
                    viral_taxonomy=viral_tax,
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

    # Library-prep artefacts (chimeras + PCR duplicates). Applied last so they
    # compose with every other upstream feature (coverage, error models,
    # similarity filter, train/test). We derive a dedicated sub-seed so the
    # injection is reproducible without perturbing the chunking RNG stream
    # consumed above. The chimera pass runs first because PCR duplication then
    # replicates the already-chimeric reads, matching real wet-lab artefact
    # ordering (fragmentation → PCR amplification).
    if chimera_rate > 0 or pcr_duplicate_rate > 0:
        artefact_seed = (seed + 40000) if seed is not None else None
        artefact_rng = random.Random(artefact_seed) if artefact_seed is not None else random.Random()
        if chimera_rate > 0:
            all_chunks, n_chim = _introduce_chimeras(all_chunks, chimera_rate, artefact_rng)
            logger.info("Chimera injection: replaced %d / %d reads (rate=%.3f)", n_chim, len(all_chunks), chimera_rate)
        if pcr_duplicate_rate > 0:
            before = len(all_chunks)
            all_chunks, n_dup = _introduce_pcr_duplicates(all_chunks, pcr_duplicate_rate, artefact_rng)
            logger.info(
                "PCR duplicates: added %d duplicate reads (rate=%.3f, base=%d -> %d)",
                n_dup, pcr_duplicate_rate, before, len(all_chunks),
            )

    def _write_abundance_file(records: list[SeqRecord], path: Path) -> None:
        """Write genome_id, read_count, proportion to a tab-separated file (ground-truth for benchmarking)."""
        from collections import Counter
        counts: Counter[str] = Counter()
        def _strip_category_label(s: str) -> str:
            # We prefix output FASTA/FASTQ ids with the category label (e.g. "bacteria_"),
            # but the abundance ground-truth should be keyed by the original genome prefix.
            for cat in ("bacteria", "virus", "archaea", "plasmid"):
                marker = f"{cat}_"
                if s.startswith(marker):
                    return s[len(marker):]
            return s

        for rec in records:
            if "_read_" in rec.id:
                prefix = _strip_category_label(rec.id.rsplit("_read_", 1)[0])
            elif "_contig_" in rec.id:
                prefix = _strip_category_label(rec.id.rsplit("_contig_", 1)[0])
            else:
                prefix = _strip_category_label(rec.id.split()[0])
            counts[prefix] += 1
        total = sum(counts.values())
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w") as f:
            f.write("genome_id\tread_count\tproportion\n")
            for prefix in sorted(counts):
                n = counts[prefix]
                f.write(f"{prefix}\t{n}\t{n / total}\n")

    if return_records:
        if output_fastq:
            for rec in all_chunks:
                add_illumina_qualities_to_record(rec, error_model=error_model)
        return len(all_chunks), all_chunks
    if output_fastq:
        for rec in all_chunks:
            add_illumina_qualities_to_record(rec, error_model=error_model)
        write_path = out_path if out_path.suffix.lower() == ".fastq" else out_path.with_suffix(".fastq")
        count = SeqIO.write(all_chunks, write_path, "fastq")
        if write_abundance:
            _write_abundance_file(all_chunks, write_path.with_stem(write_path.stem + "_abundance").with_suffix(".txt"))
        return int(count)
    count = SeqIO.write(all_chunks, out_path, "fasta")
    if write_abundance:
        _write_abundance_file(all_chunks, out_path.with_stem(out_path.stem + "_abundance").with_suffix(".txt"))
    return int(count)


def parse_multi_length_spec(spec: str) -> list[int]:
    """Parse a comma-separated multi-length spec like ``"300,500,1000,3000"``.

    Returns a list of unique positive integers preserving input order. Raises
    ``ValueError`` on empty input, non-integer tokens, or non-positive values.
    """
    if not spec or not spec.strip():
        raise ValueError("multi-length spec must be a non-empty comma-separated list, e.g. '300,500,1000,3000'")
    tokens = [t.strip() for t in spec.split(",") if t.strip()]
    if not tokens:
        raise ValueError("multi-length spec contains no lengths after splitting on commas")
    lengths: list[int] = []
    seen: set[int] = set()
    for tok in tokens:
        try:
            value = int(tok)
        except ValueError as exc:
            raise ValueError(f"multi-length token {tok!r} is not an integer") from exc
        if value <= 0:
            raise ValueError(f"multi-length value {value} must be > 0")
        if value in seen:
            continue
        seen.add(value)
        lengths.append(value)
    return lengths


def build_metagenome_multi_length(
    input_path: Path,
    out_path: Path,
    lengths: list[int],
    reads_per_organism: int | None,
    **kwargs,
) -> list[tuple[int, int, Path]]:
    """Write one metagenome output per sequence length.

    Calls :func:`build_metagenome` once per entry in ``lengths``, writing to
    ``{stem}_L{N}{suffix}`` (e.g. ``metagenome.fasta`` with lengths
    ``[300, 500]`` produces ``metagenome_L300.fasta`` and
    ``metagenome_L500.fasta``). ``kwargs`` are forwarded to each call.
    Returns a list of ``(length, written_count, written_path)`` tuples
    (path uses ``.fastq`` when ``output_fastq=True``).

    Intentionally does *not* support ``return_records=True`` or
    variable-length modes (``min_length``/``max_length``/
    ``contig_quality_profile``) — multi-length is defined as fixed-length
    reads of each requested size, and each length would otherwise produce
    its own record list that isn't meaningful to merge.
    """
    if not lengths:
        raise ValueError("build_metagenome_multi_length requires at least one length")
    if kwargs.get("return_records"):
        raise ValueError("build_metagenome_multi_length does not support return_records=True")
    for incompat in ("min_length", "max_length", "contig_quality_profile"):
        if kwargs.get(incompat) is not None:
            raise ValueError(
                f"build_metagenome_multi_length: {incompat!r} is incompatible with multi-length (fixed-length only)"
            )
    results: list[tuple[int, int, Path]] = []
    output_fastq = bool(kwargs.get("output_fastq"))
    stem = out_path.stem
    suffix = ".fastq" if output_fastq else (out_path.suffix if out_path.suffix else ".fasta")
    for length in lengths:
        per_length_name = f"{stem}_L{length}{suffix}"
        per_length_path = out_path.with_name(per_length_name)
        count = build_metagenome(
            input_path,
            per_length_path,
            length,
            reads_per_organism,
            **kwargs,
        )
        if isinstance(count, tuple):
            count = count[0]
        results.append((length, int(count), per_length_path))
        logger.info("Multi-length: wrote %d reads of length %d to %s", int(count), length, per_length_path)
    return results


def _collect_pairs_for_file(
    fp: Path,
    read_length: int,
    insert_size: int,
    insert_size_sd: float,
    n_pairs_limit: int | None,
    prefix: str,
    *,
    seed: int | None = None,
    allow_ambiguous: bool = True,
    substitution_rate: float = 0.0,
    indel_rate: float = 0.0,
    error_model: str | None = None,
    mutation_rng: random.Random | None = None,
    embed_taxonomy: bool = False,
    viral_taxonomy: dict[str, str] | None = None,
    category_override: str | None = None,
) -> list[tuple[SeqRecord, SeqRecord]]:
    """Return a list of (R1, R2) tuples for a single FASTA file.

    Per-file cap ``n_pairs_limit`` counts **pairs** (not individual reads), so
    coverage/abundance/reads-per-organism semantics stay consistent with the
    single-end case: 5 reads/organism = 5 pairs = 10 physical records
    written. Error models, mutations, embed-taxonomy, and ambiguous-base
    filtering are applied to each read independently; if either read fails
    ambiguity checking, the whole pair is dropped (paired-end concordance).
    """
    cat = category_override or _category_from_path(fp)
    rng = random.Random(seed) if seed is not None else None
    apply_platform = bool(error_model) and mutation_rng is not None
    apply_mutations = (
        not apply_platform and (substitution_rate > 0 or indel_rate > 0) and mutation_rng is not None
    )
    pairs: list[tuple[SeqRecord, SeqRecord]] = []
    for rec_idx, record in enumerate(SeqIO.parse(fp, "fasta")):
        if n_pairs_limit is not None and len(pairs) >= n_pairs_limit:
            break
        remaining = None if n_pairs_limit is None else n_pairs_limit - len(pairs)
        record_prefix = prefix if rec_idx == 0 else f"{prefix}_seg{rec_idx}"
        for r1, r2 in chunk_sequence_paired(
            record,
            record_prefix,
            read_length,
            insert_size,
            insert_size_sd,
            remaining,
            rng=rng,
        ):
            r1.id = f"{cat}_{r1.id}"
            r2.id = f"{cat}_{r2.id}"
            if not _is_allowed_sequence(r1, allow_ambiguous) or not _is_allowed_sequence(r2, allow_ambiguous):
                continue
            if apply_platform:
                r1 = _apply_error_model_to_record(r1, error_model, mutation_rng)
                r2 = _apply_error_model_to_record(r2, error_model, mutation_rng)
            elif apply_mutations:
                r1 = _apply_mutations_to_record(r1, substitution_rate, indel_rate, mutation_rng)
                r2 = _apply_mutations_to_record(r2, substitution_rate, indel_rate, mutation_rng)
            if embed_taxonomy:
                r1 = _embed_tax_in_description(r1, cat, prefix, viral_taxonomy)
                r2 = _embed_tax_in_description(r2, cat, prefix, viral_taxonomy)
            pairs.append((r1, r2))
            if n_pairs_limit is not None and len(pairs) >= n_pairs_limit:
                break
    return pairs


def _introduce_chimeras_paired(
    pairs: list[tuple[SeqRecord, SeqRecord]],
    chimera_rate: float,
    rng: random.Random,
) -> tuple[list[tuple[SeqRecord, SeqRecord]], int]:
    """Paired-end analogue of :func:`_introduce_chimeras`.

    Both mates of a chimeric pair come from the same two parent pairs so the
    pair stays physically coherent (R1 is a chimera of the two parents' R1s,
    R2 of their R2s). This matches what a wet-lab PCR chimera would look like
    when sequenced paired-end (both reads from the same hybrid fragment).
    """
    if chimera_rate <= 0 or len(pairs) < 2:
        return pairs, 0
    chimera_rate = min(1.0, chimera_rate)
    n = len(pairs)
    n_chimeras = int(round(n * chimera_rate))
    if n_chimeras == 0:
        return pairs, 0
    indices = list(range(n))
    rng.shuffle(indices)
    to_replace = set(indices[:n_chimeras])
    out = list(pairs)
    for idx in to_replace:
        a_r1, a_r2 = pairs[idx]
        candidates = [i for i in range(n) if i != idx]
        b_idx = rng.choice(candidates)
        b_r1, b_r2 = pairs[b_idx]
        def _splice(left: SeqRecord, right: SeqRecord, mate_tag: str) -> SeqRecord:
            sl = str(left.seq)
            sr = str(right.seq)
            cut_l = len(sl) // 2
            cut_r = len(sr) // 2
            new_seq = sl[:cut_l] + sr[cut_r:]
            rec = SeqRecord(
                Seq(new_seq),
                id=f"chimera_{idx}_pair/{mate_tag}",
                description=f"chimera parents={left.id}|{right.id} mate={mate_tag}",
            )
            ql = left.letter_annotations.get("phred_quality")
            qr = right.letter_annotations.get("phred_quality")
            if ql is not None and qr is not None:
                rec.letter_annotations["phred_quality"] = list(ql[:cut_l]) + list(qr[cut_r:])
            return rec

        new_r1 = _splice(a_r1, b_r1, "1")
        new_r2 = _splice(a_r2, b_r2, "2")
        out[idx] = (new_r1, new_r2)
    return out, n_chimeras


def _introduce_pcr_duplicates_paired(
    pairs: list[tuple[SeqRecord, SeqRecord]],
    duplicate_rate: float,
    rng: random.Random,
) -> tuple[list[tuple[SeqRecord, SeqRecord]], int]:
    """Paired-end analogue of :func:`_introduce_pcr_duplicates`.

    Each pair is duplicated as a unit (both mates) with probability
    ``duplicate_rate`` — real PCR duplicates in paired-end libraries always
    clone both mates together.
    """
    if duplicate_rate <= 0 or not pairs:
        return pairs, 0
    duplicate_rate = min(1.0, duplicate_rate)
    out = list(pairs)
    dup_counter = 0
    for r1, r2 in pairs:
        if rng.random() < duplicate_rate:
            dup_counter += 1
            def _clone(rec: SeqRecord) -> SeqRecord:
                clone = SeqRecord(
                    Seq(str(rec.seq)),
                    id=f"{rec.id}_dup",
                    description=(
                        f"{rec.description} pcr_duplicate=true"
                        if rec.description
                        else "pcr_duplicate=true"
                    ),
                )
                q = rec.letter_annotations.get("phred_quality")
                if q is not None:
                    clone.letter_annotations["phred_quality"] = list(q)
                return clone
            out.append((_clone(r1), _clone(r2)))
    return out, dup_counter


def build_metagenome_paired(
    input_path: Path,
    out_path: Path,
    read_length: int,
    reads_per_organism: int | None,
    *,
    insert_size: int,
    insert_size_sd: float,
    seed: int | None = None,
    allow_ambiguous: bool = True,
    substitution_rate: float = 0.0,
    indel_rate: float = 0.0,
    abundance_profile: dict[str, float] | None = None,
    abundance_distribution: str | None = None,
    viral_taxonomy_json: Path | None = None,
    balance_viral_by_taxonomy: bool = False,
    error_model: str | None = None,
    output_fastq: bool = False,
    write_abundance: bool = False,
    embed_taxonomy: bool = False,
    chimera_rate: float = 0.0,
    pcr_duplicate_rate: float = 0.0,
    coverage: float | None = None,
    coverage_cv: float = 0.0,
    extra_viral_fasta: Path | None = None,
    cap_total_reads: int | None = None,
) -> tuple[int, Path, Path]:
    """Paired-end metagenome generator. Writes ``{stem}_R1.*`` and ``{stem}_R2.*``.

    Intentionally supports a focused subset of :func:`build_metagenome`
    options (coverage, abundance, error models, taxonomy balancing/embedding,
    chimeras, PCR duplicates, extra viral FASTA, write-abundance). EVE
    filtering, similarity filtering, train/test split, multi-length and
    variable-length modes are **not** supported in paired mode and must be
    rejected by the caller. Returns ``(n_pairs, r1_path, r2_path)``; the
    total physical records written is ``2 * n_pairs``.

    ``reads_per_organism`` caps **pairs** per file (so N pairs = 2N records).
    ``cap_total_reads`` caps the **total pair count** after building (again
    so 2 * cap_total_reads records get written).
    """
    if read_length <= 0:
        raise ValueError("read_length must be > 0")
    if insert_size < read_length:
        raise ValueError("insert_size must be >= read_length")
    if insert_size_sd < 0:
        raise ValueError("insert_size_sd must be >= 0")
    if extra_viral_fasta is not None and not extra_viral_fasta.is_file():
        raise FileNotFoundError(f"extra_viral_fasta not found or not a file: {extra_viral_fasta}")
    if output_fastq and not error_model:
        error_model = "illumina"
    use_mutations = substitution_rate > 0 or indel_rate > 0 or bool(error_model)
    if (use_mutations or output_fastq) and seed is None:
        seed = 42

    mutation_rng = random.Random(seed) if use_mutations else None

    if input_path.is_file():
        prefixes_and_files = [(input_path.stem, input_path)]
    else:
        prefixes_and_files = iter_genome_fastas(input_path)
    num_files = len(prefixes_and_files)

    # Per-file pair limit derivation — mirrors build_metagenome's logic but
    # for paired mode the effective "read length" for coverage math is still
    # ``read_length`` (not 2*read_length) because each pair consumes
    # ``insert_size`` bases of genome, of which ``2*read_length`` land in
    # records and the rest is the un-sequenced middle. Using ``read_length``
    # produces the same reads-per-bp that a single-end run would at the same
    # coverage, which is the user-intuitive definition.
    if coverage is not None and coverage > 0:
        read_limits = _compute_coverage_read_limits(
            prefixes_and_files,
            coverage,
            max(1, read_length),
            coverage_cv=coverage_cv,
            abundance_profile=abundance_profile,
            seed=seed,
        )
        if reads_per_organism is not None:
            read_limits = [min(rl, reads_per_organism) for rl in read_limits]
    elif reads_per_organism is None:
        read_limits = [None] * num_files
    else:
        read_limits = _compute_read_limits(
            prefixes_and_files,
            reads_per_organism,
            abundance_profile=abundance_profile,
            abundance_distribution=abundance_distribution,
            seed=seed,
        )

    viral_tax: dict[str, str] | None = None
    if viral_taxonomy_json and viral_taxonomy_json.is_file() and (balance_viral_by_taxonomy or embed_taxonomy):
        from .viral_taxonomy import load_viral_taxonomy
        viral_tax = load_viral_taxonomy(viral_taxonomy_json) or None

    if viral_tax and balance_viral_by_taxonomy:
        stats = get_file_stats(input_path, read_length)
        prefix_to_max_reads = {p: max_r for p, _fp, _tb, max_r in stats} if stats else {}
        read_limits = _apply_viral_taxonomy_balance(
            prefixes_and_files,
            read_limits,
            viral_tax,
            prefix_to_max_reads=prefix_to_max_reads,
        )
        logger.info("Viral taxonomy balance (paired): applied for %d viral prefixes", len(viral_tax))

    all_pairs: list[tuple[SeqRecord, SeqRecord]] = []
    for i, (prefix, fp) in enumerate(prefixes_and_files):
        file_seed = (seed + i) if seed is not None else None
        file_mutation_rng = random.Random(seed + 10000 + i) if mutation_rng is not None else None
        pairs = _collect_pairs_for_file(
            fp,
            read_length,
            insert_size,
            insert_size_sd,
            read_limits[i],
            prefix,
            seed=file_seed,
            allow_ambiguous=allow_ambiguous,
            substitution_rate=substitution_rate,
            indel_rate=indel_rate,
            error_model=error_model,
            mutation_rng=file_mutation_rng,
            embed_taxonomy=embed_taxonomy,
            viral_taxonomy=viral_tax,
        )
        all_pairs.extend(pairs)

    if extra_viral_fasta is not None:
        extra_seed = seed if seed is not None else 42
        extra_pairs = _collect_pairs_for_file(
            extra_viral_fasta,
            read_length,
            insert_size,
            insert_size_sd,
            reads_per_organism,
            "extra_viral",
            seed=extra_seed,
            allow_ambiguous=allow_ambiguous,
            substitution_rate=substitution_rate,
            indel_rate=indel_rate,
            error_model=error_model,
            mutation_rng=random.Random(extra_seed + 10000) if mutation_rng is not None else None,
            embed_taxonomy=embed_taxonomy,
            viral_taxonomy=viral_tax,
            category_override="virus",
        )
        all_pairs.extend(extra_pairs)
        logger.info("Extra viral FASTA (paired): added %d pairs from %s", len(extra_pairs), extra_viral_fasta)

    if cap_total_reads is not None and len(all_pairs) > cap_total_reads:
        rng = random.Random(seed)
        all_pairs = rng.sample(all_pairs, cap_total_reads)

    if chimera_rate > 0 or pcr_duplicate_rate > 0:
        artefact_seed = (seed + 40000) if seed is not None else None
        artefact_rng = random.Random(artefact_seed) if artefact_seed is not None else random.Random()
        if chimera_rate > 0:
            all_pairs, n_chim = _introduce_chimeras_paired(all_pairs, chimera_rate, artefact_rng)
            logger.info("Chimera injection (paired): replaced %d / %d pairs (rate=%.3f)", n_chim, len(all_pairs), chimera_rate)
        if pcr_duplicate_rate > 0:
            before = len(all_pairs)
            all_pairs, n_dup = _introduce_pcr_duplicates_paired(all_pairs, pcr_duplicate_rate, artefact_rng)
            logger.info(
                "PCR duplicates (paired): added %d duplicate pairs (rate=%.3f, base=%d -> %d)",
                n_dup, pcr_duplicate_rate, before, len(all_pairs),
            )

    r1_records = [p[0] for p in all_pairs]
    r2_records = [p[1] for p in all_pairs]

    stem = out_path.stem
    ext = ".fastq" if output_fastq else (out_path.suffix if out_path.suffix else ".fasta")
    r1_path = out_path.with_name(f"{stem}_R1{ext}")
    r2_path = out_path.with_name(f"{stem}_R2{ext}")

    if output_fastq:
        for rec in r1_records:
            add_illumina_qualities_to_record(rec, error_model=error_model)
        for rec in r2_records:
            add_illumina_qualities_to_record(rec, error_model=error_model)
        write_format = "fastq"
    else:
        write_format = "fasta"

    out_path.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(r1_records, r1_path, write_format)
    SeqIO.write(r2_records, r2_path, write_format)

    if write_abundance:
        # Group by genome prefix (strip category + /1 | /2 + _dup, walk back
        # through _pair_N so we reconstruct the file-level prefix).
        from collections import Counter

        def _prefix_of(rec: SeqRecord) -> str:
            rid = rec.id
            if rid.startswith("chimera_"):
                return "chimera"
            # Strip mate suffix
            if "/" in rid:
                rid = rid.rsplit("/", 1)[0]
            if rid.endswith("_dup"):
                rid = rid[: -len("_dup")]
            if "_pair_" in rid:
                rid = rid.rsplit("_pair_", 1)[0]
            for cat in ("bacteria", "virus", "archaea", "plasmid"):
                marker = f"{cat}_"
                if rid.startswith(marker):
                    return rid[len(marker):]
            return rid

        counts: Counter[str] = Counter()
        # One pair contributes 1 to the genome's count (so values match the
        # single-end ``write_abundance`` semantics of "n reads" where a read
        # is a physical fragment, not a physical record).
        for rec in r1_records:
            counts[_prefix_of(rec)] += 1
        total = sum(counts.values()) or 1
        abundance_path = out_path.with_name(f"{stem}_abundance.txt")
        with open(abundance_path, "w") as f:
            f.write("genome_id\tread_count\tproportion\n")
            for p in sorted(counts):
                n = counts[p]
                f.write(f"{p}\t{n}\t{n / total}\n")

    return len(all_pairs), r1_path, r2_path


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
    write_fastq: bool = False,
) -> tuple[int, int]:
    """Split records into train (train_pct%%) and test; remove from test any sequence similar to train; write train and test FASTA or FASTQ.

    Similarity check: BLAST (megablast task, faster for high identity) of test vs train DB in batches;
    test sequences with a hit >= similarity_threshold over min_coverage of length are dropped.
    Returns (n_train, n_test_after_filter). If write_fastq, outputs .fastq (records must have phred_quality); else .fasta.
    """
    from .similarity_filter import filter_candidates_against_kept

    if not records:
        return 0, 0
    train_pct = normalize_train_split_percent(train_pct)
    train_pct = max(0.0, min(100.0, train_pct))
    n_train_target = max(1, int(len(records) * train_pct / 100.0))
    n_test_target = len(records) - n_train_target
    ext = "fastq" if write_fastq else "fasta"
    fmt = "fastq" if write_fastq else "fasta"
    if n_train_target <= 0 or n_test_target <= 0:
        logger.warning("Train-test split: not enough sequences (%d); writing all to train", len(records))
        train_path = output_dir / f"{output_stem}_train.{ext}"
        SeqIO.write(records, train_path, fmt)
        (output_dir / f"{output_stem}_test.{ext}").write_text("")
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
    train_path = output_dir / f"{output_stem}_train.{ext}"
    test_path = output_dir / f"{output_stem}_test.{ext}"
    SeqIO.write(train_list, train_path, fmt)
    SeqIO.write(test_filtered, test_path, fmt)
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
        "--contig-quality-profile",
        type=str,
        default=None,
        choices=list(CONTIG_QUALITY_PRESET_NAMES),
        help="Preset non-overlapping contig-length stratification (low/medium/high quality bins). Mutually exclusive with --min-contig-length/--max-contig-length.",
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
        help="Split output into train and test: accepts train percentage (e.g. 80) or fraction (e.g. 0.8). Test sequences similar to train are removed. Writes {output_stem}_train.fasta and {output_stem}_test.fasta.",
    )
    parser.add_argument(
        "--train-test-similarity-threshold",
        type=float,
        default=90.0,
        help="Max BLASTN percent identity for train-test: remove from test if similar to train above this. Default: 90",
    )
    parser.add_argument(
        "--output-fastq",
        action="store_true",
        help="Write FASTQ instead of FASTA, with per-base Phred quality scores (Illumina-like position-dependent). Use --seed for reproducibility.",
    )
    parser.add_argument(
        "--write-abundance",
        action="store_true",
        help="Write a tab-separated abundance file (genome_id, read_count, proportion) next to the output for ground-truth benchmarking.",
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
    if args.contig_quality_profile is not None and (args.min_contig_length is not None or args.max_contig_length is not None):
        parser.error("--contig-quality-profile cannot be combined with --min-contig-length/--max-contig-length")
    if use_variable and (args.min_contig_length < 1 or args.max_contig_length < args.min_contig_length):
        parser.error("--min-contig-length and --max-contig-length must be positive and min <= max")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    out_path = args.output_dir / args.output

    effective_length = args.sequence_length
    min_len, max_len = None, None
    quality_profile = args.contig_quality_profile
    if quality_profile is not None:
        effective_length = contig_quality_profile_mean_length(quality_profile)
        print(f"Contig-quality profile: {quality_profile} (mean length ~{effective_length} bp)")
    elif use_variable:
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
            contig_quality_profile=quality_profile,
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
    output_fastq_flag = getattr(args, "output_fastq", False)
    result = build_metagenome(
        args.input,
        out_path,
        args.sequence_length,
        reads_per_organism,
        min_length=min_len,
        max_length=max_len,
        contig_quality_profile=quality_profile,
        seed=args.seed,
        cap_total_reads=args.cap_total_reads,
        eve_intervals=eve_intervals,
        allow_ambiguous=allow_ambiguous,
        filter_similar=getattr(args, "filter_similar", False),
        similarity_threshold=getattr(args, "similarity_threshold", 90.0),
        similarity_min_coverage=getattr(args, "similarity_min_coverage", 0.8),
        oversample_factor=getattr(args, "oversample_factor", 2.0),
        return_records=do_train_test_split,
        output_fastq=output_fastq_flag,
        write_abundance=getattr(args, "write_abundance", False),
    )
    if do_train_test_split:
        _count, records = result
        output_stem = Path(args.output).stem
        train_split_pct = normalize_train_split_percent(args.train_test_split)
        n_train, n_test = split_train_test_and_write(
            records,
            train_split_pct,
            args.seed,
            args.output_dir,
            output_stem,
            similarity_threshold=getattr(args, "train_test_similarity_threshold", 90.0),
            similarity_min_coverage=0.8,
            write_fastq=output_fastq_flag,
        )
        ext = "fastq" if output_fastq_flag else "fasta"
        train_path = args.output_dir / f"{output_stem}_train.{ext}"
        test_path = args.output_dir / f"{output_stem}_test.{ext}"
        print(f"Train-test split ({train_split_pct}% train): wrote {n_train} to {train_path}, {n_test} to {test_path}")
    else:
        count = result
        if output_fastq_flag:
            out_path = out_path if out_path.suffix.lower() == ".fastq" else out_path.with_suffix(".fastq")
        if args.cap_total_reads is not None and count == args.cap_total_reads:
            print(f"Wrote {count} sequences (capped) to {out_path}")
        else:
            print(f"Wrote {count} sequences to {out_path}")


def main() -> None:
    _cli()


if __name__ == "__main__":
    main()
