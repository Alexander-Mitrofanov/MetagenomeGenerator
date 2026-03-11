#!/usr/bin/env python3
from __future__ import annotations
"""Unified CLI for the metagenome simulator.

Subcommands: download, snapshot, migrate-snapshot, chunk, pipeline, blastn-filter, seeker, temporal-split, temporal-split-info.
"""

import argparse
import logging
from pathlib import Path

from .accession_snapshot import get_default_snapshot_path, migrate_snapshot_to_categories, run_snapshot
from .download_genomes import download_genomes
from .chunk_genomes import build_metagenome, get_file_stats, split_train_test_and_write
from .seeker_wrapper import SEEKER_MIN_LENGTH, run_seeker
from .blastn_filter import export_eve_regions_fasta, load_eve_intervals, run_blastn_from_dirs
from .genome_layout import validate_genome_dir
from .temporal_split import run_temporal_split, run_temporal_split_info
from .viral_taxonomy import run_viral_taxonomy

# Organized output layout (pipeline): one root dir with step-based subdirs for easy navigation.
OUTPUT_DIR_DOWNLOADED = "downloaded"
OUTPUT_DIR_BLASTN = "blastn"
OUTPUT_DIR_METAGENOME = "metagenome"
OUTPUT_DIR_SEEKER = "seeker"
OUTPUT_DIR_LOGS = "logs"


def _parse_abundance_profile(s: str | None) -> dict[str, float] | None:
    """Parse KEY=VAL,KEY2=VAL2 into dict. Returns None if s is None or empty."""
    if not s or not s.strip():
        return None
    out: dict[str, float] = {}
    for part in s.split(","):
        part = part.strip()
        if "=" in part:
            k, v = part.split("=", 1)
            out[k.strip()] = float(v.strip())
    return out if out else None


def _add_download_subparser(subparsers) -> None:
    p = subparsers.add_parser(
        "download",
        help="Download bacterial and viral genomes from NCBI",
    )
    p.add_argument(
        "--num-organisms",
        type=int,
        default=10,
        help="Number of organisms to download per group (bacterial and viral). Default: 10",
    )
    p.add_argument(
        "--output-dir",
        type=Path,
        default=Path("output"),
        help="Directory for downloaded genome FASTAs. Default: output",
    )
    p.add_argument(
        "--num-archaea",
        type=int,
        default=0,
        help="Number of archaeal genomes (negative samples). Default: 0",
    )
    p.add_argument(
        "--num-plasmid",
        type=int,
        default=0,
        help="Number of plasmids (negative samples). Default: 0",
    )
    p.add_argument(
        "--accessions-file",
        type=Path,
        default=None,
        metavar="PATH",
        help="Load accession IDs from JSON (skip NCBI search) for reproducibility.",
    )
    p.add_argument(
        "--save-accessions",
        type=Path,
        default=None,
        metavar="PATH",
        help="After searching NCBI, save accession list and UTC timestamp to JSON for later --accessions-file runs.",
    )
    p.add_argument(
        "--complete-only",
        action="store_true",
        help="When searching NCBI (no --accessions-file), restrict to complete genomes only (exclude WGS/draft). For reproducible runs use a snapshot created with snapshot --complete-only.",
    )
    p.set_defaults(func=_run_download)


def _add_snapshot_subparser(subparsers) -> None:
    p = subparsers.add_parser(
        "snapshot",
        help="Catalog all NCBI accession IDs for bacteria, virus, archaea, plasmid (no downloads). Saves to snapshots/; use with --accessions-file.",
    )
    p.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output JSON path. Default: snapshots/accession_snapshot_YYYY-MM-DD.json (run date)",
    )
    p.add_argument(
        "--no-metadata",
        action="store_true",
        help="Do not fetch CreateDate and title per accession (lists only).",
    )
    p.add_argument(
        "--metadata-batch-size",
        type=int,
        default=500,
        help="Batch size for esummary when fetching per-accession metadata. Default: 500",
    )
    p.add_argument(
        "--log",
        type=Path,
        default=None,
        metavar="PATH",
        help="Write progress log to this file. Default: snapshot_YYYY-MM-DD.log next to the output JSON.",
    )
    p.add_argument(
        "--complete-only",
        action="store_true",
        help="Restrict to complete genomes only (NCBI: complete[Properties], NOT WGS[Properties]). Use this snapshot with --accessions-file for reproducible complete-only runs.",
    )
    p.set_defaults(func=_run_snapshot)


def _add_migrate_snapshot_subparser(subparsers) -> None:
    p = subparsers.add_parser(
        "migrate-snapshot",
        help="Convert a snapshot JSON to category-inline format (metadata per category; removes ncbi_db_info). No re-download.",
    )
    p.add_argument(
        "input",
        type=Path,
        metavar="SNAPSHOT_JSON",
        help="Path to snapshot JSON (legacy format with accession_metadata and/or ncbi_db_info).",
    )
    p.set_defaults(func=_run_migrate_snapshot)


def _run_migrate_snapshot(args) -> None:
    migrate_snapshot_to_categories(args.input)
    print(f"Migrated {args.input} to category-inline format.")


def _add_chunk_subparser(subparsers) -> None:
    p = subparsers.add_parser(
        "chunk",
        help="Chunk genomes into fixed-length reads and write a metagenome FASTA",
    )
    p.add_argument(
        "--input",
        type=Path,
        required=True,
        help="Input FASTA file or directory of genome FASTAs",
    )
    p.add_argument(
        "--output",
        type=str,
        required=True,
        help="Output filename for the metagenome FASTA",
    )
    p.add_argument(
        "--output-dir",
        type=Path,
        default=Path("output"),
        help="Directory to write the metagenome FASTA. Default: output",
    )
    p.add_argument(
        "--sequence-length",
        type=int,
        default=250,
        help="Read length in nucleotides. Default: 250",
    )
    p.add_argument(
        "--reads-per-organism",
        type=int,
        default=None,
        help="Max reads to extract per organism (per input file). Default: all reads",
    )
    p.add_argument(
        "--balanced",
        action="store_true",
        help="Use equal reads per file (min of max possible reads across files).",
    )
    p.add_argument(
        "--min-contig-length",
        type=int,
        default=None,
        help="Variable-length mode: min contig length (e.g. 300). Use with --max-contig-length.",
    )
    p.add_argument(
        "--max-contig-length",
        type=int,
        default=None,
        help="Variable-length mode: max contig length (e.g. 2000). Use with --min-contig-length.",
    )
    p.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed for variable-length or cap-total-reads (reproducibility).",
    )
    p.add_argument(
        "--cap-total-reads",
        type=int,
        default=None,
        help="Downsample output to this many reads (e.g. to match positive set size).",
    )
    p.add_argument(
        "--eve-intervals",
        type=Path,
        default=None,
        help="Path to eve_intervals.json (from blastn-filter); exclude chunks overlapping EVE regions.",
    )
    p.add_argument(
        "--forbid-ambiguous",
        action="store_true",
        help="Discard chunks containing ambiguous bases (non-ACGT, e.g. N). By default, such chunks are kept.",
    )
    p.add_argument(
        "--substitution-rate",
        type=float,
        default=0.0,
        metavar="R",
        help="Per-base substitution rate (0–1). Apply random base changes for robustness benchmarks. Use --seed for reproducibility. Default: 0",
    )
    p.add_argument(
        "--indel-rate",
        type=float,
        default=0.0,
        metavar="R",
        help="Per-base indel rate (0–1). 50%% insert, 50%% delete; may change chunk length. Use with --seed. Default: 0",
    )
    p.add_argument(
        "--extra-viral-fasta",
        type=Path,
        default=None,
        metavar="PATH",
        help="Extra FASTA of viral sequences (e.g. metavirome contigs) to chunk and merge with RefSeq viral chunks. Multi-record OK.",
    )
    p.add_argument(
        "--abundance-profile",
        type=str,
        default=None,
        metavar="KEY=VAL,...",
        help="Per-category read weights, e.g. bacterial=0.5,viral=2,archaea=1,plasmid=1. Scales reads per file by category.",
    )
    p.add_argument(
        "--abundance-distribution",
        type=str,
        default=None,
        choices=["exponential"],
        help="Per-genome abundance: exponential draws weights from Exp(1); use --seed for reproducibility.",
    )
    p.add_argument(
        "--viral-taxonomy",
        type=Path,
        default=None,
        metavar="PATH",
        help="JSON mapping viral_1, viral_2, ... to taxonomy group (e.g. family). Use with --balance-viral-by-taxonomy.",
    )
    p.add_argument(
        "--balance-viral-by-taxonomy",
        action="store_true",
        help="Balance viral reads by taxonomy group so each group contributes equally. Requires --viral-taxonomy.",
    )
    p.set_defaults(func=_run_chunk)


def _add_pipeline_subparser(subparsers) -> None:
    p = subparsers.add_parser(
        "pipeline",
        help="Run download + chunk (+ optional BLASTN, Seeker). Uses organized layout: output_dir/downloaded, blastn, metagenome, seeker.",
    )
    p.add_argument(
        "--num-organisms",
        type=int,
        default=10,
        help="Number of organisms to download per group (bacterial and viral). Default: 10",
    )
    p.add_argument(
        "--output-dir",
        type=Path,
        default=Path("output"),
        help="Root output directory. Pipeline creates: downloaded/, blastn/, metagenome/, seeker/ under it. Default: output",
    )
    p.add_argument(
        "--genome-dir",
        type=Path,
        default=None,
        metavar="PATH",
        help="Use existing genome data at PATH instead of downloading. PATH must contain bacteria/, virus/, archaea/, plasmid/ (virus/ non-empty; at least one of bacteria/, archaea/, plasmid/ non-empty). Skips Step 1.",
    )
    p.add_argument(
        "--output",
        type=str,
        default="metagenome.fasta",
        help="Metagenome FASTA filename (written under output-dir/metagenome/). Default: metagenome.fasta",
    )
    p.add_argument(
        "--sequence-length",
        type=int,
        default=250,
        help="Read length in nucleotides. Default: 250",
    )
    p.add_argument(
        "--reads-per-organism",
        type=int,
        default=1000,
        help="Max reads to extract per organism (per input file). Default: 1000",
    )
    p.add_argument(
        "--balanced",
        action="store_true",
        help="Use equal reads per file (min of max possible reads across files).",
    )
    p.add_argument(
        "--num-archaea",
        type=int,
        default=0,
        help="Number of archaeal genomes to download. Default: 0",
    )
    p.add_argument(
        "--num-plasmid",
        type=int,
        default=0,
        help="Number of plasmids to download. Default: 0",
    )
    p.add_argument(
        "--min-contig-length",
        type=int,
        default=None,
        help="Variable-length contigs: min length (e.g. 300). Use with --max-contig-length.",
    )
    p.add_argument(
        "--max-contig-length",
        type=int,
        default=None,
        help="Variable-length contigs: max length (e.g. 2000). Use with --min-contig-length.",
    )
    p.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed for variable-length or cap-total-reads.",
    )
    p.add_argument(
        "--cap-total-reads",
        type=int,
        default=None,
        help="Downsample metagenome to this many reads.",
    )
    p.add_argument(
        "--run-seeker",
        action="store_true",
        help="After building the metagenome FASTA, run Seeker (predict-metagenome) on it",
    )
    p.add_argument(
        "--seeker-conda-env",
        type=str,
        default="seeker",
        help="Conda env name to run Seeker in (used via `conda run -n ...`). Default: seeker",
    )
    p.add_argument(
        "--seeker-threshold",
        type=float,
        default=0.5,
        help="Score threshold for exporting predicted phage reads. Default: 0.5",
    )
    p.add_argument(
        "--run-blastn-filter",
        action="store_true",
        help="Before chunking, run BLASTN to find EVEs in non-viral genomes and exclude those regions when chunking.",
    )
    p.add_argument(
        "--blastn-out-dir",
        type=Path,
        default=None,
        help="Directory for BLASTN output (default: output-dir/blastn).",
    )
    p.add_argument(
        "--blastn-evalue",
        type=float,
        default=1e-5,
        help="BLASTN E-value threshold. Default: 1e-5",
    )
    p.add_argument(
        "--blastn-perc-identity",
        type=float,
        default=70.0,
        help="BLASTN percent identity threshold. Default: 70",
    )
    p.add_argument(
        "--blastn-export-eve-fasta",
        type=Path,
        default=None,
        metavar="PATH",
        help="Optional FASTA output path to write EVE/provirus intervals as sequences (one record per interval).",
    )
    p.add_argument(
        "--blastn-export-eve-min-length",
        type=int,
        default=200,
        metavar="N",
        help="Minimum interval length to export with --blastn-export-eve-fasta. Default: 200",
    )
    p.add_argument(
        "--filter-similar",
        action="store_true",
        help="Filter chunks by similarity: drop sequences >=90%% similar to already-kept; oversample and refill to target.",
    )
    p.add_argument(
        "--similarity-threshold",
        type=float,
        default=90.0,
        help="Max BLASTN percent identity for 'similar' (chunks above this vs kept set are removed). Default: 90",
    )
    p.add_argument(
        "--similarity-min-coverage",
        type=float,
        default=0.8,
        help="Min fraction of query length in alignment for similarity. Default: 0.8",
    )
    p.add_argument(
        "--oversample-factor",
        type=float,
        default=2.0,
        help="When --filter-similar: generate up to this many times target count before filtering. Default: 2.0",
    )
    p.add_argument(
        "--train-test-split",
        type=float,
        default=None,
        metavar="PCT",
        help="Split metagenome into train (PCT%%) and test; remove from test sequences similar to train. Writes {stem}_train.fasta and {stem}_test.fasta.",
    )
    p.add_argument(
        "--train-test-similarity-threshold",
        type=float,
        default=90.0,
        help="Max BLASTN percent identity for train-test: remove from test if similar to train above this. Default: 90",
    )
    p.add_argument(
        "--train-test-blast-threads",
        type=int,
        default=4,
        help="Number of threads for BLAST when filtering test vs train. Default: 4",
    )
    p.add_argument(
        "--train-test-blast-batch-size",
        type=int,
        default=2000,
        help="Test sequences per BLAST run (larger = fewer runs, more memory). Default: 2000",
    )
    p.add_argument(
        "--accessions-file",
        type=Path,
        default=None,
        metavar="PATH",
        help="Load accession IDs from JSON for download step (reproducibility; skip NCBI search).",
    )
    p.add_argument(
        "--save-accessions",
        type=Path,
        default=None,
        metavar="PATH",
        help="After NCBI search, save accession list and UTC timestamp to JSON for later --accessions-file runs.",
    )
    p.add_argument(
        "--complete-only",
        action="store_true",
        help="When downloading (no --accessions-file), restrict to complete genomes only (exclude WGS/draft). Use snapshot --complete-only for reproducible complete-only runs.",
    )
    p.add_argument(
        "--forbid-ambiguous",
        action="store_true",
        help="Discard chunks containing ambiguous bases (non-ACGT, e.g. N). By default, such chunks are kept.",
    )
    p.add_argument(
        "--substitution-rate",
        type=float,
        default=0.0,
        metavar="R",
        help="Per-base substitution rate (0–1) when chunking. Use --seed for reproducibility. Default: 0",
    )
    p.add_argument(
        "--indel-rate",
        type=float,
        default=0.0,
        metavar="R",
        help="Per-base indel rate (0–1) when chunking; may change chunk length. Default: 0",
    )
    p.add_argument(
        "--extra-viral-fasta",
        type=Path,
        default=None,
        metavar="PATH",
        help="Extra FASTA of viral sequences (e.g. metavirome contigs) to chunk and merge with RefSeq viral chunks.",
    )
    p.add_argument(
        "--abundance-profile",
        type=str,
        default=None,
        metavar="KEY=VAL,...",
        help="Per-category read weights, e.g. bacterial=0.5,viral=2,archaea=1,plasmid=1.",
    )
    p.add_argument(
        "--abundance-distribution",
        type=str,
        default=None,
        choices=["exponential"],
        help="Per-genome abundance: exponential weights (use --seed for reproducibility).",
    )
    p.add_argument(
        "--viral-taxonomy",
        type=Path,
        default=None,
        metavar="PATH",
        help="JSON mapping viral_1, viral_2, ... to taxonomy group. Use with --balance-viral-by-taxonomy.",
    )
    p.add_argument(
        "--balance-viral-by-taxonomy",
        action="store_true",
        help="Balance viral reads by taxonomy group (requires --viral-taxonomy).",
    )
    p.set_defaults(func=_run_pipeline)


def _add_blastn_filter_subparser(subparsers) -> None:
    p = subparsers.add_parser(
        "blastn-filter",
        help="Run BLASTN (non-viral vs viral), write eve_intervals.json for use with chunk --eve-intervals",
    )
    p.add_argument(
        "--genome-dir",
        type=Path,
        required=True,
        help="Directory containing bacterial_*, viral_*, and optionally archaea_*, plasmid_*.fasta",
    )
    p.add_argument(
        "--out-dir",
        type=Path,
        required=True,
        help="Output directory for viral_concat.fasta, blastn/, and eve_intervals.json",
    )
    p.add_argument(
        "--evalue",
        type=float,
        default=1e-5,
        help="BLASTN E-value. Default: 1e-5",
    )
    p.add_argument(
        "--perc-identity",
        type=float,
        default=70.0,
        help="BLASTN percent identity. Default: 70",
    )
    p.add_argument(
        "--export-eve-fasta",
        type=Path,
        default=None,
        metavar="PATH",
        help="Optional FASTA output path to write EVE/provirus intervals as sequences (one record per interval).",
    )
    p.add_argument(
        "--export-eve-min-length",
        type=int,
        default=200,
        metavar="N",
        help="Minimum interval length to export with --export-eve-fasta. Default: 200",
    )
    p.set_defaults(func=_run_blastn_filter)


def _add_seeker_subparser(subparsers) -> None:
    p = subparsers.add_parser(
        "seeker",
        help="Run Seeker (predict-metagenome) on a metagenome FASTA",
    )
    p.add_argument(
        "--input",
        type=Path,
        required=True,
        help="Input metagenome FASTA (e.g. output from chunk or pipeline)",
    )
    p.add_argument(
        "--output-dir",
        type=Path,
        default=Path("output"),
        help="Directory to write Seeker outputs. Default: output",
    )
    p.add_argument(
        "--threshold",
        type=float,
        default=0.5,
        help="Score threshold for exporting predicted phage reads (0-1). Default: 0.5",
    )
    p.add_argument(
        "--conda-env",
        type=str,
        default="seeker",
        help="Conda env name to run Seeker in (used via `conda run -n ...`). Default: seeker",
    )
    p.add_argument(
        "--predictions-tsv",
        type=str,
        default=None,
        help="Filename for Seeker predictions TSV (name/prediction/score). Default: derived from input",
    )
    p.add_argument(
        "--phage-fasta",
        type=str,
        default=None,
        help="Filename for FASTA of reads with score >= threshold. Default: derived from input",
    )
    p.add_argument(
        "--min-length",
        type=int,
        default=200,
        help="Filter out reads shorter than this before running Seeker (Seeker requires >=200). Default: 200",
    )
    p.set_defaults(func=_run_seeker)


def _add_temporal_split_subparser(subparsers) -> None:
    p = subparsers.add_parser(
        "temporal-split",
        help="Split accession list by NCBI CreateDate into train (before date) and test (on/after). Outputs two JSONs for use with download --accessions-file.",
    )
    p.add_argument(
        "--accessions-file",
        type=Path,
        required=True,
        metavar="PATH",
        help="Input JSON with bacterial, viral, archaea, plasmid accession lists (e.g. from snapshot).",
    )
    p.add_argument(
        "--split-date",
        type=str,
        required=True,
        metavar="YYYY-MM-DD",
        help="Cutoff date (YYYY-MM-DD). Train = CreateDate < this date, test = CreateDate >= this date (e.g. 2015-05-01 for DeepVirFinder-style split).",
    )
    p.add_argument(
        "--output-train",
        type=Path,
        default=None,
        metavar="PATH",
        help="Output path for train accessions JSON. Default: same dir as accessions-file, name train_<basename>.",
    )
    p.add_argument(
        "--output-test",
        type=Path,
        default=None,
        metavar="PATH",
        help="Output path for test accessions JSON. Default: same dir as accessions-file, name test_<basename>.",
    )
    p.add_argument(
        "--batch-size",
        type=int,
        default=200,
        help="NCBI esummary batch size (rate-limited). Default: 200",
    )
    p.set_defaults(func=_run_temporal_split)


def _add_temporal_split_info_subparser(subparsers) -> None:
    p = subparsers.add_parser(
        "temporal-split-info",
        help="Show train/test counts for a temporal split by date (no files written). Use before running temporal-split.",
    )
    p.add_argument(
        "--accessions-file",
        type=Path,
        required=True,
        metavar="PATH",
        help="Input JSON with bacterial, viral, archaea, plasmid lists (e.g. from snapshot).",
    )
    p.add_argument(
        "--split-date",
        type=str,
        required=True,
        metavar="YYYY-MM-DD",
        help="Cutoff date (YYYY-MM-DD). Train = CreateDate < this, test = CreateDate >= this.",
    )
    p.add_argument(
        "--batch-size",
        type=int,
        default=500,
        help="NCBI esummary batch size when metadata not in snapshot. Default: 500",
    )
    p.set_defaults(func=_run_temporal_split_info)


def _add_viral_taxonomy_subparser(subparsers) -> None:
    p = subparsers.add_parser(
        "viral-taxonomy",
        help="Fetch viral taxonomy from NCBI and write viral_1 -> family (or realm) JSON for --viral-taxonomy / --balance-viral-by-taxonomy.",
    )
    p.add_argument(
        "--accessions-file",
        type=Path,
        required=True,
        metavar="PATH",
        help="Input JSON with viral accession list (e.g. from snapshot).",
    )
    p.add_argument(
        "--output",
        type=Path,
        required=True,
        metavar="PATH",
        help="Output JSON path (e.g. viral_taxonomy.json).",
    )
    p.add_argument(
        "--level",
        type=str,
        default="family",
        choices=["family", "realm"],
        help="Taxonomy level to use as group. Default: family",
    )
    p.add_argument(
        "--batch-size",
        type=int,
        default=100,
        help="Batch size for NCBI elink/efetch. Default: 100",
    )
    p.set_defaults(func=_run_viral_taxonomy)


def _run_viral_taxonomy(args) -> None:
    if not args.accessions_file.exists():
        raise SystemExit(f"--accessions-file not found: {args.accessions_file}")
    run_viral_taxonomy(
        args.accessions_file,
        args.output,
        level=args.level,
        batch_size=args.batch_size,
    )


def _run_temporal_split_info(args) -> None:
    run_temporal_split_info(
        args.accessions_file,
        args.split_date,
        batch_size=getattr(args, "batch_size", 500),
        verbose=True,
    )


def _run_download(args) -> None:
    download_genomes(
        args.num_organisms,
        args.output_dir,
        num_archaea=getattr(args, "num_archaea", 0),
        num_plasmid=getattr(args, "num_plasmid", 0),
        accessions_file=getattr(args, "accessions_file", None),
        save_accessions_to=getattr(args, "save_accessions", None),
        complete_only=getattr(args, "complete_only", False),
    )


def _run_snapshot(args) -> None:
    output = getattr(args, "output", None) or get_default_snapshot_path()
    run_snapshot(
        output,
        fetch_metadata=not getattr(args, "no_metadata", False),
        metadata_batch_size=getattr(args, "metadata_batch_size", 500),
        log_path=getattr(args, "log", None),
        complete_only=getattr(args, "complete_only", False),
    )


def _run_chunk(args) -> None:
    if args.input.is_dir():
        ok, err = validate_genome_dir(args.input)
        if not ok:
            print(f"Error: {err}", file=__import__("sys").stderr)
            raise SystemExit(1)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    out_path = args.output_dir / args.output

    min_len = getattr(args, "min_contig_length", None)
    max_len = getattr(args, "max_contig_length", None)
    use_variable = min_len is not None and max_len is not None
    effective_length = (min_len + max_len) // 2 if use_variable else args.sequence_length

    reads_per_organism = args.reads_per_organism
    if args.balanced:
        stats = get_file_stats(
            args.input,
            effective_length,
            min_length=min_len,
            max_length=max_len,
        )
        if not stats:
            raise SystemExit("No FASTA files found under --input")
        for prefix, _fp, total_bases, max_reads in stats:
            print(f"  {prefix}: {total_bases} bp -> ~{max_reads} reads")
        reads_per_organism = min(max_reads for (_, _, _, max_reads) in stats)
        print(f"  Balanced: using {reads_per_organism} reads per file.")

    eve_intervals = None
    if getattr(args, "eve_intervals", None) is not None:
        eve_path = args.eve_intervals
        if not eve_path.exists():
            raise SystemExit(f"--eve-intervals not found: {eve_path}")
        eve_intervals = load_eve_intervals(eve_path)
        print(f"EVE intervals loaded for {len(eve_intervals)} sequences.")

    allow_ambiguous = not getattr(args, "forbid_ambiguous", False)
    sub_rate = getattr(args, "substitution_rate", 0.0)
    indel_r = getattr(args, "indel_rate", 0.0)
    if (sub_rate > 0 or indel_r > 0) and getattr(args, "seed", None) is None:
        args.seed = 42
    extra_viral = getattr(args, "extra_viral_fasta", None)
    if extra_viral is not None and not extra_viral.exists():
        raise SystemExit(f"--extra-viral-fasta not found: {extra_viral}")
    abundance_profile = _parse_abundance_profile(getattr(args, "abundance_profile", None))
    abundance_dist = getattr(args, "abundance_distribution", None)
    viral_tax_path = getattr(args, "viral_taxonomy", None)
    balance_viral = getattr(args, "balance_viral_by_taxonomy", False)
    if balance_viral and (viral_tax_path is None or not viral_tax_path.exists()):
        raise SystemExit("--balance-viral-by-taxonomy requires --viral-taxonomy PATH to an existing JSON file.")
    count = build_metagenome(
        args.input,
        out_path,
        args.sequence_length,
        reads_per_organism,
        min_length=min_len,
        max_length=max_len,
        seed=getattr(args, "seed", None),
        cap_total_reads=getattr(args, "cap_total_reads", None),
        eve_intervals=eve_intervals,
        allow_ambiguous=allow_ambiguous,
        substitution_rate=sub_rate,
        indel_rate=indel_r,
        extra_viral_fasta=extra_viral,
        abundance_profile=abundance_profile,
        abundance_distribution=abundance_dist,
        viral_taxonomy_json=viral_tax_path,
        balance_viral_by_taxonomy=balance_viral,
    )
    print(f"Wrote {count} sequences to {out_path}")


def _run_blastn_filter(args) -> None:
    ok, err = validate_genome_dir(args.genome_dir)
    if not ok:
        print(f"Error: {err}", file=__import__("sys").stderr)
        raise SystemExit(1)
    eve_intervals = run_blastn_from_dirs(
        args.genome_dir,
        args.out_dir,
        evalue=args.evalue,
        perc_identity=args.perc_identity,
    )
    export_path = getattr(args, "export_eve_fasta", None)
    if export_path is not None:
        from .genome_layout import get_nonviral_fasta_paths
        nonviral = get_nonviral_fasta_paths(args.genome_dir)
        export_eve_regions_fasta(
            nonviral,
            eve_intervals,
            export_path,
            min_interval_length=getattr(args, "export_eve_min_length", 200),
        )


def _run_pipeline(args) -> None:
    base = args.output_dir
    download_dir = base / OUTPUT_DIR_DOWNLOADED
    blastn_dir = base / OUTPUT_DIR_BLASTN
    metagenome_dir = base / OUTPUT_DIR_METAGENOME
    seeker_dir = base / OUTPUT_DIR_SEEKER
    log_dir = base / OUTPUT_DIR_LOGS

    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / "pipeline.log"
    file_handler = logging.FileHandler(log_file, mode="w", encoding="utf-8")
    file_handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(name)s: %(message)s"))
    file_handler.setLevel(logging.DEBUG)
    root = logging.getLogger()
    root.addHandler(file_handler)
    root.setLevel(logging.DEBUG)
    plog = logging.getLogger("pipeline")
    plog.info("Pipeline started; output root=%s", base.resolve())

    genome_dir_provided = getattr(args, "genome_dir", None) is not None
    if genome_dir_provided:
        download_dir = Path(args.genome_dir).resolve()
        ok, err = validate_genome_dir(download_dir)
        if not ok:
            plog.error("Genome dir validation failed: %s", err)
            print(f"Error: {err}", file=__import__("sys").stderr)
            raise SystemExit(1)
        plog.info("Step 1: Download — skipped (using existing genome dir: %s)", download_dir)
        print(f"Using existing genome data: {download_dir}")
    else:
        download_dir = base / OUTPUT_DIR_DOWNLOADED
        download_dir.mkdir(parents=True, exist_ok=True)
        metagenome_dir.mkdir(parents=True, exist_ok=True)
        plog.info("Step 1: Download genomes (num_organisms=%s, num_archaea=%s, num_plasmid=%s)",
                  args.num_organisms, getattr(args, "num_archaea", 0), getattr(args, "num_plasmid", 0))
        download_genomes(
            args.num_organisms,
            download_dir,
            num_archaea=getattr(args, "num_archaea", 0),
            num_plasmid=getattr(args, "num_plasmid", 0),
            accessions_file=getattr(args, "accessions_file", None),
            save_accessions_to=getattr(args, "save_accessions", None),
            complete_only=getattr(args, "complete_only", False),
        )

    metagenome_dir.mkdir(parents=True, exist_ok=True)

    eve_intervals = None
    if getattr(args, "run_blastn_filter", False):
        plog.info("Step 2: BLASTN filter (EVE removal); evalue=%s, perc_identity=%s",
                  getattr(args, "blastn_evalue", 1e-5), getattr(args, "blastn_perc_identity", 70.0))
        blastn_out = getattr(args, "blastn_out_dir", None) or blastn_dir
        blastn_out.mkdir(parents=True, exist_ok=True)
        eve_intervals = run_blastn_from_dirs(
            download_dir,
            blastn_out,
            evalue=getattr(args, "blastn_evalue", 1e-5),
            perc_identity=getattr(args, "blastn_perc_identity", 70.0),
        )
        export_path = getattr(args, "blastn_export_eve_fasta", None)
        if export_path is not None:
            from .genome_layout import get_nonviral_fasta_paths
            nonviral = get_nonviral_fasta_paths(download_dir)
            export_eve_regions_fasta(
                nonviral,
                eve_intervals,
                export_path,
                min_interval_length=getattr(args, "blastn_export_eve_min_length", 200),
            )
        plog.info("BLASTN filter: eve_intervals loaded for %d sequences", len(eve_intervals))
        print(f"EVE intervals loaded for chunk step ({len(eve_intervals)} sequences).")
    else:
        plog.info("Step 2: BLASTN filter (EVE removal) — skipped (use --run-blastn-filter to enable)")

    plog.info("Step 3: Chunk genomes -> %s", metagenome_dir / args.output)
    out_path = metagenome_dir / args.output

    min_len = getattr(args, "min_contig_length", None)
    max_len = getattr(args, "max_contig_length", None)
    use_variable = min_len is not None and max_len is not None
    effective_length = (min_len + max_len) // 2 if use_variable else args.sequence_length

    reads_per_organism = args.reads_per_organism
    if args.balanced:
        stats = get_file_stats(
            download_dir,
            effective_length,
            min_length=min_len,
            max_length=max_len,
        )
        if not stats:
            raise SystemExit("No FASTA files found in download-dir.")
        for prefix, _fp, total_bases, max_reads in stats:
            print(f"  {prefix}: {total_bases} bp -> ~{max_reads} reads")
        reads_per_organism = min(max_reads for (_, _, _, max_reads) in stats)
        print(f"  Balanced: using {reads_per_organism} reads per file.")

    filter_similar = getattr(args, "filter_similar", False)
    train_test_split = getattr(args, "train_test_split", None)
    do_train_test_split = train_test_split is not None
    if filter_similar:
        plog.info("Similarity filter enabled: threshold=%.1f%%, min_coverage=%.2f, oversample_factor=%.1f",
                  getattr(args, "similarity_threshold", 90.0), getattr(args, "similarity_min_coverage", 0.8),
                  getattr(args, "oversample_factor", 2.0))
    if do_train_test_split:
        plog.info("Train-test split: %.1f%% train; similarity threshold for test removal = %.1f%%",
                  train_test_split, getattr(args, "train_test_similarity_threshold", 90.0))
    allow_ambiguous = not getattr(args, "forbid_ambiguous", False)
    sub_rate = getattr(args, "substitution_rate", 0.0)
    indel_r = getattr(args, "indel_rate", 0.0)
    if (sub_rate > 0 or indel_r > 0) and getattr(args, "seed", None) is None:
        args.seed = 42
    extra_viral = getattr(args, "extra_viral_fasta", None)
    if extra_viral is not None and not extra_viral.exists():
        raise SystemExit(f"--extra-viral-fasta not found: {extra_viral}")
    abundance_profile = _parse_abundance_profile(getattr(args, "abundance_profile", None))
    abundance_dist = getattr(args, "abundance_distribution", None)
    viral_tax_path = getattr(args, "viral_taxonomy", None)
    balance_viral = getattr(args, "balance_viral_by_taxonomy", False)
    if balance_viral and (viral_tax_path is None or not viral_tax_path.exists()):
        raise SystemExit("--balance-viral-by-taxonomy requires --viral-taxonomy PATH to an existing JSON file.")
    result = build_metagenome(
        download_dir,
        out_path,
        args.sequence_length,
        reads_per_organism,
        min_length=min_len,
        max_length=max_len,
        seed=getattr(args, "seed", None),
        cap_total_reads=getattr(args, "cap_total_reads", None),
        eve_intervals=eve_intervals,
        filter_similar=filter_similar,
        similarity_threshold=getattr(args, "similarity_threshold", 90.0),
        similarity_min_coverage=getattr(args, "similarity_min_coverage", 0.8),
        oversample_factor=getattr(args, "oversample_factor", 2.0),
        similarity_work_dir=metagenome_dir / ".simfilter_work",
        return_records=do_train_test_split,
        allow_ambiguous=allow_ambiguous,
        substitution_rate=sub_rate,
        indel_rate=indel_r,
        extra_viral_fasta=extra_viral,
        abundance_profile=abundance_profile,
        abundance_distribution=abundance_dist,
        viral_taxonomy_json=viral_tax_path,
        balance_viral_by_taxonomy=balance_viral,
    )
    if do_train_test_split:
        _count, records = result
        output_stem = Path(args.output).stem
        n_train, n_test = split_train_test_and_write(
            records,
            train_test_split,
            getattr(args, "seed", None),
            metagenome_dir,
            output_stem,
            similarity_threshold=getattr(args, "train_test_similarity_threshold", 90.0),
            similarity_min_coverage=0.8,
            work_dir=metagenome_dir / ".train_test_sim_work",
            blast_batch_size=getattr(args, "train_test_blast_batch_size", 2000),
            blast_num_threads=getattr(args, "train_test_blast_threads", 4),
        )
        train_path = metagenome_dir / f"{output_stem}_train.fasta"
        test_path = metagenome_dir / f"{output_stem}_test.fasta"
        plog.info("Train-test split: train=%d -> %s, test=%d -> %s", n_train, train_path, n_test, test_path)
        print(f"Train-test split ({train_test_split}% train): wrote {n_train} to {train_path}, {n_test} to {test_path}")
        out_path = train_path
    else:
        count = result
        plog.info("Chunk step: wrote %d sequences to %s", count, out_path)
        print(f"Wrote {count} sequences to {out_path}")

    if args.run_seeker:
        plog.info("Step 4: Seeker (phage prediction) -> %s", seeker_dir)
        seeker_dir.mkdir(parents=True, exist_ok=True)
        run_seeker(
            out_path,
            seeker_dir,
            threshold=args.seeker_threshold,
            conda_env=args.seeker_conda_env,
            min_length=SEEKER_MIN_LENGTH,
        )
    else:
        plog.info("Step 4: Seeker (phage prediction) — skipped (use --run-seeker to enable)")
    plog.info("Pipeline finished. Log file: %s", log_file.resolve())


def _run_seeker(args) -> None:
    try:
        run_seeker(
            args.input,
            args.output_dir,
            threshold=args.threshold,
            conda_env=args.conda_env or None,
            predictions_tsv=args.predictions_tsv,
            phage_fasta=args.phage_fasta,
            min_length=args.min_length,
        )
    except (FileNotFoundError, RuntimeError) as e:
        raise SystemExit(e) from e


def _run_temporal_split(args) -> None:
    p = args.accessions_file.resolve()
    if not p.exists():
        raise SystemExit(f"--accessions-file not found: {p}")
    base = p.parent
    stem = p.stem
    out_train = args.output_train or (base / f"train_{stem}.json")
    out_test = args.output_test or (base / f"test_{stem}.json")
    try:
        run_temporal_split(
            p,
            args.split_date,
            out_train,
            out_test,
            batch_size=args.batch_size,
        )
    except ValueError as e:
        raise SystemExit(e) from e


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Metagenome simulator (download, snapshot, migrate-snapshot, chunk, pipeline, blastn-filter, viral-taxonomy, seeker, temporal-split, temporal-split-info)",
    )
    subparsers = parser.add_subparsers(dest="command")
    subparsers.required = True

    _add_download_subparser(subparsers)
    _add_snapshot_subparser(subparsers)
    _add_migrate_snapshot_subparser(subparsers)
    _add_chunk_subparser(subparsers)
    _add_pipeline_subparser(subparsers)
    _add_blastn_filter_subparser(subparsers)
    _add_seeker_subparser(subparsers)
    _add_temporal_split_subparser(subparsers)
    _add_temporal_split_info_subparser(subparsers)
    _add_viral_taxonomy_subparser(subparsers)

    args = parser.parse_args()
    args.func(args)
