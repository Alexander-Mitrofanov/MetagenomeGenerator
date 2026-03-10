#!/usr/bin/env python3
from __future__ import annotations
"""Unified CLI for the metagenome simulator.

Subcommands: download, snapshot, chunk, pipeline, blastn-filter, seeker.
"""

import argparse
import logging
from pathlib import Path

from .accession_snapshot import get_default_snapshot_path, run_snapshot
from .download_genomes import download_genomes
from .chunk_genomes import build_metagenome, get_file_stats, split_train_test_and_write
from .seeker_wrapper import SEEKER_MIN_LENGTH, run_seeker
from .blastn_filter import load_eve_intervals, run_blastn_from_dirs
from .genome_layout import validate_genome_dir

# Organized output layout (pipeline): one root dir with step-based subdirs for easy navigation.
OUTPUT_DIR_DOWNLOADED = "downloaded"
OUTPUT_DIR_BLASTN = "blastn"
OUTPUT_DIR_METAGENOME = "metagenome"
OUTPUT_DIR_SEEKER = "seeker"
OUTPUT_DIR_LOGS = "logs"


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
        "--no-db-info",
        action="store_true",
        help="Do not fetch NCBI nucleotide db metadata (einfo).",
    )
    p.set_defaults(func=_run_snapshot)


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
        "--forbid-ambiguous",
        action="store_true",
        help="Discard chunks containing ambiguous bases (non-ACGT, e.g. N). By default, such chunks are kept.",
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


def _run_download(args) -> None:
    download_genomes(
        args.num_organisms,
        args.output_dir,
        num_archaea=getattr(args, "num_archaea", 0),
        num_plasmid=getattr(args, "num_plasmid", 0),
        accessions_file=getattr(args, "accessions_file", None),
        save_accessions_to=getattr(args, "save_accessions", None),
    )


def _run_snapshot(args) -> None:
    output = getattr(args, "output", None) or get_default_snapshot_path()
    run_snapshot(output, fetch_db_info=not getattr(args, "no_db_info", False))


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
    )
    print(f"Wrote {count} sequences to {out_path}")


def _run_blastn_filter(args) -> None:
    ok, err = validate_genome_dir(args.genome_dir)
    if not ok:
        print(f"Error: {err}", file=__import__("sys").stderr)
        raise SystemExit(1)
    run_blastn_from_dirs(
        args.genome_dir,
        args.out_dir,
        evalue=args.evalue,
        perc_identity=args.perc_identity,
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
        )

    metagenome_dir.mkdir(parents=True, exist_ok=True)

    eve_intervals = None
    if getattr(args, "run_blastn_filter", False):
        plog.info("Step 2: BLASTN filter (EVE removal); evalue=%s, perc_identity=%s",
                  getattr(args, "blastn_evalue", 1e-5), getattr(args, "blastn_perc_identity", 70.0))
        blastn_out = getattr(args, "blastn_out_dir", None) or blastn_dir
        blastn_out.mkdir(parents=True, exist_ok=True)
        run_blastn_from_dirs(
            download_dir,
            blastn_out,
            evalue=getattr(args, "blastn_evalue", 1e-5),
            perc_identity=getattr(args, "blastn_perc_identity", 70.0),
        )
        eve_json = blastn_out / "eve_intervals.json"
        eve_intervals = load_eve_intervals(eve_json)
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


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Metagenome simulator (download, snapshot, chunk, pipeline, blastn-filter, seeker)",
    )
    subparsers = parser.add_subparsers(dest="command")
    subparsers.required = True

    _add_download_subparser(subparsers)
    _add_snapshot_subparser(subparsers)
    _add_chunk_subparser(subparsers)
    _add_pipeline_subparser(subparsers)
    _add_blastn_filter_subparser(subparsers)
    _add_seeker_subparser(subparsers)

    args = parser.parse_args()
    args.func(args)
