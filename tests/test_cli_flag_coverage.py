from __future__ import annotations

import subprocess
import sys
import re
import argparse
import time
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# Project root
ROOT = Path(__file__).resolve().parent.parent
CLI = [sys.executable, "-m", "metagenome_generator"]


def run_cmd(args: list[str], timeout_s: int = 60) -> tuple[int, str, str]:
    proc = subprocess.run(
        CLI + args,
        cwd=ROOT,
        capture_output=True,
        text=True,
        timeout=timeout_s,
    )
    return proc.returncode, proc.stdout or "", proc.stderr or ""


def _write_fasta(path: Path, rec_id: str, seq: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f">{rec_id}\n{seq}\n")


def _make_tiny_genome_dir(tmp_path: Path) -> Path:
    # validate_genome_dir requires virus/ to be non-empty and at least one non-viral category non-empty.
    genome_dir = tmp_path / "genomes"
    seq = "ACGT" * 6 + "A"  # length = 25
    _write_fasta(genome_dir / "bacteria" / "B1.fasta", "B1", seq)
    _write_fasta(genome_dir / "virus" / "V1.fasta", "V1", seq)
    return genome_dir


def _make_large_genome_dir(tmp_path: Path) -> Path:
    genome_dir = tmp_path / "genomes_large"
    seq = ("ACGT" * 3000)  # length 12,000
    _write_fasta(genome_dir / "bacteria" / "B1.fasta", "B1", seq)
    _write_fasta(genome_dir / "virus" / "V1.fasta", "V1", seq)
    return genome_dir


def _first_record_start_and_id(fasta_path: Path) -> tuple[int, str]:
    recs = list(SeqIO.parse(str(fasta_path), "fasta"))
    assert recs, f"No records parsed from {fasta_path}"
    rec0 = recs[0]
    m = re.search(r"start=(\d+)", rec0.description)
    assert m, f"Could not find start= in description: {rec0.description!r}"
    return int(m.group(1)), rec0.id


def test_chunk_random_chunking_changes_start_offset(tmp_path: Path) -> None:
    genome_dir = _make_tiny_genome_dir(tmp_path)
    out_no = tmp_path / "out_no"
    out_yes = tmp_path / "out_yes"
    out_no.mkdir()
    out_yes.mkdir()

    base = [
        "chunk",
        "--input",
        str(genome_dir),
        "--output",
        "mg.fasta",
        "--output-dir",
        str(out_no),
        "--sequence-length",
        "10",
        "--reads-per-organism",
        "10",
        "--seed",
        "1",
    ]
    code, out, err = run_cmd(base)
    assert code == 0, f"chunk failed: code={code}\nstdout={out}\nstderr={err}"
    start0, rid0 = _first_record_start_and_id(out_no / "mg.fasta")
    assert rid0.startswith("bacteria_"), f"Expected bacteria_* read id, got {rid0}"
    assert start0 == 0

    base[5 + 1] = str(out_yes)  # output-dir
    base_yes = base[:-8] + [
        "--output-dir",
        str(out_yes),
        "--random-chunking",
    ]
    # base above is a bit brittle; build the final args explicitly for clarity.
    base_yes = [
        "chunk",
        "--input",
        str(genome_dir),
        "--output",
        "mg.fasta",
        "--output-dir",
        str(out_yes),
        "--sequence-length",
        "10",
        "--reads-per-organism",
        "10",
        "--seed",
        "1",
        "--random-chunking",
    ]
    code, out, err = run_cmd(base_yes)
    assert code == 0, f"chunk --random-chunking failed: code={code}\nstdout={out}\nstderr={err}"
    start1, rid1 = _first_record_start_and_id(out_yes / "mg.fasta")
    assert rid1.startswith("bacteria_")
    # With seed=1 and chunk_size=10, randint(0,9) => 2, so first window start should be 2.
    assert start1 == 2


def test_pipeline_random_chunking_flag(tmp_path: Path) -> None:
    genome_dir = _make_tiny_genome_dir(tmp_path)
    outdir = tmp_path / "pipeline_out"
    outdir.mkdir()

    code, out, err = run_cmd(
        [
            "pipeline",
            "--genome-dir",
            str(genome_dir),
            "--output-dir",
            str(outdir),
            "--output",
            "mg.fasta",
            "--sequence-length",
            "10",
            "--reads-per-organism",
            "10",
            "--seed",
            "1",
            "--random-chunking",
        ]
    )
    assert code == 0, f"pipeline failed: code={code}\nstdout={out}\nstderr={err}"

    start, rid = _first_record_start_and_id(outdir / "mg.fasta")
    assert rid.startswith("bacteria_")
    assert start == 2


def test_temporal_pipeline_random_chunking_flag_is_parsed(tmp_path: Path) -> None:
    code, out, err = run_cmd(
        [
            "temporal-pipeline",
            "--accessions-file",
            "/nonexistent/snap.json",
            "--split-date",
            "2020-01-01",
            "--output-dir",
            str(tmp_path / "tp_out"),
            "--random-chunking",
        ],
        timeout_s=15,
    )
    assert code != 0
    combined = (out + "\n" + err).lower()
    assert "unrecognized arguments" not in combined
    assert "not found" in combined or "accessions-file" in combined


def test_split_metagenome_train_test_smoke_fastq(tmp_path: Path) -> None:
    in_fasta = tmp_path / "mg.fasta"
    in_fasta.write_text(">r1\n" + ("A" * 30) + "\n")

    outdir = tmp_path / "split_out"
    outdir.mkdir()

    code, out, err = run_cmd(
        [
            "split-metagenome-train-test",
            "--input",
            str(in_fasta),
            "--output-dir",
            str(outdir),
            "--output-stem",
            "mg",
            "--train-test-split",
            "80",
            "--output-fastq",
        ],
        timeout_s=60,
    )
    assert code == 0, f"split-metagenome-train-test failed: code={code}\nstdout={out}\nstderr={err}"

    train_p = outdir / "mg_train.fastq"
    test_p = outdir / "mg_test.fastq"
    assert train_p.exists()
    assert test_p.exists()

    train_recs = list(SeqIO.parse(str(train_p), "fastq"))
    test_recs = list(SeqIO.parse(str(test_p), "fastq"))
    assert len(train_recs) == 1
    assert len(test_recs) == 0


def test_subcommand_help_smoke_covers_flags() -> None:
    # `--help` exercises argparse for each subcommand, ensuring every flag is registered
    # and does not cause an argparse crash.
    subcommands = [
        "chunk",
        "pipeline",
        "temporal-pipeline",
        "filter-test-against-train",
        "split-metagenome-train-test",
        "benchmark-recipe",
        "biome-metagenome",
    ]
    for sub in subcommands:
        code, out, err = run_cmd([sub, "--help"], timeout_s=15)
        assert code == 0, f"{sub} --help failed: code={code}\nstdout={out}\nstderr={err}"
        assert ("usage:" in (out.lower() + err.lower())) or (out.strip() != ""), f"{sub} --help produced no output"


def test_biome_metagenome_with_genome_dir_smoke(tmp_path: Path) -> None:
    genome_dir = _make_tiny_genome_dir(tmp_path)
    outdir = tmp_path / "biome_out"
    outdir.mkdir()
    code, out, err = run_cmd(
        [
            "biome-metagenome",
            "--biome-profile",
            "marine",
            "--genome-dir",
            str(genome_dir),
            "--output-dir",
            str(outdir),
            "--output",
            "mg.fasta",
            "--reads-per-organism",
            "2",
            "--sequence-length",
            "10",
            "--seed",
            "3",
        ]
    )
    assert code == 0, f"biome-metagenome failed: code={code}\nstdout={out}\nstderr={err}"
    out_fa = outdir / "mg.fasta"
    assert out_fa.exists()
    recs = list(SeqIO.parse(str(out_fa), "fasta"))
    assert recs, "Expected generated metagenome reads"


def test_chunk_contig_quality_profile_generates_variable_contigs(tmp_path: Path) -> None:
    genome_dir = _make_large_genome_dir(tmp_path)
    outdir = tmp_path / "profile_out"
    outdir.mkdir()
    code, out, err = run_cmd(
        [
            "chunk",
            "--input",
            str(genome_dir),
            "--output",
            "mg.fasta",
            "--output-dir",
            str(outdir),
            "--reads-per-organism",
            "6",
            "--seed",
            "2",
            "--contig-quality-profile",
            "low-quality",
        ]
    )
    assert code == 0, f"chunk contig-quality-profile failed: code={code}\nstdout={out}\nstderr={err}"
    recs = list(SeqIO.parse(str(outdir / "mg.fasta"), "fasta"))
    assert recs
    lens = [len(r.seq) for r in recs]
    assert min(lens) >= 300
    assert max(lens) <= 6000
    assert len(set(lens)) > 1


def test_split_train_test_accepts_fraction(monkeypatch, tmp_path: Path) -> None:
    from metagenome_generator import similarity_filter
    from metagenome_generator.chunk_genomes import split_train_test_and_write

    # Keep all test reads so this test checks split parsing, not BLAST behavior.
    monkeypatch.setattr(
        similarity_filter,
        "filter_candidates_against_kept",
        lambda candidates, *_args, **_kwargs: candidates,
    )

    records = [
        SeqRecord(Seq("ACGT" * 20), id=f"r{i}", description="")
        for i in range(10)
    ]
    n_train, n_test = split_train_test_and_write(
        records,
        0.8,  # should be interpreted as 80%
        seed=7,
        output_dir=tmp_path,
        output_stem="mg",
    )
    assert n_train == 8
    assert n_test == 2


def test_split_metagenome_default_train_test_split_is_80_percent() -> None:
    from metagenome_generator.cli import _add_split_metagenome_train_test_subparser

    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command")
    _add_split_metagenome_train_test_subparser(subparsers)
    args = parser.parse_args(["split-metagenome-train-test", "--input", "x.fasta"])
    assert args.train_test_split == 80.0


def test_split_train_test_written_files_match_expected_ratio(monkeypatch, tmp_path: Path) -> None:
    from metagenome_generator import similarity_filter
    from metagenome_generator.chunk_genomes import split_train_test_and_write

    monkeypatch.setattr(
        similarity_filter,
        "filter_candidates_against_kept",
        lambda candidates, *_args, **_kwargs: candidates,
    )

    records = [
        SeqRecord(Seq(("ACGT" * 30)[:100]), id=f"r{i}", description="")
        for i in range(25)
    ]
    split_train_test_and_write(
        records,
        80.0,
        seed=11,
        output_dir=tmp_path,
        output_stem="ratio",
    )
    train_recs = list(SeqIO.parse(str(tmp_path / "ratio_train.fasta"), "fasta"))
    test_recs = list(SeqIO.parse(str(tmp_path / "ratio_test.fasta"), "fasta"))
    assert len(train_recs) == 20
    assert len(test_recs) == 5


def test_split_train_test_performance_smoke(monkeypatch, tmp_path: Path) -> None:
    from metagenome_generator import similarity_filter
    from metagenome_generator.chunk_genomes import split_train_test_and_write

    monkeypatch.setattr(
        similarity_filter,
        "filter_candidates_against_kept",
        lambda candidates, *_args, **_kwargs: candidates,
    )

    records = [
        SeqRecord(Seq(("ACGT" * 50)[:200]), id=f"perf_{i}", description="")
        for i in range(2000)
    ]
    t0 = time.perf_counter()
    n_train, n_test = split_train_test_and_write(
        records,
        80.0,
        seed=3,
        output_dir=tmp_path,
        output_stem="perf",
    )
    elapsed = time.perf_counter() - t0
    assert n_train == 1600
    assert n_test == 400
    # generous bound: catches major regressions without flakiness
    assert elapsed < 5.0


@pytest.mark.parametrize(
    ("split_value", "expected_train", "expected_test"),
    [
        (0.7, 140, 60),
        (70.0, 140, 60),
        (80.0, 160, 40),
        (0.9, 180, 20),
        (90.0, 180, 20),
    ],
)
def test_split_values_and_performance_with_min_100_reads(
    monkeypatch,
    tmp_path: Path,
    split_value: float,
    expected_train: int,
    expected_test: int,
) -> None:
    from metagenome_generator import similarity_filter
    from metagenome_generator.chunk_genomes import split_train_test_and_write

    # Keep all test reads so split percentages are directly testable.
    monkeypatch.setattr(
        similarity_filter,
        "filter_candidates_against_kept",
        lambda candidates, *_args, **_kwargs: candidates,
    )

    # 200 reads >= requested minimum of 100.
    records = [
        SeqRecord(Seq(("ACGT" * 50)[:200]), id=f"s{i}", description="")
        for i in range(200)
    ]

    t0 = time.perf_counter()
    n_train, n_test = split_train_test_and_write(
        records,
        split_value,
        seed=13,
        output_dir=tmp_path,
        output_stem=f"split_{str(split_value).replace('.', '_')}",
    )
    elapsed = time.perf_counter() - t0

    assert n_train == expected_train
    assert n_test == expected_test
    # Performance guard for this small synthetic case.
    assert elapsed < 5.0

