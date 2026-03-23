from __future__ import annotations

import subprocess
import sys
import re
from pathlib import Path

from Bio import SeqIO


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
    ]
    for sub in subcommands:
        code, out, err = run_cmd([sub, "--help"], timeout_s=15)
        assert code == 0, f"{sub} --help failed: code={code}\nstdout={out}\nstderr={err}"
        assert ("usage:" in (out.lower() + err.lower())) or (out.strip() != ""), f"{sub} --help produced no output"

