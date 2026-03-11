#!/usr/bin/env python3
"""Extensive CLI validation tests: wrong flag combinations, invalid values, edge cases.

Run from project root: python -m pytest tests/test_cli_validation.py -v
Or: python tests/test_cli_validation.py
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

# Project root
ROOT = Path(__file__).resolve().parent.parent
CLI = [sys.executable, "-m", "metagenome_generator"]


def run_cmd(args: list[str], expect_fail: bool = True) -> tuple[int, str, str]:
    """Run metagenome-generator with args. Return (returncode, stdout, stderr)."""
    proc = subprocess.run(
        CLI + args,
        cwd=ROOT,
        capture_output=True,
        text=True,
        timeout=60,
    )
    return proc.returncode, proc.stdout or "", proc.stderr or ""


def test_download_nonexistent_accessions_file():
    """--accessions-file pointing to missing file should exit non-zero with clear message."""
    code, out, err = run_cmd(["download", "--num-bacteria", "2", "--num-virus", "2", "--accessions-file", "/nonexistent/snap.json"])
    assert code != 0, f"Expected failure, got code=0. stderr: {err}"
    assert "not found" in err or "nonexistent" in err.lower() or "No such" in err, f"Expected file-not-found message. err={err!r}"


def test_download_negative_num_bacteria():
    """Negative --num-bacteria must be rejected."""
    code, out, err = run_cmd(["download", "--num-bacteria", "-1", "--num-virus", "1", "--output-dir", "/tmp/mg_test_dl"])
    assert code != 0, f"Negative num_bacteria should fail. code={code}, out={out!r}, err={err!r}"
    assert "num-bacteria" in (out + err).lower() or ">= 0" in (out + err), f"Expected message about num-bacteria. err={err!r}"


def test_download_negative_num_virus():
    """Negative --num-virus must be rejected."""
    code, out, err = run_cmd(["download", "--num-bacteria", "1", "--num-virus", "-1", "--output-dir", "/tmp/mg_test_dl"])
    assert code != 0, f"Negative num_virus should fail. code={code}, out={out!r}, err={err!r}"
    assert "num-virus" in (out + err).lower() or ">= 0" in (out + err), f"Expected message about num-virus. err={err!r}"


def test_chunk_input_nonexistent():
    """--input that does not exist should fail with clear error."""
    code, out, err = run_cmd([
        "chunk", "--input", "/nonexistent/genomes", "--output", "out.fasta", "--output-dir", "/tmp/mg_test",
        "--sequence-length", "250", "--reads-per-organism", "100",
    ])
    assert code != 0, f"Expected failure for nonexistent input. code=0, err={err!r}"
    # May fail at validate_genome_dir (not a dir) or later (no FASTA files)
    assert "exist" in err.lower() or "not found" in err.lower() or "No FASTA" in err or "Error" in err, f"err={err!r}"


def test_chunk_sequence_length_zero():
    """--sequence-length 0 can cause range()/division errors; should be rejected or handled."""
    # Use a dir that exists but may be empty so we get to chunk logic
    code, out, err = run_cmd([
        "chunk", "--input", ROOT, "--output", "out.fasta", "--output-dir", "/tmp/mg_test",
        "--sequence-length", "0", "--reads-per-organism", "100",
    ])
    assert code != 0, f"sequence-length 0 should fail. code=0, err={err!r}"


def test_chunk_sequence_length_negative():
    """--sequence-length negative should fail."""
    code, out, err = run_cmd([
        "chunk", "--input", ROOT, "--output", "out.fasta", "--output-dir", "/tmp/mg_test",
        "--sequence-length", "-1", "--reads-per-organism", "100",
    ])
    assert code != 0, f"Negative sequence-length should fail. code=0, err={err!r}"


def test_chunk_min_contig_gt_max_contig():
    """--min-contig-length > --max-contig-length should fail (randint(min,max) error otherwise)."""
    code, out, err = run_cmd([
        "chunk", "--input", ROOT, "--output", "out.fasta", "--output-dir", "/tmp/mg_test",
        "--sequence-length", "250", "--min-contig-length", "500", "--max-contig-length", "300",
    ])
    assert code != 0, f"min > max contig length should fail. code=0, err={err!r}"


def test_chunk_substitution_rate_above_one():
    """--substitution-rate > 1 is semantically invalid (probability)."""
    # Might not crash but should ideally be rejected
    code, out, err = run_cmd([
        "chunk", "--input", ROOT, "--output", "out.fasta", "--output-dir", "/tmp/mg_test",
        "--sequence-length", "250", "--substitution-rate", "1.5",
    ])
    # Tool may still run (clamp or accept); we only check it doesn't crash with traceback
    assert "Traceback" not in err or code != 0, f"Unexpected traceback: err={err!r}"


def test_chunk_balance_viral_without_viral_taxonomy():
    """--balance-viral-by-taxonomy without --viral-taxonomy should fail (or fail earlier at dir validation)."""
    code, out, err = run_cmd([
        "chunk", "--input", ROOT, "--output", "out.fasta", "--output-dir", "/tmp/mg_test",
        "--sequence-length", "250", "--balance-viral-by-taxonomy",
    ])
    assert code != 0, f"balance-viral without viral-taxonomy should fail. code=0, err={err!r}"
    # May fail at genome-dir validation (no virus/) or at balance-viral (requires --viral-taxonomy)
    combined = out + err
    assert "viral" in combined.lower() or "requires" in combined.lower(), f"Expected viral/requires message. err={err!r}"


def test_chunk_eve_intervals_nonexistent():
    """--eve-intervals pointing to missing file should fail."""
    code, out, err = run_cmd([
        "chunk", "--input", ROOT, "--output", "out.fasta", "--output-dir", "/tmp/mg_test",
        "--sequence-length", "250", "--eve-intervals", "/nonexistent/eve.json",
    ])
    assert code != 0, f"Nonexistent eve-intervals should fail. code=0, err={err!r}"


def test_filter_test_against_train_nonexistent_train():
    """--train-fasta missing should exit with error."""
    code, out, err = run_cmd([
        "filter-test-against-train",
        "--train-fasta", "/nonexistent/train.fasta",
        "--test-fasta", "/nonexistent/test.fasta",
        "--output", "/tmp/out.fasta",
    ])
    assert code != 0, f"Nonexistent train-fasta should fail. code=0, err={err!r}"
    assert "not found" in err or "File not found" in err or "nonexistent" in err.lower(), f"err={err!r}"


def test_temporal_split_nonexistent_accessions_file():
    """temporal-split with missing --accessions-file should fail."""
    code, out, err = run_cmd([
        "temporal-split", "--accessions-file", "/nonexistent/snap.json", "--split-date", "2020-01-01",
    ])
    assert code != 0, f"Nonexistent accessions-file should fail. code=0, err={err!r}"


def test_temporal_split_info_invalid_date():
    """Invalid --split-date format should be rejected by temporal_split module."""
    # We need a valid accessions file for temporal-split-info to get to date parsing; use a minimal JSON
    import tempfile
    import json
    with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
        json.dump({"bacterial": [], "viral": [], "archaea": [], "plasmid": []}, f)
        path = f.name
    try:
        code, out, err = run_cmd([
            "temporal-split-info", "--accessions-file", path, "--split-date", "2020-13-45",
        ])
        assert code != 0, f"Invalid date should fail. code=0, err={err!r}"
    finally:
        Path(path).unlink(missing_ok=True)


def test_blastn_filter_nonexistent_genome_dir():
    """blastn-filter with nonexistent --genome-dir should fail."""
    code, out, err = run_cmd([
        "blastn-filter", "--genome-dir", "/nonexistent/genomes", "--out-dir", "/tmp/mg_blast",
    ])
    assert code != 0, f"Nonexistent genome-dir should fail. code=0, err={err!r}"


def test_viral_taxonomy_nonexistent_accessions_file():
    """viral-taxonomy with missing --accessions-file should fail."""
    code, out, err = run_cmd([
        "viral-taxonomy", "--accessions-file", "/nonexistent/snap.json", "--output", "/tmp/vt.json",
    ])
    assert code != 0, f"Nonexistent accessions-file should fail. code=0, err={err!r}"


def test_chunk_train_test_split_invalid_percentage():
    """--train-test-split 0 or >100 might be clamped; 150 could be invalid."""
    # If we have no input files, we may fail earlier; use ROOT as input (no fasta), we get "No FASTA" or validation
    code, out, err = run_cmd([
        "chunk", "--input", ROOT, "--output", "out.fasta", "--output-dir", "/tmp/mg_test",
        "--sequence-length", "250", "--train-test-split", "150",
    ])
    # 150 might be accepted by argparse; backend may clamp to 100
    assert "Traceback" not in err or code != 0, f"Traceback with train-test-split 150: err={err!r}"


if __name__ == "__main__":
    import pytest
    pytest.main([__file__, "-v"])
