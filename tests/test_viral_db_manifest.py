from __future__ import annotations

import importlib
import json
from pathlib import Path

import pytest

bf = importlib.import_module("metagenome_generator.blastn_filter")


def _make_fake_db(prefix: Path) -> None:
    prefix.parent.mkdir(parents=True, exist_ok=True)
    for ext, content in (
        ("nhr", b"header"),
        ("nin", b"index"),
        ("nsq", b"sequence"),
    ):
        (prefix.parent / f"{prefix.name}.{ext}").write_bytes(content)


def test_manifest_write_and_verify_roundtrip(tmp_path: Path):
    db_prefix = tmp_path / "blastn_db" / "viral_db"
    _make_fake_db(db_prefix)
    manifest = tmp_path / "viral_db_manifest.json"

    bf.write_viral_db_manifest(
        db_prefix,
        manifest,
        snapshot_timestamp="2026-03-10T00:00:00Z",
        snapshot_date="2026-03-10",
        source_snapshot="snapshots/accession_snapshot_2026-03-10.json",
        viral_accession_count=123,
    )

    data = json.loads(manifest.read_text())
    assert data["snapshot_date"] == "2026-03-10"
    assert data["viral_accession_count"] == 123
    assert data["aggregate_sha256"]

    # Should validate successfully with matching manifest.
    bf.verify_viral_db(db_prefix, manifest_path=manifest)


def test_verify_detects_checksum_mismatch(tmp_path: Path):
    db_prefix = tmp_path / "blastn_db" / "viral_db"
    _make_fake_db(db_prefix)
    manifest = tmp_path / "viral_db_manifest.json"
    bf.write_viral_db_manifest(db_prefix, manifest)

    # Corrupt one file after manifest generation.
    (db_prefix.parent / f"{db_prefix.name}.nsq").write_bytes(b"tampered")

    with pytest.raises(ValueError):
        bf.verify_viral_db(db_prefix, manifest_path=manifest)


def test_run_blastn_from_dirs_reuses_eve_cache(tmp_path: Path, monkeypatch):
    genome_dir = tmp_path / "genomes"
    (genome_dir / "bacteria").mkdir(parents=True, exist_ok=True)
    (genome_dir / "bacteria" / "B1.fasta").write_text(">B1\n" + ("ACGT" * 50) + "\n")
    out_dir = tmp_path / "out"

    db_prefix = tmp_path / "blastn_db" / "viral_db"
    _make_fake_db(db_prefix)

    calls = {"n": 0}

    def _fake_run_blastn(query_fasta: Path, _db_prefix: Path, out_tsv: Path, **_kwargs):
        calls["n"] += 1
        out_tsv.parent.mkdir(parents=True, exist_ok=True)
        # qseqid must match the FASTA record id ("B1")
        out_tsv.write_text("B1\tV1\t99.0\t100\t0\t0\t1\t100\t1\t100\t1e-20\t200\n")

    monkeypatch.setattr(bf, "run_blastn", _fake_run_blastn)

    first = bf.run_blastn_from_dirs(
        genome_dir,
        out_dir,
        viral_db_prefix=db_prefix,
        perc_identity=70.0,
        evalue=1e-5,
        num_threads=2,
        task="dc-megablast",
        reuse_cache=True,
    )
    assert first
    assert calls["n"] == 1
    assert (out_dir / "eve_intervals.json").exists()
    assert (out_dir / "eve_cache_meta.json").exists()

    # If cache works, this replacement should never be called.
    def _should_not_run(*_args, **_kwargs):
        raise AssertionError("run_blastn should not execute on cache hit")

    monkeypatch.setattr(bf, "run_blastn", _should_not_run)
    second = bf.run_blastn_from_dirs(
        genome_dir,
        out_dir,
        viral_db_prefix=db_prefix,
        perc_identity=70.0,
        evalue=1e-5,
        num_threads=2,
        task="dc-megablast",
        reuse_cache=True,
    )
    assert second == first

