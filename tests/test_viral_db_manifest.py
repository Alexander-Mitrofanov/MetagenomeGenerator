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

