from __future__ import annotations

import importlib
from pathlib import Path

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

dg = importlib.import_module("metagenome_generator.download_genomes")


def test_download_category_reuses_existing_fastas(monkeypatch, tmp_path: Path):
    out_dir = tmp_path / "downloads" / "bacteria"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Existing cached file should be reused.
    cached = out_dir / "ACC_A.fasta"
    cached.write_text(">ACC_A\nACGT\n")

    calls: list[list[str]] = []

    def fake_fetch(ids: list[str], max_retries: int = 3):
        calls.append(list(ids))
        recs = []
        for acc in ids:
            recs.append(SeqRecord(Seq("ACGT"), id=acc, description=""))
        return recs

    monkeypatch.setattr(dg, "fetch_sequences", fake_fetch)
    monkeypatch.setattr(dg.time, "sleep", lambda *_args, **_kwargs: None)

    dg._download_category_batched(
        ["ACC_A", "ACC_B"],
        out_dir,
        "bacteria",
        "bacteria",
    )

    # Only missing accession should be fetched.
    assert calls == [["ACC_B"]]
    assert (out_dir / "ACC_A.fasta").exists()
    assert (out_dir / "ACC_B.fasta").exists()

