from __future__ import annotations

import json
from pathlib import Path

import pytest

from metagenome_generator import benchmark_recipe as br


def test_count_similar_queries_in_blast_tsv_threshold_and_coverage(tmp_path: Path):
    tsv = tmp_path / "blast.tsv"
    # outfmt6: qseqid sseqid pident length
    tsv.write_text(
        "\n".join(
            [
                "q1\ts1\t95.0\t900",  # pass (>=90, 900 >= 0.8*1000)
                "q2\ts1\t80.0\t900",  # fail identity
                "q3\ts1\t95.0\t500",  # fail coverage (500 < 0.8*1000=800)
                "q1\ts1\t92.0\t800",  # redundant hit; still counts once
            ]
        )
        + "\n"
    )
    qlens = {"q1": 1000, "q2": 1000, "q3": 1000}

    penalty = br._count_similar_queries_in_blast_tsv(
        tsv,
        qlens,
        perc_identity=90.0,
        min_coverage=0.8,
    )
    assert penalty == 1


def test_sample_ids_with_optional_overlap_prefers_exclusion():
    ids = [f"id{i}" for i in range(10)]
    rng = br.random.Random(123)
    exclude = {"id0", "id1", "id2", "id3", "id4", "id5", "id6", "id7"}

    # Only id8/id9 are available without overlap; request 2 -> should return both.
    picked = br._sample_ids_with_optional_overlap(rng, ids, 2, exclude=exclude)
    assert set(picked).issubset({"id8", "id9"})
    assert len(picked) == 2

    # Request 3: available without overlap=2, so overlap becomes unavoidable; length stays 3.
    rng2 = br.random.Random(124)
    picked2 = br._sample_ids_with_optional_overlap(rng2, ids, 3, exclude=exclude)
    assert len(picked2) == 3
    assert set(picked2).issuperset({"id8", "id9"})


def test_ensure_cache_has_accessions_downloads_only_missing(monkeypatch, tmp_path: Path):
    cache_root = tmp_path / "cache"

    # Create dummy cached files for some accessions so they should be treated as present.
    (cache_root / "bacteria").mkdir(parents=True, exist_ok=True)
    (cache_root / "virus").mkdir(parents=True, exist_ok=True)

    # Minimal valid FASTA is needed because downstream code might read; we only check file existence.
    (cache_root / "bacteria" / "B1.fasta").write_text(">B1\nACGT\n")
    (cache_root / "virus" / "V1.fasta").write_text(">V1\nACGT\n")

    seen = {}

    def fake_download_genomes(num_bacteria, num_virus, output_dir, *, accessions_file=None, **kwargs):
        assert output_dir == cache_root
        assert accessions_file is not None
        data = json.loads(Path(accessions_file).read_text())
        seen["data"] = data
        return None

    import importlib
    dg = importlib.import_module("metagenome_generator.download_genomes")

    monkeypatch.setattr(dg, "download_genomes", fake_download_genomes)

    # B1 and V1 exist; B2 and V2 should be missing and thus requested to download.
    br._ensure_cache_has_accessions(
        cache_root,
        bacterial=["B1", "B2"],
        viral=["V1", "V2"],
        archaea=[],
        plasmid=[],
        sample_seed=42,
    )

    assert "data" in seen
    d = seen["data"]
    assert set(d["bacterial"]) == {"B2"}
    assert set(d["viral"]) == {"V2"}
    assert d.get("archaea", []) == []
    assert d.get("plasmid", []) == []

