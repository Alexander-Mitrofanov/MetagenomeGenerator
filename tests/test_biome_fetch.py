from __future__ import annotations

import json
from pathlib import Path

from metagenome_generator.biome_fetch import run_biome_dataset_pipeline, run_fetch_biome_data


def test_fetch_biome_data_dry_run_fraction(tmp_path: Path) -> None:
    manifest = run_fetch_biome_data(
        output_dir=tmp_path,
        biome="marine",
        source="all",
        level="metadata",
        fraction=0.25,  # with 8 samples -> 2 selected
        seed=7,
        dry_run=True,
    )
    assert manifest.exists()
    data = json.loads(manifest.read_text())
    selected = data["biomes"]["marine"]["selected_samples"]
    assert len(selected) == 2
    # dry-run records resolved targets but does not download files.
    for p in data["downloaded_files"]:
        assert not Path(p).exists()


def test_fetch_biome_data_download_uses_stub(monkeypatch, tmp_path: Path) -> None:
    # Replace network fetch with a local stub write.
    import metagenome_generator.biome_fetch as bf

    def _stub(url: str, target_path: Path, timeout_s: int = 45) -> None:  # noqa: ARG001
        target_path.parent.mkdir(parents=True, exist_ok=True)
        target_path.write_text(f"stub from {url}\n")

    monkeypatch.setattr(bf, "_download_url", _stub)
    manifest = run_fetch_biome_data(
        output_dir=tmp_path,
        biome="soil",
        source="sra",
        level="metadata",
        fraction=0.5,  # with 8 samples -> 4 selected
        seed=3,
        dry_run=False,
    )
    data = json.loads(manifest.read_text())
    selected = data["biomes"]["soil"]["selected_samples"]
    assert len(selected) == 4
    downloaded = data["downloaded_files"]
    assert downloaded, "Expected at least one downloaded metadata file"
    for p in downloaded:
        assert Path(p).exists()


def test_biome_dataset_pipeline_generates_metagenome(tmp_path: Path) -> None:
    # Prepare local FASTAs and a custom manifest with sample-level contig_urls + category.
    virus_fa = tmp_path / "virus_src.fasta"
    bact_fa = tmp_path / "bacteria_src.fasta"
    virus_fa.write_text(">v1\n" + ("ACGT" * 100) + "\n")
    bact_fa.write_text(">b1\n" + ("TGCA" * 100) + "\n")

    manifest_data = {
        "schema_version": 1,
        "biomes": {
            "marine": {
                "description": "test biome",
                "samples": [
                    {
                        "id": "marine_pair_01",
                        "contig_urls": [
                            {"url": virus_fa.resolve().as_uri(), "source": "zenodo", "category": "virus"},
                            {"url": bact_fa.resolve().as_uri(), "source": "zenodo", "category": "bacteria"},
                        ],
                    }
                ],
            }
        },
    }
    manifest_file = tmp_path / "manifest.json"
    manifest_file.write_text(json.dumps(manifest_data))

    selection_manifest, metagenome_out = run_biome_dataset_pipeline(
        output_dir=tmp_path / "out",
        biome="marine",
        source="zenodo",
        fraction=1.0,
        seed=9,
        manifest_file=manifest_file,
        sequence_length=100,
        reads_per_organism=2,
        metagenome_name="mg.fasta",
    )
    assert selection_manifest.exists()
    assert metagenome_out.exists()
    content = metagenome_out.read_text()
    assert ">virus_" in content
    assert ">bacteria_" in content
