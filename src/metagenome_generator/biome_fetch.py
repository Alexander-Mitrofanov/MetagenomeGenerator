#!/usr/bin/env python3
from __future__ import annotations
"""Fetch a reproducible subset of biome benchmark resources.

This helper is designed for profile-building workflows where users may want to
download only a fraction of a benchmark corpus first.
"""

import json
import math
import random
import shutil
from pathlib import Path
from typing import Any
from urllib.parse import urlparse
from urllib.request import urlopen

from .chunk_genomes import build_metagenome
from .genome_layout import ARCHAEA_DIR, BACTERIA_DIR, PLASMID_DIR, VIRUS_DIR


DEFAULT_BIOME_SOURCES: dict[str, Any] = {
    "schema_version": 1,
    "biomes": {
        "marine": {
            "description": "Antarctic seawater benchmark cohort (8 paired samples).",
            "metadata_urls": [
                {"url": "https://www.ebi.ac.uk/ena/browser/view/PRJEB71789", "source": "sra"},
                {"url": "https://doi.org/10.5281/zenodo.10886947", "source": "zenodo"},
            ],
            "samples": [{"id": f"marine_pair_{i:02d}"} for i in range(1, 9)],
        },
        "soil": {
            "description": "Agricultural soil benchmark cohort (8 paired samples).",
            "metadata_urls": [
                {"url": "https://www.ncbi.nlm.nih.gov/bioproject/PRJNA646773", "source": "sra"},
                {"url": "https://doi.org/10.5281/zenodo.10886947", "source": "zenodo"},
            ],
            "samples": [{"id": f"soil_pair_{i:02d}"} for i in range(1, 9)],
        },
        "gut": {
            "description": "Human gut benchmark cohort (8 paired samples).",
            "metadata_urls": [
                {"url": "https://www.ncbi.nlm.nih.gov/bioproject/PRJNA389927", "source": "sra"},
                {"url": "https://doi.org/10.5281/zenodo.10886947", "source": "zenodo"},
            ],
            "samples": [{"id": f"gut_pair_{i:02d}"} for i in range(1, 9)],
        },
    },
}


def _basename_from_url(url: str) -> str:
    parsed = urlparse(url)
    name = Path(parsed.path).name
    if not name:
        name = parsed.netloc.replace(".", "_")
    if "." not in name:
        name = f"{name}.html"
    return name


def _download_url(url: str, target_path: Path, timeout_s: int = 45) -> None:
    target_path.parent.mkdir(parents=True, exist_ok=True)
    with urlopen(url, timeout=timeout_s) as r:
        data = r.read()
    target_path.write_bytes(data)


def _load_sources(manifest_file: Path | None) -> dict[str, Any]:
    if manifest_file is None:
        return DEFAULT_BIOME_SOURCES
    with manifest_file.open() as f:
        data = json.load(f)
    if "biomes" not in data or not isinstance(data["biomes"], dict):
        raise ValueError("Manifest must contain a top-level 'biomes' object.")
    return data


def _sample_subset(sample_ids: list[str], fraction: float, max_samples: int | None, seed: int) -> list[str]:
    if not 0 < fraction <= 1:
        raise ValueError("--fraction must be in (0, 1].")
    if max_samples is not None and max_samples < 1:
        raise ValueError("--max-samples must be >= 1 when provided.")
    n_total = len(sample_ids)
    if n_total == 0:
        return []
    n = max(1, int(math.floor(n_total * fraction)))
    if max_samples is not None:
        n = min(n, max_samples)
    n = min(n, n_total)
    rng = random.Random(seed)
    chosen = sorted(rng.sample(sample_ids, n))
    return chosen


def run_fetch_biome_data(
    *,
    output_dir: Path,
    biome: str = "all",
    source: str = "all",
    level: str = "metadata",
    fraction: float = 1.0,
    max_samples: int | None = None,
    seed: int = 42,
    manifest_file: Path | None = None,
    dry_run: bool = False,
) -> Path:
    sources = _load_sources(manifest_file)
    biomes: dict[str, Any] = sources["biomes"]
    if biome == "all":
        selected_biomes = sorted(biomes.keys())
    else:
        if biome not in biomes:
            raise ValueError(f"Unknown biome: {biome}")
        selected_biomes = [biome]

    output_dir.mkdir(parents=True, exist_ok=True)
    downloaded: list[str] = []
    downloaded_records: list[dict[str, str]] = []
    selection: dict[str, Any] = {
        "schema_version": 1,
        "selection": {
            "biome": biome,
            "source": source,
            "level": level,
            "fraction": fraction,
            "max_samples": max_samples,
            "seed": seed,
            "dry_run": dry_run,
        },
        "biomes": {},
        "downloaded_files": downloaded,
        "downloaded_records": downloaded_records,
    }

    url_key = {
        "metadata": "metadata_urls",
        "contigs": "contig_urls",
        "reads": "read_urls",
    }[level]

    for b in selected_biomes:
        bdata = biomes[b]
        sample_ids = [str(s.get("id")) for s in bdata.get("samples", []) if s.get("id")]
        chosen = _sample_subset(sample_ids, fraction, max_samples, seed)
        selection["biomes"][b] = {
            "description": bdata.get("description", ""),
            "total_samples": len(sample_ids),
            "selected_samples": chosen,
        }

        urls: list[dict[str, str]] = []
        # Use sample-level urls if present; otherwise for metadata use biome-level urls.
        for s in bdata.get("samples", []):
            sid = str(s.get("id", ""))
            if sid in chosen:
                for e in s.get(url_key, []):
                    entry = dict(e)
                    entry["sample_id"] = sid
                    urls.append(entry)
        if not urls and level == "metadata":
            for e in bdata.get("metadata_urls", []):
                entry = dict(e)
                entry["sample_id"] = ""
                urls.append(entry)

        # Deduplicate by URL + source + category to preserve category-specific files.
        seen: set[tuple[str, str, str]] = set()
        unique_urls: list[dict[str, str]] = []
        for entry in urls:
            u = str(entry.get("url", "")).strip()
            src = str(entry.get("source", "")).strip() or "unknown"
            category = str(entry.get("category", "")).strip().lower()
            if not u:
                continue
            key = (u, src, category)
            if key in seen:
                continue
            seen.add(key)
            unique_urls.append({
                "url": u,
                "source": src,
                "category": category,
                "sample_id": str(entry.get("sample_id", "")),
            })

        for entry in unique_urls:
            if source != "all" and entry["source"] != source:
                continue
            fname = _basename_from_url(entry["url"])
            target = output_dir / b / level / entry["source"] / fname
            rec = {
                "biome": b,
                "level": level,
                "sample_id": entry.get("sample_id", ""),
                "source": entry["source"],
                "category": entry.get("category", ""),
                "url": entry["url"],
                "path": str(target),
            }
            if dry_run:
                downloaded.append(str(target))
                downloaded_records.append(rec)
                continue
            try:
                _download_url(entry["url"], target)
                downloaded.append(str(target))
                downloaded_records.append(rec)
            except Exception as e:  # pragma: no cover - network variability
                # Keep manifest informative even if some URLs fail.
                err_file = target.with_suffix(target.suffix + ".error.txt")
                err_file.parent.mkdir(parents=True, exist_ok=True)
                err_file.write_text(f"Failed URL: {entry['url']}\nError: {e}\n")

    manifest_out = output_dir / "selection_manifest.json"
    with manifest_out.open("w") as f:
        json.dump(selection, f, indent=2)
    return manifest_out


def run_biome_dataset_pipeline(
    *,
    output_dir: Path,
    biome: str = "all",
    source: str = "all",
    fraction: float = 1.0,
    max_samples: int | None = None,
    seed: int = 42,
    manifest_file: Path | None = None,
    sequence_length: int = 250,
    reads_per_organism: int = 100,
    metagenome_name: str = "metagenome_from_biome.fasta",
) -> tuple[Path, Path]:
    """Fetch a sampled subset of biome contig FASTAs and build a metagenome FASTA.

    Requires sample-level `contig_urls` entries in the manifest with a `category`
    field (virus, bacteria, archaea, plasmid). Default built-in manifest provides
    metadata links only, so pass --manifest-file for dataset generation.
    """
    selection_manifest = run_fetch_biome_data(
        output_dir=output_dir / "fetched",
        biome=biome,
        source=source,
        level="contigs",
        fraction=fraction,
        max_samples=max_samples,
        seed=seed,
        manifest_file=manifest_file,
        dry_run=False,
    )
    with selection_manifest.open() as f:
        selection = json.load(f)

    records: list[dict[str, str]] = selection.get("downloaded_records", [])
    if not records:
        raise ValueError(
            "No contig files were fetched. Provide a manifest with sample-level "
            "'contig_urls' entries."
        )

    genome_dir = output_dir / "genome_dir"
    by_cat = {
        "bacteria": genome_dir / BACTERIA_DIR,
        "virus": genome_dir / VIRUS_DIR,
        "archaea": genome_dir / ARCHAEA_DIR,
        "plasmid": genome_dir / PLASMID_DIR,
    }
    for d in by_cat.values():
        d.mkdir(parents=True, exist_ok=True)

    copied = 0
    for rec in records:
        cat = str(rec.get("category", "")).lower().strip()
        if cat not in by_cat:
            continue
        src = Path(rec["path"])
        if not src.exists():
            continue
        dst_name = src.name if src.suffix.lower() == ".fasta" else f"{src.stem}.fasta"
        dst = by_cat[cat] / dst_name
        if dst.resolve() != src.resolve():
            shutil.copy2(src, dst)
        copied += 1

    if copied == 0:
        raise ValueError(
            "Fetched contig files did not include categorized entries. "
            "Add `category` fields in manifest contig_urls (virus/bacteria/archaea/plasmid)."
        )
    nonviral_has = any(any(True for _ in d.glob("*.fasta")) for d in [by_cat["bacteria"], by_cat["archaea"], by_cat["plasmid"]])
    if not nonviral_has:
        raise ValueError("Need at least one non-viral FASTA (bacteria/archaea/plasmid) in fetched data.")
    viral_has = any(True for _ in by_cat["virus"].glob("*.fasta"))
    if not viral_has:
        raise ValueError("Need at least one viral FASTA in fetched data.")

    metagenome_out = output_dir / metagenome_name
    n = build_metagenome(
        genome_dir,
        metagenome_out,
        sequence_length=sequence_length,
        reads_per_organism=reads_per_organism,
        seed=seed,
    )
    if n <= 0:
        raise ValueError("Generated metagenome has zero reads; check input contig lengths and sequence_length.")
    return selection_manifest, metagenome_out

