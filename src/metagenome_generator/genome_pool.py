#!/usr/bin/env python3
from __future__ import annotations

"""Shared genome pool: download a superset once, then materialize seed-specific subsets via symlinks.

Use when multiple runs differ only by ``sample_seed`` but should reuse the same underlying FASTAs
and incremental EVE/BLAST caches (see ``blastn_filter.run_blastn_from_dirs`` + ``--eve-query-store``).
"""

import json
import logging
import random
import shutil
from pathlib import Path

from .download_genomes import (
    download_accession_lists,
    get_accession_lists_from_data,
    load_accessions,
)
from .genome_layout import ARCHAEA_DIR, BACTERIA_DIR, PLASMID_DIR, VIRUS_DIR

logger = logging.getLogger(__name__)

POOL_MANIFEST_NAME = "pool_manifest.json"


def _sample_list(lst: list[str], max_n: int | None, seed: int) -> list[str]:
    if max_n is None or len(lst) <= max_n:
        return list(lst)
    rng = random.Random(seed)
    return rng.sample(lst, max_n)


def prepare_pool(
    accessions_file: Path,
    pool_dir: Path,
    *,
    max_bacteria: int | None,
    max_virus: int | None,
    max_archaea: int | None,
    max_plasmid: int | None,
    pool_seed: int,
) -> None:
    """Sample up to max_* IDs per category (reproducible with pool_seed) and download into pool_dir."""
    data = load_accessions(accessions_file)
    bacterial_ids, viral_ids, archaea_ids, plasmid_ids = get_accession_lists_from_data(data)
    ts = data.get("timestamp", "unknown")
    print(f"Using accessions from {accessions_file} (timestamp: {ts})")
    print(
        f"  Pool sampling (pool_seed={pool_seed}): max bacterial={max_bacteria}, viral={max_virus}, "
        f"archaea={max_archaea}, plasmid={max_plasmid}"
    )
    seed = pool_seed
    bacterial_ids = _sample_list(bacterial_ids, max_bacteria, seed)
    viral_ids = _sample_list(viral_ids, max_virus, seed)
    archaea_ids = _sample_list(archaea_ids, max_archaea, seed)
    plasmid_ids = _sample_list(plasmid_ids, max_plasmid, seed)
    print(
        f"  Pool lists: bacterial={len(bacterial_ids)}, viral={len(viral_ids)}, "
        f"archaea={len(archaea_ids)}, plasmid={len(plasmid_ids)}"
    )
    pool_dir.mkdir(parents=True, exist_ok=True)
    download_accession_lists(pool_dir, bacterial_ids, viral_ids, archaea_ids, plasmid_ids)
    manifest = {
        "schema": 1,
        "accessions_file": str(accessions_file.resolve()),
        "pool_seed": int(pool_seed),
        "max_bacteria": max_bacteria,
        "max_virus": max_virus,
        "max_archaea": max_archaea,
        "max_plasmid": max_plasmid,
        "bacterial": bacterial_ids,
        "viral": viral_ids,
        "archaea": archaea_ids,
        "plasmid": plasmid_ids,
    }
    mp = pool_dir / POOL_MANIFEST_NAME
    with mp.open("w") as f:
        json.dump(manifest, f, indent=2)
    logger.info("Wrote pool manifest %s", mp)
    print(f"Pool manifest written: {mp}")


def _subset_from_pool(pool_list: list[str], n: int | None, sample_seed: int) -> list[str]:
    if n is None:
        return list(pool_list)
    if len(pool_list) <= n:
        return list(pool_list)
    rng = random.Random(sample_seed)
    return rng.sample(pool_list, n)


def _link_or_copy(src: Path, dest: Path, *, use_symlinks: bool) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() or dest.is_symlink():
        dest.unlink()
    if use_symlinks:
        dest.symlink_to(src.resolve())
    else:
        shutil.copy2(src, dest)


def materialize_from_pool(
    pool_dir: Path,
    dest_dir: Path,
    *,
    sample_seed: int,
    max_bacteria: int | None,
    max_virus: int | None,
    max_archaea: int | None,
    max_plasmid: int | None,
    use_symlinks: bool = True,
    clean_dest: bool = True,
) -> None:
    """Build bacteria/, virus/, … under dest_dir containing only the subset for this sample_seed."""
    mp = pool_dir / POOL_MANIFEST_NAME
    if not mp.is_file():
        raise FileNotFoundError(f"Pool manifest not found: {mp} (run genome-pool prepare first)")
    with mp.open() as f:
        manifest = json.load(f)
    pool_b = manifest.get("bacterial") or []
    pool_v = manifest.get("viral") or []
    pool_a = manifest.get("archaea") or []
    pool_p = manifest.get("plasmid") or []

    sub_b = _subset_from_pool(pool_b, max_bacteria, sample_seed)
    sub_v = _subset_from_pool(pool_v, max_virus, sample_seed)
    sub_a = _subset_from_pool(pool_a, max_archaea, sample_seed)
    sub_p = _subset_from_pool(pool_p, max_plasmid, sample_seed)

    print(
        f"Materialize (sample_seed={sample_seed}): bacterial={len(sub_b)}, viral={len(sub_v)}, "
        f"archaea={len(sub_a)}, plasmid={len(sub_p)} (symlinks={use_symlinks})"
    )

    if clean_dest and dest_dir.exists():
        shutil.rmtree(dest_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)

    def _materialize_category(ids: list[str], sub: str, label: str) -> None:
        if not ids:
            return
        src_d = pool_dir / sub
        dst_d = dest_dir / sub
        dst_d.mkdir(parents=True, exist_ok=True)
        for acc in ids:
            src = src_d / f"{acc}.fasta"
            if not src.is_file():
                raise FileNotFoundError(f"Pool missing {label} FASTA for {acc}: {src}")
            _link_or_copy(src, dst_d / f"{acc}.fasta", use_symlinks=use_symlinks)

    _materialize_category(sub_b, BACTERIA_DIR, "bacteria")
    _materialize_category(sub_v, VIRUS_DIR, "virus")
    _materialize_category(sub_a, ARCHAEA_DIR, "archaea")
    _materialize_category(sub_p, PLASMID_DIR, "plasmid")

    run_manifest = {
        "schema": 1,
        "pool_dir": str(pool_dir.resolve()),
        "sample_seed": int(sample_seed),
        "max_bacteria": max_bacteria,
        "max_virus": max_virus,
        "max_archaea": max_archaea,
        "max_plasmid": max_plasmid,
        "bacterial": sub_b,
        "viral": sub_v,
        "archaea": sub_a,
        "plasmid": sub_p,
        "use_symlinks": use_symlinks,
    }
    with (dest_dir / "materialized_manifest.json").open("w") as f:
        json.dump(run_manifest, f, indent=2)
    print(f"Materialized genomes -> {dest_dir.resolve()}")
