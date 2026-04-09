from __future__ import annotations

import json
from pathlib import Path

from metagenome_generator.genome_pool import (
    POOL_MANIFEST_NAME,
    materialize_from_pool,
    prepare_pool,
)


def test_materialize_from_pool_symlinks(tmp_path: Path) -> None:
    pool = tmp_path / "pool"
    (pool / "bacteria").mkdir(parents=True)
    (pool / "virus").mkdir(parents=True)
    (pool / "bacteria" / "A1.fasta").write_text(">A1\nATGC\n")
    (pool / "bacteria" / "A2.fasta").write_text(">A2\nGGC\n")
    (pool / "virus" / "V1.fasta").write_text(">V1\nATGCATGC\n")
    manifest = {
        "bacterial": ["A1", "A2"],
        "viral": ["V1"],
        "archaea": [],
        "plasmid": [],
    }
    (pool / POOL_MANIFEST_NAME).write_text(json.dumps(manifest))

    dest = tmp_path / "seed1"
    materialize_from_pool(
        pool,
        dest,
        sample_seed=1,
        max_bacteria=1,
        max_virus=1,
        max_archaea=None,
        max_plasmid=None,
        use_symlinks=True,
        clean_dest=True,
    )
    bac = list((dest / "bacteria").glob("*.fasta"))
    vir = list((dest / "virus").glob("*.fasta"))
    assert len(bac) == 1 and len(vir) == 1
    assert bac[0].is_symlink()
    assert bac[0].resolve() == (pool / "bacteria" / bac[0].name).resolve()
