#!/usr/bin/env python3
from __future__ import annotations
"""Genome folder layout: bacteria/, virus/, archaea/, plasmid/ under the download root.

Download step writes into these subdirs; chunk and BLASTN discover files from them
or from a flat directory (backward compatible).
"""

from pathlib import Path

# Category subdir names under the download root (nouns: bacteria, virus, archaea, plasmid)
BACTERIA_DIR = "bacteria"
VIRUS_DIR = "virus"
ARCHAEA_DIR = "archaea"
PLASMID_DIR = "plasmid"

# Legacy flat-file prefix (for backward compatibility when reading old layouts)
BACTERIA_PREFIX = "bacterial_"
VIRUS_PREFIX = "viral_"
ARCHAEA_PREFIX = "archaea_"
PLASMID_PREFIX = "plasmid_"


def _category_dirs() -> list[tuple[str, str]]:
    """(subdir_name, file_prefix) for each category. Prefix used only for flat fallback."""
    return [
        (BACTERIA_DIR, BACTERIA_PREFIX),
        (VIRUS_DIR, VIRUS_PREFIX),
        (ARCHAEA_DIR, ARCHAEA_PREFIX),
        (PLASMID_DIR, PLASMID_PREFIX),
    ]


def iter_genome_fastas(root: Path) -> list[tuple[str, Path]]:
    """List all genome FASTAs under root: (stem, path). Prefer subdirs bacteria/, virus/, archaea/, plasmid/; fallback to flat *.fasta."""
    out: list[tuple[str, Path]] = []
    for subdir_name, _prefix in _category_dirs():
        d = root / subdir_name
        if d.is_dir():
            for fp in sorted(d.glob("*.fasta")):
                out.append((fp.stem, fp))
    if out:
        return sorted(out, key=lambda x: (x[0], str(x[1])))
    # Backward compatibility: flat layout
    for fp in sorted(root.glob("*.fasta")):
        out.append((fp.stem, fp))
    return out


def get_viral_fasta_paths(root: Path) -> list[Path]:
    """Paths to all viral genome FASTAs (virus/*.fasta or flat viral_*.fasta for legacy)."""
    virus_dir = root / VIRUS_DIR
    if virus_dir.is_dir():
        return sorted(virus_dir.glob("*.fasta"))
    return sorted(root.glob(f"{VIRUS_PREFIX}*.fasta"))


def get_nonviral_fasta_paths(root: Path) -> list[Path]:
    """Paths to all non-viral genome FASTAs (bacteria/, archaea/, plasmid/ or flat patterns)."""
    paths: list[Path] = []
    for subdir_name, prefix in _category_dirs():
        if subdir_name == VIRUS_DIR:
            continue
        d = root / subdir_name
        if d.is_dir():
            paths.extend(sorted(d.glob("*.fasta")))
        else:
            paths.extend(sorted(root.glob(f"{prefix}*.fasta")))
    return paths


def validate_genome_dir(root: Path) -> tuple[bool, str | None]:
    """Check that genome root has the expected layout: virus/ must be non-empty; at least one of bacteria/, archaea/, plasmid/ must be non-empty.

    Expects 4 folders (bacteria, virus, archaea, plasmid). Returns (True, None) if valid,
    (False, error_message) otherwise. Caller should print and log the message and exit.
    """
    root = Path(root)
    if not root.is_dir():
        return False, f"Genome directory does not exist or is not a directory: {root}"

    viral = get_viral_fasta_paths(root)
    if not viral:
        return (
            False,
            "Viral folder (virus/) must not be empty. Place viral genome FASTAs (from download or your own in-house files) in virus/.",
        )

    nonviral = get_nonviral_fasta_paths(root)
    if not nonviral:
        return (
            False,
            "At least one of bacteria/, archaea/, or plasmid/ must be non-empty. Place genome FASTAs (from download or your own in-house files) in the corresponding folders.",
        )

    return True, None
