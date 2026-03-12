#!/usr/bin/env python3
from __future__ import annotations
"""Download bacterial and viral genomes from NCBI.

Saves into category subfolders: bacteria/, virus/, archaea/, plasmid/ under output_dir.
Exposes `download_genomes` for programmatic use and a CLI when run directly.

Reproducibility: NCBI search results change over time. Use --save-accessions to write
the accession list (and timestamp) to a JSON file after searching; use --accessions-file
to skip the search and download exactly those accessions on a later run.
"""

import argparse
import json
import logging
import random
import time
from pathlib import Path

from Bio import Entrez, SeqIO

from .genome_layout import (
    ARCHAEA_DIR,
    BACTERIA_DIR,
    PLASMID_DIR,
    VIRUS_DIR,
)
from .ncbi_search import get_queries, search_genomes

logger = logging.getLogger(__name__)


# Keys used in the accessions JSON file (reproducibility)
ACCESSIONS_KEY_TIMESTAMP = "timestamp"
ACCESSIONS_KEY_BACTERIAL = "bacterial"
ACCESSIONS_KEY_VIRAL = "viral"
ACCESSIONS_KEY_ARCHAEA = "archaea"
ACCESSIONS_KEY_PLASMID = "plasmid"
# Optional: per-accession metadata (create_date, title) from NCBI esummary
ACCESSION_METADATA_KEY = "accession_metadata"

# Category keys in order (bacterial, viral, archaea, plasmid)
ACCESSIONS_CATEGORY_KEYS = (
    ACCESSIONS_KEY_BACTERIAL,
    ACCESSIONS_KEY_VIRAL,
    ACCESSIONS_KEY_ARCHAEA,
    ACCESSIONS_KEY_PLASMID,
)


def _category_value_to_id_list(value: list) -> list[str]:
    """Normalize a category value to a list of accession ID strings.

    Supports: list of strings (legacy) or list of dicts with 'accession' key (new format).
    """
    if not value:
        return []
    first = value[0]
    if isinstance(first, str):
        return list(value)
    if isinstance(first, dict) and "accession" in first:
        return [item["accession"] for item in value if isinstance(item, dict) and item.get("accession")]
    return []


def get_accession_lists_from_data(data: dict) -> tuple[list[str], list[str], list[str], list[str]]:
    """Extract (bacterial, viral, archaea, plasmid) ID lists from loaded snapshot/accessions JSON.

    Handles both legacy format (category = list of ID strings) and new format
    (category = list of {accession, create_date, title}).
    """
    return (
        _category_value_to_id_list(data.get(ACCESSIONS_KEY_BACTERIAL, [])),
        _category_value_to_id_list(data.get(ACCESSIONS_KEY_VIRAL, [])),
        _category_value_to_id_list(data.get(ACCESSIONS_KEY_ARCHAEA, [])),
        _category_value_to_id_list(data.get(ACCESSIONS_KEY_PLASMID, [])),
    )


def get_accession_metadata_from_data(data: dict) -> dict[str, dict]:
    """Extract accession -> {create_date, title} from loaded JSON.

    Uses top-level accession_metadata if present (legacy), otherwise builds from
    per-category list of objects (new format: each item has accession, create_date, title).
    """
    legacy = data.get(ACCESSION_METADATA_KEY)
    if isinstance(legacy, dict) and legacy:
        return dict(legacy)
    result: dict[str, dict] = {}
    for key in ACCESSIONS_CATEGORY_KEYS:
        value = data.get(key, [])
        if not value:
            continue
        first = value[0] if value else None
        if isinstance(first, dict) and "accession" in first:
            for item in value:
                if isinstance(item, dict) and item.get("accession"):
                    acc = item["accession"]
                    result[acc] = {
                        "create_date": item.get("create_date") or "",
                        "title": item.get("title") or "",
                    }
    return result


def load_accessions(path: Path) -> dict:
    """Load accession lists and optional metadata from a JSON file.

    Expected keys: bacterial, viral, archaea, plasmid (lists of accession IDs).
    Optional: timestamp (ISO8601 string); accession_metadata (dict of accession -> {create_date, title}).
    Returns the full JSON dict; missing list keys default to [] when used.
    """
    with path.open() as f:
        return json.load(f)


def save_accessions(
    path: Path,
    bacterial: list[str],
    viral: list[str],
    archaea: list[str],
    plasmid: list[str],
) -> None:
    """Write accession lists and current UTC timestamp to a JSON file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    now = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    data = {
        ACCESSIONS_KEY_TIMESTAMP: now,
        ACCESSIONS_KEY_BACTERIAL: bacterial,
        ACCESSIONS_KEY_VIRAL: viral,
        ACCESSIONS_KEY_ARCHAEA: archaea,
        ACCESSIONS_KEY_PLASMID: plasmid,
    }
    with path.open("w") as f:
        json.dump(data, f, indent=2)
    logger.info("Saved accessions (timestamp=%s) to %s", now, path)


# Batch size for efetch (NCBI accepts multiple IDs; reduces API calls)
EFETCH_BATCH_SIZE = 20


def fetch_sequences(ids: list[str], max_retries: int = 3) -> list:
    """Fetch FASTA sequences for given IDs from NCBI. Retries on transient failures.
    Pass multiple IDs to reduce API calls; records are returned in the same order as ids.
    """
    if not ids:
        return []
    id_str = ",".join(ids)
    last_error = None
    for attempt in range(max_retries):
        try:
            handle = Entrez.efetch(
                db="nucleotide",
                id=id_str,
                rettype="fasta",
                retmode="text",
            )
            records = list(SeqIO.parse(handle, "fasta"))
            handle.close()
            return records
        except Exception as e:
            last_error = e
            if attempt < max_retries - 1:
                time.sleep(1.0 * (attempt + 1))
            continue
    raise last_error


def _download_category_batched(
    ids: list[str],
    out_dir: Path,
    _prefix: str,
    category_label: str,
) -> None:
    """Download genomes in batches; write one FASTA per accession, named by accession (e.g. NC_000001.1.fasta). On batch failure, fall back to per-ID fetch."""
    for start in range(0, len(ids), EFETCH_BATCH_SIZE):
        batch = ids[start : start + EFETCH_BATCH_SIZE]
        try:
            time.sleep(0.4)
            records = fetch_sequences(batch)
            if len(records) != len(batch):
                raise RuntimeError(f"Batch returned {len(records)} records for {len(batch)} IDs")
            for j, rec in enumerate(records):
                acc = rec.id.split()[0]
                path = out_dir / f"{acc}.fasta"
                print(f"Fetching {category_label} {start + j + 1}/{len(ids)}: {batch[j]} -> {path.name}")
                logger.info("Download: accession=%s origin=%s path=%s", batch[j], category_label, path)
                SeqIO.write([rec], path, "fasta")
        except Exception as e:
            for k, gid in enumerate(batch):
                print(f"Fetching {category_label} {start + k + 1}/{len(ids)}: {gid} -> {gid}.fasta")
                logger.info("Download: accession=%s origin=%s", gid, category_label)
                time.sleep(0.4)
                try:
                    records = fetch_sequences([gid])
                    acc = records[0].id.split()[0]
                    path = out_dir / f"{acc}.fasta"
                    SeqIO.write(records, path, "fasta")
                except Exception as e2:
                    print(f"  Warning: failed to fetch {gid}: {e2}")
                    logger.warning("Download failed: accession=%s error=%s", gid, e2)


def download_genomes(
    num_bacteria: int,
    num_virus: int,
    output_dir: Path,
    *,
    num_archaea: int = 0,
    num_plasmid: int = 0,
    accessions_file: Path | None = None,
    save_accessions_to: Path | None = None,
    complete_only: bool = False,
    max_bacteria: int | None = None,
    max_virus: int | None = None,
    max_archaea: int | None = None,
    max_plasmid: int | None = None,
    sample_seed: int | None = None,
) -> None:
    """Download bacterial and viral genomes (and optionally archaea, plasmid) into `output_dir`.

    Uses category subfolders: bacteria/, virus/, archaea/, plasmid/. You specify how many
    bacterial and how many viral genomes separately (num_bacteria, num_virus). Negative
    samples (bacteria, archaea, plasmid) support training classifiers to distinguish viruses.

    If accessions_file is set, load accession IDs from that JSON (and optional timestamp)
    instead of searching NCBI, so the run is reproducible. Optionally limit how many to
    use per category with max_bacteria, max_virus, max_archaea, max_plasmid (sample from
    the file with sample_seed for reproducibility). If save_accessions_to is set,
    write the accession list and current UTC timestamp to that file after searching.
    When not using accessions_file, complete_only=True restricts NCBI search to complete
    genomes only (excludes WGS/draft via complete[Properties] and NOT WGS[Properties]).
    For reproducible complete-only runs, create a snapshot with --complete-only then use
    that JSON as --accessions-file.
    """
    if num_bacteria < 0:
        raise ValueError("num_bacteria must be >= 0")
    if num_virus < 0:
        raise ValueError("num_virus must be >= 0")
    for name, val in [("max_bacteria", max_bacteria), ("max_virus", max_virus), ("max_archaea", max_archaea), ("max_plasmid", max_plasmid)]:
        if val is not None and val < 0:
            raise ValueError(f"{name} must be >= 0 when set")
    output_dir.mkdir(parents=True, exist_ok=True)
    bacteria_dir = output_dir / BACTERIA_DIR
    virus_dir = output_dir / VIRUS_DIR
    bacteria_dir.mkdir(parents=True, exist_ok=True)
    virus_dir.mkdir(parents=True, exist_ok=True)

    if accessions_file is not None:
        data = load_accessions(accessions_file)
        bacterial_ids, viral_ids, archaea_ids, plasmid_ids = get_accession_lists_from_data(data)
        ts = data.get(ACCESSIONS_KEY_TIMESTAMP, "unknown")
        print(f"Using accessions from {accessions_file} (timestamp: {ts})")
        print(f"  Bacterial: {len(bacterial_ids)}, viral: {len(viral_ids)}, archaea: {len(archaea_ids)}, plasmid: {len(plasmid_ids)}")

        def _sample(lst: list[str], max_n: int | None, seed: int | None) -> list[str]:
            if max_n is None or len(lst) <= max_n:
                return lst
            rng = random.Random(seed)
            return list(rng.sample(lst, max_n))

        if max_bacteria is not None or max_virus is not None or max_archaea is not None or max_plasmid is not None:
            seed = sample_seed if sample_seed is not None else 42
            bacterial_ids = _sample(bacterial_ids, max_bacteria, seed)
            viral_ids = _sample(viral_ids, max_virus, seed)
            archaea_ids = _sample(archaea_ids, max_archaea, seed)
            plasmid_ids = _sample(plasmid_ids, max_plasmid, seed)
            print(f"  Sampled (seed={seed}): bacterial={len(bacterial_ids)}, viral={len(viral_ids)}, archaea={len(archaea_ids)}, plasmid={len(plasmid_ids)}")
    else:
        queries = get_queries(complete_only=complete_only)
        if complete_only:
            print("Complete-only: restricting to complete genomes (excluding WGS/draft).")
        print("Searching for bacterial genomes...")
        bacterial_ids = search_genomes(queries["bacterial"], num_bacteria)
        print(f"  Found {len(bacterial_ids)} IDs")

        print("Searching for viral genomes...")
        viral_ids = search_genomes(queries["viral"], num_virus)
        print(f"  Found {len(viral_ids)} IDs")

        archaea_ids = []
        if num_archaea > 0:
            print("Searching for archaeal genomes...")
            archaea_ids = search_genomes(queries["archaea"], num_archaea)
            print(f"  Found {len(archaea_ids)} IDs")

        plasmid_ids = []
        if num_plasmid > 0:
            print("Searching for plasmids...")
            plasmid_ids = search_genomes(queries["plasmid"], num_plasmid)
            print(f"  Found {len(plasmid_ids)} IDs")

        if save_accessions_to is not None:
            save_accessions(
                save_accessions_to,
                bacterial_ids,
                viral_ids,
                archaea_ids,
                plasmid_ids,
            )
            print(f"Saved accession list to {save_accessions_to}")

    _download_category_batched(bacterial_ids, bacteria_dir, BACTERIA_DIR, "bacteria")
    _download_category_batched(viral_ids, virus_dir, VIRUS_DIR, "virus")

    if archaea_ids:
        archaea_dir = output_dir / ARCHAEA_DIR
        archaea_dir.mkdir(parents=True, exist_ok=True)
        _download_category_batched(archaea_ids, archaea_dir, ARCHAEA_DIR, "archaea")

    if plasmid_ids:
        plasmid_dir = output_dir / PLASMID_DIR
        plasmid_dir.mkdir(parents=True, exist_ok=True)
        _download_category_batched(plasmid_ids, plasmid_dir, PLASMID_DIR, "plasmid")

    logger.info("Download step complete: genomes saved in %s", output_dir.resolve())
    print("Done. Genome FASTAs saved in", output_dir.resolve())


def _cli(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="Download bacterial and viral genomes (and optionally archaea, plasmid) from NCBI."
    )
    parser.add_argument(
        "--num-bacteria",
        type=int,
        default=10,
        help="Number of bacterial genomes to download. Default: 10",
    )
    parser.add_argument(
        "--num-virus",
        type=int,
        default=10,
        help="Number of viral genomes to download. Default: 10",
    )
    parser.add_argument(
        "--num-archaea",
        type=int,
        default=0,
        help="Number of archaeal genomes to download (negative samples). Default: 0",
    )
    parser.add_argument(
        "--num-plasmid",
        type=int,
        default=0,
        help="Number of plasmids to download (negative samples). Default: 0",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("output"),
        help="Directory for downloaded genome FASTAs. Default: output",
    )
    parser.add_argument(
        "--accessions-file",
        type=Path,
        default=None,
        metavar="PATH",
        help="Load accession IDs from JSON (skip NCBI search). File must contain bacterial, viral (and optionally archaea, plasmid) lists and may include a timestamp.",
    )
    parser.add_argument(
        "--save-accessions",
        type=Path,
        default=None,
        metavar="PATH",
        help="After searching NCBI, save accession list and UTC timestamp to this JSON for reproducible runs (use with --accessions-file later).",
    )
    parser.add_argument(
        "--max-bacteria",
        type=int,
        default=None,
        metavar="N",
        help="When using --accessions-file: use at most N bacterial accessions (random sample). Omit to use all.",
    )
    parser.add_argument(
        "--max-virus",
        type=int,
        default=None,
        metavar="N",
        help="When using --accessions-file: use at most N viral accessions (random sample). Omit to use all.",
    )
    parser.add_argument(
        "--max-archaea",
        type=int,
        default=None,
        metavar="N",
        help="When using --accessions-file: use at most N archaeal accessions (random sample). Omit to use all.",
    )
    parser.add_argument(
        "--max-plasmid",
        type=int,
        default=None,
        metavar="N",
        help="When using --accessions-file: use at most N plasmid accessions (random sample). Omit to use all.",
    )
    parser.add_argument(
        "--sample-seed",
        type=int,
        default=None,
        metavar="SEED",
        help="When using --max-* with --accessions-file: random seed for sampling (default 42).",
    )
    args = parser.parse_args(argv)

    if args.accessions_file is not None and not args.accessions_file.exists():
        parser.error(f"--accessions-file not found: {args.accessions_file}")

    download_genomes(
        args.num_bacteria,
        args.num_virus,
        args.output_dir,
        num_archaea=args.num_archaea,
        num_plasmid=args.num_plasmid,
        accessions_file=args.accessions_file,
        save_accessions_to=args.save_accessions,
        max_bacteria=getattr(args, "max_bacteria", None),
        max_virus=getattr(args, "max_virus", None),
        max_archaea=getattr(args, "max_archaea", None),
        max_plasmid=getattr(args, "max_plasmid", None),
        sample_seed=getattr(args, "sample_seed", None),
    )


def main() -> None:
    _cli()


if __name__ == "__main__":
    main()
