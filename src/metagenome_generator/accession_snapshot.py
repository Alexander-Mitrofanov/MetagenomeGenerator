#!/usr/bin/env python3
from __future__ import annotations
"""Accession snapshot: catalog all matching NCBI accession IDs without downloading sequences.

Queries NCBI for every bacterial, viral, archaeal, and plasmid genome matching the
project's RefSeq + complete-genome criteria and saves the lists to the snapshots/
folder. No sequences are downloaded—only presence in the database is recorded.
The JSON can be used with download/pipeline --accessions-file (e.g. after trimming
to a subset) and can be updated manually or re-run to refresh the catalog.

Uses NCBI History server and paging (retmax 10,000) to retrieve full result sets.
Optionally stores nucleotide DB metadata (einfo). A local UTC timestamp is recorded
since NCBI does not expose a database "as of" date.
"""

import argparse
import json
import logging
import time
from pathlib import Path

from Bio import Entrez

from .download_genomes import (
    ACCESSION_METADATA_KEY,
    ACCESSIONS_KEY_ARCHAEA,
    ACCESSIONS_KEY_BACTERIAL,
    ACCESSIONS_KEY_PLASMID,
    ACCESSIONS_KEY_TIMESTAMP,
    ACCESSIONS_KEY_VIRAL,
)
from .ncbi_search import DEFAULT_QUERIES, search_genomes_all
from .temporal_split import fetch_accession_metadata

logger = logging.getLogger(__name__)

# Optional key for NCBI db metadata (einfo result)
NCBI_DB_INFO_KEY = "ncbi_db_info"

# Default path: snapshots folder in project; filename includes run date (YYYY-MM-DD)
SNAPSHOTS_DIR = Path("snapshots")


def get_default_snapshot_path() -> Path:
    """Return snapshots/accession_snapshot_YYYY-MM-DD.json using current UTC date."""
    date_str = time.strftime("%Y-%m-%d", time.gmtime())
    return SNAPSHOTS_DIR / f"accession_snapshot_{date_str}.json"


def fetch_nucleotide_db_info() -> dict | None:
    """Fetch nucleotide database info from NCBI via einfo. Returns a small dict or None on failure.

    NCBI EInfo does not provide a 'last updated' timestamp; we may get DbName, Count, etc.
    Useful for documentation (e.g. total record count at snapshot time).
    """
    try:
        handle = Entrez.einfo(db="nucleotide")
        record = Entrez.read(handle)
        handle.close()
        info = record.get("DbInfo", {})
        # Keep only simple, JSON-serializable fields
        return {
            "DbName": str(info.get("DbName", "")),
            "Count": int(info.get("Count", 0)),
            "MenuName": str(info.get("MenuName", "")),
        }
    except Exception as e:
        logger.debug("einfo for nucleotide failed: %s", e)
        return None


def run_snapshot(
    output_path: Path,
    *,
    fetch_db_info: bool = True,
    fetch_metadata: bool = True,
    metadata_batch_size: int = 500,
) -> None:
    """Query NCBI for all matching accession IDs (no downloads), write JSON to snapshots.

    Fetches every bacterial, viral, archaeal, and plasmid genome matching
    DEFAULT_QUERIES (RefSeq, complete genome, length filters). Output is compatible
    with download_genomes --accessions-file. Optionally includes ncbi_db_info.
    When fetch_metadata is True, also fetches CreateDate and Title (genome header)
    per accession via esummary and stores them under accession_metadata; temporal-split
    can then use this to avoid a second NCBI round-trip.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    now = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    ncbi_info = fetch_nucleotide_db_info() if fetch_db_info else None

    print("Querying NCBI nucleotide for all matching accessions (no sequences downloaded)...")
    print("  Bacterial (all RefSeq complete genomes in range)...")
    bacterial = search_genomes_all(DEFAULT_QUERIES["bacterial"])
    print(f"    Found {len(bacterial)} IDs")
    print("  Viral (all RefSeq complete genomes in range)...")
    viral = search_genomes_all(DEFAULT_QUERIES["viral"])
    print(f"    Found {len(viral)} IDs")
    print("  Archaea (all RefSeq complete genomes in range)...")
    archaea = search_genomes_all(DEFAULT_QUERIES["archaea"])
    print(f"    Found {len(archaea)} IDs")
    print("  Plasmid (all RefSeq in range)...")
    plasmid = search_genomes_all(DEFAULT_QUERIES["plasmid"])
    print(f"    Found {len(plasmid)} IDs")

    data = {
        ACCESSIONS_KEY_TIMESTAMP: now,
        ACCESSIONS_KEY_BACTERIAL: bacterial,
        ACCESSIONS_KEY_VIRAL: viral,
        ACCESSIONS_KEY_ARCHAEA: archaea,
        ACCESSIONS_KEY_PLASMID: plasmid,
    }
    if ncbi_info:
        data[NCBI_DB_INFO_KEY] = ncbi_info
        print(f"  NCBI nucleotide db: {ncbi_info.get('Count', '?')} records (einfo)")

    # Write snapshot immediately so we have a valid file even if metadata fetch fails
    with output_path.open("w") as f:
        json.dump(data, f, indent=2)
    print(f"Snapshot (accession lists) saved to {output_path}")

    if fetch_metadata:
        all_ids = list(dict.fromkeys(bacterial + viral + archaea + plasmid))
        if all_ids:
            total_batches = (len(all_ids) + metadata_batch_size - 1) // metadata_batch_size
            print(f"  Fetching CreateDate and title for {len(all_ids)} accessions (batches of {metadata_batch_size}, ~{total_batches} batches)...")

            def _progress(batch_num: int, total: int, fetched: int) -> None:
                print(f"    Metadata: batch {batch_num}/{total}, {fetched} accessions so far")

            metadata = fetch_accession_metadata(
                all_ids,
                batch_size=metadata_batch_size,
                progress_callback=_progress,
            )
            data[ACCESSION_METADATA_KEY] = metadata
            with output_path.open("w") as f:
                json.dump(data, f, indent=2)
            print(f"    Stored metadata for {len(metadata)} accessions; snapshot updated.")

    logger.info("Accession snapshot written to %s (timestamp=%s)", output_path, now)
    print(f"Snapshot saved to {output_path} (timestamp {now})")
    print("Use with: download --accessions-file ... or pipeline --accessions-file ...")
    print("Edit the JSON to subset accessions, then run download with --accessions-file.")


def _cli(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="Catalog all NCBI accession IDs for bacteria, virus, archaea, plasmid (no downloads). "
        "Saves to snapshots/ for use with --accessions-file; file can be edited to subset."
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output JSON path. Default: snapshots/accession_snapshot_YYYY-MM-DD.json (run date)",
    )
    parser.add_argument(
        "--no-db-info",
        action="store_true",
        help="Do not fetch NCBI nucleotide db metadata (einfo).",
    )
    parser.add_argument(
        "--no-metadata",
        action="store_true",
        help="Do not fetch CreateDate and title per accession (faster snapshot; temporal-split will fetch dates when needed).",
    )
    parser.add_argument(
        "--metadata-batch-size",
        type=int,
        default=500,
        help="Batch size for esummary when fetching metadata (default 500).",
    )
    args = parser.parse_args(argv)
    output_path = args.output or get_default_snapshot_path()
    run_snapshot(
        output_path,
        fetch_db_info=not args.no_db_info,
        fetch_metadata=not args.no_metadata,
        metadata_batch_size=args.metadata_batch_size,
    )


def main() -> None:
    _cli()


if __name__ == "__main__":
    main()
