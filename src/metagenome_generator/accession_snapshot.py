#!/usr/bin/env python3
from __future__ import annotations
"""Accession snapshot: catalog all matching NCBI accession IDs without downloading sequences.

Queries NCBI for every bacterial, viral, archaeal, and plasmid genome matching the
project's RefSeq + complete-genome criteria and saves the lists to the snapshots/
folder. No sequences are downloaded—only presence in the database is recorded.
The JSON can be used with download/pipeline --accessions-file (e.g. after trimming
to a subset) and can be updated manually or re-run to refresh the catalog.

Uses NCBI History server and paging (retmax 10,000) to retrieve full result sets.
A local UTC timestamp is recorded since NCBI does not expose a database "as of" date.
"""

import argparse
import json
import logging
import sys
import time
from pathlib import Path

from .download_genomes import (
    get_accession_lists_from_data,
    get_accession_metadata_from_data,
    ACCESSIONS_KEY_ARCHAEA,
    ACCESSIONS_KEY_BACTERIAL,
    ACCESSIONS_KEY_PLASMID,
    ACCESSIONS_KEY_TIMESTAMP,
    ACCESSIONS_KEY_VIRAL,
)
from .ncbi_search import get_queries, search_genomes_all
from .temporal_split import fetch_accession_metadata

logger = logging.getLogger(__name__)

# Default path: snapshots folder in project; filename includes run date (YYYY-MM-DD)
SNAPSHOTS_DIR = Path("snapshots")


def get_default_snapshot_path() -> Path:
    """Return snapshots/accession_snapshot_YYYY-MM-DD.json using current UTC date."""
    date_str = time.strftime("%Y-%m-%d", time.gmtime())
    return SNAPSHOTS_DIR / f"accession_snapshot_{date_str}.json"


def _default_log_path(output_path: Path) -> Path:
    """Default log path next to the snapshot JSON: snapshot_YYYY-MM-DD.log."""
    date_str = time.strftime("%Y-%m-%d", time.gmtime())
    return output_path.parent / f"snapshot_{date_str}.log"


class _Tee:
    """Write to both stdout and a log file."""

    def __init__(self, log_path: Path):
        self._file = log_path.open("w", encoding="utf-8")
        self._stdout = sys.stdout

    def write(self, data: str) -> int:
        self._file.write(data)
        self._file.flush()
        return self._stdout.write(data)

    def flush(self) -> None:
        self._file.flush()
        self._stdout.flush()

    def close(self) -> None:
        self._file.close()


def run_snapshot(
    output_path: Path,
    *,
    fetch_metadata: bool = True,
    metadata_batch_size: int = 500,
    log_path: Path | None = None,
    complete_only: bool = False,
) -> None:
    """Query NCBI for all matching accession IDs (no downloads), write JSON to snapshots.

    Fetches every bacterial, viral, archaeal, and plasmid genome matching
    the RefSeq queries (complete genome, length filters). If complete_only is True,
    uses stricter NCBI filters (complete[Properties], NOT WGS[Properties]). Output
    is compatible with download_genomes --accessions-file. By default, CreateDate
    and Title are stored per category. Use fetch_metadata=False for ID lists only.
    If log_path is set, progress is written to that file (default: snapshot_YYYY-MM-DD.log).
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    log_file = log_path if log_path is not None else _default_log_path(output_path)
    tee = _Tee(log_file)
    try:
        sys.stdout = tee
        _run_snapshot_impl(
            output_path,
            fetch_metadata=fetch_metadata,
            metadata_batch_size=metadata_batch_size,
            complete_only=complete_only,
        )
    finally:
        sys.stdout = tee._stdout
        tee.close()
        print(f"Log saved to {log_file}", file=sys.stdout)


def _category_with_metadata(
    ids: list[str],
    metadata: dict[str, dict],
) -> list[dict]:
    """Build list of {accession, create_date, title} for a category from ID list and flat metadata."""
    out = []
    for acc in ids:
        m = metadata.get(acc) or {}
        out.append({
            "accession": acc,
            "create_date": m.get("create_date") or "",
            "title": m.get("title") or "",
        })
    return out


def _run_snapshot_impl(
    output_path: Path,
    *,
    fetch_metadata: bool = True,
    metadata_batch_size: int = 500,
    complete_only: bool = False,
) -> None:
    """Implementation of run_snapshot (stdout may be redirected to Tee)."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    now = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    queries = get_queries(complete_only=complete_only)

    print("Querying NCBI nucleotide for all matching accessions (no sequences downloaded)...")
    if complete_only:
        print("  (complete-only: excluding WGS/draft; requiring complete[Properties])")
    search_start = time.time()

    def _search_progress(category: str):
        def _cb(page: int, total: int, ids_so_far: int) -> None:
            print(f"    {category}: page {page}/{total} ({ids_so_far:,} IDs)...", flush=True)
        return _cb

    print("  Bacterial (all RefSeq complete genomes in range)...")
    bacterial = search_genomes_all(
        queries["bacterial"],
        progress_callback=_search_progress("Bacterial"),
    )
    print(f"    Found {len(bacterial):,} IDs")
    print("  Viral (all RefSeq complete genomes in range)...")
    viral = search_genomes_all(
        queries["viral"],
        progress_callback=_search_progress("Viral"),
    )
    print(f"    Found {len(viral):,} IDs")
    print("  Archaea (all RefSeq complete genomes in range)...")
    archaea = search_genomes_all(
        queries["archaea"],
        progress_callback=_search_progress("Archaea"),
    )
    print(f"    Found {len(archaea):,} IDs")
    print("  Plasmid (all RefSeq in range)...")
    plasmid = search_genomes_all(
        queries["plasmid"],
        progress_callback=_search_progress("Plasmid"),
    )
    print(f"    Found {len(plasmid):,} IDs")
    print(f"  Search phase done in {time.time() - search_start:.0f}s")

    if fetch_metadata:
        all_ids = list(dict.fromkeys(bacterial + viral + archaea + plasmid))
        if all_ids:
            total_batches = (len(all_ids) + metadata_batch_size - 1) // metadata_batch_size
            print(f"  Fetching CreateDate and title for {len(all_ids):,} accessions (~{total_batches} batches of {metadata_batch_size})...")
            meta_start = time.time()

            def _progress(batch_num: int, total: int, fetched: int) -> None:
                pct = 100 * batch_num / total if total else 0
                elapsed = time.time() - meta_start
                print(f"    Metadata: batch {batch_num}/{total} ({pct:.0f}%) | {fetched:,} accessions | {elapsed:.0f}s elapsed", flush=True)

            metadata = fetch_accession_metadata(
                all_ids,
                batch_size=metadata_batch_size,
                progress_callback=_progress,
            )
            data = {
                ACCESSIONS_KEY_TIMESTAMP: now,
                ACCESSIONS_KEY_BACTERIAL: _category_with_metadata(bacterial, metadata),
                ACCESSIONS_KEY_VIRAL: _category_with_metadata(viral, metadata),
                ACCESSIONS_KEY_ARCHAEA: _category_with_metadata(archaea, metadata),
                ACCESSIONS_KEY_PLASMID: _category_with_metadata(plasmid, metadata),
            }
            with output_path.open("w", encoding="utf-8") as f:
                json.dump(data, f, indent=2, ensure_ascii=False)
            print(f"    Stored metadata for {len(metadata):,} accessions in {time.time() - meta_start:.0f}s; snapshot saved.")
        else:
            data = {
                ACCESSIONS_KEY_TIMESTAMP: now,
                ACCESSIONS_KEY_BACTERIAL: [],
                ACCESSIONS_KEY_VIRAL: [],
                ACCESSIONS_KEY_ARCHAEA: [],
                ACCESSIONS_KEY_PLASMID: [],
            }
            with output_path.open("w", encoding="utf-8") as f:
                json.dump(data, f, indent=2)
    else:
        data = {
            ACCESSIONS_KEY_TIMESTAMP: now,
            ACCESSIONS_KEY_BACTERIAL: bacterial,
            ACCESSIONS_KEY_VIRAL: viral,
            ACCESSIONS_KEY_ARCHAEA: archaea,
            ACCESSIONS_KEY_PLASMID: plasmid,
        }
        with output_path.open("w", encoding="utf-8") as f:
            json.dump(data, f, indent=2)
        print(f"Snapshot (accession lists) saved to {output_path}")

    logger.info("Accession snapshot written to %s (timestamp=%s)", output_path, now)
    print(f"Snapshot saved to {output_path} (timestamp {now})")
    print("Use with: download --accessions-file ... or pipeline --accessions-file ...")
    print("Edit the JSON to subset accessions, then run download with --accessions-file.")


def migrate_snapshot_to_categories(path: Path) -> None:
    """Convert a snapshot from legacy format to category-inline format (no re-download).

    Reads the JSON file: removes ncbi_db_info; moves accession_metadata into each
    category so each accession is stored as {accession, create_date, title} within
    bacterial/viral/archaea/plasmid. Overwrites the file in place.
    """
    path = Path(path)
    with path.open(encoding="utf-8") as f:
        data = json.load(f)
    # Remove legacy top-level keys we no longer want
    data.pop("ncbi_db_info", None)
    metadata = get_accession_metadata_from_data(data)
    b_ids, v_ids, a_ids, p_ids = get_accession_lists_from_data(data)
    data[ACCESSIONS_KEY_BACTERIAL] = _category_with_metadata(b_ids, metadata)
    data[ACCESSIONS_KEY_VIRAL] = _category_with_metadata(v_ids, metadata)
    data[ACCESSIONS_KEY_ARCHAEA] = _category_with_metadata(a_ids, metadata)
    data[ACCESSIONS_KEY_PLASMID] = _category_with_metadata(p_ids, metadata)
    data.pop("accession_metadata", None)
    with path.open("w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
    total = len(b_ids) + len(v_ids) + len(a_ids) + len(p_ids)
    logger.info("Migrated snapshot to category-inline format: %s (%s accessions)", path, total)


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
        "--no-metadata",
        action="store_true",
        help="Do not fetch CreateDate and title per accession (lists only; temporal-split will fetch dates when needed).",
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
        fetch_metadata=not args.no_metadata,
        metadata_batch_size=args.metadata_batch_size,
    )


def main() -> None:
    _cli()


if __name__ == "__main__":
    main()
