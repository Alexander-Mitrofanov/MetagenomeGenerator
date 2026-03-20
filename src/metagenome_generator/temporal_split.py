#!/usr/bin/env python3
"""Temporal train/test split of accession lists by NCBI CreateDate.

Enables DeepVirFinder/VirFinder-style evaluation: train on sequences "discovered"
before a cutoff date, test on sequences added on or after that date. Uses
NCBI esummary (CreateDate) to partition accessions. Outputs two JSON files
compatible with download --accessions-file.
"""

from __future__ import annotations

import json
import logging
import os
import time
from pathlib import Path
from typing import Callable

from Bio import Entrez

Entrez.email = os.environ.get("ENTREZ_EMAIL", "your_email@example.com")
Entrez.api_key = os.environ.get("ENTREZ_API_KEY")

from .download_genomes import (
    ACCESSIONS_KEY_ARCHAEA,
    ACCESSIONS_KEY_BACTERIAL,
    ACCESSIONS_KEY_PLASMID,
    ACCESSIONS_KEY_TIMESTAMP,
    ACCESSIONS_KEY_VIRAL,
    get_accession_lists_from_data,
    get_accession_metadata_from_data,
    load_accessions,
)

logger = logging.getLogger(__name__)

# NCBI rate limiting; esummary supports batches (500 is typically accepted)
ESUMMARY_BATCH_SIZE = 500
ESUMMARY_SLEEP = 0.34  # ~3 req/s with key
ESUMMARY_MAX_RETRIES = 3
ESUMMARY_RETRY_SLEEP = 1.0


def _parse_create_date(record: dict) -> str | None:
    """Get CreateDate from esummary record. NCBI format: YYYY/MM/DD."""
    raw = record.get("CreateDate")
    if not raw:
        return None
    return str(raw).strip()


def _parse_title(record: dict) -> str:
    """Get Title from esummary record (genome description / FASTA header text)."""
    raw = record.get("Title")
    return str(raw).strip() if raw else ""


def fetch_accession_metadata(
    accessions: list[str],
    *,
    batch_size: int = ESUMMARY_BATCH_SIZE,
    max_retries: int = ESUMMARY_MAX_RETRIES,
    progress_callback: Callable[[int, int, int], None] | None = None,
) -> dict[str, dict]:
    """Fetch NCBI nucleotide CreateDate and Title for each accession via esummary.

    Returns dict accession_id -> {"create_date": "YYYY/MM/DD", "title": "..."}.
    Title is the genome description (equivalent to FASTA header after the accession).
    Missing or failed IDs are absent from the result. Retries each batch on failure.
    progress_callback(batch_index, total_batches, fetched_count) is called periodically.
    """
    result: dict[str, dict] = {}
    total_batches = (len(accessions) + batch_size - 1) // batch_size
    for i in range(0, len(accessions), batch_size):
        batch = accessions[i : i + batch_size]
        batch_idx = i // batch_size
        id_str = ",".join(batch)
        time.sleep(ESUMMARY_SLEEP)
        for attempt in range(max_retries):
            try:
                handle = Entrez.esummary(db="nucleotide", id=id_str)
                records = Entrez.read(handle)
                handle.close()
                break
            except Exception as e:
                logger.warning("esummary batch %d/%d attempt %d failed: %s", batch_idx + 1, total_batches, attempt + 1, e)
                if attempt < max_retries - 1:
                    time.sleep(ESUMMARY_RETRY_SLEEP * (attempt + 1))
                else:
                    continue
        else:
            continue
        for rec in records:
            acc = rec.get("AccessionVersion") or rec.get("Caption") or str(rec.get("Id", ""))
            if acc:
                acc = str(acc)
                create_date = _parse_create_date(rec)
                title = _parse_title(rec)
                result[acc] = {"create_date": create_date or "", "title": title}
        if progress_callback and (
            (batch_idx + 1) % 5 == 0 or (batch_idx + 1) == total_batches
        ):
            progress_callback(batch_idx + 1, total_batches, len(result))
    return result


def fetch_accession_dates(
    accessions: list[str],
    *,
    batch_size: int = ESUMMARY_BATCH_SIZE,
) -> dict[str, str]:
    """Fetch NCBI nucleotide CreateDate for each accession via esummary.

    Returns dict accession_id -> "YYYY/MM/DD". Missing or failed IDs are absent.
    """
    result: dict[str, str] = {}
    for i in range(0, len(accessions), batch_size):
        batch = accessions[i : i + batch_size]
        id_str = ",".join(batch)
        time.sleep(ESUMMARY_SLEEP)
        try:
            handle = Entrez.esummary(db="nucleotide", id=id_str)
            records = Entrez.read(handle)
            handle.close()
        except Exception as e:
            logger.warning("esummary batch failed for %d IDs: %s", len(batch), e)
            continue
        # esummary returns list of dicts; use AccessionVersion or Caption to match our IDs
        for rec in records:
            acc = rec.get("AccessionVersion") or rec.get("Caption") or str(rec.get("Id", ""))
            if acc:
                date_str = _parse_create_date(rec)
                if date_str:
                    result[str(acc)] = date_str
    return result


def _date_to_comparable(s: str) -> tuple[int, int, int]:
    """Parse YYYY/MM/DD or YYYY-MM-DD to (y, m, d) for comparison."""
    s = s.replace("-", "/")
    parts = s.split("/")
    if len(parts) != 3:
        return (0, 0, 0)
    try:
        return (int(parts[0]), int(parts[1]), int(parts[2]))
    except ValueError:
        return (0, 0, 0)


def _comparable_to_iso(d: tuple[int, int, int]) -> str:
    """Format (y, m, d) as YYYY-MM-DD."""
    return f"{d[0]:04d}-{d[1]:02d}-{d[2]:02d}"


def split_ids_by_date(
    ids: list[str],
    date_by_id: dict[str, str],
    cutoff_iso: str,
) -> tuple[list[str], list[str]]:
    """Split ID list into train (CreateDate < cutoff) and test (CreateDate >= cutoff).

    cutoff_iso: YYYY-MM-DD. IDs with no date are assigned to train (conservative).
    """
    cutoff = _date_to_comparable(cutoff_iso)
    train, test = [], []
    for acc in ids:
        date_str = date_by_id.get(acc)
        if not date_str:
            train.append(acc)
            continue
        acc_date = _date_to_comparable(date_str)
        if acc_date < cutoff:
            train.append(acc)
        else:
            test.append(acc)
    return train, test


def _validate_split_date(s: str) -> None:
    """Raise ValueError if s is not YYYY-MM-DD."""
    if len(s) != 10 or s[4] != "-" or s[7] != "-":
        raise ValueError(f"split_date must be YYYY-MM-DD, got: {s!r}")
    try:
        y, m, d = int(s[:4]), int(s[5:7]), int(s[8:10])
        if not (1 <= m <= 12 and 1 <= d <= 31):
            raise ValueError("invalid month or day")
    except (ValueError, TypeError) as e:
        raise ValueError(f"split_date must be YYYY-MM-DD, got: {s!r}") from e


def run_temporal_split_info(
    accessions_file: Path,
    split_date: str,
    *,
    batch_size: int = ESUMMARY_BATCH_SIZE,
    verbose: bool = True,
) -> dict:
    """Compute and optionally print train/test counts for a temporal split (no files written).

    Returns a dict with: split_date, source, timestamp, categories (bacterial, viral, archaea, plasmid),
    each with total, train, test, train_pct, test_pct; and totals (total, train, test).
    """
    _validate_split_date(split_date)
    path = Path(accessions_file)
    data = load_accessions(path)
    bacterial, viral, archaea, plasmid = get_accession_lists_from_data(data)

    all_ids = list(dict.fromkeys(bacterial + viral + archaea + plasmid))
    if not all_ids:
        raise ValueError(f"No accessions found in {path}")

    metadata = get_accession_metadata_from_data(data)
    if metadata:
        date_by_id = {
            acc: m["create_date"]
            for acc, m in metadata.items()
            if isinstance(m, dict) and m.get("create_date")
        }
        if verbose:
            print(f"Using CreateDate from snapshot metadata for {len(date_by_id):,} accessions.")
    else:
        if verbose:
            print(f"Fetching CreateDate for {len(all_ids):,} accessions from NCBI (batches of {batch_size})...")
        date_by_id = fetch_accession_dates(all_ids, batch_size=batch_size)

    missing = len(all_ids) - len(date_by_id)
    if missing and verbose:
        print(f"  Note: {missing} accessions had no CreateDate; counted as train.")

    def split_list(ids: list[str]) -> tuple[list[str], list[str]]:
        return split_ids_by_date(ids, date_by_id, split_date)

    categories = (
        ("bacterial", bacterial),
        ("viral", viral),
        ("archaea", archaea),
        ("plasmid", plasmid),
    )
    rows = []
    total_all, train_all, test_all = 0, 0, 0
    for name, ids in categories:
        train_ids, test_ids = split_list(ids)
        total = len(ids)
        train_n, test_n = len(train_ids), len(test_ids)
        total_all += total
        train_all += train_n
        test_all += test_n
        train_pct = (100.0 * train_n / total) if total else 0
        test_pct = (100.0 * test_n / total) if total else 0
        rows.append({
            "category": name,
            "total": total,
            "train": train_n,
            "test": test_n,
            "train_pct": train_pct,
            "test_pct": test_pct,
        })
    ts = data.get(ACCESSIONS_KEY_TIMESTAMP, "unknown")
    result = {
        "split_date": split_date,
        "source": str(path.resolve()),
        "timestamp": ts,
        "categories": {r["category"]: r for r in rows},
        "totals": {"total": total_all, "train": train_all, "test": test_all},
    }

    if verbose:
        print()
        print(f"Temporal split preview (split date: {split_date})")
        print(f"  Source: {path}")
        print(f"  Snapshot timestamp: {ts}")
        print()
        print("  Category   |   Total |   Train |    Test | Train % | Test %")
        print("  -----------+---------+---------+---------+---------+--------")
        for r in rows:
            print(f"  {r['category']:11} | {r['total']:>7,} | {r['train']:>7,} | {r['test']:>7,} | {r['train_pct']:>6.1f}% | {r['test_pct']:>5.1f}%")
        print("  -----------+---------+---------+---------+---------+--------")
        print(f"  {'Total':11} | {total_all:>7,} | {train_all:>7,} | {test_all:>7,} | {100.0 * train_all / total_all if total_all else 0:>6.1f}% | {100.0 * test_all / total_all if total_all else 0:>5.1f}%")
        print()
        print("  Train = CreateDate < split date (older entries)")
        print("  Test  = CreateDate >= split date (newer entries)")
        print()
        print("  Run 'temporal-split' with this --split-date to write train/test JSONs.")

    return result


def run_temporal_split_search(
    accessions_file: Path,
    min_train: int,
    min_test: int,
    *,
    batch_size: int = ESUMMARY_BATCH_SIZE,
    # Optional per-category minima.
    # For test: if minima for bacteria/virus are not explicitly provided, they default
    # to min_test (so suggested test splits don't end up with 0 viral reads).
    min_train_bacteria: int | None = None,
    min_train_viral: int | None = None,
    min_train_archaea: int | None = None,
    min_train_plasmid: int | None = None,
    min_test_bacteria: int | None = None,
    min_test_viral: int | None = None,
    min_test_archaea: int | None = None,
    min_test_plasmid: int | None = None,
    verbose: bool = True,
) -> dict:
    """Find a split date such that train set has at least min_train and test set at least min_test.

    Uses CreateDate from snapshot metadata or NCBI. Returns the latest (most recent) split date
    that satisfies both bounds, so the test set is as large as possible while keeping enough train.
    """
    path = Path(accessions_file)
    data = load_accessions(path)
    bacterial, viral, archaea, plasmid = get_accession_lists_from_data(data)
    all_ids = list(dict.fromkeys(bacterial + viral + archaea + plasmid))
    if not all_ids:
        raise ValueError(f"No accessions found in {path}")

    metadata = get_accession_metadata_from_data(data)
    if metadata:
        date_by_id = {
            acc: m["create_date"]
            for acc, m in metadata.items()
            if isinstance(m, dict) and m.get("create_date")
        }
        if verbose:
            print(f"Using CreateDate from snapshot metadata for {len(date_by_id):,} accessions.")
    else:
        if verbose:
            print(f"Fetching CreateDate for {len(all_ids):,} accessions from NCBI (batches of {batch_size})...")
        date_by_id = fetch_accession_dates(all_ids, batch_size=batch_size)

    missing = len(all_ids) - len(date_by_id)
    if missing and verbose:
        print(f"  Note: {missing} accessions had no CreateDate; counted as train.")

    acc_with_dates = [
        (acc, _date_to_comparable(date_by_id[acc]))
        for acc in all_ids
        if date_by_id.get(acc) and _date_to_comparable(date_by_id[acc]) != (0, 0, 0)
    ]
    acc_with_dates.sort(key=lambda x: x[1])
    unique_dates = []
    seen = set()
    for _acc, d in acc_with_dates:
        if d not in seen:
            seen.add(d)
            unique_dates.append(d)

    best_date_iso = None
    best_train = 0
    best_test = 0
    best_counts: dict[str, dict[str, int]] | None = None

    min_train_bacteria_eff = 0 if min_train_bacteria is None else min_train_bacteria
    min_train_viral_eff = 0 if min_train_viral is None else min_train_viral
    min_train_archaea_eff = 0 if min_train_archaea is None else min_train_archaea
    min_train_plasmid_eff = 0 if min_train_plasmid is None else min_train_plasmid

    min_test_bacteria_eff = min_test if min_test_bacteria is None else min_test_bacteria
    min_test_viral_eff = min_test if min_test_viral is None else min_test_viral
    min_test_archaea_eff = 0 if min_test_archaea is None else min_test_archaea
    min_test_plasmid_eff = 0 if min_test_plasmid is None else min_test_plasmid

    # Precompute category for each accession and comparable CreateDate tuples.
    acc_cat: dict[str, str] = {}
    for acc in bacterial:
        acc_cat[acc] = "bacterial"
    for acc in viral:
        acc_cat[acc] = "viral"
    for acc in archaea:
        acc_cat[acc] = "archaea"
    for acc in plasmid:
        acc_cat[acc] = "plasmid"

    date_comp_by_acc: dict[str, tuple[int, int, int]] = {}
    for acc, date_str in date_by_id.items():
        if not date_str:
            continue
        comp = _date_to_comparable(date_str)
        if comp != (0, 0, 0):
            date_comp_by_acc[acc] = comp

    for d in unique_dates:
        # Count train/test sizes in one pass (avoid materializing ID lists).
        b_train = b_test = 0
        v_train = v_test = 0
        a_train = a_test = 0
        p_train = p_test = 0
        train_total = 0
        test_total = 0

        for acc in all_ids:
            cat = acc_cat.get(acc)
            if cat is None:
                # Should not happen, but keep conservative behavior.
                train_total += 1
                continue

            acc_date = date_comp_by_acc.get(acc)
            # Missing or invalid CreateDate is counted as train (conservative).
            is_test = acc_date is not None and acc_date >= d

            if is_test:
                test_total += 1
                if cat == "bacterial":
                    b_test += 1
                elif cat == "viral":
                    v_test += 1
                elif cat == "archaea":
                    a_test += 1
                elif cat == "plasmid":
                    p_test += 1
            else:
                train_total += 1
                if cat == "bacterial":
                    b_train += 1
                elif cat == "viral":
                    v_train += 1
                elif cat == "archaea":
                    a_train += 1
                elif cat == "plasmid":
                    p_train += 1

        if (
            train_total >= min_train
            and test_total >= min_test
            and b_test >= min_test_bacteria_eff
            and v_test >= min_test_viral_eff
            and a_test >= min_test_archaea_eff
            and p_test >= min_test_plasmid_eff
            and b_train >= min_train_bacteria_eff
            and v_train >= min_train_viral_eff
            and a_train >= min_train_archaea_eff
            and p_train >= min_train_plasmid_eff
        ):
            split_iso = _comparable_to_iso(d)
            best_date_iso = split_iso
            best_train = train_total
            best_test = test_total
            best_counts = {
                "bacterial": {"train": b_train, "test": b_test},
                "viral": {"train": v_train, "test": v_test},
                "archaea": {"train": a_train, "test": a_test},
                "plasmid": {"train": p_train, "test": p_test},
            }

    if best_date_iso is None:
        raise ValueError(
            "No split date found that satisfies the requested minima. "
            f"Total train>={min_train}, total test>={min_test}, "
            f"test bacterial>={min_test_bacteria_eff}, test viral>={min_test_viral_eff}, "
            f"test archaea>={min_test_archaea_eff}, test plasmid>={min_test_plasmid_eff}. "
            "Try lowering min_test (or set explicit per-category minima), or use a snapshot with more accessions."
        )
    assert best_counts is not None

    result = {
        "suggested_date": best_date_iso,
        "train_count": best_train,
        "test_count": best_test,
        "per_category": {
            "bacterial": best_counts["bacterial"],
            "viral": best_counts["viral"],
            "archaea": best_counts["archaea"],
            "plasmid": best_counts["plasmid"],
        },
        "minima_used": {
            "min_train": min_train,
            "min_test": min_test,
            "min_train_bacteria": min_train_bacteria_eff,
            "min_train_viral": min_train_viral_eff,
            "min_train_archaea": min_train_archaea_eff,
            "min_train_plasmid": min_train_plasmid_eff,
            "min_test_bacteria": min_test_bacteria_eff,
            "min_test_viral": min_test_viral_eff,
            "min_test_archaea": min_test_archaea_eff,
            "min_test_plasmid": min_test_plasmid_eff,
        },
    }

    if verbose:
        print()
        print(f"Suggested --split-date: {best_date_iso}  (train={best_train:,}, test={best_test:,})")
        print()
        print("  Requested minima (effective):")
        print("  ------------------------------")
        print(f"  Total train         >= {min_train:,}")
        print(f"  Total test          >= {min_test:,}")
        print(f"  Test bacterial      >= {min_test_bacteria_eff:,}")
        print(f"  Test viral          >= {min_test_viral_eff:,}")
        print(f"  Test archaea        >= {min_test_archaea_eff:,}")
        print(f"  Test plasmid        >= {min_test_plasmid_eff:,}")
        print()
        print("  Category   |   Train |    Test")
        print("  -----------+---------+--------")
        for cat, key in [("bacterial", "bacterial"), ("viral", "viral"), ("archaea", "archaea"), ("plasmid", "plasmid")]:
            r = result["per_category"][key]
            print(f"  {cat:11} | {r['train']:>7,} | {r['test']:>7,}")
        print("  -----------+---------+--------")
        print(f"  {'Total':11} | {best_train:>7,} | {best_test:>7,}")
        print()
        print("  Run: metagenome-generator temporal-split --accessions-file <path> --split-date", best_date_iso)

    return result


def run_temporal_split(
    accessions_file: Path,
    split_date: str,
    output_train: Path,
    output_test: Path,
    *,
    batch_size: int = ESUMMARY_BATCH_SIZE,
) -> None:
    """Load accessions JSON, fetch CreateDates, split by split_date, write train/test JSONs.

    split_date: YYYY-MM-DD (e.g. 2015-05-01). Train = CreateDate < split_date, Test = >=.
    Output files have the same structure as the input (bacterial, viral, archaea, plasmid,
    timestamp) and can be used with download --accessions-file.
    """
    _validate_split_date(split_date)
    data = load_accessions(accessions_file)
    bacterial, viral, archaea, plasmid = get_accession_lists_from_data(data)

    all_ids = list(dict.fromkeys(bacterial + viral + archaea + plasmid))
    if not all_ids:
        raise ValueError(f"No accessions found in {accessions_file}")

    # Use metadata from snapshot if present (legacy accession_metadata or per-category objects)
    metadata = get_accession_metadata_from_data(data)
    if metadata:
        date_by_id = {
            acc: m["create_date"]
            for acc, m in metadata.items()
            if isinstance(m, dict) and m.get("create_date")
        }
        print(f"Using CreateDate from snapshot metadata for {len(date_by_id)} accessions (skipping NCBI fetch).")
    else:
        print(f"Fetching CreateDate for {len(all_ids)} accessions from NCBI (batches of {batch_size})...")
        date_by_id = fetch_accession_dates(all_ids, batch_size=batch_size)
    missing = len(all_ids) - len(date_by_id)
    if missing:
        print(f"  Note: {missing} accessions had no CreateDate; assigned to train set.")

    def split_list(ids: list[str]) -> tuple[list[str], list[str]]:
        return split_ids_by_date(ids, date_by_id, split_date)

    b_train, b_test = split_list(bacterial)
    v_train, v_test = split_list(viral)
    a_train, a_test = split_list(archaea)
    p_train, p_test = split_list(plasmid)

    now = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    train_ts = f"temporal split train (CreateDate < {split_date}); generated {now}"
    test_ts = f"temporal split test (CreateDate >= {split_date}); generated {now}"

    output_train.parent.mkdir(parents=True, exist_ok=True)
    output_test.parent.mkdir(parents=True, exist_ok=True)

    train_data = {
        ACCESSIONS_KEY_TIMESTAMP: train_ts,
        ACCESSIONS_KEY_BACTERIAL: b_train,
        ACCESSIONS_KEY_VIRAL: v_train,
        ACCESSIONS_KEY_ARCHAEA: a_train,
        ACCESSIONS_KEY_PLASMID: p_train,
    }
    test_data = {
        ACCESSIONS_KEY_TIMESTAMP: test_ts,
        ACCESSIONS_KEY_BACTERIAL: b_test,
        ACCESSIONS_KEY_VIRAL: v_test,
        ACCESSIONS_KEY_ARCHAEA: a_test,
        ACCESSIONS_KEY_PLASMID: p_test,
    }

    with output_train.open("w") as f:
        json.dump(train_data, f, indent=2)
    with output_test.open("w") as f:
        json.dump(test_data, f, indent=2)

    print(f"Train: {len(b_train)} bacterial, {len(v_train)} viral, {len(a_train)} archaea, {len(p_train)} plasmid -> {output_train}")
    print(f"Test:  {len(b_test)} bacterial, {len(v_test)} viral, {len(a_test)} archaea, {len(p_test)} plasmid -> {output_test}")
    logger.info("Temporal split: train=%s test=%s", output_train, output_test)
