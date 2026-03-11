#!/usr/bin/env python3
from __future__ import annotations
"""Shared NCBI Entrez search logic for nucleotide genome queries.

Used by download_genomes and accession_snapshot. Single place for query definitions
and search so both modules stay in sync and avoid duplication.
"""

import os
import time
from pathlib import Path
from typing import Callable

from Bio import Entrez

Entrez.email = os.environ.get("ENTREZ_EMAIL", "your_email@example.com")
Entrez.api_key = os.environ.get("ENTREZ_API_KEY")

# Default search queries (RefSeq, complete genomes, length-filtered)
DEFAULT_QUERIES = {
    "bacterial": (
        "bacteria[Organism] AND complete genome[Title] "
        "AND refseq[filter] AND (10000:5000000[Sequence Length])"
    ),
    "viral": (
        "viruses[Organism] AND complete genome[Title] "
        "AND refseq[filter] AND (5000:500000[Sequence Length])"
    ),
    "archaea": (
        "archaea[Organism] AND complete genome[Title] "
        "AND refseq[filter] AND (10000:5000000[Sequence Length])"
    ),
    "plasmid": (
        "plasmid[Title] AND refseq[filter] AND (5000:500000[Sequence Length])"
    ),
}

# Stricter completeness: exclude WGS/draft; require complete [Properties] (NCBI nucleotide)
COMPLETE_ONLY_SUFFIX = " AND complete[Properties] AND NOT WGS[Properties]"


def get_queries(*, complete_only: bool = False) -> dict[str, str]:
    """Return category -> query dict. If complete_only, append NCBI completeness filters."""
    if not complete_only:
        return dict(DEFAULT_QUERIES)
    return {k: (q + COMPLETE_ONLY_SUFFIX) for k, q in DEFAULT_QUERIES.items()}


# NCBI limits esearch to 10,000 UIDs per request; we page with this batch size.
ESEARCH_BATCH_SIZE = 10_000


def search_genomes(query: str, count: int) -> list[str]:
    """Search NCBI nucleotide database and return up to `count` accession IDs."""
    time.sleep(0.4)
    handle = Entrez.esearch(
        db="nucleotide",
        term=query,
        retmax=max(count * 2, 100),
        sort="relevance",
        idtype="acc",
    )
    record = Entrez.read(handle)
    handle.close()
    return record["IdList"][:count]


def search_genomes_all(
    query: str,
    *,
    progress_callback: Callable[[int, int, int], None] | None = None,
) -> list[str]:
    """Retrieve all accession IDs matching the query using History server and paging.

    NCBI returns at most 10,000 UIDs per esearch request. We use usehistory=y to
    store the result on the server, then page with retstart/retmax to collect
    every ID. Rate-limiting sleep is applied between requests.
    progress_callback(page_num, total_pages, ids_so_far) is called after each page.
    """
    max_initial_retries = 3
    record = None
    for init_attempt in range(max_initial_retries):
        handle = None
        try:
            time.sleep(0.4)
            handle = Entrez.esearch(
                db="nucleotide",
                term=query,
                usehistory="y",
                retmax=0,
                sort="relevance",
                idtype="acc",
            )
            record = Entrez.read(handle)
            break
        except (RuntimeError, Exception):
            if handle is not None:
                try:
                    handle.close()
                except Exception:
                    pass
            if init_attempt < max_initial_retries - 1:
                time.sleep(2.0 * (init_attempt + 1))
                continue
            raise
        finally:
            if handle is not None:
                try:
                    handle.close()
                except Exception:
                    pass
    if record is None:
        return []
    count = int(record["Count"])
    webenv = record["WebEnv"]
    query_key = record["QueryKey"]
    if count == 0:
        return []

    total_pages = (count + ESEARCH_BATCH_SIZE - 1) // ESEARCH_BATCH_SIZE
    all_ids: list[str] = []
    max_page_retries = 3
    for page_zero in range(total_pages):
        start = page_zero * ESEARCH_BATCH_SIZE
        time.sleep(0.4)
        for attempt in range(max_page_retries):
            handle = None
            try:
                handle = Entrez.esearch(
                    db="nucleotide",
                    term=query,
                    webenv=webenv,
                    query_key=query_key,
                    retstart=start,
                    retmax=ESEARCH_BATCH_SIZE,
                    idtype="acc",
                )
                batch = Entrez.read(handle)
                all_ids.extend(batch["IdList"])
                break
            except (RuntimeError, Exception):
                if attempt < max_page_retries - 1:
                    time.sleep(2.0 * (attempt + 1))
                    continue
                raise
            finally:
                if handle is not None:
                    try:
                        handle.close()
                    except Exception:
                        pass
        if progress_callback is not None:
            progress_callback(page_zero + 1, total_pages, len(all_ids))
    return all_ids
