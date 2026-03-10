#!/usr/bin/env python3
from __future__ import annotations
"""Shared NCBI Entrez search logic for nucleotide genome queries.

Used by download_genomes and accession_snapshot. Single place for query definitions
and search so both modules stay in sync and avoid duplication.
"""

import os
import time
from pathlib import Path

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


def search_genomes_all(query: str) -> list[str]:
    """Retrieve all accession IDs matching the query using History server and paging.

    NCBI returns at most 10,000 UIDs per esearch request. We use usehistory=y to
    store the result on the server, then page with retstart/retmax to collect
    every ID. Rate-limiting sleep is applied between requests.
    """
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
    handle.close()
    count = int(record["Count"])
    webenv = record["WebEnv"]
    query_key = record["QueryKey"]
    if count == 0:
        return []

    all_ids: list[str] = []
    for start in range(0, count, ESEARCH_BATCH_SIZE):
        time.sleep(0.4)
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
        handle.close()
        all_ids.extend(batch["IdList"])
    return all_ids
