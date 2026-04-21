#!/usr/bin/env python3
from __future__ import annotations
"""Fetch viral taxonomy from NCBI and produce accession->group mapping for taxonomy-aware balancing.

Uses elink (nuccore -> taxonomy) to get taxid per accession, then efetch (taxonomy)
to get lineage; extracts family or realm for each. Output JSON maps each viral accession
to taxonomy group (e.g. NC_001234.1 -> Herpesviridae).
"""

import json
import logging
import time
from pathlib import Path

from Bio import Entrez
import os

Entrez.email = os.environ.get("ENTREZ_EMAIL", "your_email@example.com")
Entrez.api_key = os.environ.get("ENTREZ_API_KEY")

from .download_genomes import get_accession_lists_from_data, load_accessions

logger = logging.getLogger(__name__)

BATCH_SIZE = 100
SLEEP = 0.34
MAX_RETRIES = 3


def _taxid_from_elink(accessions: list[str]) -> dict[str, str]:
    """Return accession -> taxid (from nuccore elink to taxonomy). Missing/failed -> absent."""
    if not accessions:
        return {}
    id_str = ",".join(accessions)
    result: dict[str, str] = {}
    time.sleep(SLEEP)
    for attempt in range(MAX_RETRIES):
        try:
            handle = Entrez.elink(dbfrom="nuccore", db="taxonomy", id=id_str)
            data = Entrez.read(handle)
            handle.close()
            break
        except Exception as e:
            logger.warning("elink attempt %d failed: %s", attempt + 1, e)
            if attempt < MAX_RETRIES - 1:
                time.sleep(1.0 * (attempt + 1))
            else:
                return {}
    else:
        return {}
    # Entrez.read elink: with multiple ids returns list of LinkSet (one per id), order matches.
    link_sets = data if isinstance(data, list) else [data]
    for i, ls in enumerate(link_sets):
        if i >= len(accessions):
            break
        if not isinstance(ls, dict):
            continue
        for ldb in ls.get("LinkSetDb") or []:
            if ldb.get("DbTo") != "taxonomy":
                continue
            links = ldb.get("Link") or []
            if links:
                link0 = links[0]
                tid = link0.get("Id") if isinstance(link0, dict) else getattr(link0, "Id", None)
                if tid:
                    result[accessions[i]] = str(tid)
            break
    # Single LinkSet with multiple Link entries (one per id in IdList order)
    if not result and link_sets and len(link_sets) == 1:
        single = link_sets[0]
        if isinstance(single, dict):
            for ldb in single.get("LinkSetDb") or []:
                if ldb.get("DbTo") == "taxonomy":
                    for j, link in enumerate(ldb.get("Link") or []):
                        if j >= len(accessions):
                            break
                        tid = link.get("Id") if isinstance(link, dict) else getattr(link, "Id", None)
                        if tid:
                            result[accessions[j]] = str(tid)
                    break
    return result


def _lineage_from_efetch(taxids: list[str], level: str) -> dict[str, str]:
    """Fetch taxonomy records and return taxid -> group name (family or realm)."""
    if not taxids:
        return {}
    id_str = ",".join(taxids)
    result: dict[str, str] = {}
    time.sleep(SLEEP)
    for attempt in range(MAX_RETRIES):
        try:
            handle = Entrez.efetch(db="taxonomy", id=id_str, retmode="xml")
            # Parse XML: each TaxaSet has Taxonomys; each has LineageEx with Rank/ScientificName
            import xml.etree.ElementTree as ET
            tree = ET.parse(handle)
            handle.close()
            root = tree.getroot()
            break
        except Exception as e:
            logger.warning("efetch taxonomy attempt %d failed: %s", attempt + 1, e)
            if attempt < MAX_RETRIES - 1:
                time.sleep(1.0 * (attempt + 1))
            else:
                return {}
    else:
        return {}
    # Taxonomy XML: <TaxaSet><Taxon><TaxId>10310</TaxId><LineageEx><Taxon><Rank>...</Rank><ScientificName>...
    # Handle optional XML namespace
    def local_tag(tag: str) -> str:
        return tag.split("}")[-1] if "}" in tag else tag

    level_lower = level.lower()
    for taxon in root.iter():
        if local_tag(taxon.tag) != "Taxon":
            continue
        # ElementTree elements with no children are *falsy*, so the old
        # ``a or b`` pattern silently preferred the namespaced find whenever the
        # plain one matched a leaf element (TaxId/Rank/ScientificName). Use
        # explicit ``is None`` chaining so we always take the first real match.
        taxid_el = taxon.find("TaxId")
        if taxid_el is None:
            taxid_el = taxon.find("{*}TaxId")
        taxid = taxid_el.text if taxid_el is not None and taxid_el.text else None
        if not taxid:
            continue
        # Only process top-level Taxon (has LineageEx); skip nested Taxon inside LineageEx
        lineage = taxon.find("LineageEx")
        if lineage is None:
            lineage = taxon.find("{*}LineageEx")
        if lineage is None:
            continue
        name = None
        fallback = None
        lex_children = lineage.findall("Taxon")
        if not lex_children:
            lex_children = lineage.findall("{*}Taxon")
        for lex in lex_children:
            rank_el = lex.find("Rank")
            if rank_el is None:
                rank_el = lex.find("{*}Rank")
            name_el = lex.find("ScientificName")
            if name_el is None:
                name_el = lex.find("{*}ScientificName")
            if rank_el is None or name_el is None or not rank_el.text or not name_el.text:
                continue
            r = rank_el.text.lower()
            if r == level_lower:
                name = name_el.text
                break
            if fallback is None and r in ("realm", "kingdom", "order", "family", "subfamily"):
                fallback = name_el.text
        result[taxid] = name or fallback or "unknown"
    return result


UNKNOWN_TAXONOMY_WARN_FRACTION = 0.5


def fetch_viral_taxonomy_groups(
    viral_accessions: list[str],
    level: str = "family",
    *,
    batch_size: int = BATCH_SIZE,
    progress_callback=None,
) -> list[str]:
    """Return list of taxonomy group names (family or realm) for each viral accession in order.

    Uses elink nuccore->taxonomy then efetch taxonomy. Missing/failed -> "unknown".

    Emits a WARNING when a large fraction of accessions (>= 50%) resolve to
    ``"unknown"``. This commonly indicates an accession-format mismatch
    (e.g. versioned vs. unversioned IDs between the snapshot and the viral
    taxonomy JSON) or an NCBI connectivity issue, and would otherwise silently
    funnel every sequence into a single ``unknown`` bucket — defeating
    ``--balance-viral-by-taxonomy``.
    """
    groups: list[str] = []
    total_batches = (len(viral_accessions) + batch_size - 1) // batch_size
    for start in range(0, len(viral_accessions), batch_size):
        batch = viral_accessions[start : start + batch_size]
        batch_idx = start // batch_size
        acc_to_taxid = _taxid_from_elink(batch)
        if not acc_to_taxid:
            groups.extend(["unknown"] * len(batch))
            if progress_callback:
                progress_callback(batch_idx + 1, total_batches, len(groups))
            continue
        taxids = [acc_to_taxid.get(acc, "") for acc in batch]
        taxids = [t for t in taxids if t]
        if not taxids:
            groups.extend(["unknown"] * len(batch))
            if progress_callback:
                progress_callback(batch_idx + 1, total_batches, len(groups))
            continue
        taxid_to_group = _lineage_from_efetch(taxids, level)
        for acc in batch:
            tid = acc_to_taxid.get(acc, "")
            groups.append(taxid_to_group.get(tid, "unknown"))
        if progress_callback:
            progress_callback(batch_idx + 1, total_batches, len(groups))

    total = len(groups)
    if total > 0:
        unknown_count = sum(1 for g in groups if g == "unknown")
        unknown_frac = unknown_count / total
        if unknown_frac >= UNKNOWN_TAXONOMY_WARN_FRACTION:
            logger.warning(
                "Viral taxonomy lookup returned 'unknown' for %d/%d accessions (%.1f%%) at level=%s. "
                "This typically indicates an accession format mismatch (versioned vs. unversioned IDs, "
                "or snapshot accessions not present in NCBI) or a transient NCBI failure. "
                "Balancing by taxonomy will be ineffective if most accessions share the 'unknown' bucket.",
                unknown_count, total, 100.0 * unknown_frac, level,
            )
    return groups


def build_prefix_to_taxonomy(
    viral_accessions: list[str],
    level: str = "family",
    *,
    batch_size: int = BATCH_SIZE,
    progress_callback=None,
) -> dict[str, str]:
    """Return dict accession -> taxonomy group (e.g. NC_001234.1 -> Herpesviridae)."""
    groups = fetch_viral_taxonomy_groups(
        viral_accessions,
        level=level,
        batch_size=batch_size,
        progress_callback=progress_callback,
    )
    return dict(zip(viral_accessions, groups))


def run_viral_taxonomy(
    accessions_file: Path,
    output_path: Path,
    *,
    level: str = "family",
    batch_size: int = BATCH_SIZE,
) -> dict[str, str]:
    """Load viral list from accessions file, fetch taxonomy, write accession->group JSON. Returns the mapping."""
    data = load_accessions(accessions_file)
    _, viral_ids, _, _ = get_accession_lists_from_data(data)
    if not viral_ids:
        logger.warning("No viral accessions in %s", accessions_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with output_path.open("w") as f:
            json.dump({}, f)
        return {}
    total = len(viral_ids)
    total_batches = (total + batch_size - 1) // batch_size

    def progress(batch_num: int, total_batches: int, fetched: int) -> None:
        pct = 100 * batch_num / total_batches if total_batches else 0
        print(f"  Viral taxonomy: batch {batch_num}/{total_batches} ({pct:.0f}%) | {fetched}/{total} accessions", flush=True)

    mapping = build_prefix_to_taxonomy(
        viral_ids,
        level=level,
        batch_size=batch_size,
        progress_callback=progress,
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w") as f:
        json.dump(mapping, f, indent=2)
    n_unknown = sum(1 for v in mapping.values() if v == "unknown")
    logger.info("Viral taxonomy: wrote %s (%d viral, %d unknown)", output_path, len(mapping), n_unknown)
    print(f"Wrote {output_path} ({len(mapping)} viral accessions, {n_unknown} unknown)")
    return mapping


def load_viral_taxonomy(path: Path) -> dict[str, str]:
    """Load accession -> taxonomy group from JSON (e.g. NC_001234.1 -> Herpesviridae)."""
    with path.open() as f:
        return json.load(f)
