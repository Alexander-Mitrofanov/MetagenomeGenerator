#!/usr/bin/env python3
"""One-off experiment: split snapshot viral accessions into prok vs euk heuristic bins.

See PLAN.md. Does not modify CHIMERA CLI or package metadata.
"""
from __future__ import annotations

import argparse
import json
import os
import random
import sys
import time
import xml.etree.ElementTree as ET
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_REPO_ROOT / "src"))

from Bio import Entrez  # noqa: E402

from metagenome_generator.download_genomes import (  # noqa: E402
    ACCESSIONS_KEY_ARCHAEA,
    ACCESSIONS_KEY_BACTERIAL,
    ACCESSIONS_KEY_PLASMID,
    ACCESSIONS_KEY_VIRAL,
    load_accessions,
)

Entrez.email = os.environ.get("ENTREZ_EMAIL", "your_email@example.com")
Entrez.api_key = os.environ.get("ENTREZ_API_KEY")

SLEEP = 0.34
ELINK_BATCH = 25
EFETCH_TAX_BATCH = 100
MAX_RETRIES = 3


def _local_tag(tag: str) -> str:
    return tag.split("}")[-1] if "}" in tag else tag


def taxid_from_elink(accessions: list[str]) -> dict[str, str]:
    """accession -> virus taxid (nuccore elink to taxonomy)."""
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
        except Exception:
            if attempt < MAX_RETRIES - 1:
                time.sleep(1.0 * (attempt + 1))
            else:
                return {}
    else:
        return {}
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


def lineage_string_for_taxids(taxids: list[str]) -> dict[str, str]:
    """taxid -> 'name1 | name2 | ...' from LineageEx + terminal taxon name."""
    out: dict[str, str] = {}
    if not taxids:
        return out
    id_str = ",".join(taxids)
    time.sleep(SLEEP)
    for attempt in range(MAX_RETRIES):
        try:
            handle = Entrez.efetch(db="taxonomy", id=id_str, retmode="xml")
            tree = ET.parse(handle)
            handle.close()
            break
        except Exception:
            if attempt < MAX_RETRIES - 1:
                time.sleep(1.0 * (attempt + 1))
            else:
                return out
    else:
        return out
    root = tree.getroot()
    for taxon in root.iter():
        if _local_tag(taxon.tag) != "Taxon":
            continue
        taxid_el = taxon.find("TaxId")
        if taxid_el is None:
            taxid_el = taxon.find("{*}TaxId")
        if taxid_el is None or not taxid_el.text:
            continue
        tid = taxid_el.text.strip()
        names: list[str] = []
        lineage = taxon.find("LineageEx")
        if lineage is None:
            lineage = taxon.find("{*}LineageEx")
        if lineage is not None:
            for lex in lineage.findall("Taxon") or lineage.findall("{*}Taxon") or []:
                name_el = lex.find("ScientificName")
                if name_el is None:
                    name_el = lex.find("{*}ScientificName")
                if name_el is not None and name_el.text:
                    names.append(name_el.text.strip())
        sci_el = taxon.find("ScientificName")
        if sci_el is None:
            sci_el = taxon.find("{*}ScientificName")
        if sci_el is not None and sci_el.text:
            term = sci_el.text.strip()
            if not names or names[-1] != term:
                names.append(term)
        out[tid] = " | ".join(names)
    return out


def classify_lineage(lineage: str, prok_markers: list[str], euk_markers: list[str]) -> str:
    low = lineage.lower()
    for m in prok_markers:
        if m.lower() in low:
            return "prok"
    for m in euk_markers:
        if m.lower() in low:
            return "euk"
    return "unknown"


def viral_entries(data: dict) -> list[dict]:
    raw = data.get(ACCESSIONS_KEY_VIRAL, [])
    if not raw:
        return []
    if isinstance(raw[0], str):
        return [{"accession": s, "create_date": "", "title": ""} for s in raw]
    return [dict(x) for x in raw if isinstance(x, dict) and x.get("accession")]


def build_output_snapshot(
    source: dict,
    viral_entries_filtered: list,
    *,
    experiment_note: str,
) -> dict:
    out = {
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "experiment_note": experiment_note,
        "source_timestamp": source.get("timestamp"),
        ACCESSIONS_KEY_BACTERIAL: source.get(ACCESSIONS_KEY_BACTERIAL, []),
        ACCESSIONS_KEY_VIRAL: viral_entries_filtered,
        ACCESSIONS_KEY_ARCHAEA: source.get(ACCESSIONS_KEY_ARCHAEA, []),
        ACCESSIONS_KEY_PLASMID: source.get(ACCESSIONS_KEY_PLASMID, []),
    }
    if source.get("accession_metadata"):
        out["accession_metadata"] = source["accession_metadata"]
    return out


def main() -> None:
    p = argparse.ArgumentParser(description="Build prok-only and euk-only viral snapshot JSONs (experiment).")
    p.add_argument(
        "--accessions-file",
        type=Path,
        required=True,
        help="Input CHIMERA snapshot JSON (e.g. snapshots/accession_snapshot_2026-03-10.json).",
    )
    p.add_argument(
        "--markers",
        type=Path,
        default=Path(__file__).resolve().parent / "heuristic_markers.json",
        help="JSON with prokaryotic_infecting_lineage_markers and eukaryotic_infecting_lineage_markers lists.",
    )
    p.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parent / "out",
        help="Directory for TSV + two snapshot JSON files.",
    )
    p.add_argument("--max-viral", type=int, default=None, help="Subsample this many viral records (reproducible).")
    p.add_argument("--sample-seed", type=int, default=42, help="Used with --max-viral.")
    args = p.parse_args()

    markers_path = args.markers
    with markers_path.open() as f:
        markers = json.load(f)
    prok_m = markers["prokaryotic_infecting_lineage_markers"]
    euk_m = markers["eukaryotic_infecting_lineage_markers"]

    data = load_accessions(args.accessions_file)
    entries = viral_entries(data)
    if args.max_viral is not None and args.max_viral < len(entries):
        rng = random.Random(args.sample_seed)
        entries = rng.sample(entries, args.max_viral)
        print(f"Subsampled viral: {len(entries)} (seed={args.sample_seed})", flush=True)
    accessions = [e["accession"] for e in entries]

    args.output_dir.mkdir(parents=True, exist_ok=True)
    tsv_path = args.output_dir / "viral_host_classification.tsv"

    acc_to_taxid: dict[str, str] = {}
    for i in range(0, len(accessions), ELINK_BATCH):
        batch = accessions[i : i + ELINK_BATCH]
        acc_to_taxid.update(taxid_from_elink(batch))
        print(f"  elink {min(i + ELINK_BATCH, len(accessions))}/{len(accessions)}", flush=True)
    missing = [a for a in accessions if not acc_to_taxid.get(a)]
    if missing:
        print(f"  elink fallback (single id) for {len(missing)} accessions", flush=True)
        for a in missing:
            acc_to_taxid.update(taxid_from_elink([a]))

    unique_taxids = sorted({t for t in acc_to_taxid.values() if t})
    taxid_to_lineage: dict[str, str] = {}
    for i in range(0, len(unique_taxids), EFETCH_TAX_BATCH):
        batch = unique_taxids[i : i + EFETCH_TAX_BATCH]
        taxid_to_lineage.update(lineage_string_for_taxids(batch))
        print(f"  taxonomy efetch {min(i + EFETCH_TAX_BATCH, len(unique_taxids))}/{len(unique_taxids)}", flush=True)

    rows: list[tuple[str, str, str, str]] = []
    acc_to_class: dict[str, str] = {}
    for e in entries:
        acc = e["accession"]
        tid = acc_to_taxid.get(acc, "")
        lin = taxid_to_lineage.get(tid, "") if tid else ""
        cls = classify_lineage(lin, prok_m, euk_m) if lin else "unknown"
        acc_to_class[acc] = cls
        rows.append((acc, tid, cls, lin[:500]))

    with tsv_path.open("w") as f:
        f.write("accession\ttaxid\theuristic_host_bin\tlineage_truncated\n")
        for acc, tid, cls, lin in rows:
            f.write(f"{acc}\t{tid}\t{cls}\t{lin.replace(chr(9), ' ')}\n")
    print(f"Wrote {tsv_path}", flush=True)

    prok_entries = [e for e in entries if acc_to_class.get(e["accession"]) == "prok"]
    euk_entries = [e for e in entries if acc_to_class.get(e["accession"]) == "euk"]
    unk = sum(1 for e in entries if acc_to_class.get(e["accession"]) == "unknown")
    print(
        f"Counts: prok={len(prok_entries)} euk={len(euk_entries)} unknown={unk} (of {len(entries)} classified)",
        flush=True,
    )

    prok_path = args.output_dir / "accession_snapshot_viral_prok.json"
    euk_path = args.output_dir / "accession_snapshot_viral_euk.json"
    note = (
        "Heuristic viral stratification from build_host_stratified_snapshots.py; "
        "see experiments/viral_host_stratify/PLAN.md. Not ICTV-validated host labels."
    )
    with prok_path.open("w") as f:
        json.dump(
            build_output_snapshot(
                data,
                prok_entries,
                experiment_note=note + " viral=prok-inferring subset.",
            ),
            f,
            indent=2,
        )
    with euk_path.open("w") as f:
        json.dump(
            build_output_snapshot(
                data,
                euk_entries,
                experiment_note=note + " viral=euk-inferring subset.",
            ),
            f,
            indent=2,
        )
    print(f"Wrote {prok_path}", flush=True)
    print(f"Wrote {euk_path}", flush=True)


if __name__ == "__main__":
    main()
