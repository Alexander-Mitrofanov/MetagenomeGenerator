# Experiment: stratify viral accessions by inferred host domain (euk vs prok)

## Goal

Produce **two CHIMERA-compatible accession JSON files** from an existing snapshot:

1. **Prokaryote-infecting viruses** (bacterial + archaeal viruses / “phage-heavy” lineages).
2. **Eukaryote-infecting viruses** (animal, plant, fungal, protist viruses).

Non-viral categories (`bacterial`, `archaea`, `plasmid`) are **copied unchanged** so `download` / `genome-pool` behavior stays the same; only the `viral` list is filtered.

## Why this is experimental

- **ICTV / community best practice**: host-based bins are **operational**, not phylogeny ([ICTV four principles](https://ictv.global/news/four_principles)). Many taxa are unambiguous; edge cases go to `unknown`.
- We use **NCBI taxonomy lineage** (nuccore → taxid → `LineageEx`) plus a **frozen substring allowlist** (`heuristic_markers.json`). This is **not** curated per-isolate host metadata (BioSample / lab host).
- Outputs include a **TSV audit trail** (`viral_host_classification.tsv`) for manual QC.

## Implementation steps (done in repo)

1. **`heuristic_markers.json`** — Ordered markers: if lineage contains a **prok** marker first (substring match), label `prok`; else if any **euk** marker, label `euk`; else `unknown`.
2. **`build_host_stratified_snapshots.py`** — Loads snapshot → optional `--max-viral` subsample → batched `elink` + `efetch` → classify → writes:
   - `viral_host_classification.tsv`
   - `accession_snapshot_viral_prok.json`
   - `accession_snapshot_viral_euk.json`
3. **`README.md`** — How to run (no changes to main package CLI / PyPI description).

## How you use the outputs

```bash
# Example: download only prok-infecting viral slice (same bacteria/archaea/plasmid lists)
metagenome-generator download \
  --accessions-file experiments/viral_host_stratify/out/accession_snapshot_viral_prok.json \
  --max-bacteria 50 --max-virus 200 --sample-seed 1 \
  --output-dir output/prok_virus_slice
```

Repeat with `accession_snapshot_viral_euk.json` for the euk slice.

## Limits

- **Unknown** fraction should be reported; tune markers after inspecting TSV.
- Full viral list (~8k for 2026-03-10 snapshot) needs many Entrez calls; respect `ENTREZ_EMAIL` / `ENTREZ_API_KEY`.
