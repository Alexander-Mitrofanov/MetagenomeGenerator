# Viral host stratification (experiment)

See **PLAN.md** for rationale and limits.

## Run

From repository root (requires `ENTREZ_EMAIL`, optional `ENTREZ_API_KEY`, and `pip install -e .` for `metagenome_generator`):

```bash
cd /path/to/MetagenomeGenerator
python experiments/viral_host_stratify/build_host_stratified_snapshots.py \
  --accessions-file snapshots/accession_snapshot_2026-03-10.json \
  --output-dir experiments/viral_host_stratify/out
```

Smoke test on a subsample (faster):

```bash
python experiments/viral_host_stratify/build_host_stratified_snapshots.py \
  --accessions-file snapshots/accession_snapshot_2026-03-10.json \
  --output-dir experiments/viral_host_stratify/out \
  --max-viral 400 \
  --sample-seed 42
```

## Outputs

| File | Purpose |
|------|--------|
| `out/viral_host_classification.tsv` | accession, taxid, `prok` / `euk` / `unknown`, truncated lineage |
| `out/accession_snapshot_viral_prok.json` | Same non-viral lists as input; `viral` = prok-inferred subset only |
| `out/accession_snapshot_viral_euk.json` | Same non-viral lists; `viral` = euk-inferred subset only |

Use either JSON with `metagenome-generator download --accessions-file ...` (optionally `--max-virus` to sample further).

## Tuning

Edit **heuristic_markers.json** after inspecting the TSV, then rerun.
