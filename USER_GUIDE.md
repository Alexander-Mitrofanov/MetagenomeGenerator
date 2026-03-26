# CHIMERA User Guide (Detailed)

CHIMERA (Configurable Hybrid In-silico Metagenome Emulator for Read Analysis) is a dataset generation tool for building synthetic metagenomes from curated genomes. It is designed for users who need reproducible train/test datasets for sequence classification, benchmarking, and method comparison.

This document is intentionally detailed. It explains not only "which command to run," but also:

- what each workflow is for,
- why each option matters,
- how to choose parameter values,
- and how to keep results reproducible across time, systems, and collaborators.

---

## Contents

1. What CHIMERA is for
2. Mental model of the pipeline
3. Installation and prerequisites
4. Reproducibility model (snapshot + viral DB manifest)
5. Quick-start workflow
6. Detailed workflows
   - A. Single dataset with `pipeline`
   - B. Reproducible snapshot-driven generation
   - C. Temporal train/test generation
   - D. Structured benchmark generation
   - E. Splitting existing metagenomes
7. EVE/prophage filtering and viral DB pinning
8. Header labels and output semantics
9. Output folder layout reference
10. Performance tuning
11. Troubleshooting and common mistakes
12. Recommended release workflow for teams

---

## 1) What CHIMERA Is For

CHIMERA is built for users who need controlled synthetic sequence data where generation choices are explicit, traceable, and repeatable.

Typical use cases:

- **Train classifiers** (e.g., viral vs non-viral, phage vs prokaryotic) on synthetic data with known provenance.
- **Evaluate generalization** by separating train/test with either random split + similarity filtering or temporal split.
- **Run benchmark suites** with fixed counts per class and multi-replicate runs.
- **Reduce confounders** such as endogenous viral elements (EVEs) in non-viral genomes by BLAST-based interval exclusion.
- **Use public and in-house data** (RefSeq snapshots and/or your own FASTA collections).

In short: CHIMERA is not just a sequence generator; it is a reproducible benchmarking framework.

---

## 2) Mental Model of the Pipeline

CHIMERA workflows are composed of four conceptual stages:

1. **Genome selection**
   - From NCBI query or from a pinned snapshot.
2. **Genome materialization**
   - Download category FASTAs (`bacteria`, `virus`, `archaea`, `plasmid`).
3. **Optional contamination control**
   - BLAST non-viral sequences against a viral DB to detect EVE/prophage-like regions.
4. **Read generation and split logic**
   - Chunk genomes into reads/contigs.
   - Optionally split into train/test and remove test reads too similar to train.

You can run all stages at once (`pipeline`, `temporal-pipeline`) or run modular commands separately (`download`, `blastn-filter`, `chunk`, `filter-test-against-train`).

---

## 3) Installation and Prerequisites

## 3.1 Install CHIMERA

From repository root:

```bash
pip install -e .
```

Verify:

```bash
metagenome-generator --help
```

## 3.2 Required tools

- Python 3.10+
- BLAST+ tools:
  - `makeblastdb`
  - `blastn`

No additional optional tool is required for the current CHIMERA pipeline workflows.

## 3.3 Optional NCBI config

If you download from NCBI, set:

```bash
export NCBI_EMAIL="you@example.com"
export NCBI_API_KEY="..."
```

This improves rate behavior and reliability.

---

## 4) Reproducibility Model (Important)

CHIMERA reproducibility has two pillars:

1. **Accession snapshot** (which genomes are eligible)
2. **Viral DB identity** (which BLAST reference defines EVE filtering)

If either changes, your final datasets can change.

## 4.1 Snapshot pinning

Use a fixed snapshot JSON (e.g. `accession_snapshot_2026-03-10.json`) for reproducible accession universe.

## 4.2 Viral DB pinning

Use a fixed viral DB plus manifest:

- DB prefix:
  - `.../blastn_db/viral_db`
- Manifest:
  - `.../viral_db_manifest.json`
- Optional strict pin:
  - `--require-viral-db-sha256 <aggregate_sha256>`

This prevents silent drift as viral reference assets evolve.

## 4.3 Download cache behavior

`download` now defaults to `downloads/` and reuses existing accession FASTAs.
If `<accession>.fasta` is already present and non-empty, CHIMERA skips re-download.

This improves speed and keeps repeated runs stable.

---

## 5) Quick Start (Practical First Run)

Goal: create one small train/test dataset quickly using a pinned snapshot.

```bash
metagenome-generator pipeline \
  --accessions-file snapshots/accession_snapshot_2026-03-10.json \
  --max-bacteria 5 --max-virus 5 --max-archaea 2 --max-plasmid 2 \
  --sample-seed 1 \
  --output-dir output/quickstart \
  --output metagenome.fasta \
  --sequence-length 250 \
  --reads-per-organism 200 \
  --train-test-split 80 \
  --balanced \
  --seed 1
```

What this does:

- selects a small deterministic subset from the snapshot,
- downloads genomes,
- chunks into 250-nt reads,
- enforces equal per-genome read contribution (`--balanced`),
- splits train/test at 80/20,
- removes test reads too similar to train.

Expected outputs:

- `output/quickstart/metagenome_train.fasta`
- `output/quickstart/metagenome_test.fasta`
- `output/quickstart/downloaded/`
- `output/quickstart/logs/pipeline.log`

---

## 6) Detailed Workflows

## 6A) Single dataset generation with `pipeline`

Use when:

- you want one standalone dataset,
- you do not need temporal split,
- and you prefer one command.

### Example

```bash
metagenome-generator pipeline \
  --accessions-file snapshots/accession_snapshot_2026-03-10.json \
  --max-bacteria 100 --max-virus 100 --max-archaea 20 --max-plasmid 20 \
  --sample-seed 42 \
  --output-dir output/run1 \
  --output metagenome.fasta \
  --sequence-length 250 \
  --reads-per-organism 1000 \
  --balanced \
  --train-test-split 80 \
  --seed 42
```

### Why these options exist

- `--max-*` controls class complexity and runtime.
- `--sample-seed` makes subset choice deterministic.
- `--reads-per-organism` controls dataset size.
- `--balanced` avoids long-genome dominance.
- `--train-test-split 80` creates model-ready split directly.
- `--seed` controls chunking randomness and split shuffle.

### Typical modifications

- Increase `--sequence-length` to 1000 for long-fragment tasks.
- Remove `--train-test-split` if you want one combined metagenome.
- Add `--run-blastn-filter` + pinned viral DB for EVE exclusion.

---

## 6B) Reproducible snapshot-driven modular workflow

Use when:

- you want separate control over download, filtering, chunking,
- or you reuse one genome set across many read-generation runs.

### Step 1: Create or pick snapshot

Create:

```bash
metagenome-generator snapshot
```

Or use release-provided snapshot asset.

### Step 2: Download from snapshot

```bash
metagenome-generator download \
  --accessions-file snapshots/accession_snapshot_2026-03-10.json \
  --max-bacteria 200 --max-virus 200 --max-archaea 50 --max-plasmid 50 \
  --sample-seed 7 \
  --output-dir downloads/my_set
```

What this does:

- pins accession universe to snapshot,
- selects deterministic subset,
- caches FASTAs under `downloads/my_set`.

### Step 3: Optional EVE detection

```bash
metagenome-generator blastn-filter \
  --genome-dir downloads/my_set \
  --out-dir output/my_set/blastn \
  --viral-db /path/to/viral_db_YYYY-MM-DD/blastn_db/viral_db \
  --viral-db-manifest /path/to/viral_db_YYYY-MM-DD/viral_db_manifest.json
```

### Step 4: Chunk and generate final metagenome

```bash
metagenome-generator chunk \
  --input downloads/my_set \
  --output my_metagenome.fasta \
  --output-dir output/my_set \
  --sequence-length 250 \
  --reads-per-organism 1000 \
  --balanced \
  --eve-intervals output/my_set/blastn/eve_intervals.json \
  --seed 7
```

This modular approach is best for repeated experimentation.

---

## 6C) Temporal train/test generation

Use when:

- you want train on older genomes and test on newer genomes,
- and you care about realistic generalization over time.

### Step 1: Find viable split date

```bash
metagenome-generator temporal-split-search \
  --accessions-file snapshots/accession_snapshot_2026-03-10.json \
  --min-train 400 --min-test 100
```

This helps avoid unusable splits (e.g. too few viruses in test).

### Step 2: Run one-shot temporal pipeline

```bash
metagenome-generator temporal-pipeline \
  --accessions-file snapshots/accession_snapshot_2026-03-10.json \
  --split-date 2025-01-01 \
  --output-dir working_directory/temporal_run \
  --max-bacteria-train 100 --max-virus-train 100 --max-archaea-train 100 --max-plasmid-train 100 \
  --max-bacteria-test 25 --max-virus-test 25 --max-archaea-test 25 --max-plasmid-test 25 \
  --sequence-length 1000 \
  --reads-per-organism 30 \
  --sample-seed 42 --train-seed 42 --test-seed 43 \
  --viral-db /path/to/viral_db_YYYY-MM-DD/blastn_db/viral_db \
  --viral-db-manifest /path/to/viral_db_YYYY-MM-DD/viral_db_manifest.json
```

What this command orchestrates:

1. Temporal split on snapshot metadata.
2. Download train and test genomes.
3. Optional EVE filtering on both sets.
4. Chunk train and test.
5. Remove test reads similar to train.

Final files:

- `train_metagenome.fasta`
- `test_metagenome.fasta` (filtered against train)

---

## 6D) Structured benchmark replicates (`benchmark-recipe`)

Use when:

- you need multiple runs for statistical comparison,
- and want controlled composition and diversity across replicates.

### Example

```bash
metagenome-generator benchmark-recipe \
  --accessions-file snapshots/accession_snapshot_2026-03-10.json \
  --output-dir benchmarks/run1 \
  --per-category 50 \
  --archaea 10 --plasmid 10 \
  --replicates 5 \
  --train-test-split 80 \
  --seed 42 \
  --sequence-length 250 \
  --reads-per-organism 1000
```

What it provides:

- fixed counts per category per replicate,
- diversity-aware replicate selection,
- explicit train/test files in each replicate directory.

Per replicate outputs:

- `replicate_XXX/downloaded/`
- `replicate_XXX/metagenome_train.fasta`
- `replicate_XXX/metagenome_test.fasta`

---

## 6E) Split existing metagenome files

Use when:

- you already generated a metagenome without train/test split,
- and now need a comparable train/test partition with leakage control.

```bash
metagenome-generator split-metagenome-train-test \
  --input output/existing_metagenome.fasta \
  --output-dir output/existing_split \
  --output-stem metagenome \
  --train-test-split 80 \
  --similarity-threshold 90 \
  --min-coverage 0.8 \
  --seed 42
```

---

## 6F) Worked example: temporal 100/25 benchmark set

This is the most commonly requested evaluation setup and replaces the old standalone
`Detailed_Use_Cases.md` scenario.

### Objective

Build a temporal dataset with:

- Train: 100 genomes per category (bacteria, virus, archaea, plasmid) from genomes before 2025.
- Test: 25 genomes per category from genomes in 2025 and later.
- Read length: 1000 nt.
- Minimum train read budget: at least 10,000 reads.
- Leakage control: remove test reads highly similar to train.
- Single output root folder under `working_directory/`.

### Recommended one-shot command

Use `temporal-pipeline` unless you explicitly need every intermediate command:

```bash
metagenome-generator temporal-pipeline \
  --accessions-file snapshots/accession_snapshot_2026-03-10.json \
  --split-date 2025-01-01 \
  --output-dir working_directory/temporal_100_25 \
  --max-bacteria-train 100 --max-virus-train 100 --max-archaea-train 100 --max-plasmid-train 100 \
  --max-bacteria-test 25 --max-virus-test 25 --max-archaea-test 25 --max-plasmid-test 25 \
  --sequence-length 1000 \
  --reads-per-organism 30 \
  --sample-seed 42 --train-seed 42 --test-seed 43 \
  --viral-db /path/to/viral_db_YYYY-MM-DD/blastn_db/viral_db \
  --viral-db-manifest /path/to/viral_db_YYYY-MM-DD/viral_db_manifest.json
```

### Why this satisfies the objective

- `100 + 100 + 100 + 100 = 400` train genomes.
- `reads-per-organism 30` gives up to ~12,000 reads before any exclusions, so it generally meets the `>=10,000` target.
- Temporal split at `2025-01-01` enforces “train on older / test on newer”.
- Built-in train/test similarity filtering removes near-duplicate leakage.

### Expected final outputs

- `working_directory/temporal_100_25/train_downloaded/`
- `working_directory/temporal_100_25/test_downloaded/`
- `working_directory/temporal_100_25/blastn/`
- `working_directory/temporal_100_25/train_metagenome.fasta`
- `working_directory/temporal_100_25/test_metagenome.fasta`

### Equivalent modular command chain (advanced users)

Only use this if you want manual control per step:

1. `temporal-split`
2. `download` train (max 100/category)
3. `download` test (max 25/category)
4. `blastn-filter` train and test (optional but recommended)
5. `chunk` train and test (test to temporary unfiltered FASTA)
6. `filter-test-against-train`

---

## 7) EVE/Prophage Filtering and Viral DB Pinning

EVE filtering is crucial when non-viral genomes may contain viral-like regions that can bias classifiers.

### Recommended mode

- Use release-pinned DB + manifest:

```bash
metagenome-generator blastn-filter \
  --genome-dir output/downloaded \
  --out-dir output/blastn \
  --viral-db /path/to/viral_db_YYYY-MM-DD/blastn_db/viral_db \
  --viral-db-manifest /path/to/viral_db_YYYY-MM-DD/viral_db_manifest.json
```

### Strict mode

```bash
metagenome-generator blastn-filter \
  --genome-dir output/downloaded \
  --out-dir output/blastn \
  --viral-db /path/to/viral_db_YYYY-MM-DD/blastn_db/viral_db \
  --require-viral-db-sha256 <aggregate_sha256_from_manifest>
```

Use strict mode for publication-grade reproducibility.

---

## 8) Header Labels and Output Semantics

CHIMERA writes class-prefixed read IDs in generated outputs:

- `bacteria_...`
- `virus_...`
- `archaea_...`
- `plasmid_...`

This is useful for:

- quick sanity checks,
- downstream parsing,
- class-level auditing.

`--train-test-split` values are percentages:

- correct: `80`
- incorrect: `0.8`

---

## 9) Output Folder Layout Reference

## 9.1 `pipeline`

- `output_dir/downloaded/`
- `output_dir/blastn/` (if enabled)
- `output_dir/logs/pipeline.log`
- final metagenome file(s) in `output_dir/`

## 9.2 `temporal-pipeline`

- `output_dir/train_downloaded/`
- `output_dir/test_downloaded/`
- `output_dir/blastn/`
- `output_dir/train_metagenome.fasta`
- `output_dir/test_metagenome.fasta`

## 9.3 `build-viral-db`

- `output_dir/viral_db_YYYY-MM-DD/blastn_db/viral_db*`
- `output_dir/viral_db_YYYY-MM-DD/viral_db_manifest.json`

---

## 10) Performance Tuning

Most runtime usually comes from:

1. NCBI downloads
2. BLAST-based filtering

Tips:

- Reduce `--max-*` while prototyping.
- Reuse download folders (cache skips already-downloaded FASTAs).
- Reuse viral DB from releases instead of rebuilding.
- For many seeds, reuse one downloaded set and vary chunk seed only.
- Keep BLAST filtering for final benchmark runs if runtime is tight.

---

## 11) Troubleshooting and Common Mistakes

### Missing labels in headers

Check installation path and version:

```bash
which metagenome-generator
metagenome-generator --help
```

### Viral DB not found

If `--viral-db` prefix is provided, corresponding `.nhr` etc. must exist.

### Manifest validation error

DB files and manifest do not match (changed/corrupted/wrong pair). Re-extract matching release asset.

### Slow run

Use smaller subset and confirm workflow first; then scale up.

### Train/test split unexpectedly tiny

Check that you used percent (`80`) instead of fraction (`0.8`).

---

## 12) Recommended Team Workflow with Releases

For robust collaboration:

1. Pin CHIMERA release version (e.g. `v0.7`).
2. Use the matching release snapshot asset.
3. Use the matching release viral DB asset and manifest.
4. Store all command lines and seeds in run metadata.
5. Archive outputs with pipeline logs.

This allows another user to reconstruct your dataset exactly.

---

## Appendix: Command Help

For full flag details:

```bash
metagenome-generator <command> --help
```

Most commonly consulted:

```bash
metagenome-generator pipeline --help
metagenome-generator temporal-pipeline --help
metagenome-generator blastn-filter --help
metagenome-generator benchmark-recipe --help
metagenome-generator split-metagenome-train-test --help
```

