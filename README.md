# Metagenome Generator

**Simulated metagenome FASTA datasets from NCBI RefSeq for training and evaluating sequence classifiers (e.g. viral vs. prokaryotic, phage vs. bacteria).**

---

## Table of contents

- [What it does](#what-it-does)
- [When to use it](#when-to-use-it)
- [Use cases at a glance](#use-cases-at-a-glance)
- [Installation](#installation)
- [Requirements](#requirements)
- [Quick start](#quick-start)
- [Extended usage](#extended-usage)
- [Command reference](#command-reference)
- [Project structure](#project-structure)
- [Notes](#notes)

---

## What it does

Metagenome Generator downloads bacterial, viral, archaeal, and plasmid genomes from NCBI, splits them into fixed- or variable-length reads, and writes a single metagenome FASTA. It supports:

- **Reproducible runs** — snapshot accession lists to JSON and re-download or re-chunk the same set anytime.
- **Rigorous train/test evaluation** — temporal split by NCBI submission date, or percentage split with BLAST-based similarity filtering so test reads similar to train are removed.
- **EVE handling** — BLAST non-viral vs viral to exclude (and optionally export) endogenous viral element regions.
- **Extra options** — mutation simulation, user-provided viral FASTA, genome completeness filter, abundance models, taxonomy-aware viral balancing.

---

## When to use it

- You need **synthetic training data** for viral/prokaryotic or phage/bacteria classifiers.
- You want **reproducible** train/test sets (fixed accession lists, optional temporal split).
- You want to **avoid data leakage**: remove test reads that are highly similar to train (e.g. different strains).
- You need **EVE-aware** non-viral genomes (exclude or export provirus regions).

---

## Use cases at a glance

| Use case | What you want | Command or flow |
|----------|----------------|------------------|
| **Single metagenome** | One FASTA for training or benchmarking | `pipeline --num-organisms N --output-dir out --output metagenome.fasta --sequence-length 250 --reads-per-organism 1000` |
| **Reproducible run** | Same genomes every time | `snapshot` → save JSON; then `download --accessions-file <json>` (and chunk) or use that file in `pipeline` |
| **Temporal train/test** | Train on “old” genomes, test on “new” (e.g. for generalization) | `temporal-split-info` → `temporal-split` → build train and test metagenomes → **`filter-test-against-train`** (important: removes test reads similar to train) |
| **Single metagenome + train/test** | One dataset split 80/20 with similarity filter | `chunk` or `pipeline` with `--train-test-split 80` (similarity filter applied automatically) |

For all runs, use a dedicated output directory (e.g. `working_directory/`).

---

## Installation

**From PyPI** (when published):

```bash
pip install metagenome-generator
```

**From source:**

```bash
cd MetagenomeGenerator
pip install -e .
```

**With BLAST+** (needed for EVE removal and similarity filtering):

```bash
conda env create -f environment.yml
conda activate metagenome-simulator
pip install -e .
```

---

## Requirements

- **Python** 3.8+
- **Biopython** ≥ 1.83
- **BLAST+** (optional; for EVE removal and train/test similarity filtering)

Set NCBI Entrez credentials (required for download and snapshot):

```bash
export ENTREZ_EMAIL="your_email@example.com"
export ENTREZ_API_KEY="your_ncbi_api_key"   # optional, for higher rate limits
```

---

## Quick start

**Get one metagenome in one command:**

```bash
metagenome-generator pipeline \
  --num-organisms 10 \
  --output-dir output \
  --output metagenome.fasta \
  --sequence-length 250 \
  --reads-per-organism 1000
```

Result: genomes in `output/downloaded/`, metagenome FASTA in `output/metagenome/metagenome.fasta`.

**Do it in two steps (download, then chunk):**

```bash
metagenome-generator download --num-organisms 10 --output-dir output
metagenome-generator chunk \
  --input output/downloaded \
  --output metagenome.fasta \
  --output-dir output/metagenome \
  --sequence-length 250 \
  --reads-per-organism 1000
```

**Equal reads per genome (balanced):** add `--balanced` to `chunk`.

---

## Extended usage

Sections below describe each feature and how to use it.

---

### Download genomes

Obtain RefSeq genomes by category (bacteria, virus, archaea, plasmid). Use `--accessions-file` with a snapshot JSON for reproducible runs.

```bash
metagenome-generator download --num-organisms 10 --output-dir output
```

| Option | Use |
|--------|-----|
| `--num-organisms` | Number of bacterial and viral genomes. |
| `--num-archaea`, `--num-plasmid` | Extra archaeal/plasmid genomes. |
| `--accessions-file` | Load accession IDs from JSON (no NCBI search); use with snapshot for reproducibility. |
| `--save-accessions` | Save the chosen accession list to JSON for later `--accessions-file` runs. |
| `--complete-only` | Restrict to complete genomes (exclude WGS/draft) when searching; for reproducibility, create a snapshot with `snapshot --complete-only` and use that as `--accessions-file`. |

---

### Chunk genomes into reads

Turn genome FASTAs into one metagenome FASTA. Input is typically the download directory (`bacteria/`, `virus/`, `archaea/`, `plasmid/`).

```bash
metagenome-generator chunk \
  --input output/downloaded \
  --output metagenome.fasta \
  --output-dir output/metagenome \
  --sequence-length 250 \
  --reads-per-organism 1000
```

| Option | Use |
|--------|-----|
| `--sequence-length` | **Fixed read length (nt).** Each chunk is exactly this many nucleotides. Typical values: 250–500 for short-read style; match your classifier’s expected input (e.g. 250 for some viral classifiers). Required unless you use variable-length mode. |
| `--reads-per-organism` | **Max reads per genome.** Upper limit on how many non-overlapping (or sampled) reads are taken from each genome file. Omit to use all possible reads from every genome (can produce very large outputs). Use a fixed value (e.g. 1000) for controlled dataset size and balance across genomes when combined with other options. |
| `--balanced` | **Same number of reads per genome.** Each genome contributes the same count of reads (the minimum across all genomes). Use when you want to avoid one category (e.g. bacteria) dominating simply because genomes are longer; good for training classifiers with even representation per genome. |
| `--cap-total-reads` | **Cap total reads.** Downsample the whole metagenome to at most N reads. Use to match a target size (e.g. cap to the size of your positive set) or to keep evaluation sets manageable. Applied after per-genome limits and balancing. |
| `--min-contig-length`, `--max-contig-length` | **Variable-length contigs.** Instead of fixed-length reads, sample contigs with lengths uniformly between these two values (nt). Use for long-read or contig-level benchmarks (e.g. 300–2000 bp). Omit both to use fixed `--sequence-length`. |
| `--seed` | **Random seed.** Fixes randomness for variable-length sampling, cap, mutation, and train/test split. Use the same seed to reproduce the exact same metagenome; change the seed to get a different sample. |
| `--eve-intervals` | **EVE exclusion.** Path to `eve_intervals.json` produced by `blastn-filter`. Chunks that overlap these endogenous viral element (EVE) intervals on non-viral genomes are excluded from the metagenome. Use to avoid bacterial/archaeal regions that look viral and would confound viral vs. non-viral classifiers. |
| `--forbid-ambiguous` | **Exclude ambiguous bases.** Discard any read that contains non-ACGT characters (e.g. N, R, Y). Use when your pipeline or classifier assumes strict ACGT-only sequence, or to simulate cleaner sequencing. |
| `--substitution-rate`, `--indel-rate` | **Mutation simulation.** Introduce substitutions and/or indels at the given per-base rate (0–1). Use to test classifier robustness to sequencing error or divergence (e.g. 0.01 for 1% substitution rate). Combine with `--seed` for reproducible mutated datasets. |
| `--extra-viral-fasta` | **Merge user viral sequences.** Path to a FASTA of additional viral sequences (e.g. metavirome contigs, custom viral set). They are chunked like RefSeq viral genomes and merged into the viral pool. Use to combine public RefSeq viral data with your own viral contigs in one metagenome. |
| `--abundance-profile` | **Per-category read weights.** Comma-separated `category=weight`, e.g. `bacterial=0.5,viral=2,archaea=1,plasmid=1`. Scales how many reads are taken from each category relative to the base limit. Use to simulate uneven community composition (e.g. more viral, less bacterial) without changing genome lists. |
| `--abundance-distribution` | **Per-genome abundance model.** Set to `exponential` to assign each genome a weight from an exponential distribution (then normalized). Produces a few “abundant” and many “rare” genomes, similar to real communities. Use `--seed` for reproducibility. |
| `--viral-taxonomy`, `--balance-viral-by-taxonomy` | **Taxonomy-aware viral balancing.** `--viral-taxonomy` is the path to the JSON from the `viral-taxonomy` command (viral prefix → taxonomy group). With `--balance-viral-by-taxonomy`, viral read limits are set so each taxonomy group (e.g. family) contributes equally. Use to avoid a few viral families dominating and to better train on under-represented groups. |
| `--filter-similar` | **Within-metagenome similarity filter.** Remove any read that is ≥90% similar (identity and coverage) to a read already kept. The tool oversamples and refills to try to reach the target count. Use to reduce near-duplicate sequences in a single metagenome. |
| `--train-test-split` | **Train/test split with similarity filter.** Percentage of reads for training (e.g. 80). Outputs `*_train.fasta` and `*_test.fasta`. Any test read that is ≥ similarity threshold (default 90% identity over 80% length) to a train read is removed. Use for quick evaluation from one metagenome while avoiding inflated metrics from near-duplicate train/test pairs. |

---

### Pipeline (download + chunk)

One command to download and chunk; optionally run BLASTN (EVE) and Seeker. Layout: `output-dir/downloaded/`, `metagenome/`, `blastn/`, `seeker/`, `logs/`.

```bash
metagenome-generator pipeline \
  --num-organisms 10 \
  --output-dir output \
  --output metagenome.fasta \
  --sequence-length 250 \
  --reads-per-organism 1000
```

Pipeline accepts the same chunk options (e.g. `--train-test-split`, `--balanced`, `--eve-intervals`) plus `--run-blastn-filter`, `--run-seeker`, `--accessions-file`, `--complete-only`. See `metagenome-generator pipeline --help`.

---

### Accession snapshot (reproducible runs)

Save a full catalog of matching NCBI accessions to a date-stamped JSON. No sequences are downloaded. Use the JSON as `--accessions-file` for reproducible download/chunk.

```bash
metagenome-generator snapshot
```

By default, each category stores `create_date` and `title` per accession (for temporal split and auditing). Use `--no-metadata` for lists only. Use `--complete-only` to restrict to complete genomes; then use that snapshot with `--accessions-file` for reproducible complete-only runs. Legacy snapshots: `metagenome-generator migrate-snapshot <path>`.

---

### Train/test split and similarity filtering

Two workflows:

- **Temporal split** — Split accessions by NCBI submission date; build train and test metagenomes separately; then run **`filter-test-against-train`** to remove test reads that are highly similar to train. Use when you want “train on past, test on future” (e.g. generalization to novel viruses).
- **Percentage split** — Build one metagenome and split reads (e.g. 80% train, 20% test); the tool automatically removes from test any read ≥ threshold similar to a train read. Use for quick train/test from a single dataset.

**Why similarity filtering matters:** Different strains or closely related genomes can appear in both train and test and inflate metrics. Removing test reads similar to train is standard practice (e.g. DeepVirFinder, VirFinder: temporal holdout + avoiding near-duplicates). The tool supports both workflows and applies or offers this filter in each.

---

#### Temporal split (by NCBI CreateDate)

1. **Preview** (no files written):

   ```bash
   metagenome-generator temporal-split-info \
     --accessions-file snapshots/accession_snapshot_2026-03-10.json \
     --split-date 2019-06-01
   ```

2. **Write train/test JSONs:**

   ```bash
   metagenome-generator temporal-split \
     --accessions-file snapshots/accession_snapshot_2026-03-10.json \
     --split-date 2019-06-01
   ```

3. **Build train and test metagenomes:** Run `download` (or `pipeline`) twice with `--accessions-file train_<basename>.json` and `--accessions-file test_<basename>.json` into separate dirs; chunk each to get `train_metagenome.fasta` and `test_metagenome.fasta`.

4. **Filter test against train (important for rigorous evaluation):** Remove test reads that are ≥ threshold similar to train so different strains or near-duplicates are not counted as “novel” test.

   ```bash
   metagenome-generator filter-test-against-train \
     --train-fasta output_train/metagenome/train_metagenome.fasta \
     --test-fasta output_test/metagenome/test_metagenome.fasta \
     --output output_test/metagenome/test_metagenome_filtered.fasta \
     --similarity-threshold 90
   ```

   Use `test_metagenome_filtered.fasta` as the test set. Options: `--min-coverage` (default 0.8), `--threads`, `--batch-size`. Requires BLAST+.

---

#### Percentage split with similarity check (single metagenome)

Build one metagenome; the tool splits reads and removes from test any read ≥ similarity threshold (default 90% identity over 80% length) to train.

```bash
metagenome-generator chunk --input output/downloaded --output metagenome.fasta --output-dir output/metagenome \
  --sequence-length 250 --reads-per-organism 1000 \
  --train-test-split 80 --seed 42
```

Output: `metagenome_train.fasta` and `metagenome_test.fasta`. Options: `--train-test-similarity-threshold`, `--train-test-blast-threads`, `--train-test-blast-batch-size`. Same behavior when using `pipeline` with `--train-test-split 80`.

---

### BLASTN filtering (EVE removal)

EVEs in non-viral genomes can be misclassified as viral. BLAST non-viral vs viral; exclude chunks overlapping hits when building the metagenome.

**Standalone:**

```bash
metagenome-generator blastn-filter --genome-dir output/downloaded --out-dir output/blastn --evalue 1e-5 --perc-identity 70 \
  --export-eve-fasta output/blastn/eve_intervals.fasta --export-eve-min-length 200
metagenome-generator chunk --input output/downloaded --output metagenome.fasta --output-dir output/metagenome \
  --balanced --eve-intervals output/blastn/eve_intervals.json
```

**In pipeline:** add `--run-blastn-filter`; optional `--blastn-evalue`, `--blastn-perc-identity`, `--blastn-export-eve-fasta`, `--blastn-export-eve-min-length`. Requires BLAST+.

---

### Viral taxonomy (taxonomy-aware balancing)

Fetch viral taxonomy from NCBI and write a JSON mapping viral prefixes to group (e.g. family). Use with `chunk` or `pipeline` and `--balance-viral-by-taxonomy` so each viral taxonomy group contributes equally.

```bash
metagenome-generator viral-taxonomy \
  --accessions-file snapshots/accession_snapshot_2026-03-10.json \
  --output output/viral_taxonomy.json \
  --level family
```

Then chunk with balancing:

```bash
metagenome-generator chunk --input output/downloaded --output metagenome.fasta --output-dir output/metagenome \
  --sequence-length 250 --reads-per-organism 1000 \
  --viral-taxonomy output/viral_taxonomy.json --balance-viral-by-taxonomy
```

---

### Seeker (phage/bacteria prediction)

Run [Seeker](https://github.com/gussow/seeker) on the generated metagenome. Seeker must be installed in a separate conda env.

```bash
metagenome-generator seeker --input output/metagenome/metagenome.fasta --output-dir output/seeker --conda-env seeker
```

Or add `--run-seeker` to the pipeline.

---

## Command reference

| Command | Purpose |
|---------|--------|
| `download` | Download genomes from NCBI into category folders. |
| `chunk` | Turn genome FASTAs into one metagenome FASTA. |
| `pipeline` | Download + chunk (+ optional BLASTN, Seeker). |
| `snapshot` | Save full accession catalog to JSON (no downloads). |
| `temporal-split-info` | Show train/test counts for a split date (no files written). |
| `temporal-split` | Write train and test accession JSONs by CreateDate. |
| `filter-test-against-train` | Remove from test FASTA reads similar to train (BLAST). **Use after temporal split** for rigorous evaluation. |
| `migrate-snapshot` | Convert legacy snapshot to per-category metadata format. |
| `blastn-filter` | BLAST non-viral vs viral; EVE intervals for chunk and optional provirus/EVE FASTA. |
| `viral-taxonomy` | Fetch viral taxonomy; write viral_1→group JSON for `--balance-viral-by-taxonomy`. |
| `seeker` | Run Seeker on a metagenome FASTA. |

Full options: `metagenome-generator <command> --help`.

---

## Capabilities summary

| Capability | Description |
|------------|-------------|
| Download by category | RefSeq bacteria, virus, archaea, plasmid from NCBI. |
| Chunk into reads | Fixed or variable length; balanced or weighted. |
| Reproducible datasets | Snapshot accessions to JSON; re-download/re-chunk same set. |
| Temporal train/test | Split by CreateDate; **filter-test-against-train** for similarity. |
| EVE removal / export | BLAST non-viral vs viral; exclude/export provirus regions. |
| Similarity filtering | Drop test reads similar to train (both temporal and percentage split). |
| Mutation simulation | Substitution/indel rates for robustness tests. |
| Extra viral FASTA | Merge user viral sequences with RefSeq viral. |
| Genome completeness | `--complete-only` in snapshot and download. |
| Abundance model | Per-category or per-genome (e.g. exponential) weights. |
| Taxonomy-aware viral balance | Balance viral reads by family/realm. |
| Ambiguous-base filter | Exclude reads with N (e.g. `--forbid-ambiguous`). |
| Seeker integration | Run Seeker from the same workflow. |

---

## Project structure

```
MetagenomeGenerator/
├── pyproject.toml
├── README.md
├── LICENSE
├── environment.yml
├── main.py
├── src/metagenome_generator/
│   ├── cli.py
│   ├── download_genomes.py
│   ├── ncbi_search.py
│   ├── accession_snapshot.py
│   ├── chunk_genomes.py
│   ├── genome_layout.py
│   ├── blastn_filter.py
│   ├── similarity_filter.py
│   ├── temporal_split.py
│   ├── viral_taxonomy.py
│   └── seeker_wrapper.py
├── snapshots/
└── working_directory/
```

Programmatic use: `from metagenome_generator import build_metagenome, download_genomes, load_accessions, validate_genome_dir`.

---

## Notes

- NCBI rate limits apply; the tool uses delays and retries (up to 3 with backoff).
- Genome selection uses RefSeq and length filters; see `DEFAULT_QUERIES` in `ncbi_search.py` to change criteria.
- Prefer a dedicated working directory for runs (e.g. `working_directory/`).
