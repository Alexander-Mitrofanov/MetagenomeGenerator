# Metagenome Generator

**Simulated metagenome datasets from NCBI RefSeq for training and evaluating sequence classifiers.**

Metagenome Generator produces FASTA training data for viral vs. prokaryotic (or phage vs. bacteria) classifiers. It downloads bacterial, viral, archaeal, and plasmid genomes from NCBI, segments them into fixed- or variable-length reads, and writes a single metagenome FASTA. The tool supports reproducible runs, temporal train/test splits, EVE removal, similarity filtering, mutation simulation, and merging of user-provided viral contigs—aligning with data-preparation practices used in tools such as DeepVirFinder and VirFinder.

---

## Table of contents

- [Capabilities](#capabilities)
- [Installation](#installation)
- [Requirements](#requirements)
- [Quick start](#quick-start)
- [Extended usage](#extended-usage)
- [Command reference](#command-reference)
- [Project structure](#project-structure)
- [Notes](#notes)

---

## Capabilities

| Capability | Purpose |
|------------|--------|
| **Download by category** | Fetch RefSeq genomes (bacteria, virus, archaea, plasmid) from NCBI Nucleotide into standard folder layout. |
| **Chunk into reads** | Convert genome FASTAs into fixed-length (e.g. 250 nt) or variable-length contigs; control reads per genome or balance across genomes. |
| **Single metagenome FASTA** | One output file with reads from all categories, suitable for training or benchmarking. |
| **Reproducible datasets** | Snapshot full accession lists (with optional CreateDate/title); re-download or re-chunk the same set anytime. |
| **Temporal train/test split** | Split by NCBI CreateDate (e.g. train on pre-2019, test on 2019+); evaluate generalization to “novel” sequences. |
| **EVE removal / export** | BLAST non-viral vs viral; exclude chunks overlapping endogenous viral elements and optionally export provirus/EVE intervals as a separate FASTA. |
| **Similarity filtering** | Drop reads ≥90% similar to already-kept; optional train/test split with test sequences removed if similar to train. |
| **Mutation simulation** | Per-base substitution and optional indel rates; test classifier robustness to sequencing error or divergence (e.g. 1% substitutions). |
| **Extra viral FASTA** | Merge user-provided viral sequences (e.g. metavirome contigs) with RefSeq viral chunks in one run. |
| **Genome quality / completeness** | Restrict to complete genomes only (exclude WGS/draft) via `--complete-only` in snapshot and download; uses NCBI `complete[Properties]` and `NOT WGS[Properties]`. |
| **Abundance / coverage model** | Per-category weights (`--abundance-profile`) or per-genome exponential weights (`--abundance-distribution exponential`); use `--seed` for reproducibility. |
| **Taxonomy-aware viral balancing** | Balance viral reads by taxonomy group (e.g. family/realm) so under-represented groups contribute equally; `viral-taxonomy` fetches NCBI taxonomy and writes prefix→group JSON. |
| **Ambiguous-base filter** | Exclude reads containing non-ACGT characters (e.g. N). |
| **Seeker integration** | Run Seeker (phage/bacteria prediction) on the generated metagenome from the same workflow. |

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

**With BLAST+ (for EVE removal and similarity filtering):**

```bash
conda env create -f environment.yml
conda activate metagenome-simulator
pip install -e .
```

---

## Requirements

- **Python** 3.8+
- **Biopython** ≥ 1.83
- **BLAST+** (optional; for EVE removal and similarity filtering)

Set NCBI Entrez credentials (required for download/snapshot):

```bash
export ENTREZ_EMAIL="your_email@example.com"
export ENTREZ_API_KEY="your_ncbi_api_key"   # optional, for higher rate limits
```

---

## Quick start

Use the `metagenome-generator` command (or `python -m metagenome_generator`). For all outputs, use a dedicated directory (e.g. `working_directory/`); see [working_directory/README.md](working_directory/README.md).

**1. Download and build one metagenome (simplest):**

```bash
metagenome-generator pipeline \
  --num-organisms 10 \
  --output-dir output \
  --output metagenome.fasta \
  --sequence-length 250 \
  --reads-per-organism 1000
```

Result: `output/downloaded/` (genomes) and `output/metagenome/metagenome.fasta` (reads).

**2. Download only:**

```bash
metagenome-generator download --num-organisms 10 --output-dir output
```

**3. Chunk existing genomes:**

```bash
metagenome-generator chunk \
  --input output/downloaded \
  --output metagenome.fasta \
  --output-dir output/metagenome \
  --sequence-length 250 \
  --reads-per-organism 1000
```

**4. Balanced reads (equal contribution per genome):**

```bash
metagenome-generator chunk \
  --input output/downloaded \
  --output metagenome.fasta \
  --output-dir output/metagenome \
  --sequence-length 250 \
  --balanced
```

---

## Extended usage

Each section below describes **what a feature is for**, then **how to use it**.

---

### Download genomes

**Purpose:** Obtain RefSeq genomes by category (bacteria, virus, archaea, plasmid) for chunking. Use `--accessions-file` for reproducible runs (same set on every run).

```bash
metagenome-generator download \
  --num-organisms 10 \
  --output-dir output
```

| Option | Use |
|--------|-----|
| `--num-organisms` | Number of genomes per group (bacterial and viral). |
| `--num-archaea`, `--num-plasmid` | Include archaeal/plasmid genomes as additional negatives. |
| `--accessions-file` | Load accession IDs from JSON (skip NCBI search); use with snapshot for reproducibility. |
| `--save-accessions` | After searching, save accession list and timestamp to JSON for later `--accessions-file` runs. |
| `--complete-only` | When searching NCBI (no `--accessions-file`), restrict to complete genomes only (exclude WGS/draft). For reproducible complete-only runs, create a snapshot with `snapshot --complete-only` then use that JSON as `--accessions-file`. |

---

### Chunk genomes into reads

**Purpose:** Turn genome FASTAs into a single metagenome FASTA of fixed- or variable-length reads. Input is typically the download directory (`bacteria/`, `virus/`, `archaea/`, `plasmid/`).

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
| `--sequence-length` | Fixed read length (nt). |
| `--reads-per-organism` | Max reads per genome file; omit for all possible reads. |
| `--balanced` | Same number of reads per genome (min across files); avoids bacterial dominance. |
| `--cap-total-reads` | Downsample to at most N reads (e.g. match positive set size for balanced train data). |
| `--min-contig-length`, `--max-contig-length` | Variable-length contigs (e.g. 300–2000 bp) instead of fixed length. |
| `--seed` | Random seed for variable-length sampling, cap, and mutation; use for reproducibility. |
| `--eve-intervals` | Path to `eve_intervals.json` from `blastn-filter`; exclude chunks overlapping EVE regions. |
| `--forbid-ambiguous` | Discard reads containing non-ACGT (e.g. N). |
| `--substitution-rate` | Per-base substitution rate (0–1) for robustness benchmarks (e.g. 0.01). |
| `--indel-rate` | Per-base indel rate (0–1); may change read length. |
| `--extra-viral-fasta` | FASTA of additional viral sequences (e.g. metavirome); chunked and merged with RefSeq viral. |
| `--abundance-profile` | Per-category read weights, e.g. `bacterial=0.5,viral=2,archaea=1,plasmid=1`. |
| `--abundance-distribution` | `exponential`: per-genome weights from Exp(1); use `--seed` for reproducibility. |
| `--viral-taxonomy` | JSON mapping viral_1, viral_2, ... to taxonomy group (from `viral-taxonomy` command). |
| `--balance-viral-by-taxonomy` | Balance viral reads so each taxonomy group contributes equally (requires `--viral-taxonomy`). |
| `--filter-similar` | Drop reads ≥90% similar to already-kept; oversample and refill to target. |
| `--train-test-split` | Write train/test FASTAs (e.g. 80% train); remove from test reads similar to train. |

---

### Full pipeline (download + chunk)

**Purpose:** One command to download and chunk; optionally run BLASTN (EVE removal) and Seeker. Uses a fixed layout: `output-dir/downloaded/`, `metagenome/`, `blastn/`, `seeker/`, `logs/`.

```bash
metagenome-generator pipeline \
  --num-organisms 10 \
  --output-dir output \
  --output metagenome.fasta \
  --sequence-length 250 \
  --reads-per-organism 1000
```

Pipeline options include all chunk options plus: `--run-blastn-filter`, `--run-seeker`, `--accessions-file`, `--complete-only`, `--forbid-ambiguous`, `--substitution-rate`, `--indel-rate`, `--extra-viral-fasta`, `--abundance-profile`, `--abundance-distribution`. See `metagenome-generator pipeline --help`.

---

### Accession snapshot (no downloads)

**Purpose:** Save a full catalog of matching NCBI accessions (bacterial, viral, archaeal, plasmid) to a date-stamped JSON. Use to subset or re-download the same set later; no sequences are fetched.

```bash
metagenome-generator snapshot
```

With metadata (default), each category stores `{accession, create_date, title}` per entry for temporal split and auditing. Use `--no-metadata` for lists only. Use `--complete-only` to restrict the catalog to complete genomes only (exclude WGS/draft); then use that snapshot with `--accessions-file` for reproducible complete-only downloads. Convert legacy snapshots with `metagenome-generator migrate-snapshot <path>`.

---

### Viral taxonomy (for taxonomy-aware balancing)

**Purpose:** Fetch viral taxonomy from NCBI (family or realm) and write a JSON mapping `viral_1`, `viral_2`, ... to group name. Use this file with `chunk` or `pipeline` and `--balance-viral-by-taxonomy` so each viral taxonomy group contributes equally to the metagenome.

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

### Train/test split and similarity

The tool supports two evaluation patterns. In both, a **similarity check** (BLAST: remove from test any read ≥ threshold similar to train) is important. For temporal split you run it as a dedicated step after building train and test metagenomes; for percentage split it is applied automatically.

**Best practices (references):** Temporal holdout is recommended for viral/metagenomic classifiers so the test set reflects “novel” sequences not seen at training time. DeepVirFinder (Ren et al., *Bioinformatics* 2020; [PubMed 34084563](https://pubmed.ncbi.nlm.nih.gov/34084563/)) trained on viral RefSeq discovered **before May 2015** and evaluated on sequences **after May 2015**; VirFinder (Ren et al., *Microbiome* 2017) used a similar temporal split (before/after Jan 2014). Splitting by date alone is not enough: **different strains or closely related genomes** (e.g. same species, different isolate) can appear in both train and test and inflate metrics. Removing test reads that are highly similar to train is therefore a best practice; the tool provides this for both workflows below.

---

#### Option A — Temporal split (by NCBI CreateDate)

Split **accessions** by submission date so train and test are disjoint in time. Build two separate metagenomes, then **run the similarity filter** so that test reads that are highly similar to train (e.g. different strains of the same species) are removed. **This step is important:** without it, closely related sequences can appear in both sets and inflate evaluation metrics.

**Step 1 — Preview counts (no files written):**

```bash
metagenome-generator temporal-split-info \
  --accessions-file snapshots/accession_snapshot_2026-03-10.json \
  --split-date 2019-06-01
```

**Step 2 — Write train and test JSONs:**

```bash
metagenome-generator temporal-split \
  --accessions-file snapshots/accession_snapshot_2026-03-10.json \
  --split-date 2019-06-01
```

**Step 3 — Build train and test datasets:**

Run `download` (or `pipeline`) twice: once with `--accessions-file train_<basename>.json`, once with `--accessions-file test_<basename>.json`, into separate output dirs. Chunk each to get `train_metagenome.fasta` and `test_metagenome.fasta`.

**Step 4 — Filter test against train (similarity check; required for rigorous evaluation):**

Remove from the test set any read that is ≥ threshold similar to a train read (BLAST). This avoids counting different strains or near-duplicates as “novel” test sequences.

```bash
metagenome-generator filter-test-against-train \
  --train-fasta output_train/metagenome/train_metagenome.fasta \
  --test-fasta output_test/metagenome/test_metagenome.fasta \
  --output output_test/metagenome/test_metagenome_filtered.fasta \
  --similarity-threshold 90
```

Use `test_metagenome_filtered.fasta` as your test set for evaluation. Options: `--min-coverage` (default 0.8), `--threads`, `--batch-size`. Requires BLAST+.

---

#### Option B — Percentage split with similarity check (single metagenome)

Build **one** metagenome, then split the resulting reads by percentage (e.g. 80% train, 20% test). The tool **removes from the test set any read that is ≥ similarity threshold** (default 90% identity over 80% of length) similar to a train read (BLAST of test vs train). This avoids inflated metrics from near-duplicate train/test sequences.

Use with **chunk** or **pipeline**:

```bash
metagenome-generator chunk --input output/downloaded --output metagenome.fasta --output-dir output/metagenome \
  --sequence-length 250 --reads-per-organism 1000 \
  --train-test-split 80 --seed 42
```

Output: `metagenome_train.fasta` and `metagenome_test.fasta`; test reads similar to train are dropped (see log). Options: `--train-test-similarity-threshold` (default 90), `--train-test-blast-threads`, `--train-test-blast-batch-size`.

In the **pipeline**, add `--train-test-split 80` (and optionally `--seed 42`); the same similarity check runs after the single metagenome is built.

---

### BLASTN filtering (EVE removal)

**Purpose:** Endogenous viral elements (EVEs) in non-viral genomes can be misclassified as viral. BLAST non-viral vs viral; exclude chunks overlapping hits when building the metagenome.

**Standalone:**

```bash
metagenome-generator blastn-filter --genome-dir output/downloaded --out-dir output/blastn --evalue 1e-5 --perc-identity 70 \
  --export-eve-fasta output/blastn/eve_intervals.fasta --export-eve-min-length 200
metagenome-generator chunk --input output/downloaded --output metagenome.fasta --output-dir output/metagenome \
  --balanced --eve-intervals output/blastn/eve_intervals.json
```

**In pipeline:** add `--run-blastn-filter` and optionally:
- `--blastn-evalue`, `--blastn-perc-identity`
- `--blastn-export-eve-fasta PATH` (with `--blastn-export-eve-min-length N`)

Requires BLAST+.

---

### Seeker (phage/bacteria prediction)

**Purpose:** Run [Seeker](https://github.com/gussow/seeker) on the generated metagenome to label reads and export predicted phage reads. Seeker must be installed in a separate conda env.

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
| `blastn-filter` | BLAST non-viral vs viral; write EVE intervals for chunk and optional provirus/EVE FASTA. |
| `viral-taxonomy` | Fetch viral taxonomy from NCBI; write viral_1→group JSON for `--balance-viral-by-taxonomy`. |
| `seeker` | Run Seeker on a metagenome FASTA. |

Full options: `metagenome-generator <command> --help`.

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

- NCBI rate limits apply; the tool delays requests and retries (up to 3 times with backoff).
- Genome selection uses RefSeq and length filters; see `DEFAULT_QUERIES` in `ncbi_search.py` to change criteria.
- Prefer a dedicated working directory for runs; see `working_directory/README.md`.
