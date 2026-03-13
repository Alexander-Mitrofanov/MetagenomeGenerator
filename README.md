# CHIMERA

**CHIMERA:** **C**onfigurable **H**ybrid **I**n-silico **M**etagenome **E**mulator for **R**ead **A**nalysis.

CHIMERA is a pipeline for generating synthetic metagenomic datasets in FASTA or FASTQ format for training and benchmarking sequence classifiers, such as viral vs. prokaryotic or phage vs. bacterial classifiers. It constructs artificial metagenomes from reference genomes obtained from sources such as NCBI RefSeq or from user-provided FASTA files organized into category-specific directories (e.g., bacteria, virus, archaea, plasmid).

Input genomes are fragmented into reads or contigs of fixed or variable length, with optional sequencing error and mutation models applied to simulate realistic data. The resulting dataset is written as a single metagenome file in which read identifiers remain traceable to their source genomes, with optional ground-truth metadata describing dataset composition and abundance.

CHIMERA is designed to support reproducible benchmarking. Genome sets can be defined using date-stamped accession snapshots, and fixed random seeds ensure deterministic read sampling across runs. The pipeline also supports strategies to reduce train–test contamination, including temporal splits based on genome release dates and similarity-based filtering of test sequences.

Additional features include optional filtering of endogenous viral elements (EVEs) in microbial genomes, standardized dataset generation with fixed genome counts per category, and support for custom genome collections by simply placing FASTA files into the corresponding category folders.

---

## Table of contents

- [Pre-built snapshots and viral reference](#pre-built-snapshots-and-viral-reference)
- [Accession snapshot](#accession-snapshot)
- [Use cases at a glance](#use-cases-at-a-glance)
- [Installation](#installation)
- [Requirements](#requirements)
- [Quick start](#quick-start)
- [Extended usage](#extended-usage)
- [Command reference](#command-reference)
- [Capabilities summary](#capabilities-summary)
- [Project structure](#project-structure)
- [Notes](#notes)

---

## Pre-built snapshots and viral reference

You do **not** need to build your own accession snapshot or viral reference database. This repository provides both, and they are updated regularly.

- **Accession snapshots** are in the [`snapshots/`](snapshots/) directory. Use any `accession_snapshot_YYYY-MM-DD.json` with `--accessions-file` for reproducible downloads and pipelines. New snapshots are added as the RefSeq catalog is refreshed.
- **Viral reference BLAST databases** (for EVE/prophage detection) are published in the [Releases](https://github.com/Alexander-Mitrofanov/MetagenomeGenerator/releases) section. Each release includes a date-stamped viral DB (e.g. `viral_db_YYYY-MM-DD`) built from the same snapshot date. Download the latest release asset (e.g. `viral_db_2026-03-10.tar.gz`), extract it, and pass the BLAST DB path to `blastn-filter --viral-db`. That way you can run EVE filtering without running `build-viral-db` (which downloads thousands of viral genomes and can take over an hour).

If you need a custom snapshot date or a DB built from a different snapshot, use the `snapshot` and `build-viral-db` commands as described below.

*Maintainers:* To publish a new viral DB, run `build-viral-db`, then create a GitHub Release and attach the dated folder as a tarball (e.g. `tar czvf viral_db_YYYY-MM-DD.tar.gz -C viral_reference viral_db_YYYY-MM-DD`).

---

## Accession snapshot

NCBI's RefSeq catalog changes over time (new submissions, retractions, taxonomy updates). If you search and download "N bacterial and N viral genomes" on different dates, you get different accession sets, so experiments are not reproducible. The **snapshot** command records the **current catalog** of matching accessions to a date-stamped JSON **without downloading any sequences**. That JSON is a frozen accession list: later, `download --accessions-file <json>` (or `pipeline` with the same file) fetches exactly those accessions, so the same snapshot always yields the same genome set. **Snapshot = reproducible genome selection**; **download with accessions-file = same genomes every time.**

**Using a pre-built snapshot from this repo:** Clone or download this repository and use any snapshot from the [`snapshots/`](snapshots/) directory with `--accessions-file snapshots/accession_snapshot_YYYY-MM-DD.json`. Snapshots are updated regularly, so you can use the latest without running `snapshot` yourself.

**Creating your own snapshot (optional):**

```bash
metagenome-generator snapshot
```

The output path is **automatic**: the file is written to `snapshots/accession_snapshot_YYYY-MM-DD.json` (current date). The `snapshots/` directory is created if needed. To use a custom path, pass `--output <path>`.

The command queries NCBI for all RefSeq accessions that match the chosen categories and writes them to the JSON. **It can take a long time** (tens of minutes to over an hour, depending on catalog size and rate limits). Use a [pre-built snapshot](#pre-built-snapshots-and-viral-reference) from this repo when possible; run `snapshot` only when you need a custom date or a fresh catalog.

**Snapshot contents.** The JSON contains four category lists: **bacterial**, **viral**, **archaeal**, and **plasmid**. Each list holds accession IDs (e.g. `NC_000001.1`). By default, each accession also has **create_date** (NCBI submission date, for temporal train/test splits) and **title** (genome description, for auditing). Use `--no-metadata` to store only the ID lists. Use `--complete-only` to restrict the catalog to complete genomes only; then use that snapshot with `--accessions-file` for reproducible complete-only runs.

**Using a snapshot.** Pass the JSON to `download` or `pipeline` via `--accessions-file <path>`. To download only a subset (snapshots can contain tens or hundreds of thousands of IDs), use **`--max-bacteria N --max-virus M`** (and optionally `--max-archaea`, `--max-plasmid`) with **`--sample-seed`** so the subset is reproducible.

**Replicating past experiments.** To reproduce results from a paper or an earlier run, use the **same snapshot file** that was used then. Keep old snapshots (e.g. `accession_snapshot_2025-01-15.json`) in version control or an archive; with that file and the same seeds, you can regenerate the same genome set and metagenome.

---

## Use cases at a glance

(For snapshot-based flows, see [Accession snapshot](#accession-snapshot) above.)

| Use case | Objective | Command or flow |
|----------|-----------|------------------|
| **Single metagenome** | Generate one synthetic metagenome FASTA (fixed or variable read length) for classifier training or method benchmarking, with controlled genome counts and read parameters. | `pipeline --num-bacteria N --num-virus N --output-dir out --output metagenome.fasta --sequence-length 250 --reads-per-organism 1000` |
| **In-house genome set** | Use your own genome FASTA files (isolates, assemblies, phages) instead of NCBI: place them in `bacteria/`, `virus/`, etc., then run `chunk` with `--input` pointing to that directory. | Create folder layout → put FASTAs in the right category folders → `chunk --input my_genomes --output metagenome.fasta --output-dir out --sequence-length 250 --reads-per-organism 1000` |
| **Reproducible genome set** | Fix the set of genomes used across runs and machines: record the catalog once with `snapshot`, then download and generate reads from that list so that results depend only on the snapshot and seeds, not on NCBI’s current state. Optionally limit to a subset with `--max-bacteria`, `--max-virus`, and `--sample-seed`. | `snapshot` → save JSON; then `download --accessions-file <json>` (and run read generation) or use that file in `pipeline`. To use a subset: add `--max-bacteria N --max-virus M --sample-seed 42`. |
| **Temporal train/test** | Evaluate generalization to “future” genomes: train on accessions submitted before a cutoff date and test on accessions on or after that date, with BLAST-based removal of test reads that are highly similar to train (avoids inflated metrics from near-identical strains). | Optional: `temporal-split-search --min-train N --min-test M` to get a split date. Then `temporal-split-info` (preview) → `temporal-split` → build train and test metagenomes → **`filter-test-against-train`**. |
| **Single-dataset train/test** | Split one synthetic metagenome into train and test fractions (e.g. 80/20) with automatic removal of test reads that are ≥ threshold similar to train, for quick evaluation without temporal split. | `chunk` or `pipeline` with `--train-test-split 80` (similarity filter applied automatically). |
| **Structured benchmark** | Produce multiple replicate datasets with fixed N genomes per category (e.g. 50 bacterial, 50 viral per replicate) sampled from a snapshot, for method comparison and reporting mean ± std across replicates in a standardized way. | `snapshot` → `benchmark-recipe --accessions-file <snap> --output-dir out --per-category 50 --replicates 5` |

For a detailed walkthrough (temporal train/test with fixed genome counts, read budget, and similarity filter), see [Detailed Use Cases](Detailed_Use_Cases.md).

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
  --num-bacteria 10 \
  --num-virus 10 \
  --output-dir output \
  --output metagenome.fasta \
  --sequence-length 250 \
  --reads-per-organism 1000
```

Result: genomes in `output/downloaded/`, metagenome FASTA in `output/metagenome/metagenome.fasta`.

**Do it in two steps (download, then generate reads):**

```bash
metagenome-generator download --num-bacteria 10 --num-virus 10 --output-dir output
metagenome-generator chunk \
  --input output/downloaded \
  --output metagenome.fasta \
  --output-dir output/metagenome \
  --sequence-length 250 \
  --reads-per-organism 1000
```

**Equal reads per genome (balanced):** add `--balanced` to the read-generation command (`chunk`).

---

## Extended usage

### Download genomes

Obtain RefSeq genomes by category (bacteria, virus, archaea, plasmid) from NCBI Nucleotide. You specify **how many bacteria** and **how many virus** genomes separately; optionally add archaea and plasmid as extra negative samples. Each genome is saved as **`{accession}.fasta`** (e.g. `NC_000001.1.fasta`) in the corresponding category folder. Output layout: `bacteria/`, `virus/`, `archaea/`, `plasmid/` under the output directory.

```bash
metagenome-generator download --num-bacteria 10 --num-virus 10 --output-dir output
```

To use a **reproducible subset** from an existing snapshot (e.g. 50 bacterial + 50 viral) instead of downloading the whole file:

```bash
metagenome-generator download \
  --accessions-file snapshots/accession_snapshot_2026-03-10.json \
  --max-bacteria 50 --max-virus 50 \
  --sample-seed 42 \
  --output-dir output/downloaded
```

| Option | Use |
|--------|-----|
| `--num-bacteria` | **Number of bacteria genomes.** How many RefSeq bacterial genomes to fetch. Use for negative (non-viral) samples in viral vs. prokaryotic classifiers. Set to 0 only when using `--accessions-file` with a snapshot that has no bacteria list. |
| `--num-virus` | **Number of virus genomes.** How many RefSeq viral genomes to fetch. Use for positive (viral) samples. Set to 0 only when using `--accessions-file` with a snapshot that has no virus list. |
| `--num-archaea` | **Number of archaea genomes.** Optional; default 0. Archaea are additional negative samples (non-viral). Use to broaden the diversity of non-viral sequences (e.g. for phage vs. bacteria + archaea). |
| `--num-plasmid` | **Number of plasmid sequences.** Optional; default 0. Plasmids are additional negative samples. Use when you want to avoid classifying plasmid-derived reads as viral. |
| `--output-dir` | **Output directory.** All category folders (`bacteria/`, `virus/`, etc.) are created under this path. Use a dedicated directory (e.g. `working_directory/downloaded/`) to keep runs organized. |
| `--accessions-file` | **Reproducible run.** Path to a JSON file containing accession IDs (e.g. from `snapshot` or a previous `--save-accessions` run). NCBI search is skipped; by default **all** accessions in the file are downloaded. Use when you need the same genome set on every run (e.g. for benchmarks or paper reproducibility). |
| `--max-bacteria`, `--max-virus`, `--max-archaea`, `--max-plasmid` | **Limit how many to use from the snapshot.** When using `--accessions-file`, these set an upper bound per category: the tool takes a **random sample** of that many accessions (or all if the file has fewer). Omit to download the full snapshot. Example: `--accessions-file snap.json --max-bacteria 50 --max-virus 50` downloads 50 bacterial + 50 viral from the file. |
| `--sample-seed` | **Reproducible subset.** When using `--max-*` with `--accessions-file`, seed for the random sample (default 42). Use the same seed to get the same subset on every run. |
| `--save-accessions` | **Save chosen accessions.** After searching NCBI, write the selected accession lists and a UTC timestamp to this JSON path. Use this file later as `--accessions-file` to re-download the same set. Ignored when `--accessions-file` is set. |
| `--complete-only` | **Complete genomes only.** When searching NCBI (no `--accessions-file`), restrict results to complete genomes and exclude WGS/draft (uses NCBI `complete[Properties]` and `NOT WGS[Properties]`). For reproducible complete-only runs, create a snapshot with `snapshot --complete-only` and use that JSON as `--accessions-file`. Ignored when using `--accessions-file`. |

For large snapshots, use `--max-*` and `--sample-seed` (see [Accession snapshot](#accession-snapshot)).

---

### Using your own (in-house) genome set

You can skip the download step and use **your own genome FASTA files**. Use the same folder layout the tool expects:

| Folder       | Contents |
|-------------|----------|
| `bacteria/` | One or more FASTA files (e.g. your bacterial isolates or assemblies). |
| `virus/`    | One or more FASTA files (e.g. your viral sequences or phages). |
| `archaea/`  | Optional. FASTA files for archaeal genomes. |
| `plasmid/`  | Optional. FASTA files for plasmid sequences. |

**Requirements:** At least one file in **`virus/`** and at least one file in **one of** `bacteria/`, `archaea/`, or `plasmid/` (so that both viral and non-viral categories are present). Empty folders are ignored.

**File naming:** Use any filename (e.g. `isolate_001.fasta`, `NC_12345.fasta`). The **file stem** (filename without `.fasta`) becomes the genome identifier in the output (e.g. `isolate_001_read_0` with description `start=0 end=250`). Multi-record FASTA files are supported: all records in a file share the same prefix.

**Workflow:** Create the directory, place your FASTA files in the correct category folders, then run the read-generation command **`chunk`** with `--input` pointing to that directory. You do **not** need to run `download`.

```bash
# Example: in-house data in my_genomes/
# my_genomes/bacteria/isolate_A.fasta  my_genomes/bacteria/isolate_B.fasta
# my_genomes/virus/phage_1.fasta

metagenome-generator chunk \
  --input my_genomes \
  --output metagenome.fasta \
  --output-dir output/metagenome \
  --sequence-length 250 \
  --reads-per-organism 1000
```

You can also **mix** NCBI-downloaded and in-house data: run `download` into a directory, then copy or symlink your own FASTA files into the same `bacteria/`, `virus/`, etc. folders before running the `chunk` (read-generation) step.

---

### Generate reads from genomes (`chunk`)

The `chunk` subcommand turns genome FASTAs into one metagenome FASTA (or FASTQ) by splitting each genome into fixed-length simulated reads or variable-length contigs. Input is either the **download** output directory or **your own directory** with the same layout (`bacteria/`, `virus/`, `archaea/`, `plasmid/`). See [Using your own (in-house) genome set](#using-your-own-in-house-genome-set) above.

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
| `--sequence-length` | **Fixed read length (nt).** Each simulated read is exactly this many nucleotides. Typical values: 250–500 for short-read style; match your classifier’s expected input (e.g. 250 for some viral classifiers). Required unless you use variable-length mode. |
| `--reads-per-organism` | **Max reads per genome.** Upper limit on how many non-overlapping (or sampled) reads are taken from each genome file. Omit to use all possible reads from every genome (can produce very large outputs). Use a fixed value (e.g. 1000) for controlled dataset size and balance across genomes when combined with other options. |
| `--balanced` | **Same number of reads per genome.** Each genome contributes the same count of reads (the minimum across all genomes). Use when you want to avoid one category (e.g. bacteria) dominating simply because genomes are longer; good for training classifiers with even representation per genome. |
| `--cap-total-reads` | **Cap total reads.** Downsample the whole metagenome to at most N reads. Use to match a target size (e.g. cap to the size of your positive set) or to keep evaluation sets manageable. Applied after per-genome limits and balancing. |
| `--min-contig-length`, `--max-contig-length` | **Variable-length contigs.** Instead of fixed-length reads, sample contigs with lengths uniformly between these two values (nt). Use for long-read or contig-level benchmarks (e.g. 300–2000 bp). Omit both to use fixed `--sequence-length`. |
| `--seed` | **Random seed.** Fixes randomness for variable-length sampling, cap, mutation, and train/test split. Use the same seed to reproduce the exact same metagenome; change the seed to get a different sample. |
| `--eve-intervals` | **EVE exclusion.** Path to `eve_intervals.json` produced by `blastn-filter`. Reads/contigs that overlap these endogenous viral element (EVE) intervals on non-viral genomes are excluded from the metagenome. Use to avoid bacterial/archaeal regions that look viral and would confound viral vs. non-viral classifiers. |
| `--forbid-ambiguous` | **Exclude ambiguous bases.** Discard any read that contains non-ACGT characters (e.g. N, R, Y). Use when your pipeline or classifier assumes strict ACGT-only sequence, or to simulate cleaner sequencing. |
| `--substitution-rate`, `--indel-rate` | **Mutation simulation.** Introduce substitutions and/or indels at the given per-base rate (0–1). Use to test classifier robustness to sequencing error or divergence (e.g. 0.01 for 1% substitution rate). Combine with `--seed` for reproducible mutated datasets. |
| `--error-model` | **Platform-specific sequencing errors.** Set to `illumina` for position-dependent substitution. See [Error model, FASTQ, and abundance file](#error-model-fastq-and-abundance-file) below. |
| `--output-fastq` | **FASTQ output.** Write single-end FASTQ with per-base Phred qualities. See [Error model, FASTQ, and abundance file](#error-model-fastq-and-abundance-file) below. |
| `--write-abundance` | **Ground-truth abundance file.** Write `{output_stem}_abundance.txt` (genome_id, read_count, proportion). See [Error model, FASTQ, and abundance file](#error-model-fastq-and-abundance-file) below. |
| `--extra-viral-fasta` | **Merge user viral sequences.** Path to a FASTA of additional viral sequences (e.g. metavirome contigs, custom viral set). Reads are generated from them as for RefSeq viral genomes and merged into the viral pool. Use to combine public RefSeq viral data with your own viral contigs in one metagenome. |
| `--abundance-profile` | **Per-category read weights.** Comma-separated `category=weight`, e.g. `bacteria=0.5,virus=2,archaea=1,plasmid=1`. Scales how many reads are taken from each category relative to the base limit. Use to simulate uneven community composition (e.g. more virus, less bacteria) without changing genome lists. |
| `--abundance-distribution` | **Per-genome abundance model.** Set to `exponential` to assign each genome a weight from an exponential distribution (then normalized). Produces a few “abundant” and many “rare” genomes, similar to real communities. Use `--seed` for reproducibility. |
| `--viral-taxonomy`, `--balance-viral-by-taxonomy` | **Taxonomy-aware viral balancing.** `--viral-taxonomy` is the path to the JSON from the `viral-taxonomy` command (viral accession → taxonomy group). With `--balance-viral-by-taxonomy`, viral read limits are set so each taxonomy group (e.g. family) contributes equally. Use to avoid a few viral families dominating and to better train on under-represented groups. |
| `--filter-similar` | **Within-metagenome similarity filter.** Remove any read that is ≥90% similar (identity and coverage) to a read already kept. The tool oversamples and refills to try to reach the target count. Use to reduce near-duplicate sequences in a single metagenome. |
| `--train-test-split` | **Train/test split with similarity filter.** Percentage of reads for training (e.g. 80). Outputs `*_train.fasta` and `*_test.fasta`. Any test read that is ≥ similarity threshold (default 90% identity over 80% length) to a train read is removed. Use for quick evaluation from one metagenome while avoiding inflated metrics from near-duplicate train/test pairs. |

**Read and contig IDs; traceability.** Fixed-length segments are named **reads** (`{accession}_read_{idx}`); variable-length segments are **contigs** (`{accession}_contig_{idx}`). The FASTA/FASTQ description includes **`start=` and `end=`** (0-based positions on the source genome) so you can trace each read or contig back to its origin. With accession-named genome files (e.g. `NC_000001.1.fasta`), the prefix in the ID is the accession.

---

### Error model, FASTQ, and abundance file

- **`--error-model illumina`** — Position-dependent substitution (low at 5′, higher at 3′); use for realistic short-read benchmarking. Use `--seed` for reproducibility.
- **`--output-fastq`** — Write FASTQ with per-base Phred qualities (Illumina-like). Use when downstream tools expect FASTQ (e.g. aligners).
- **`--write-abundance`** — Write `{stem}_abundance.txt` next to the metagenome (columns: genome_id, read_count, proportion). Use as ground truth for abundance estimators or method papers.

Example: add `--error-model illumina --output-fastq --write-abundance --seed 42` to `chunk` or `pipeline`.

---

### Pipeline (download + read generation)

One command to download genomes and generate reads; optionally run BLASTN (EVE) and Seeker. Layout: `output-dir/downloaded/`, `metagenome/`, `blastn/`, `seeker/`, `logs/`.

```bash
metagenome-generator pipeline \
  --num-bacteria 10 \
  --num-virus 10 \
  --output-dir output \
  --output metagenome.fasta \
  --sequence-length 250 \
  --reads-per-organism 1000
```

Pipeline accepts the same read-generation options (e.g. `--train-test-split`, `--balanced`, `--eve-intervals`) plus `--run-blastn-filter`, `--run-seeker`, `--accessions-file`, `--complete-only`. See `metagenome-generator pipeline --help`.

---

### Accession snapshot (reference)

See [Accession snapshot](#accession-snapshot) above. Command: `metagenome-generator snapshot` (writes to `snapshots/accession_snapshot_YYYY-MM-DD.json`); optional `--output <path>`. Options: `--no-metadata` (IDs only), `--complete-only`. Legacy: `metagenome-generator migrate-snapshot <path>`.

---

### Structured benchmark recipe

Fixed N per category (e.g. 50 bacterial, 50 viral), optional replicates; one command, reproducible and comparable to published benchmarks. No NCBI search at recipe time—samples from your snapshot.

**Example:**

```bash
# 1. Create a snapshot once (or use an existing one); output is snapshots/accession_snapshot_YYYY-MM-DD.json
metagenome-generator snapshot

# 2. Run the recipe: 50 bacterial + 50 viral per replicate, 5 replicates, seed 42
metagenome-generator benchmark-recipe \
  --accessions-file snapshots/accession_snapshot_$(date +%Y-%m-%d).json \
  --output-dir benchmarks/run1 \
  --per-category 50 \
  --replicates 5 \
  --seed 42 \
  --sequence-length 250 \
  --reads-per-organism 1000
```

Output: `benchmarks/run1/replicate_001/downloaded/`, `benchmarks/run1/replicate_001/metagenome/metagenome.fasta`, and so on for `replicate_002` … `replicate_005`. Optional: `--archaea 50`, `--plasmid 50` to include archaea/plasmid in each replicate; `--output metagenome.fasta` to set the FASTA filename.

---

### Train/test split and similarity filtering

Two workflows:

- **Temporal split** — Split accessions by NCBI submission date; build train and test metagenomes separately; then run **`filter-test-against-train`** to remove test reads that are highly similar to train. Use when you want “train on past, test on future” (e.g. generalization to novel viruses).
- **Percentage split** — Build one metagenome and split reads (e.g. 80% train, 20% test); the tool automatically removes from test any read ≥ threshold similar to a train read. Use for quick train/test from a single dataset.

Removing test reads similar to train avoids inflated metrics from near-identical strains; CHIMERA supports this for both temporal and percentage split.

---

#### Temporal split (by NCBI CreateDate)

1. **Find a split date** that gives at least N train and M test genomes (optional; use when you want fixed sizes):

   ```bash
   metagenome-generator temporal-split-search \
     --accessions-file snapshots/accession_snapshot_YYYY-MM-DD.json \
     --min-train 100 --min-test 20
   ```
   Prints a suggested `--split-date` and counts. Then use that date in the steps below.

2. **Preview** counts for a chosen date (no files written):

   ```bash
   metagenome-generator temporal-split-info \
     --accessions-file snapshots/accession_snapshot_2026-03-10.json \
     --split-date 2019-06-01
   ```

3. **Write train/test JSONs:**

   ```bash
   metagenome-generator temporal-split \
     --accessions-file snapshots/accession_snapshot_2026-03-10.json \
     --split-date 2019-06-01
   ```

4. **Build train and test metagenomes:** Run `download` (or `pipeline`) twice with `--accessions-file train_<basename>.json` and `--accessions-file test_<basename>.json` into separate dirs; generate reads from each to get `train_metagenome.fasta` and `test_metagenome.fasta`.

5. **Filter test against train (important for rigorous evaluation):** Remove test reads that are ≥ threshold similar to train so different strains or near-duplicates are not counted as “novel” test.

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

EVEs in non-viral genomes can be misclassified as viral. BLAST non-viral vs virus; exclude reads/contigs overlapping hits when building the metagenome.

**Viral reference for proper prophage/EVE detection:** By default, the viral BLAST DB is built from the `virus/` folder in `--genome-dir` (i.e. only the viral genomes you downloaded for that run). Prophage/EVE regions that match viruses *not* in that set are missed. To check against the full viral catalog, use a pre-built viral DB or build one yourself:

- **Pre-built (recommended):** [Pre-built viral reference DBs](#pre-built-snapshots-and-viral-reference) are available from this repository’s [Releases](https://github.com/Alexander-Mitrofanov/MetagenomeGenerator/releases) page. Download the latest `viral_db_YYYY-MM-DD.tar.gz`, extract it, then run:

  ```bash
  # After extracting the release asset (e.g. viral_db_2026-03-10/)
  metagenome-generator blastn-filter --genome-dir output/downloaded --out-dir output/blastn \
    --viral-db /path/to/viral_db_YYYY-MM-DD/blastn_db/viral_db
  ```

- **Build your own:** If you need a DB for a snapshot date not yet in Releases, run `build-viral-db` once (creates `viral_reference/viral_db_YYYY-MM-DD/` using the snapshot date), then pass the printed DB path to `blastn-filter --viral-db`:

  ```bash
  metagenome-generator build-viral-db --accessions-file snapshots/accession_snapshot_YYYY-MM-DD.json --output-dir viral_reference
  metagenome-generator blastn-filter --genome-dir output/downloaded --out-dir output/blastn --viral-db viral_reference/viral_db_YYYY-MM-DD/blastn_db/viral_db
  ```

You can instead pass a FASTA of viral sequences with `--viral-reference-fasta` (the tool will run `makeblastdb` on it).

**Standalone (default: viral DB from genome-dir):**

```bash
metagenome-generator blastn-filter --genome-dir output/downloaded --out-dir output/blastn --evalue 1e-5 --perc-identity 70 \
  --export-eve-fasta output/blastn/eve_intervals.fasta --export-eve-min-length 200
metagenome-generator chunk --input output/downloaded --output metagenome.fasta --output-dir output/metagenome \
  --balanced --eve-intervals output/blastn/eve_intervals.json
```

**In pipeline:** add `--run-blastn-filter`; optional `--blastn-evalue`, `--blastn-perc-identity`, `--blastn-export-eve-fasta`, `--blastn-export-eve-min-length`. Requires BLAST+.

---

### Viral taxonomy (taxonomy-aware balancing)

Fetch viral taxonomy from NCBI and write a JSON mapping **viral accession → taxonomy group** (e.g. `NC_001234.1` → Herpesviridae). Use with the read-generation step (`chunk`) or `pipeline` and `--balance-viral-by-taxonomy` so each viral taxonomy group contributes equally.

```bash
metagenome-generator viral-taxonomy \
  --accessions-file snapshots/accession_snapshot_2026-03-10.json \
  --output output/viral_taxonomy.json \
  --level family
```

Then run read generation with balancing:

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
| `chunk` | Generate reads from genome FASTAs and write one metagenome FASTA (or FASTQ). |
| `pipeline` | Download + read generation (+ optional BLASTN, Seeker). |
| `snapshot` | Save full accession catalog to JSON (no downloads). |
| `temporal-split-info` | Show train/test counts for a split date (no files written). |
| `temporal-split-search` | Find a split date so train set has at least N and test set at least M genomes. |
| `temporal-split` | Write train and test accession JSONs by CreateDate. |
| `filter-test-against-train` | Remove from test FASTA reads similar to train (BLAST). **Use after temporal split** for rigorous evaluation. |
| `migrate-snapshot` | Convert legacy snapshot to per-category metadata format. |
| `blastn-filter` | BLAST non-viral vs viral; EVE intervals for read generation. Use `--viral-db` or `--viral-reference-fasta` for full viral catalog. |
| `build-viral-db` | Download all viral genomes from a snapshot and build a BLAST DB for use with `blastn-filter --viral-db` (proper prophage/EVE detection). |
| `viral-taxonomy` | Fetch viral taxonomy; write accession→group JSON for `--balance-viral-by-taxonomy`. |
| `benchmark-recipe` | **Structured benchmark:** fixed N per category, R replicates; samples from snapshot, no NCBI search. |
| `seeker` | Run Seeker on a metagenome FASTA. |

Full options: `metagenome-generator <command> --help`.

---

## Capabilities summary

| Area | Features |
|------|----------|
| Genomes | Download by category (RefSeq); in-house FASTAs; snapshot for reproducibility; `--complete-only`. |
| Read generation | Fixed/variable length; balanced or weighted; `--forbid-ambiguous`; mutation rates; Illumina-like errors; FASTQ + abundance file. |
| Train/test | Temporal split by CreateDate or percentage split; **filter-test-against-train** / similarity filtering. |
| EVE | BLAST non-viral vs viral; exclude or export provirus regions. |
| Benchmark | `benchmark-recipe`: fixed N per category, R replicates; viral taxonomy balancing; extra viral FASTA; Seeker. |

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
│   ├── benchmark_recipe.py
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
