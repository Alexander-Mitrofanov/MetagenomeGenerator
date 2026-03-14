# Detailed Use Cases

This document describes concrete workflows with step-by-step commands and explanations. All paths assume you run from the project root and use the existing snapshot in `snapshots/`.

---

## Use case: Temporal train/test with fixed genome counts and minimum read budget

### Objective

Build a temporal train/test dataset for classifier evaluation with:

- **Train set:** 100 genomes per category (bacteria, virus, archaea, plasmid) from accessions with **CreateDate before 2025** (older genomes).
- **Test set:** 25 genomes per category from accessions with **CreateDate in 2025 or later** (newer genomes).
- **Reads:** At least 10,000 reads in the train set, each of length **1000 nucleotides**.
- **Similarity:** Remove from the test set any read that is highly similar to a train read (avoid train–test contamination).
- **Location:** All outputs in a single folder under the work directory, e.g. `working_directory/temporal_100_25/`.

This mimics a “train on the past, test on the future” evaluation: the model never sees 2025+ genomes during training, and the test set is checked for similarity so that near-duplicate strains do not inflate metrics.

### Why these choices

- **Temporal split (2025-01-01):** Splitting by NCBI submission date gives a realistic evaluation setting (generalization to newly submitted genomes). The cutoff 2025-01-01 separates “before 2025” (train) from “2025 and later” (test).
- **100 train / 25 test per category:** Fixed counts make the benchmark reproducible and balanced across the four categories. You can subsample from the temporal split JSONs with `--max-bacteria`, `--max-virus`, etc., and `--sample-seed` for a deterministic subset.
- **10,000 reads of length 1000:** Ensures a minimum useful size for training. With 400 train genomes (100×4 categories), 10,000 reads implies at least 25 reads per genome on average; using `--reads-per-organism 30` (or higher) guarantees ≥10,000 reads.
- **Similarity filter:** `filter-test-against-train` removes test reads that are too similar to train reads (e.g. same species, different strain), so reported performance reflects generalization rather than memorization.

### Prerequisites

- Existing snapshot, e.g. `snapshots/accession_snapshot_2026-03-10.json` (from `metagenome-generator snapshot`).
- NCBI Entrez credentials set: `ENTREZ_EMAIL` and optionally `ENTREZ_API_KEY`.
- Working directory created: `working_directory/temporal_100_25/` (or your chosen name).

### Step-by-step commands

All commands are run from the **project root** (`MetagenomeGenerator/`). The snapshot used here is `snapshots/accession_snapshot_2026-03-10.json`; replace the date with your snapshot’s date if different.

---

#### Step 1: Create train and test accession lists by date

Split the snapshot into two JSONs: one for accessions with CreateDate **before** 2025-01-01 (train) and one for 2025-01-01 **and later** (test). No genomes are downloaded yet; only accession lists are written.

```bash
metagenome-generator temporal-split \
  --accessions-file snapshots/accession_snapshot_2026-03-10.json \
  --split-date 2025-01-01 \
  --output-train snapshots/train_accession_snapshot_2026-03-10.json \
  --output-test snapshots/test_accession_snapshot_2026-03-10.json
```

**What it does:** Reads the snapshot, uses the stored `create_date` per accession (no NCBI call if metadata is present), and writes:

- `train_...json`: accessions with CreateDate &lt; 2025-01-01  
- `test_...json`: accessions with CreateDate ≥ 2025-01-01  

These files have the same structure as the snapshot (bacterial, viral, archaea, plasmid lists) and are used as `--accessions-file` in the next steps.

---

#### Step 2: Download train genomes (100 per category)

Download exactly 100 bacterial, 100 viral, 100 archaeal, and 100 plasmid genomes from the **train** list. The tool takes a reproducible random sample when the list is larger than 100 per category.

```bash
metagenome-generator download \
  --accessions-file snapshots/train_accession_snapshot_2026-03-10.json \
  --max-bacteria 100 \
  --max-virus 100 \
  --max-archaea 100 \
  --max-plasmid 100 \
  --sample-seed 42 \
  --output-dir working_directory/temporal_100_25/train_downloaded
```

**What it does:** Does not search NCBI; only downloads the accessions listed in the train JSON. For each category, it randomly samples up to 100 accessions (or fewer if the list has fewer). `--sample-seed 42` makes this subset identical on every run. Genomes are saved as `{accession}.fasta` under `train_downloaded/bacteria/`, `train_downloaded/virus/`, etc.

---

#### Step 3: Download test genomes (25 per category)

Download 25 genomes per category from the **test** list (accessions from 2025 onward).

```bash
metagenome-generator download \
  --accessions-file snapshots/test_accession_snapshot_2026-03-10.json \
  --max-bacteria 25 \
  --max-virus 25 \
  --max-archaea 25 \
  --max-plasmid 25 \
  --sample-seed 42 \
  --output-dir working_directory/temporal_100_25/test_downloaded
```

**What it does:** Same as Step 2 but for the test JSON and with a cap of 25 per category. Output directory is `working_directory/temporal_100_25/test_downloaded/`.

---

#### Step 4: Generate train metagenome (≥10,000 reads, length 1000)

Turn the train genomes into a single FASTA of simulated reads: length 1000 nt and at least 10,000 reads. With 400 train genomes, `--reads-per-organism 30` gives at least 12,000 reads. Write directly into the output base so the final folder has no redundant subfolders.

```bash
metagenome-generator chunk \
  --input working_directory/temporal_100_25/train_downloaded \
  --output train_metagenome.fasta \
  --output-dir working_directory/temporal_100_25 \
  --sequence-length 1000 \
  --reads-per-organism 30 \
  --seed 42
```

**What it does:** Scans all FASTA files under `train_downloaded/`, splits each genome into non-overlapping 1000 nt reads, and takes up to 30 reads per genome. Writes `working_directory/temporal_100_25/train_metagenome.fasta`. Read headers include source accession and position for traceability.

---

#### Step 5: Generate test metagenome (same length)

Generate the test metagenome with the same read length so train and test are comparable. Write into the same output base.

```bash
metagenome-generator chunk \
  --input working_directory/temporal_100_25/test_downloaded \
  --output test_unfiltered.fasta \
  --output-dir working_directory/temporal_100_25 \
  --sequence-length 1000 \
  --reads-per-organism 30 \
  --seed 43
```

**What it does:** Same as Step 4 but for the test genomes. A different seed (43) keeps read sampling independent of the train set. Output: `working_directory/temporal_100_25/test_unfiltered.fasta`.

---

#### Step 6: Remove test reads similar to train (similarity filter)

Remove from the test FASTA any read that is too similar to a train read (e.g. ≥90% identity over most of the read length). Write the final test set as `test_metagenome.fasta` in the same folder.

```bash
metagenome-generator filter-test-against-train \
  --train-fasta working_directory/temporal_100_25/train_metagenome.fasta \
  --test-fasta working_directory/temporal_100_25/test_unfiltered.fasta \
  --output working_directory/temporal_100_25/test_metagenome.fasta \
  --similarity-threshold 90
```

**What it does:** Runs BLAST to find test reads that match train reads above the given identity threshold and excludes them. The result is written to `test_metagenome.fasta` (use this for evaluation). Requires BLAST+ installed.

---

### Summary of outputs

The final output folder contains only the minimal set of items: two downloaded folders, one blast results folder (if EVE was run), and two FASTA files.

| Path | Description |
|------|-------------|
| `working_directory/temporal_100_25/train_downloaded/` | Train genomes (100 per category, CreateDate &lt; 2025). |
| `working_directory/temporal_100_25/test_downloaded/` | Test genomes (25 per category, CreateDate ≥ 2025). |
| `working_directory/temporal_100_25/blastn/` | BLAST/EVE results (only if you ran blastn-filter; contains `train/`, `test/`, and accession JSONs when using temporal-pipeline). |
| `working_directory/temporal_100_25/train_metagenome.fasta` | Train reads (≥10k, length 1000). |
| `working_directory/temporal_100_25/test_metagenome.fasta` | Filtered test reads (use for evaluation). |

### One-shot command (includes similarity filter)

The **`temporal-pipeline`** command runs the full workflow (split → download train/test → optional EVE with `--viral-db` → chunk both → **similarity filter**) in one go. The output dir is kept minimal: `train_downloaded/`, `test_downloaded/`, `blastn/` (with `train/`, `test/`, and accession JSONs), `train_metagenome.fasta`, and `test_metagenome.fasta`.

```bash
metagenome-generator temporal-pipeline \
  --accessions-file snapshots/accession_snapshot_2026-03-10.json \
  --split-date 2025-01-01 \
  --output-dir working_directory/temporal_100_25 \
  --viral-db viral_reference/viral_db_2026-03-10/blastn_db/viral_db
```

Defaults: 100 train and 25 test per category, 1000 nt reads, 30 reads per organism, 90% similarity threshold. Omit `--viral-db` to skip EVE detection (then `blastn/` is still created and holds only the train/test accession JSONs).

### Example run (summary)

A full run was executed with the existing snapshot `snapshots/accession_snapshot_2026-03-10.json`. Results:

- **Train:** 10,669 reads (length 1000 nt) from 400 genomes (100 per category, CreateDate &lt; 2025-01-01).
- **Test (raw):** 2,741 reads from 100 genomes (25 per category, CreateDate ≥ 2025-01-01).
- **Test (filtered):** 2,381 reads after removing test reads similar to train (≥90% identity). Use `working_directory/temporal_100_25/test_metagenome.fasta` for evaluation.
