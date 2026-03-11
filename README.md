# Metagenome Generator

A Python tool to build **simulated metagenome datasets** from NCBI RefSeq genomes for training and evaluating sequence classifiers (e.g. viral vs. prokaryotic).

---

## What it does

**Metagenome Generator** downloads bacterial, viral, archaeal, and plasmid genomes from NCBI, cuts them into fixed- or variable-length reads, and writes FASTA files you can use as training or test data. It supports reproducible runs, temporal train/test splits, and common data-prep steps used in tools like DeepVirFinder and VirFinder.

### Core functionality

- **Download genomes** — Fetch genomes from NCBI Nucleotide by category (bacteria, virus, archaea, plasmid) with RefSeq and completeness filters.
- **Chunk into reads** — Turn genome FASTAs into fixed-length (or variable-length) “reads” and optionally balance or cap reads per genome.
- **Build metagenomes** — Single FASTA with reads from all categories; optional train/test split and similarity filtering.

### What you can use it for

- **Training classifiers** — Generate positive (viral) and negative (bacterial/archaeal/plasmid) reads for virus vs. non-virus (or phage vs. bacteria) models.
- **Reproducible datasets** — Snapshot accession lists (with optional CreateDate/title) and re-download or re-chunk the same set later.
- **Temporal evaluation** — Split by NCBI CreateDate: train on “older” genomes, test on “newer” ones to mimic real-world generalization.
- **Cleaner negatives** — Optionally remove endogenous viral elements (EVEs) via BLASTN and drop near-duplicate reads with similarity filtering.
- **Downstream tools** — Output is standard FASTA; you can run e.g. Seeker (phage/bacteria prediction) on the generated metagenome.

### Commands at a glance

| Command | Purpose |
|--------|--------|
| `download` | Download genomes from NCBI into category folders. |
| `chunk` | Turn genome FASTAs into a single metagenome FASTA of reads. |
| `pipeline` | Download + chunk in one go (optional BLASTN, Seeker). |
| `snapshot` | Save full catalog of matching accessions (no downloads). |
| `temporal-split-info` | Show train/test counts for a date (no files written). |
| `temporal-split` | Write train/test accession JSONs by CreateDate. |
| `migrate-snapshot` | Convert old snapshot format to per-category metadata. |
| `blastn-filter` | BLAST non-viral vs viral, output EVE intervals for chunk step. |
| `seeker` | Run Seeker on a metagenome FASTA (phage/bacteria labels). |

---

## Installation

**From PyPI** (when published):

```bash
pip install metagenome-generator
```

**From source** (development or local install):

```bash
cd MetagenomeGenerator
pip install -e .
```

**Optional (Conda)** for BLAST+ and a dedicated env:

```bash
conda env create -f environment.yml
conda activate metagenome-simulator
pip install -e .
```

## Requirements

- Python 3.8+
- Biopython ≥ 1.83
- BLAST+ (for EVE removal and similarity filtering; install via conda or system)

NCBI Entrez requires an email (and optionally an API key for higher rate limits):

```bash
export ENTREZ_EMAIL="your_email@example.com"
export ENTREZ_API_KEY="your_ncbi_api_key"   # optional
```

---

## Usage

Use the `metagenome-generator` command (or `python -m metagenome_generator`). From the repo root you can also run `python main.py`.

**Tip:** Use a dedicated **working directory** for outputs (e.g. `working_directory/`); it is gitignored. See [working_directory/README.md](working_directory/README.md) for a suggested layout.

---

### Quick start: simple examples

**1. Download 10 bacterial + 10 viral genomes and build one metagenome (250 nt reads, 1000 per genome):**

```bash
metagenome-generator pipeline \
  --num-organisms 10 \
  --output-dir output \
  --output metagenome.fasta \
  --sequence-length 250 \
  --reads-per-organism 1000
```

Result: `output/downloaded/` (genomes) and `output/metagenome/metagenome.fasta` (reads).

**2. Download only (no chunking):**

```bash
metagenome-generator download --num-organisms 10 --output-dir output
```

**3. Chunk existing genomes (e.g. from a previous download):**

```bash
metagenome-generator chunk \
  --input output/downloaded \
  --output metagenome.fasta \
  --output-dir output/metagenome \
  --sequence-length 250 \
  --reads-per-organism 1000
```

**4. Balanced reads (same number of reads per genome, so viral and bacterial contribute equally):**

```bash
metagenome-generator chunk \
  --input output/downloaded \
  --output metagenome.fasta \
  --output-dir output/metagenome \
  --sequence-length 250 \
  --balanced
```

---

### Download genomes

Saves one FASTA per genome in category subfolders: `bacteria/`, `virus/`, `archaea/`, `plasmid/`.

```bash
metagenome-generator download \
  --num-organisms 10 \
  --output-dir output
```

| Argument | Description | Default |
|----------|-------------|---------|
| `--num-organisms` | Number per group (bacterial and viral) | 10 |
| `--output-dir` | Root for `bacteria/`, `virus/`, `archaea/`, `plasmid/` | `output` |
| `--num-archaea` | Archaeal genomes (negatives) | 0 |
| `--num-plasmid` | Plasmids (negatives) | 0 |
| `--accessions-file` | JSON of accession IDs (skip NCBI search; reproducible) | — |
| `--save-accessions` | Save accession list + timestamp to JSON after search | — |

**Reproducibility:** Use `snapshot` to write a full catalog to JSON, or `--save-accessions` during download; then use `--accessions-file` to re-download the same set.

---

### Chunk genomes into reads

Reads genome FASTAs, cuts them into fixed-length reads, writes one metagenome FASTA.

```bash
metagenome-generator chunk \
  --input output/downloaded \
  --output metagenome_250nt.fasta \
  --output-dir output/metagenome \
  --sequence-length 250 \
  --reads-per-organism 1000
```

| Argument | Description | Default |
|----------|-------------|---------|
| `--input` | FASTA file or **directory** of genome FASTAs | *(required)* |
| `--output` | **Filename** for the metagenome FASTA | *(required)* |
| `--output-dir` | Directory to write the output file | `output` |
| `--sequence-length` | Read length (nt) | 250 |
| `--reads-per-organism` | Max reads per genome file | all |
| `--balanced` | Equal reads per file (min across files) | off |
| `--eve-intervals` | JSON from `blastn-filter`; exclude EVE regions | — |
| `--filter-similar` | Drop reads ≥90% similar (BLASTN) to already-kept | off |
| `--train-test-split` | Train % (e.g. 80 → 80% train, 20% test) | off |
| `--forbid-ambiguous` | Exclude reads with non-ACGT (e.g. N) | off |

---

### Full pipeline (download + chunk)

One command for download and chunk; optional BLASTN (EVE removal) and Seeker.

```bash
metagenome-generator pipeline \
  --num-organisms 10 \
  --output-dir output \
  --output metagenome.fasta \
  --sequence-length 250 \
  --reads-per-organism 1000
```

Output layout under `output-dir/`:

| Subfolder | Contents |
|-----------|----------|
| `downloaded/` | Genome FASTAs in `bacteria/`, `virus/`, `archaea/`, `plasmid/` |
| `blastn/` | BLASTN outputs and `eve_intervals.json` (if `--run-blastn-filter`) |
| `metagenome/` | Metagenome FASTA |
| `seeker/` | Seeker predictions (if `--run-seeker`) |
| `logs/` | `pipeline.log` |

Options: `--balanced`, `--filter-similar`, `--train-test-split`, `--run-blastn-filter`, `--run-seeker`, `--accessions-file`, `--forbid-ambiguous`. See `metagenome-generator pipeline --help`.

---

### Accession snapshot (no downloads)

Catalog **all** matching NCBI accessions (bacterial, viral, archaeal, plasmid) without downloading sequences. Saves a date-stamped JSON you can subset and pass to `download` or `pipeline` via `--accessions-file`.

```bash
metagenome-generator snapshot
# or: metagenome-generator snapshot --output snapshots/accession_snapshot_2026-03-10.json
# Lists only (no CreateDate/title): metagenome-generator snapshot --no-metadata
```

**Format:** `timestamp` plus `bacterial`, `viral`, `archaea`, `plasmid`. With metadata (default), each category is a list of `{accession, create_date, title}`. With `--no-metadata`, each is a list of accession ID strings.

| Argument | Description | Default |
|----------|-------------|---------|
| `--output` | Output JSON path | `snapshots/accession_snapshot_YYYY-MM-DD.json` |
| `--no-metadata` | Do not fetch CreateDate/title per accession | off |
| `--metadata-batch-size` | esummary batch size | 500 |

**Convert old snapshot format:** `metagenome-generator migrate-snapshot snapshots/accession_snapshot_YYYY-MM-DD.json`

---

### Temporal train/test split (by NCBI CreateDate)

Split accessions by **CreateDate**: train = before a date, test = on or after (e.g. for DeepVirFinder-style evaluation). Use **`temporal-split-info`** first to see counts, then **`temporal-split`** to write JSONs.

**1. Preview (no files written):**

```bash
metagenome-generator temporal-split-info \
  --accessions-file snapshots/accession_snapshot_2026-03-10.json \
  --split-date 2019-06-01
```

Example output:

```
Using CreateDate from snapshot metadata for 133,445 accessions.

Temporal split preview (split date: 2019-06-01)
  Source: snapshots/accession_snapshot_2026-03-10.json
  Snapshot timestamp: 2026-03-10T20:25:51Z

  Category   |   Total |   Train |    Test | Train % | Test %
  -----------+---------+---------+---------+---------+--------
  bacterial   |  33,089 |   6,942 |  26,147 |   21.0% |  79.0%
  viral       |   8,324 |   4,963 |   3,361 |   59.6% |  40.4%
  archaea     |     460 |     135 |     325 |   29.3% |  70.7%
  plasmid     |  91,581 |  20,322 |  71,259 |   22.2% |  77.8%
  -----------+---------+---------+---------+---------+--------
  Total       | 133,454 |  32,362 | 101,092 |   24.2% |  75.8%

  Train = CreateDate < split date (older entries)
  Test  = CreateDate >= split date (newer entries)

  Run 'temporal-split' with this --split-date to write train/test JSONs.
```

**2. Write train and test JSONs:**

```bash
metagenome-generator temporal-split \
  --accessions-file snapshots/accession_snapshot_2026-03-10.json \
  --split-date 2019-06-01
```

Then run `download` or `pipeline` with `--accessions-file train_<basename>.json` and `--accessions-file test_<basename>.json` to build separate datasets.

| Command / argument | Description |
|--------------------|-------------|
| **temporal-split-info** | Show train/test counts (no files written). |
| `--accessions-file` | Snapshot or accessions JSON |
| `--split-date` | Cutoff **YYYY-MM-DD** (train &lt; date, test ≥ date) |
| `--output-train` / `--output-test` | Paths for output JSONs |

---

### BLASTN filtering (EVE removal)

Remove **endogenous viral elements** from non-viral genomes: BLAST non-viral vs viral, then exclude chunks overlapping hits when building the metagenome.

**Standalone:**

```bash
metagenome-generator blastn-filter --genome-dir output/downloaded --out-dir output/blastn --evalue 1e-5 --perc-identity 70
metagenome-generator chunk --input output/downloaded --output negative.fasta --output-dir output/metagenome \
  --balanced --eve-intervals output/blastn/eve_intervals.json
```

**In pipeline:**

```bash
metagenome-generator pipeline --num-organisms 5 --output-dir output \
  --balanced --run-blastn-filter --blastn-evalue 1e-5 --blastn-perc-identity 70
```

Requires BLAST+ (`conda install -c bioconda blast` or use `environment.yml`).

---

### Seeker (phage/bacteria prediction)

Run [Seeker](https://github.com/gussow/seeker) on a metagenome FASTA to label reads and export predicted phage reads. Seeker must be in a separate conda env (e.g. `seeker`).

```bash
metagenome-generator seeker --input output/metagenome/metagenome.fasta --output-dir output/seeker --conda-env seeker
```

Or add `--run-seeker` to the pipeline.

---

### More options (reference)

- **Variable-length contigs:** `--min-contig-length 300 --max-contig-length 2000` (with `--balanced`, `--seed`).
- **Down-sample negatives:** Chunk viral first, then chunk non-viral with `--cap-total-reads <viral_count>`.
- **Similarity filtering:** `--filter-similar`, `--similarity-threshold 90`, `--oversample-factor 2`.
- **Train/test split (no similarity overlap):** `--train-test-split 80`, `--train-test-similarity-threshold 90`; outputs `*_train.fasta` and `*_test.fasta`.
- **Balanced mode:** `--balanced` uses the same read count per genome (min across files).
- **Logging:** Pipeline writes `output-dir/logs/pipeline.log` (accessions, BLASTN stats, similarity filtering).

Full argument lists: `metagenome-generator <command> --help`.

---

## Project structure

```
MetagenomeGenerator/
├── pyproject.toml
├── README.md
├── LICENSE
├── environment.yml        # Optional conda env with BLAST+
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
│   └── seeker_wrapper.py
├── snapshots/              # Default snapshot output (date-stamped JSON)
└── working_directory/      # Suggested location for runs (gitignored)
```

After `pip install -e .`, the `metagenome-generator` console script is available. Programmatic use: `from metagenome_generator import build_metagenome, download_genomes, load_accessions, validate_genome_dir`.

---

## Notes

- NCBI rate limits apply; requests are delayed and retried (up to 3 times with backoff).
- Genome selection uses RefSeq and length filters; edit `DEFAULT_QUERIES` in `ncbi_search.py` to change criteria.
- For local runs and outputs, use a dedicated working directory; see `working_directory/README.md`.
