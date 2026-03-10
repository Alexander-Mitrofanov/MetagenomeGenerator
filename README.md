# Metagenome Generator

Download bacterial and viral genomes from NCBI and chunk them into fixed-length reads to build simulated metagenome FASTAs for training classifiers.

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

NCBI Entrez requires an email for the download/snapshot steps:

```bash
export ENTREZ_EMAIL="your_email@example.com"
# Optional, for higher rate limits:
export ENTREZ_API_KEY="your_ncbi_api_key"
```

## Usage

After installation, use the `metagenome-generator` command (or `python -m metagenome_generator`). From the repo root you can also run `python main.py`. Subcommands: `download`, `snapshot`, `chunk`, `pipeline`, `blastn-filter`, `seeker`, `temporal-split`.

**Keep the project clean:** use a dedicated **working directory** for all outputs and run data (e.g. `working_directory/`). It is gitignored. See [working_directory/README.md](working_directory/README.md) for layout and examples.

```bash
mkdir -p working_directory/snapshots
metagenome-generator --help
metagenome-generator download --num-organisms 10 --output-dir working_directory
metagenome-generator snapshot --output working_directory/snapshots/accession_snapshot_$(date +%Y-%m-%d).json
metagenome-generator chunk --input working_directory/downloaded --output metagenome.fasta --output-dir working_directory/metagenome
```

### Accession snapshot (no downloads)

Catalog **all** bacterial, viral, archaeal, and plasmid accession IDs that match the project’s RefSeq + complete-genome criteria—**without downloading any sequences**. Only presence in the NCBI database is recorded. Output is written to the **`snapshots/`** folder. By default the filename includes the run date: `snapshots/accession_snapshot_YYYY-MM-DD.json`. You can then subset or edit that file and pass it to `download` or `pipeline` via `--accessions-file`.

The tool uses NCBI’s History server and paging (10,000 UIDs per request) to retrieve full result sets. A **local UTC timestamp** is stored; NCBI does not provide a database “as of” date. Optionally, nucleotide DB metadata (e.g. record count) is fetched via `einfo`. By default, **per-accession metadata** is also fetched via NCBI `esummary`: **CreateDate** (YYYY/MM/DD) and **Title** (genome description, i.e. the FASTA header text). These are stored under the `accession_metadata` key in the snapshot JSON. `temporal-split` can then use this data and skip a second NCBI round-trip. Use `--no-metadata` for a faster snapshot if you do not need dates or headers.

```bash
metagenome-generator snapshot
# or: metagenome-generator snapshot --output snapshots/my_snapshot.json
# Skip per-accession metadata (faster): metagenome-generator snapshot --no-metadata
```

**Snapshot JSON format:** In addition to `timestamp`, `bacterial`, `viral`, `archaea`, `plasmid`, the file may contain `accession_metadata`: a dict mapping each accession to `{"create_date": "YYYY/MM/DD", "title": "..."}` (NCBI CreateDate and genome title/header).

| Argument | Description | Default |
|----------|-------------|---------|
| `--output` | Output JSON path | `snapshots/accession_snapshot_YYYY-MM-DD.json` (run date) |
| `--no-db-info` | Do not fetch NCBI nucleotide db metadata (einfo) | off |
| `--no-metadata` | Do not fetch CreateDate and title per accession (faster snapshot) | off |
| `--metadata-batch-size` | esummary batch size when fetching metadata | 500 |

**Typical workflow:** Run `snapshot` once (or periodically). Each run creates a date-stamped file (e.g. `snapshots/accession_snapshot_2026-03-10.json`). Edit the JSON to keep a subset of accessions if needed, then run `download --accessions-file snapshots/accession_snapshot_2026-03-10.json` (or the same with `pipeline`) to build from that catalog.

### Temporal train/test split (by NCBI CreateDate)

To align with **DeepVirFinder**- or **VirFinder**-style evaluation (train on sequences “discovered” before a date, test on those added on or after), use the `temporal-split` subcommand. It reads an accessions JSON (e.g. from `snapshot`). If the file contains **accession_metadata** (CreateDate and title from `snapshot`), it uses that and skips NCBI; otherwise it fetches CreateDate via `esummary`. It writes two JSONs: **train** (CreateDate &lt; split-date) and **test** (CreateDate ≥ split-date). Both files are valid for `download --accessions-file`.

```bash
metagenome-generator temporal-split \
  --accessions-file snapshots/accession_snapshot_2026-03-10.json \
  --split-date 2015-05-01
```

By default this writes `train_<basename>.json` and `test_<basename>.json` next to the input file. Use `--output-train` and `--output-test` to set paths. Then run download (or pipeline) twice with the train and test JSONs to build separate train/test genome sets.

| Argument | Description |
|----------|-------------|
| `--accessions-file` | Input JSON with bacterial, viral, archaea, plasmid lists (e.g. from `snapshot`) |
| `--split-date` | Cutoff date **YYYY-MM-DD** (e.g. `2015-05-01`). Train = CreateDate &lt; date, test = ≥ date |
| `--output-train` | Path for train accessions JSON (default: same dir as input, `train_<basename>.json`) |
| `--output-test` | Path for test accessions JSON (default: same dir as input, `test_<basename>.json`) |
| `--batch-size` | NCBI esummary batch size (default: 200) |

Requires `ENTREZ_EMAIL` (and optionally `ENTREZ_API_KEY`) for NCBI rate limits.

### Download genomes only

Saves one FASTA per genome in the chosen output directory.

```bash
metagenome-generator download \
  --num-organisms 10 \
  --output-dir output
```

| Argument | Description | Default |
|----------|-------------|---------|
| `--num-organisms` | Number of organisms to download **per group** (bacterial and viral) | 10 |
| `--output-dir` | Root for downloads; genomes are written into `bacteria/`, `virus/`, `archaea/`, `plasmid/` subfolders | `output` |
| `--num-archaea` | Number of archaeal genomes to download (negative samples) | 0 |
| `--num-plasmid` | Number of plasmids to download (negative samples) | 0 |
| `--accessions-file` | Load accession IDs from a JSON file (skip NCBI search). Use for reproducible runs; file may include a `timestamp` of when the list was obtained | off |
| `--save-accessions` | After searching NCBI, save the accession list and current UTC **timestamp** to this JSON path. Use this file later with `--accessions-file` to re-download the same set | off |

**Reproducibility (NCBI dataset):** NCBI search results change over time. You can fix the dataset in two ways:

1. **Snapshot only (no downloads)** — use the `snapshot` subcommand to write a full catalog to `snapshots/accession_snapshot_YYYY-MM-DD.json`. Edit to subset if needed, then run download/pipeline with `--accessions-file snapshots/accession_snapshot_YYYY-MM-DD.json`.
2. **Save during download** — run download with `--save-accessions out/accessions.json`; later use `--accessions-file out/accessions.json` to re-download the same set. The tool prints the timestamp from the file when using `--accessions-file`.

### Chunk genomes into metagenome reads

Reads the genome FASTAs, cuts them into fixed-length reads, and writes one metagenome FASTA.

```bash
metagenome-generator chunk \
  --input output \
  --output metagenome_250nt.fasta \
  --output-dir output \
  --sequence-length 250 \
  --reads-per-organism 1000
```

| Argument | Description | Default |
|----------|-------------|---------|
| `--input` | Input FASTA file or **directory** of genome FASTAs (e.g. from download step) | *(required)* |
| `--output` | **Filename** for the metagenome FASTA | *(required)* |
| `--output-dir` | Directory to write the output file | `output` |
| `--sequence-length` | Length of each read in nucleotides | 250 |
| `--reads-per-organism` | Max reads to take per organism (per input file). Omit for all reads | all |
| `--balanced` | Use equal reads per file: set reads-per-organism to the **minimum** max reads across all files (so bacteria and viruses contribute the same number of reads) | off |
| `--eve-intervals` | Path to `eve_intervals.json` from `blastn-filter`; chunks overlapping EVE regions are **excluded** from output | off |
| `--filter-similar` | Drop chunks that are ≥90% similar (BLASTN) to any already-kept chunk; oversample then refill to reach target; log warnings if target not reached | off |
| `--train-test-split` | Split into train (PCT%%) and test; remove from test sequences similar to train | off |
| `--train-test-similarity-threshold` | BLASTN percent identity above which test sequences are removed (vs train) | 90 |
| `--similarity-threshold` | Max BLASTN percent identity for "similar" (default 90) | 90 |
| `--similarity-min-coverage` | Min fraction of query length in alignment to count as similar | 0.8 |
| `--oversample-factor` | When `--filter-similar`: generate up to this many times target before filtering | 2.0 |
| `--forbid-ambiguous` | Exclude chunks that contain non-ACGT characters (e.g. N) | off |

### BLASTN filtering (EVE removal)

**Goal:** Endogenous viral elements (EVEs) in non-viral genomes (bacteria, archaea, plasmid) can be misclassified as viral. We align non-viral sequences against viral sequences with BLASTN, then **exclude** any chunk that overlaps a hit when building the training metagenome. Viral genomes are not BLAST’ed; only non-viral inputs are filtered.

**When it runs:** After download, before or during chunking. Either run the `blastn-filter` subcommand once, then pass `--eve-intervals` to `chunk`, or use `pipeline --run-blastn-filter` so the pipeline runs BLASTN after download and passes EVE data into the chunk step.

**How results are used:** BLASTN tabular output is parsed to get (query id, query start, query end) per hit. Intervals are merged per sequence. When chunking, any chunk whose [start, end) overlaps one of these intervals is **skipped** (not written). Sequences with no hits are chunked as usual.

**Requirements:** BLAST+ must be installed (`conda install -c bioconda blast` or use the project env from `environment.yml`).

**Standalone step (then chunk with EVE exclusion):**

```bash
metagenome-generator blastn-filter --genome-dir output/downloaded --out-dir output/blastn --evalue 1e-5 --perc-identity 70
# Writes output/blastn/eve_intervals.json and output/blastn/blastn/*.blastn.tsv

metagenome-generator chunk --input output/downloaded --output negative.fasta --output-dir output/metagenome \
  --balanced --eve-intervals output/blastn/eve_intervals.json
```

**Integrated in pipeline:** (BLASTN output goes to `output-dir/blastn/`)

```bash
metagenome-generator pipeline --num-organisms 5 --output-dir output \
  --balanced --run-blastn-filter --blastn-evalue 1e-5 --blastn-perc-identity 70
```

| blastn-filter argument | Description | Default |
|------------------------|-------------|---------|
| `--genome-dir` | Directory with `bacterial_*.fasta`, `viral_*.fasta`, and optionally `archaea_*.fasta`, `plasmid_*.fasta` | *(required)* |
| `--out-dir` | Where to write `viral_concat.fasta`, `blastn/`, and `eve_intervals.json` | *(required)* |
| `--evalue` | BLASTN E-value threshold | 1e-5 |
| `--perc-identity` | BLASTN percent identity threshold | 70 |

### Full pipeline: download + chunk (+ optional BLASTN, Seeker)

The pipeline uses an **organized output layout** so results are easy to navigate. One root directory (e.g. `output`) contains step-based subfolders:

| Subfolder | Contents |
|-----------|----------|
| `downloaded/` | Raw genome FASTAs from NCBI in category subfolders: `bacteria/`, `virus/`, `archaea/`, `plasmid/` |
| `blastn/` | BLASTN outputs when `--run-blastn-filter`: viral DB, TSVs, `eve_intervals.json` |
| `metagenome/` | Final chunked metagenome FASTA (e.g. `metagenome.fasta`) |
| `seeker/` | Seeker predictions and phage FASTA when `--run-seeker` is used |
| `logs/` | Pipeline log file (`pipeline.log`): download accessions/origins, BLASTN EVE stats, similarity filtering, warnings |

Example:

```bash
metagenome-generator pipeline \
  --num-organisms 10 \
  --output-dir output \
  --output metagenome.fasta \
  --sequence-length 250 \
  --reads-per-organism 1000
```

This creates `output/downloaded/`, `output/metagenome/metagenome.fasta`, and (if used) `output/blastn/`, `output/seeker/`.

| Argument | Description | Default |
|----------|-------------|---------|
| `--num-organisms` | Number of organisms to download **per group** (bacterial and viral) | 10 |
| `--output-dir` | **Root** output directory; pipeline creates `downloaded/`, `blastn/`, `metagenome/`, `seeker/` under it | `output` |
| `--output` | **Filename** for the metagenome FASTA (written under `output-dir/metagenome/`) | `metagenome.fasta` |
| `--sequence-length` | Length of each read in nucleotides | 250 |
| `--reads-per-organism` | Max reads to take per organism (per input file) | 1000 |
| `--balanced` | Use equal reads per file (min of max possible reads across files) | off |
| `--filter-similar` | Filter out chunks ≥90% similar to already-kept; oversample and refill; log warnings if target not reached | off |
| `--similarity-threshold` | Max BLASTN percent identity for "similar" | 90 |
| `--oversample-factor` | When `--filter-similar`: generate up to this many times target before filtering | 2.0 |
| `--train-test-split` | Train percentage (e.g. 80 → 80% train, 20% test); test sequences similar to train are removed | off |
| `--train-test-similarity-threshold` | Remove from test if ≥ this % identity to train | 90 |
| `--accessions-file` | For download step: load accession IDs from JSON (reproducibility; skip NCBI search) | off |
| `--save-accessions` | For download step: save accession list and UTC timestamp to JSON after searching | off |
| `--run-blastn-filter` | Run BLASTN (non-viral vs viral) and exclude EVE regions when chunking | off |
| `--run-seeker` | Run Seeker (phage/bacteria prediction) on the metagenome FASTA after chunking | off |
| `--forbid-ambiguous` | Exclude chunks containing non-ACGT characters (e.g. N) | off |

### Balanced mode

Bacterial genomes are much longer than viral ones, so without a cap you get many more bacterial than viral reads. Use **`--balanced`** so each genome file contributes the same number of reads: the tool scans each file’s length, computes max possible reads per file (`total_bases // sequence_length`), and uses the **minimum** of those as the cap for every file.

```bash
metagenome-generator chunk --input output/downloaded --output metagenome.fasta --output-dir output/metagenome --sequence-length 250 --balanced
```

You’ll see a short report (e.g. `bacterial_1: 3743357 bp -> max 14973 reads (250 nt)`) and the chosen cap (e.g. `Balanced: using 34 reads per file`).

### Example with custom paths

```bash
metagenome-generator pipeline \
  --num-organisms 5 \
  --output-dir data \
  --output metagenome.fasta \
  --sequence-length 250 \
  --reads-per-organism 500
```

Results go under `data/downloaded/`, `data/metagenome/metagenome.fasta`, etc.

Output headers look like `bacterial_1_chunk_0`, `viral_2_chunk_3`, etc.

### Seeker: phage/bacteria prediction

After building a metagenome you can run [Seeker](https://github.com/gussow/seeker) to label each read as phage or bacteria and export predicted phage reads. Seeker must be installed in a separate conda env (e.g. `seeker`, Python 3.7). Use the `seeker` subcommand or add `--run-seeker` to the pipeline.

```bash
metagenome-generator seeker --input output/metagenome/metagenome.fasta --output-dir output/seeker --conda-env seeker
```

| Argument | Description | Default |
|----------|-------------|---------|
| `--input` | Metagenome FASTA (e.g. from chunk or pipeline) | *(required)* |
| `--output-dir` | Directory for outputs | `output` |
| `--threshold` | Score threshold for phage (0–1) | 0.5 |
| `--conda-env` | Conda env for `predict-metagenome` | `seeker` |
| `--min-length` | Drop reads shorter than this (Seeker needs ≥200 nt) | 200 |

**Output files**

| File | When created | Meaning |
|------|----------------|--------|
| `*.seeker_predictions.tsv` | Always | Per-read predictions: `name`, `prediction` (Phage/Bacteria), `score` (0–1). |
| `*.seeker_phage.fasta` | Always | Reads with score ≥ threshold (predicted phage). |
| `*.seeker_filtered_min200.fasta` | Only if some reads were &lt;200 nt | Length-filtered input actually sent to Seeker. Omitted when all reads are already ≥200 nt. |

### Variable-length contigs (training-style)

To simulate variable-length metagenomic contigs (e.g. 300–2000 bp, uniform distribution), use `--min-contig-length` and `--max-contig-length` instead of a fixed `--sequence-length`:

```bash
metagenome-generator chunk --input output/downloaded --output metagenome_var.fasta --output-dir output/metagenome \
  --min-contig-length 300 --max-contig-length 2000 --balanced --seed 42
```

Use `--seed` for reproducible sampling. Balanced mode uses the mean contig length to estimate max reads per file.

### Down-sampling negative to match positive

For training data, non-viral reads are often down-sampled to match the virus set size. Chunk the positive (viral) set first, then chunk the negative set with `--cap-total-reads` set to the viral read count:

```bash
metagenome-generator chunk --input output/downloaded --output viral.fasta --output-dir output/metagenome --balanced --seed 42
# Get viral read count: grep -c ">" output/metagenome/viral.fasta
metagenome-generator chunk --input output/downloaded --output negative.fasta --output-dir output/metagenome \
  --balanced --cap-total-reads 500 --seed 42
```

Use the viral read count (e.g. `grep -c ">" output/metagenome/viral.fasta`) as `--cap-total-reads` for the negative run.

### Similarity filtering within the dataset

To avoid redundant training data, chunks that are too similar (e.g. ≥90% identity over most of the sequence) can be removed. The tool generates more chunks than needed (oversample), then keeps only sequences that are not ≥90% similar to any already-kept sequence (BLASTN), and refills by generating more chunks (different seed) until the target count is reached or max rounds.

```bash
metagenome-generator chunk --input output/downloaded --output metagenome.fasta --output-dir output/metagenome \
  --balanced --filter-similar --similarity-threshold 90 --oversample-factor 2 --seed 42
```

Warnings are written to the log (and stderr) if the dataset cannot be fully created after filtering (e.g. not enough unique sequences).

### Train-test split

The final metagenome can be split into train and test sets. The split percentage is the **train** share (e.g. 80 → 80% train, 20% test). After the split, any test sequence that is too similar to any train sequence (BLASTN, threshold set by `--train-test-similarity-threshold`) is **removed from the test set** so train and test stay non-overlapping by similarity.

```bash
metagenome-generator chunk --input output/downloaded --output metagenome.fasta --output-dir output/metagenome \
  --balanced --train-test-split 80 --train-test-similarity-threshold 90 --seed 42
```

Output: `metagenome_train.fasta` and `metagenome_test.fasta` (test may have fewer sequences after similarity removal).

Pipeline:

```bash
metagenome-generator pipeline --output-dir output --balanced --train-test-split 80 --train-test-similarity-threshold 90 --seed 42
```

| Flag | Description | Default |
|------|-------------|---------|
| `--train-test-split` | Train percentage (0–100); rest is test | off |
| `--train-test-similarity-threshold` | Remove from test if ≥ this % identity to any train sequence | 90 |
| `--train-test-blast-threads` | BLAST threads for test-vs-train similarity (speed) | 4 |
| `--train-test-blast-batch-size` | Test sequences per BLAST run (larger = fewer runs) | 2000 |

**What the train-test BLAST step does:** After splitting into train and test, the tool builds a BLAST database from the train set and runs BLAST (query = test sequences) to find any test sequence that is ≥ threshold similar to any train sequence; those are removed from the test set so train and test do not share near-identical sequences. **Optimizations** to keep this fast: (1) **megablast** task instead of default blastn (designed for high-identity, much faster), (2) **batched** test sequences (e.g. 2000 per run to reduce the number of BLAST invocations), (3) **multi-threaded** BLAST (`-num_threads`). You can tune `--train-test-blast-threads` and `--train-test-blast-batch-size` if needed.

### Logging

When running the **pipeline**, a log file is written to `output-dir/logs/pipeline.log`. It records:

- **Download step:** Accession numbers and origin (bacterial, viral, archaea, plasmid) for each downloaded genome; failures.
- **BLASTN filter:** Which files were BLASTed; number of sequences with EVE hits (excluded from chunking).
- **Chunk / similarity filter:** Generated vs removed vs kept counts; warnings when the target count could not be reached after similarity filtering.

Use this file to audit runs and debug issues (e.g. dataset too small after similarity filtering).

### Training-data design (methods alignment)

Task-specific training sets often combine:

- **Positive (virus):** RefSeq viral genomes; optionally by realm (Adnaviria, Duplodnaviria, Monodnaviria, Riboviria, etc.) or external sets (e.g. VirSorter2).
- **Negative (non-viral):** Prokaryotes (bacteria, archaea), plasmids, and optionally eukaryotes (fungi, protozoa, Insecta) as challenging negatives.

This pipeline supports:

| Method step | Supported how |
|-------------|----------------|
| Multiple negative groups | Download with `--num-archaea`, `--num-plasmid`; chunk from one dir (all prefixes) or run chunk per group. |
| Equal reads per group | `--balanced` (per-file cap). |
| Variable-length contigs (300–2000 bp) | `--min-contig-length 300 --max-contig-length 2000`. |
| Down-sample negative to viral size | Chunk viral, then chunk non-viral with `--cap-total-reads <viral_count> --seed`. |
| Train/test split (e.g. 9:1 per taxon) | Not built-in; split genome lists or output FASTAs by group (e.g. with a small script) before or after chunking. |
| BLASTN to remove endogenous viral elements | **Built-in:** `blastn-filter` subcommand (or `pipeline --run-blastn-filter`). Non-viral vs viral BLASTN; chunks overlapping hits are excluded via `chunk --eve-intervals`. |
| Similarity filtering (e.g. 90% max) | **Built-in:** `chunk --filter-similar` (or `pipeline --filter-similar`). Oversample, drop chunks similar to already-kept (BLASTN), refill; warnings in log if target not reached. |
| Train-test split (no similarity overlap) | **Built-in:** `chunk --train-test-split PCT` (or `pipeline --train-test-split PCT`). Split by percentage; remove from test any sequence ≥ threshold similar to train (`--train-test-similarity-threshold`). |
| Temporal train/test split (by NCBI CreateDate) | **Built-in:** `temporal-split` subcommand. Split accessions into train (CreateDate &lt; date) and test (≥ date); output two JSONs for use with `download --accessions-file`. Uses `accession_metadata` from snapshot when present. |

Realm-level viral taxonomy (RefSeq) can be approximated by refining the viral query in `ncbi_search.py` (e.g. by taxon name or using NCBI’s taxonomy filters).

## Project structure

Standard Python package layout for publishing:

```
MetagenomeGenerator/
├── pyproject.toml          # Build and metadata (PEP 517/621)
├── README.md
├── LICENSE
├── environment.yml        # Optional conda env with BLAST+
├── main.py                # Launcher (calls metagenome_generator.cli)
├── src/
│   └── metagenome_generator/
│       ├── __init__.py
│       ├── __main__.py     # python -m metagenome_generator
│       ├── cli.py          # CLI entry point
│       ├── download_genomes.py
│       ├── ncbi_search.py
│       ├── accession_snapshot.py
│       ├── chunk_genomes.py
│       ├── genome_layout.py
│       ├── blastn_filter.py
│       ├── similarity_filter.py
│       ├── temporal_split.py
│       └── seeker_wrapper.py
├── snapshots/              # Default output for snapshot (date-stamped JSON)
└── data/                   # Optional; user data
```

After `pip install -e .`, the `metagenome-generator` console script is available. Programmatic use: `from metagenome_generator import build_metagenome, download_genomes, load_accessions, validate_genome_dir`.

## Notes

- NCBI rate limits apply; the download and snapshot steps add short delays between requests. Failed fetches are retried up to 3 times with backoff.
- Genome selection uses RefSeq and length filters; edit `DEFAULT_QUERIES` in `ncbi_search.py` to change criteria (used by both `download` and `snapshot`).
- For local runs and outputs, use a dedicated working directory; see `working_directory/README.md`.
