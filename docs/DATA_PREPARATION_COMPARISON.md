# Data preparation comparison: MetagenomeGenerator vs. published methods

This document compares the **data preparation** (first-step) practices used in key publications for viral/metagenomic classifier training with what the current **MetagenomeGenerator** tool supports. It is based on the methods sections and benchmarks of the cited papers and the tool’s README and codebase.

---

## 1. Reference publications and their data-prep steps

### 1.1 DeepVirFinder (Ren et al., 2020) — [PubMed 34084563](https://pubmed.ncbi.nlm.nih.gov/34084563/)

**Title:** *Identifying viruses from metagenomic data using deep learning*

| Step | What they did | MetagenomeGenerator |
|------|----------------|----------------------|
| **Reference source** | Viral RefSeq (NCBI) | ✅ Same: NCBI RefSeq via `download` / `snapshot` |
| **Temporal train/test split** | Train: viral RefSeq **discovered before May 2015**; test: **after May 2015** (temporal holdout for “novel” viruses) | ❌ **Not built-in.** Snapshot stores a run date; no filtering of accessions by **sequence release/date** to define train vs test |
| **Negative class** | Prokaryotic (host) genomes | ✅ Bacteria (+ optional archaea, plasmid) as negatives |
| **Extra viral training data** | Enlarged with **millions of purified viral sequences from metavirome samples** to improve under-represented virus groups | ❌ **Not supported.** Only RefSeq viral genomes; no import of external metavirome/contig FASTA |
| **Sequence length evaluation** | Evaluated at **300, 500, 1000, 3000 bp** (stratified by length) | ⚠️ **Partial.** Fixed length (e.g. 250) or variable range (e.g. 300–2000); no built-in **multi-length benchmark sets** (e.g. separate 300/500/1000/3000 bp outputs) |
| **Mutation robustness** | Tested with **mutations at different rates** (Fig. 2C) on sequences | ❌ **Not supported.** No option to introduce substitutions/errors into chunks for robustness checks |
| **Ambiguous bases** | Not explicitly discussed | ✅ **Supported:** `--forbid-ambiguous` to exclude sequences with N (and other non-ACGT) |

---

### 1.2 VIBRANT (Kieft et al., 2020) — [PubMed 32522236](https://pubmed.ncbi.nlm.nih.gov/32522236/)

**Title:** *VIBRANT: automated recovery, annotation and curation of microbial viruses*

| Step | What they did | MetagenomeGenerator |
|------|----------------|----------------------|
| **Reference source** | Reference virus datasets + **microbiome/virome data** (not only RefSeq) | ⚠️ **RefSeq only.** No built-in use of virome/metavirome assemblies or external viral contig sets |
| **Negative class** | Bacterial/archaeal genomes, **plasmids** (to separate from true viruses) | ✅ Bacteria, archaea, plasmid (all optional) |
| **Identification method** | Protein-signature ML (KEGG, Pfam, VOG) + neural networks; **annotation-based** | ❌ **Out of scope.** Tool is **sequence-only** (no gene finding, Pfam/KEGG/VOG, or annotation pipeline) |
| **Provirus handling** | **Extract integrated provirus regions** from host scaffolds (v-score, trimming to viral annotations) | ⚠️ **Different strategy.** We **exclude** EVE regions (BLASTN non-viral vs viral → drop chunks overlapping hits); we do **not** extract/produce provirus sequences |
| **Genome quality** | Viral scaffold **quality/completeness** categories (complete, high/medium/low draft) | ❌ **Not implemented.** No filtering or labeling by completeness/quality |
| **Benchmark construction** | Compared to VirSorter, VirFinder, MARVEL on **artificial scaffolds 3–15 kb**; reference viruses, plasmids, bacterial/archaeal genomes | ✅ Compatible: we can produce fixed or variable-length chunks in that range; no built-in “benchmark recipe” (e.g. 50 viral + 50 non-viral per group) |

---

### 1.3 VirFinder (Ren et al., 2017) — [Microbiome 5:69](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-017-0283-5)

| Step | What they did | MetagenomeGenerator |
|------|----------------|----------------------|
| **Temporal split** | Train: host + viral genomes **before 1 Jan 2014**; test: **after** that date | ❌ **Not built-in.** Same gap as DeepVirFinder: no date-based train/test partition of accessions |
| **Contig lengths** | **1 kb, 3 kb, 5 kb** (subsampled from complete genomes) | ✅ Fixed length (e.g. 1000, 3000, 5000) or variable range |
| **Simulated metagenomes** | **Simulated human gut metagenomes** (assembled) for evaluation | ❌ **Not built-in.** We simulate **reads/contigs** from reference genomes, not full metagenome assembly simulation |
| **Negatives** | Host (prokaryotic) genomes | ✅ Bacteria (and optionally archaea, plasmid) |

---

### 1.4 VirSorter2 (Guo et al., 2021) — [Microbiome; PubMed 33522966](https://pubmed.ncbi.nlm.nih.gov/33522966/)

| Step | What they did | MetagenomeGenerator |
|------|----------------|----------------------|
| **Benchmark datasets** | **dsDNA phage** benchmark; **non-dsDNA phage**; **non-RefSeq viral**; **plasmids**; Pfam HMM sets (Bacteria, Archaea, Eukaryota, Viruses) | ⚠️ **Partial.** We have bacterial, viral, archaeal, plasmid RefSeq; no eukaryote/fungal/protozoan negatives; no explicit “dsDNA-only” or “non-RefSeq” viral sets |
| **Validation design** | **Equal numbers** (e.g. 50 viral, 50 non-viral) from archaea, bacteria, **fungi**, **protozoa**, plasmids; **multiple fragment lengths**; **5 replicates × 100 sequences** | ⚠️ **Partial.** `--balanced` gives equal reads **per genome file**; no built-in “50 per category” or “N replicates” benchmark generator |
| **Proviruses** | Proviruses from **microbial RefSeq** in benchmark | ⚠️ We exclude EVE regions from non-viral chunks; we do not output provirus sequences as a separate class |
| **Eukaryote negatives** | Fungi, protozoa used as negative controls | ❌ **Not supported.** Only bacteria, archaea, plasmid; no NCBI query or folder for fungi/protozoa/eukaryotes |

---

### 1.5 Other relevant practices (VirHunter, Virus2Vec, etc.)

| Practice | Description | MetagenomeGenerator |
|----------|-------------|----------------------|
| **Leave-out / novel virus validation** | Train on one set of viruses, test on “novel” (e.g. different time slice or source) | ⚠️ Possible only if user provides **two accession sets** (e.g. two snapshots or two JSON files) and runs download/chunk separately; no **temporal split by accession date** built in |
| **Same-length train/test** | Using same fragment length for train and test can improve generalisation (e.g. Virus2Vec / host prediction) | ✅ User can set one `--sequence-length` or one variable range for both |
| **Under-represented virus groups** | DeepVirFinder added metavirome data to improve under-represented groups (e.g. Crenarchaeota, Bacteroidetes viruses) | ❌ No taxonomy-aware balancing or metavirome import |

---

## 2. Summary: what the current tool already covers

- **RefSeq-based download:** Bacteria, virus, archaea, plasmid with configurable counts and length filters.
- **Reproducibility:** Accession **snapshot** (date-stamped JSON), `--accessions-file`, `--save-accessions` for fixed dataset reruns.
- **Chunking:** Fixed length or variable-length contigs (e.g. 300–2000 bp), balanced or capped reads per genome.
- **EVE handling:** BLASTN non-viral vs viral → **exclude** chunks overlapping hits (`blastn-filter` + `--eve-intervals`).
- **Similarity control:** `--filter-similar` (e.g. 90% identity) and **train/test split** with removal of test sequences similar to train.
- **Ambiguous bases:** `--forbid-ambiguous` to exclude N (and other non-ACGT).
- **Seeker integration:** Run classifier on built metagenome and export predicted phage reads.

---

## 3. Gaps vs. literature (possible extensions)

| Gap | Papers | Possible extension |
|-----|--------|--------------------|
| **Temporal train/test split** | DeepVirFinder, VirFinder | Add **accession date** (or “before/after date”) from NCBI (e.g. via esummary or release dates) and a mode like `--temporal-split 2015-05-01` to assign accessions to train vs test by date. |
| **Mutation / error simulation** | DeepVirFinder (mutation rates) | Option to **introduce substitutions** (and optionally indels) at a given rate into chunks (e.g. for robustness or error-tolerant benchmarks). |
| **Eukaryote/fungal/protozoan negatives** | VirSorter2 | Add NCBI queries (and optional folder layout) for **eukaryotes/fungi/protozoa** as negative classes. |
| **Metavirome / external viral contigs** | DeepVirFinder | Allow **user-provided FASTA** (e.g. metavirome contigs) as extra viral positives, merged with RefSeq viral chunks. |
| **Multi-length benchmark sets** | DeepVirFinder, VirFinder | Mode to output **multiple FASTA** (e.g. 300, 500, 1000, 3000 bp) in one run for length-stratified evaluation. |
| **Structured benchmark recipe** | VirSorter2 | Preset or script: **equal counts per category** (e.g. 50 viral, 50 bacterial, 50 plasmid) and **replicates** (e.g. 5× 100 sequences) for tool comparison. |
| **Provirus extraction** | VIBRANT | Optional step to **output** provirus (EVE) sequences as a separate set instead of only excluding them from non-viral chunks. |
| **Genome quality / completeness** | VIBRANT | Use NCBI/completeness metadata or heuristics to **filter or label** genomes (e.g. “complete” vs “draft”) and optionally restrict to high-quality only. |
| **Taxonomy-aware balancing** | DeepVirFinder | Option to **balance** viral training data by taxonomy (e.g. realm or family) or to oversample under-represented groups. |

---

## 4. What to do next: prioritized improvement plan

Based on the gaps above, here is a concrete, prioritized list of improvements for the current tool.

### High priority (strong alignment with DeepVirFinder / VirFinder / VirSorter2)

| # | Improvement | What to do |
|---|-------------|------------|
| 1 | **Temporal train/test split** | Fetch **release/update date** per accession (NCBI `esummary` or `efetch` with `datetype=edat`). Add `--temporal-split DATE` to `download` or a post-process: assign accessions to train vs test by date; optionally output two accession JSONs or two FASTA sets. Document that snapshot + date filter enables “train before DATE, test after DATE” as in the papers. |
| 2 | **Mutation / error simulation** | In `chunk_genomes`: add optional **substitution rate** (e.g. `--substitution-rate 0.01`) and optionally **indel rate**; for each chunk, randomly substitute bases (and optionally insert/delete) with a given rate. Use a seed for reproducibility. Enables robustness benchmarks like DeepVirFinder. |
| 3 | **External viral/metavirome FASTA** | Allow **user-provided viral FASTA** as extra positives: e.g. `--extra-viral-fasta path/to/metavirome.fasta`. Concatenate or merge with RefSeq viral chunks in the same run (same length/balance logic). No new NCBI logic; just an extra input path. |
| 4 | **Eukaryote/fungal/protozoan negatives** | Add NCBI queries in `ncbi_search.py` for **eukaryotes** (or fungi, protozoa) with RefSeq + length filters; add `eukaryote` (or `fungi`, `protozoa`) to download/snapshot layout (`genome_layout.py`: dir + prefix). Optional `--num-eukaryote` (or similar) in CLI. Matches VirSorter2-style negatives. |

### Medium priority (benchmarks and usability)

| # | Improvement | What to do |
|---|-------------|------------|
| 5 | **Multi-length benchmark output** | Add a mode (e.g. `chunk --multi-length 300,500,1000,3000`) that writes **one FASTA per length** (e.g. `metagenome_300.fasta`, `metagenome_500.fasta`, …) from the same genome set, so users get length-stratified evaluation without multiple runs. |
| 6 | **Structured benchmark recipe** | Add a subcommand or preset (e.g. `benchmark`) that: fixed number of sequences **per category** (e.g. 50 viral, 50 bacterial, 50 plasmid), optional **replicates** (e.g. 5), and optional fragment lengths. Output: e.g. `benchmark_replicate1_500bp.fasta`. Can be implemented as a small script that calls existing `download` + `chunk` with computed `--reads-per-organism` / `--cap-total-reads`. |
| 7 | **Provirus (EVE) extraction output** | Reuse BLASTN hits from `blastn-filter`: instead of only excluding EVE regions from non-viral chunks, **write** the overlapping intervals as separate sequences (e.g. one FASTA of “provirus” contigs). Optional `--export-eve-fasta` to `blastn-filter` or chunk. |

### Lower priority (optional enhancements)

| # | Improvement | What to do |
|---|-------------|------------|
| 8 | **Genome quality / completeness** | If NCBI provides completeness or “complete genome” in metadata (e.g. via `esummary`), add an optional filter: e.g. `--complete-only` or `--min-completeness`. Otherwise document that users can pre-filter accessions in the snapshot JSON. |
| 9 | **Taxonomy-aware balancing** | For viral RefSeq: fetch taxonomy (e.g. family or realm) per accession; add option to **balance** or **oversample** by taxon (e.g. equal reads per family). Requires taxonomy API or parsing from NCBI. |

### Out of scope (for this tool)

- **Full metagenome assembly simulation** (VirFinder-style “simulated human gut metagenome”): heavy and different goal; document as out of scope.
- **Protein/annotation-based features** (VIBRANT’s KEGG/Pfam/VOG): different pipeline; keep this tool sequence-only.

### Suggested implementation order

1. **#3 External viral FASTA** — Small change, no NCBI date API; unblocks metavirome use.
2. **#2 Mutation simulation** — Self-contained in `chunk_genomes`; high value for benchmarks.
3. **#1 Temporal split** — Requires NCBI date fetch and split logic; high impact for method alignment.
4. **#4 Eukaryote negatives** — New query + layout; moderate effort.
5. **#5 Multi-length output** — Convenience in chunk step.
6. **#6 Benchmark recipe** — Script or subcommand on top of existing commands.
7. **#7 Provirus extraction** — Extend `blastn_filter` and possibly chunk output.

---

## 5. References

1. Ren J, Song K, Deng C, et al. Identifying viruses from metagenomic data using deep learning. *Quant Biol*. 2020;8(1):64-77. [PubMed 34084563](https://pubmed.ncbi.nlm.nih.gov/34084563/); [PMC8172088](https://pmc.ncbi.nlm.nih.gov/articles/PMC8172088/).
2. Kieft K, Zhou Z, Anantharaman K. VIBRANT: automated recovery, annotation and curation of microbial viruses, and evaluation of viral community function from genomic sequences. *Microbiome*. 2020;8(1):90. [PubMed 32522236](https://pubmed.ncbi.nlm.nih.gov/32522236/).
3. Ren J, Ahlgren NA, Lu YY, Fuhrman JA, Sun F. VirFinder: a novel *k*-mer based tool for identifying viral sequences from assembled metagenomic data. *Microbiome*. 2017;5:69.
4. Guo J, Bolduc B, Zayed AA, et al. VirSorter2: a multi-classifier, expert-guided approach to detect diverse DNA and RNA viruses. *Microbiome*. 2021. [PubMed 33522966](https://pubmed.ncbi.nlm.nih.gov/33522966/).

---

*This comparison focuses on **data preparation** (sources, splitting, chunking, negatives, EVE handling). It does not compare classifier architectures or downstream analysis (e.g. VIBRANT’s AMG/functional annotation).*
