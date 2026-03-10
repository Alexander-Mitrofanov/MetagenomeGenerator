# Improvement ideas (backlog)

Ideas from [DATA_PREPARATION_COMPARISON.md](DATA_PREPARATION_COMPARISON.md). Check off when done.

---

## High priority

- [x] **#1 Temporal train/test split** — Fetch NCBI CreateDate per accession (esummary); split accessions by cutoff date into train (before) / test (on or after); output two accession JSONs. Enables DeepVirFinder/VirFinder-style “train before DATE, test after DATE”. *(Done: `temporal-split` subcommand; see README.)*
- [ ] **#2 Mutation / error simulation** — In chunk: optional `--substitution-rate` (and optionally indel rate); randomly mutate bases in chunks (with seed). For robustness benchmarks.
- [ ] **#3 External viral/metavirome FASTA** — `--extra-viral-fasta PATH`: merge user-provided viral FASTA with RefSeq viral chunks in same run.
- [ ] **#4 Eukaryote/fungal/protozoan negatives** — New NCBI query + category (e.g. eukaryote/fungi); folder + prefix in layout; `--num-eukaryote` (or similar) in CLI.

## Medium priority

- [ ] **#5 Multi-length benchmark output** — `chunk --multi-length 300,500,1000,3000`: write one FASTA per length from same genome set.
- [ ] **#6 Structured benchmark recipe** — Subcommand or preset: fixed N per category (e.g. 50 viral, 50 bacterial), optional replicates (e.g. 5×).
- [ ] **#7 Provirus (EVE) extraction output** — From blastn-filter: optionally write provirus intervals as a separate FASTA (in addition to excluding from non-viral chunks).

## Lower priority

- [ ] **#8 Genome quality / completeness** — Optional filter from NCBI metadata (e.g. `--complete-only`) if available.
- [ ] **#9 Taxonomy-aware balancing** — Balance or oversample viral training data by taxon (e.g. family).

## Out of scope

- Full metagenome assembly simulation (VirFinder-style).
- Protein/annotation-based features (VIBRANT’s KEGG/Pfam/VOG); keep tool sequence-only.
