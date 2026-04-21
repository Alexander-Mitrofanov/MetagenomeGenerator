# Changelog

All notable changes to this project are documented in this file.

The format is loosely based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [0.9.0] — 2026-04-21

Combined **audit bugfix release + realism feature push**. Rolls up:

- the four A / B / C / D audit pushes (behaviour fixes, documentation
  clean-up, design-quirk hardening, and polish);
- the A1 realism feature push (paired-end, multi-length, coverage-depth
  model, 3rd-gen error models, library-prep artefacts, gold-standard
  taxonomy labels).

All changes are backward compatible: no existing CLI flag, argument name,
or default value was removed or repurposed, and every new behaviour is
opt-in.

**Full test suite:** 90 / 90 passing (69 pre-release + 21 A1 regression
tests). No pre-existing test was modified.

The five blocks below detail each push individually (format preserved
from the staged PR drafts for reviewer convenience).

### A1 feature push (realism & throughput)

Six additive features that close the biggest realism gaps identified in
`docs/improvements.md` ("A1 feature push" backlog). All features ship on
both the `chunk` and `pipeline` subcommands, are opt-in, and default to
pre-existing behaviour.

#### Added

- **(F1) Multi-length output — `--multi-length 300,500,1000,3000`.**
  One invocation now emits one FASTA/FASTQ per requested length, using
  `{stem}_L{N}{suffix}` naming (e.g. `meta_L300.fasta`, `meta_L500.fasta`).
  New module API: `parse_multi_length_spec` and
  `build_metagenome_multi_length` in `chunk_genomes.py`. Incompatible with
  variable-length contig flags (`--min-contig-length` /
  `--max-contig-length` / `--contig-quality-profile`) and with
  `--train-test-split` (each length would otherwise need its own split);
  CLI validation rejects those combinations with a clear error.
  Rationale: benchmarking models across short-read / long-read regimes
  currently required running the tool once per target length.

- **(F2) Paired-end output — `--paired`, `--insert-size`,
  `--insert-size-sd`.** Generates matched `{stem}_R1{suffix}` /
  `{stem}_R2{suffix}` files (FASTA or FASTQ), with R2 produced as the
  reverse complement of the 3' end of each fragment. Read IDs use the
  Illumina `/1` / `/2` mate tags so downstream tools (BWA, bowtie2,
  seqkit) treat the files as proper pairs. Fragment length is drawn from
  a truncated normal distribution; defaults are `insert_size = 3 × read
  length` and `insert_size_sd = 0.1 × insert_size`. New module API:
  `chunk_sequence_paired`, `_collect_pairs_for_file`,
  `_introduce_chimeras_paired`, `_introduce_pcr_duplicates_paired`,
  and the top-level `build_metagenome_paired`. Paired mode supports
  abundance profiles, error models, taxonomy embedding, chimeras, and
  PCR duplicates; it explicitly rejects `--multi-length`,
  `--train-test-split`, `--filter-similar`, `--eve-intervals`,
  `--run-blastn-filter`, and variable-length options (each would need
  pair-aware implementations that aren't in scope here).

- **(F3) Coverage-depth model — `--coverage`, `--coverage-cv`.** Derives
  per-genome read counts from a target depth
  (`ceil(total_bp × coverage / read_length)`) instead of the flat
  `--reads-per-organism` budget. When `--coverage-cv > 0`, per-genome
  coverage is drawn from a log-normal with mean ≈ `coverage` and
  coefficient of variation `coverage_cv` — so different genomes get
  different depths, matching real metagenomic sequencing where coverage
  is uneven. `--reads-per-organism` still acts as a hard cap. New helper:
  `_compute_coverage_read_limits`. Rationale: training on a flat
  reads-per-genome budget biases models toward spurious "every organism
  has the same depth" assumptions.

- **(F4) 3rd-generation error models — `--error-model
  {nanopore,pacbio-hifi,pacbio-clr}`.** Long-read patterns with
  independent substitution / insertion / deletion rates. Nanopore
  additionally inflates error rates inside homopolymer runs
  (`NANOPORE_HOMOPOLYMER_MULT = 2.5` for runs ≥ 3 bases), reproducing
  the well-known basecaller weakness. FASTQ output uses a flat Phred
  score derived from each model's total error rate
  (`_phred_from_rate`, `_ERROR_MODEL_FLAT_RATE`), in contrast to
  Illumina's position-dependent quality curve. Default rates are
  calibrated to public platform benchmarks:
  nanopore ≈ 5.5% total (2.5% sub / 1.5% ins / 1.5% del),
  pacbio-hifi ≈ 0.3% total (CCS-grade),
  pacbio-clr ≈ 12.3% total. New helper: `_apply_third_gen_errors`.

- **(F5) Library-prep artefacts — `--chimera-rate`,
  `--pcr-duplicate-rate`.** Two independent artefact injectors applied
  after similarity filtering / capping and before quality annotation.
  `_introduce_chimeras` replaces a fraction of records with two-parent
  chimeras (5' half of parent A + 3' half of parent B, id
  `chimera_{idx}`, header `chimera parents=A|B`); read length is
  preserved and parent Phred qualities are spliced so FASTQ remains
  valid. `_introduce_pcr_duplicates` appends bit-identical duplicates
  (id suffix `_dup`, header tag `pcr_duplicate=true`) with per-record
  probability `duplicate_rate`; if the source has Phred qualities the
  duplicate inherits them. Both are gated on `--seed` so runs remain
  reproducible. Rationale: real metagenomic libraries routinely carry
  5–20 % duplicates and a few-percent chimera rate; models that never
  see either overfit to clean synthetic reads.

- **(F6) Gold-standard taxonomy in headers — `--embed-taxonomy`.**
  Appends `tax=<group>` to every record description. For viral reads,
  `<group>` is looked up in the JSON passed to `--viral-taxonomy`
  (keyed by accession / file stem) and falls back to `unknown`; for
  bacteria / archaea / plasmid reads, `<group>` is the category name.
  New helper: `_embed_tax_in_description`. The viral-taxonomy JSON is
  loaded once per run and reused for both the existing
  `--balance-viral-by-taxonomy` sampler and the new label embedding.
  Rationale: downstream supervised training can now read the label
  directly from the FASTA/FASTQ header instead of re-joining against a
  separate manifest.

#### Tests

- `tests/test_a1_feature_push.py` — 21 new tests covering each feature:
  - F1: `parse_multi_length_spec` happy-path + rejection cases,
    `build_metagenome_multi_length` writes one file per length, rejects
    variable-length kwargs.
  - F2: `chunk_sequence_paired` yields matching pairs with correct
    reverse-complement R2; `build_metagenome_paired` writes `_R1`/`_R2`
    files with `/1`/`/2` mate tags; FASTQ output carries Phred
    qualities on both mates.
  - F3: Coverage scales read counts proportionally to genome size
    (large:small ≈ 2:1 on a 4000:2000 bp pair); `--coverage-cv > 0`
    produces per-genome read-count variability (not all equal).
  - F4: Homopolymer-heavy sequences accumulate more indels than
    random DNA under the nanopore model; PacBio HiFi produces <1% length
    delta across 10 kb; FASTQ Phred ordering hifi > illumina > nanopore.
    `_phred_from_rate` is monotonic.
  - F5: `_introduce_chimeras` replaces the exact expected count,
    zero-rate is a no-op, chimera records carry `chimera_` IDs;
    `_introduce_pcr_duplicates` grows the list by ≈ expected count,
    duplicates have `_dup` IDs and `pcr_duplicate=true` headers;
    end-to-end `build_metagenome` run with both rates produces both
    artefact kinds.
  - F6: Viral reads get their taxonomy group, unknown accessions fall
    back to `unknown`, non-viral reads get `tax=<category>`;
    `build_metagenome` end-to-end with `--embed-taxonomy` writes the
    tag on every record.

#### Verification

- Full test suite: **90/90 passing** (69 pre-push + 21 new).
- No regressions in the A/B/C/D bugfix suites.
- `tests/test_a1_feature_push.py` is lint-clean (`ruff check`).

### A-category audit bugfix push

#### Fixed

- **(A1) Multi-record FASTA truncation in `chunk` / `pipeline` / `benchmark-recipe`.**
  `_collect_chunks_for_file` contained a stray `break` at the end of its
  `for record in SeqIO.parse(...)` loop that caused it to consume only the
  very first record of any input FASTA. Multi-segment viral genomes (e.g.
  influenza A with 8 segments), multi-contig drafts, and user-supplied
  multi-record FASTAs in category directories silently lost every record
  past the first. The loop now iterates every record, shares the
  `reads_per_organism` budget across records ("per input file" semantics,
  as documented in the CLI help), and disambiguates chunk IDs for records
  past the first with a `_seg{N}` prefix so no two chunks share the same
  header. Single-record files keep their historical
  `{accession}_read_{idx}` / `{accession}_contig_{idx}` naming.
  Files: `src/metagenome_generator/chunk_genomes.py`.

- **(A2) First batch bypassed deduplication in `filter_by_similarity`.**
  When seeding the initially empty `kept` list, the first batch of
  candidate records was appended wholesale with no BLAST comparison,
  allowing near-duplicates from within the same first batch to both
  survive the filter. This violated the function's contract that no two
  kept sequences are ≥`similarity_threshold` similar. The first batch is
  now BLAST'd against itself (all-vs-all, `max_target_seqs` sized to the
  batch), a symmetric neighbor adjacency is built from the tabular hits,
  and a greedy pass keeps each record only if it is not similar to any
  earlier-kept sibling. A new helper `_parse_similar_neighbors` parses
  the pairwise output.
  Files: `src/metagenome_generator/similarity_filter.py`.

- **(A3) `benchmark-recipe` crashed on multi-record viral FASTAs.**
  `benchmark_recipe._write_concat_fasta_and_lengths` used `SeqIO.read`,
  which raises when a FASTA holds more than one record. Any cached
  NCBI viral accession with multiple segments (common for RNA viruses)
  would crash the recipe before the BLAST step. Switched to
  `SeqIO.parse` with a per-record length accumulator, and added an
  explicit error when a cached file contains zero records.
  Files: `src/metagenome_generator/benchmark_recipe.py`.

- **(A4) `filter-test-against-train` silently mis-parsed FASTQ inputs.**
  Both train and test inputs were hard-coded as `"fasta"`, so passing
  the FASTQ outputs produced by `--output-fastq` / `split-metagenome-train-test`
  yielded empty or corrupted records with no error. Input format is
  now auto-detected from the `.fastq`/`.fq` suffix (matching the
  convention used by `split-metagenome-train-test`), the output is
  written in the test input's format so Phred qualities round-trip,
  and the CLI default output name inherits the test input's suffix.
  Files: `src/metagenome_generator/similarity_filter.py`,
  `src/metagenome_generator/cli.py`.

- **(A5) `genome-pool materialize` could delete the shared pool.**
  `materialize_from_pool` called `shutil.rmtree(dest_dir)` whenever
  `clean_dest=True` (the default). If a user passed
  `--output-dir` equal to (or an ancestor of) `--pool-dir`, the
  entire pool — including `pool_manifest.json` — would be wiped
  without warning. Resolved paths are now compared: the call raises
  `ValueError` if `dest_dir` and `pool_dir` resolve to the same
  location, or if `pool_dir` lives inside `dest_dir` while
  `clean_dest=True`.
  Files: `src/metagenome_generator/genome_pool.py`.

- **(A6) `fetch_accession_dates` lacked retry logic (transient NCBI
  failures dropped whole batches).** The twin function
  `fetch_accession_metadata` already retried with back-off, but
  `fetch_accession_dates` would silently skip an entire batch on the
  first `esummary` exception, losing up to `batch_size` accessions
  from a temporal split. Both functions now share a helper
  `_esummary_batch_with_retries` that retries up to `max_retries` times
  with back-off, simplifying the previous `for/else: continue`
  construct and producing consistent logging.
  Files: `src/metagenome_generator/temporal_split.py`.

#### Tests

- Added `tests/test_a_category_bugfixes.py` with nine regression
  tests covering A1 (multi-record fixed/variable chunking and per-file
  cap semantics), A2 (first-batch dedup), A3 (multi-record concat),
  A4 (FASTQ round-trip through `filter_test_against_train`), A5
  (pool-clobber guard, both the `dest == pool` and
  `pool inside dest` cases), and A6 (retry path for
  `fetch_accession_dates`). Tests that depend on BLAST are skipped
  automatically when `blastn`/`makeblastdb` are not on the PATH.

#### Summary

- 6 files changed, 9 new tests, full suite passes (52/52, up from 43).
- No public API signatures were removed; `fetch_accession_dates` gained
  an optional `max_retries` kwarg (default matches the existing
  `ESUMMARY_MAX_RETRIES`), and chunk IDs for multi-record inputs now
  contain a `_seg{N}` suffix — single-record NCBI downloads are
  unaffected.

### B-category audit documentation push

#### Fixed

- **(B1) Package version mismatch.** `__init__.py` hard-coded
  `__version__ = "0.1.0"` while `pyproject.toml` was at `0.8.0`, so
  `import metagenome_generator; print(mg.__version__)` disagreed with
  `pip show`. The module now reads the installed version via
  `importlib.metadata.version("chimera-metagenome-generator")` and
  falls back to `"0.0.0+unknown"` when imported from an uninstalled
  source checkout. `pyproject.toml` becomes the single source of truth.
  Files: `src/metagenome_generator/__init__.py`.

- **(B2) `USER_GUIDE.md` referenced the wrong env vars and Python
  version.** The guide instructed users to export `NCBI_EMAIL` /
  `NCBI_API_KEY`, but `ncbi_search.py`, `temporal_split.py`, and
  `viral_taxonomy.py` all read `ENTREZ_EMAIL` / `ENTREZ_API_KEY`.
  The guide also claimed "Python 3.10+", contradicting `pyproject.toml`
  and `README.md` ("Python 3.8+"). Both are now corrected and the
  guide notes that the Python version tracks `pyproject.toml`.
  Files: `USER_GUIDE.md`.

- **(B3) Stale feature matrices in
  `docs/TOOLS_AND_FEATURES.md` and
  `docs/DATA_PREPARATION_COMPARISON.md`.** Both documents still listed
  features like abundance/coverage models, sequencing error
  simulation, external viral FASTA import, EVE export, genome
  completeness filter, taxonomy-aware balancing, and the structured
  benchmark recipe as "No" or "Not supported", contradicting
  `docs/improvements.md` and the actual CLI (all shipped). Both files
  were rewritten to reflect current functionality; the "Gaps" sections
  now list only genuinely open items (multi-length-per-run, 3rd-gen
  sequencing models, eukaryote negatives [out-of-scope], chimeras /
  paired-end, full metagenome assembly simulation, protein/annotation
  features). Added a "Status note" cross-reference so the three docs
  must be updated together.
  Files: `docs/TOOLS_AND_FEATURES.md`,
  `docs/DATA_PREPARATION_COMPARISON.md`.

- **(B4) Misleading `--num-bacteria` / `--num-virus` docs and help.**
  These flags control NCBI search counts and are **ignored** when
  `--accessions-file` is set (in which case `--max-*` caps the sample
  from the snapshot). The CLI help said only "Number of bacterial
  genomes to download" with no mention of the ignored-mode, and the
  README advised users to "set to 0 only when using --accessions-file
  with a snapshot that has no bacteria list" — which is both wrong and
  unnecessary. Help strings for `download` / `pipeline` `--num-*` now
  explicitly state they are ignored with `--accessions-file` and point
  to `--max-*`; the `--accessions-file` help itself now spells out that
  per-category counts come from `--max-*`. The README table was
  rewritten accordingly.
  Files: `src/metagenome_generator/cli.py`, `README.md`.

- **(B5) `genome-pool prepare` / `materialize` options had no help
  text.** Only `--pool-seed` had a description; every other flag
  (`--accessions-file`, `--pool-dir`, `--output-dir`, `--sample-seed`,
  all `--max-*`, `--copy`, `--no-clean`) rendered with just a name and
  metavar. Each flag now has a full-sentence help entry explaining its
  effect, its default, and cross-referencing related flags (e.g.
  `--output-dir` now mentions the A5 overlap guard).
  Files: `src/metagenome_generator/cli.py`.

- **(B6) Typo in `split-metagenome-train-test --input` help.**
  `Input metagenome file (FAST A/FASTQ)` → `Input metagenome file
  (FASTA/FASTQ)`.
  Files: `src/metagenome_generator/cli.py`.

- **(B7) Incomplete read-naming description in README.** The README
  documented read IDs as `{accession}_read_{idx}` /
  `viral_1_chunk_0`, but the code prefixes every ID with the
  category label (`bacteria`/`virus`/`archaea`/`plasmid`) and, after
  the A1 fix, now also adds a `_seg{N}` infix for records past the
  first in multi-record FASTAs. The naming paragraph was rewritten
  to show the real format (`bacteria_NC_000913.3_read_17`,
  `virus_NC_001798.2_seg1_read_0`, etc.) and kept consistent with the
  accurate description already in `USER_GUIDE.md`.
  Files: `README.md`.

#### Polish

- Added `migrate-snapshot` to the top-level `metagenome-generator
  --help` description (it was previously shipped but missing from the
  command listing).
  Files: `src/metagenome_generator/cli.py`.

#### Summary

- 6 files changed (plus one doc rewrite for each comparison doc);
  full test suite still passes (52/52) with no behavioural changes
  beyond the CLI help text and README/USER_GUIDE/docs refresh.
- `__version__` now tracks `pyproject.toml` automatically — bumping
  the package version no longer requires touching two files.

### C-category audit design-quirk push

#### Fixed

- **(C1) Dead outer loop in `filter_by_similarity`.** The function
  advertised a `max_refill_rounds` parameter and a `while` loop that
  looked like it would re-process candidates across rounds, but the
  very last line of the loop body unconditionally did
  `candidates = []`, so the loop always exited after a single pass.
  This was pure cognitive overhead: the warning branch about "max
  refill rounds reached" was unreachable, and the stats field
  `rounds` was always `1`. The outer loop has been removed and the
  function is now an explicit single-pass batch filter. The
  `max_refill_rounds` keyword argument is kept as a deprecated no-op
  for backward compatibility (the real refill loop has always lived
  in `chunk_genomes.build_metagenome`, which calls
  `filter_candidates_against_kept` on freshly generated chunks).
  Files: `src/metagenome_generator/similarity_filter.py`.

- **(C2) `normalize_train_split_percent(1.0)` is ambiguous.**
  The helper accepts both percentages (`80`) and fractions (`0.8`).
  The exact value `1.0` falls into the fraction branch and is
  interpreted as 100% — which is the sane default, but a user who
  actually meant "1%" would get a silent 99-point swing in the
  opposite direction. The function now emits a WARNING when a
  `float` equal to `1.0` is passed, telling the user to write
  `0.01` for "1%" or `100` for "100%". Numeric behaviour is
  unchanged, and integer `1` (always unambiguously percent) does
  not trigger the warning.
  Files: `src/metagenome_generator/chunk_genomes.py`.

- **(C3) Hard-coded BLAST DB prefix caused silent clobbering.**
  `_build_blast_db` always wrote its DB files to
  `<db_dir>/simfilter_db.*`, so two calls pointed at the same
  `db_dir` overwrote each other. Inside a single
  `filter_by_similarity` run, the first-batch self-BLAST DB was
  immediately overwritten by the kept-set DB (by design), but any
  caller sharing a `work_dir` across processes (e.g. concurrent
  train/test filtering jobs) could race on the same prefix and
  corrupt each other's databases. The prefix now defaults to
  `<fasta stem>_db` (so `kept.fasta` → `kept_db`, `batch.fasta` →
  `batch_db`), and an optional `prefix=` kwarg lets callers override
  it entirely.
  Files: `src/metagenome_generator/similarity_filter.py`.

- **(C4) `chunk_sequence` ignored `reads_per_organism`.** The
  variable-length (`chunk_sequence_variable`) and quality-profile
  (`chunk_sequence_quality_profile`) chunkers both honoured a
  `reads_per_organism` cap internally, but the fixed-length
  `chunk_sequence` had no such parameter — only the caller
  `_collect_chunks_for_file` applied the cap. Third-party callers
  of the public `chunk_sequence` API that expected symmetric
  behaviour silently got every possible read out of a genome.
  A new keyword-only `reads_per_organism: int | None = None`
  parameter was added to `chunk_sequence` and makes the three
  chunkers behave consistently. The caller path is unchanged
  because `_collect_chunks_for_file` still enforces its own
  per-file cap on top.
  Files: `src/metagenome_generator/chunk_genomes.py`.

- **(C5) `_lineage_from_efetch` relied on ElementTree truthiness.**
  XML parsing used the pattern `taxon.find("TaxId") or
  taxon.find("{*}TaxId")` to handle optional default namespaces.
  ElementTree `Element` objects are *falsy* when they have no
  children, so for leaf tags (`TaxId`, `Rank`, `ScientificName`)
  the first match would be discarded in favour of the namespaced
  lookup — and if both missed, the fallback chain silently
  short-circuited. Every such site now uses explicit
  `x = find(a); if x is None: x = find(b)` chaining. This also
  avoids the `DeprecationWarning: Testing an element's truth
  value will raise an exception in future versions` emitted by
  Python ≥ 3.12.
  Files: `src/metagenome_generator/viral_taxonomy.py`.

- **(C6) Silent viral-taxonomy bucket collapse.** When the
  `--viral-taxonomy` JSON is keyed in a format that doesn't match
  the file-stem prefix (e.g. versioned vs. unversioned accessions,
  or with/without `.fasta`), every viral prefix falls into the
  `"unknown"` bucket, `--balance-viral-by-taxonomy` becomes a
  no-op, and the user has no indication anything went wrong.
  Two new warnings surface this failure mode: `fetch_viral_taxonomy_groups`
  emits a WARNING if ≥ 50% of accessions resolve to `"unknown"`
  after an NCBI roundtrip, and `_apply_viral_taxonomy_balance`
  emits a WARNING (with sample unmapped prefixes + sample JSON
  keys) when most viral file prefixes are missing from the
  mapping.
  Files: `src/metagenome_generator/viral_taxonomy.py`,
  `src/metagenome_generator/chunk_genomes.py`.

- **(C7) `_validate_split_date` accepted impossible calendar
  dates.** The hand-rolled validation only checked that month
  was in 1..12 and day in 1..31, so `2025-02-31`, `2025-04-31`,
  or `2023-02-29` (non-leap) were all accepted. The downstream
  comparison against NCBI CreateDate strings is lexical, so an
  invalid date like `2025-02-31` would still sort between
  `2025-02-28` and `2025-03-01` and produce an off-by-month
  train/test split. The validator now delegates to
  `datetime.date.fromisoformat`, rejecting every non-existent
  date.
  Files: `src/metagenome_generator/temporal_split.py`.

#### Tests

- `tests/test_c_category_bugfixes.py` (new): 17 regression tests
  covering every C-fix, including parametrized leap-year /
  month-length cases for C7, BLAST-free C3 coverage via a
  `subprocess.run` stub, and a fake NCBI XML response for C5.

#### Summary

- 5 source files changed, 1 new test file (17 tests).
- Full test suite: 69/69 passing (52 pre-C + 17 new).
- No public CLI or file-format changes; the only new runtime
  visible behaviour is the added WARNING logs (C2, C6) and the
  stricter error for malformed split dates (C7).
- `filter_by_similarity` keeps its signature (`max_refill_rounds`
  is now a documented no-op kwarg), so downstream scripts that
  call it directly do not need changes.

### D-category audit polish push

#### Fixed

- **(D1) CLI / README command listing.** Verified during the B push
  that the top-level `metagenome-generator --help` description and
  the README "Command reference" table both list every shipped
  subcommand including `migrate-snapshot` and `genome-pool`. The
  earlier gaps (description was missing `migrate-snapshot`, polish
  landed during B) are now closed and the lists agree with
  `cli.py` subparsers. No code change was required in the D pass;
  this is documented here so future additions touch both sites.

- **(D2) `_download_category_batched` fallback logged a fake
  filename before fetching.** In the per-ID fallback branch
  (hit when a batch fetch raises) the loop printed
  `Fetching {category}: {gid} -> {gid}.fasta` *before* calling
  `fetch_sequences([gid])`, so the message claimed the final
  filename would be `{gid}.fasta` even though the batch-success
  branch derives the filename from the returned record's normalized
  accession (`rec.id.split()[0]`). If NCBI returned a
  differently-versioned accession the fallback log and the actual
  file disagreed. The fallback now (1) emits a single `WARNING`
  about the batch fallback with its cause, (2) performs the fetch,
  and (3) logs the same `Fetching ... -> {path.name}` / structured
  `logger.info` with the real on-disk path, matching the
  batch-success branch exactly.
  Files: `src/metagenome_generator/download_genomes.py`.

- **(D3) `run_benchmark_recipe` return-type docstring.** The
  docstring said "Returns the list of test FASTA/FASTQ paths (one
  per replicate)" without mentioning that the `.fasta` vs
  `.fastq` choice is driven by the `output_fastq` flag and that
  the corresponding metagenome + split files inherit the same
  suffix. Rewritten to describe the suffix contract explicitly so
  downstream callers don't have to read the implementation to
  know which file extension they will get.
  Files: `src/metagenome_generator/benchmark_recipe.py`.

- **(D4) `_collect_chunks_for_file` docstring was the exact place
  that hid A1.** The previous one-liner only mentioned EVE overlap
  and error-model selection; the per-file vs per-record semantics
  of `reads_per_organism`, the `_seg{N}` chunk-ID disambiguation
  added in A1, and the split of cap enforcement between the three
  chunker helpers were all undocumented. Rewritten to make the
  contract explicit: `reads_per_organism` is a per-file cap shared
  across all records in a multi-record FASTA, all three chunker
  helpers (fixed, variable, quality-profile) honour the cap
  consistently after C4, and records past the first are prefixed
  with `_seg{N}` so no two chunks from the same file collide.
  Files: `src/metagenome_generator/chunk_genomes.py`.

#### Summary

- 3 source files changed (docstrings + one log-statement reshuffle);
  no behavioural change to any public API, no new tests needed —
  the existing A-category multi-record coverage already exercises
  the D4 contract and D2 is exercised implicitly by
  `tests/test_download_cache.py`.
- Full test suite: 69/69 passing.

### Bugfix push summary (A + B + C + D)

Rolled up from the four pushes above. This is the shape of the
combined bugfix release off of the audit.

#### Code (A + C + D2)

- **A1**  `chunk_genomes.py` — remove stray `break` so every record
  of a multi-record FASTA is chunked; add `_seg{N}` disambiguation.
- **A2**  `similarity_filter.py` — first batch is now all-vs-all
  BLAST'd and de-duplicated before seeding `kept`.
- **A3**  `benchmark_recipe.py` — `SeqIO.read` → `SeqIO.parse`;
  multi-record FASTAs no longer crash diversity scoring.
- **A4**  `similarity_filter.py` / `cli.py` —
  `filter-test-against-train` detects FASTQ by extension and
  preserves Phred qualities; default output suffix tracks input.
- **A5**  `genome_pool.py` — `materialize_from_pool` refuses to
  delete the pool when `--output-dir` overlaps `--pool-dir`.
- **A6**  `temporal_split.py` — shared
  `_esummary_batch_with_retries`; `fetch_accession_dates` now
  retries transient NCBI failures.
- **C1**  `similarity_filter.py` — remove dead outer loop in
  `filter_by_similarity`; `max_refill_rounds` becomes a
  deprecated no-op kwarg.
- **C2**  `chunk_genomes.py` — warn on ambiguous
  `normalize_train_split_percent(1.0)`.
- **C3**  `similarity_filter.py` — per-fasta-stem BLAST DB prefix
  (`kept_db`, `batch_db`) instead of the hard-coded
  `simfilter_db`.
- **C4**  `chunk_genomes.py` — `chunk_sequence` now accepts
  `reads_per_organism` for symmetry with the other chunkers.
- **C5**  `viral_taxonomy.py` — explicit `is None` chaining
  replaces ElementTree-truthy bugs and silences the 3.12
  DeprecationWarning.
- **C6**  `viral_taxonomy.py` + `chunk_genomes.py` — warn when
  ≥ 50% of viral accessions / prefixes map to the `"unknown"`
  taxonomy bucket (the classic versioned-vs-unversioned mismatch).
- **C7**  `temporal_split.py` — `_validate_split_date` delegates to
  `datetime.date.fromisoformat`, rejecting impossible calendar
  dates.
- **D2**  `download_genomes.py` — per-ID fallback log prints the
  real on-disk filename and emits a single `WARNING` at the
  batch-failure transition.

#### Docs (B + D1 + D3 + D4)

- **B1**  `__init__.py` — `__version__` reads from installed
  metadata; `pyproject.toml` is the single source of truth.
- **B2**  `USER_GUIDE.md` — `ENTREZ_EMAIL`/`ENTREZ_API_KEY` and
  Python 3.8+.
- **B3**  `docs/TOOLS_AND_FEATURES.md` +
  `docs/DATA_PREPARATION_COMPARISON.md` — rewritten feature
  matrices; only genuine gaps remain.
- **B4**  `cli.py` + `README.md` — `--num-bacteria`/`--num-virus`
  help strings + README table spell out the ignored-with-
  `--accessions-file` behaviour.
- **B5**  `cli.py` — full help text for every `genome-pool prepare`
  / `materialize` flag.
- **B6**  `cli.py` — "FAST A/FASTQ" → "FASTA/FASTQ".
- **B7**  `README.md` — read/contig naming paragraph rewritten
  to match the actual
  `{category}_{stem}[_seg{N}]_{read|contig}_{idx}` format.
- **D1**  `cli.py` top-level description + `README.md` command
  reference both list every subcommand.
- **D3**  `benchmark_recipe.py` — `run_benchmark_recipe`
  docstring documents the `output_fastq` → `.fastq` suffix
  contract.
- **D4**  `chunk_genomes.py` — `_collect_chunks_for_file`
  docstring makes the per-file `reads_per_organism` budget,
  the `_seg{N}` disambiguation, and the chunker-cap contract
  explicit.

#### Tests

- `tests/test_a_category_bugfixes.py` — 9 tests (A1–A6).
- `tests/test_c_category_bugfixes.py` — 17 tests (C1–C7 with
  parametrized calendar-date coverage for C7).

#### Verification

- Full test suite: **69/69 passing** (52 pre-audit + 17 new).
- No lint errors on any modified file.
- No public API / signature changes: `filter_by_similarity`
  still accepts `max_refill_rounds` (deprecated no-op);
  `chunk_sequence` gains a keyword-only default-`None`
  argument; every CLI flag keeps its name and semantics.
- New user-visible runtime behaviour is additive only:
  extra `WARNING` log lines (A2, C2, C6), stricter calendar-
  date validation (C7), and an overlap guard at
  `genome-pool materialize` (A5).
