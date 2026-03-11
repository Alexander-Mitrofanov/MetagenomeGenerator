# CLI validation tests

## Purpose

`test_cli_validation.py` runs **extensive negative tests** on the CLI: wrong flag combinations, invalid values, missing files, and values outside allowed ranges. The goal is to ensure the tool fails fast with clear messages instead of crashing or behaving unexpectedly.

## Test plan (categories)

| Category | What can go wrong | Tests |
|----------|-------------------|--------|
| **Download** | `--accessions-file` missing; negative `--num-bacteria` / `--num-virus` | `test_download_nonexistent_accessions_file`, `test_download_negative_num_bacteria`, `test_download_negative_num_virus` |
| **Chunk** | `--input` missing; `--sequence-length` 0 or negative; `--min-contig-length` > `--max-contig-length`; `--eve-intervals` missing; `--balance-viral-by-taxonomy` without `--viral-taxonomy`; substitution rate > 1 | `test_chunk_input_nonexistent`, `test_chunk_sequence_length_*`, `test_chunk_min_contig_gt_max_contig`, `test_chunk_eve_intervals_nonexistent`, `test_chunk_balance_viral_*`, `test_chunk_substitution_rate_*`, `test_chunk_train_test_split_*` |
| **filter-test-against-train** | Missing train or test FASTA | `test_filter_test_against_train_nonexistent_train` |
| **temporal-split** | Missing `--accessions-file`; invalid `--split-date` | `test_temporal_split_*`, `test_temporal_split_info_invalid_date` |
| **blastn-filter** | Missing `--genome-dir` | `test_blastn_filter_nonexistent_genome_dir` |
| **viral-taxonomy** | Missing `--accessions-file` | `test_viral_taxonomy_nonexistent_accessions_file` |

## Bugs found and fixed (during test implementation)

1. **Negative `--num-bacteria` / `--num-virus`**  
   NCBI `search_genomes` uses `record["IdList"][:count]`. For `count == -1` this becomes `[:-1]`, so the tool downloaded 99 genomes instead of 0.  
   **Fix:** Validate `num_bacteria >= 0` and `num_virus >= 0` in CLI (`_run_download`, pipeline) and in `download_genomes()`.

2. **Nonexistent `--input` in chunk**  
   When `--input` did not exist, the code continued; `iter_genome_fastas` returned `[]`, and the tool wrote 0 sequences and exited 0.  
   **Fix:** In `_run_chunk`, require `args.input.exists()` and exit with a clear message.

3. **`--sequence-length` 0 or negative**  
   `chunk_sequence(..., chunk_size=0)` leads to `range(0, N+1, 0)` → `ValueError`.  
   **Fix:** In `_run_chunk`, require `--sequence-length >= 1`.

4. **`--min-contig-length` > `--max-contig-length`**  
   `chunk_sequence_variable` uses `rng.randint(min_len, max_len)`; if `min_len > max_len`, `randint` raises.  
   **Fix:** In `_run_chunk`, require `min_contig_length <= max_contig_length` when both are set.

5. **`--accessions-file` not checked in unified CLI**  
   The standalone download script checked existence; the unified CLI did not, so `FileNotFoundError` was raised from `load_accessions`.  
   **Fix:** In `_run_download` (and pipeline when using download), check `accessions_file.exists()` before calling `download_genomes`.

6. **Test runner used wrong module**  
   Tests invoked `python -m metagenome_generator.cli`, which does not call `main()`, so all commands exited 0.  
   **Fix:** Use `python -m metagenome_generator` so `__main__.py` runs and `main()` is executed.

## How to run

From the project root:

```bash
python -m pytest tests/test_cli_validation.py -v
```

Optional: `--timeout=90` to cap long-running invocations (e.g. if a test mistakenly hits NCBI).

All tests are designed to **avoid network calls** where possible (nonexistent paths, invalid values). The negative `--num-bacteria` / `--num-virus` tests now exit immediately due to validation.
