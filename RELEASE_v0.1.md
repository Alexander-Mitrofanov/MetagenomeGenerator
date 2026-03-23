# Release v0.1

First release of CHIMERA with pre-built viral BLAST database for EVE/prophage detection.

## Assets

- **viral_db_2026-03-10.tar.gz** — BLAST database built from RefSeq viral genomes (snapshot 2026-03-10). Extract and use with `blastn-filter --viral-db path/to/viral_db_2026-03-10/blastn_db/viral_db`.

## How to create this release on GitHub

1. Go to https://github.com/Alexander-Mitrofanov/MetagenomeGenerator/releases/new
2. **Tag:** Choose "Find or create a new tag" and type `v0.1`, then click "Create new tag: v0.1 on publish".
3. **Release title:** `v0.1`
4. **Description:** Copy the text above (or use this file).
5. **Attach:** Upload `viral_db_2026-03-10.tar.gz` from the project root.
6. Click **Publish release**.

The tarball is at: `MetagenomeGenerator/viral_db_2026-03-10.tar.gz`
