#!/usr/bin/env python3
from __future__ import annotations
"""Download bacterial and viral genomes from NCBI.

Saves into category subfolders: bacteria/, virus/, archaea/, plasmid/ under output_dir.
Exposes `download_genomes` for programmatic use and a CLI when run directly.

Reproducibility: NCBI search results change over time. Use --save-accessions to write
the accession list (and timestamp) to a JSON file after searching; use --accessions-file
to skip the search and download exactly those accessions on a later run.
"""

import argparse
import json
import logging
import time
from pathlib import Path

from Bio import Entrez, SeqIO

from .genome_layout import BACTERIA_DIR, VIRUS_DIR, ARCHAEA_DIR, PLASMID_DIR
from .ncbi_search import DEFAULT_QUERIES, search_genomes

logger = logging.getLogger(__name__)


# Keys used in the accessions JSON file (reproducibility)
ACCESSIONS_KEY_TIMESTAMP = "timestamp"
ACCESSIONS_KEY_BACTERIAL = "bacterial"
ACCESSIONS_KEY_VIRAL = "viral"
ACCESSIONS_KEY_ARCHAEA = "archaea"
ACCESSIONS_KEY_PLASMID = "plasmid"


def load_accessions(path: Path) -> dict:
    """Load accession lists and optional timestamp from a JSON file.

    Expected keys: bacterial, viral, archaea, plasmid (lists of accession IDs).
    Optional key: timestamp (ISO8601 string, when the list was obtained).
    Returns the full JSON dict; missing list keys default to [] when used.
    """
    with path.open() as f:
        return json.load(f)


def save_accessions(
    path: Path,
    bacterial: list[str],
    viral: list[str],
    archaea: list[str],
    plasmid: list[str],
) -> None:
    """Write accession lists and current UTC timestamp to a JSON file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    now = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    data = {
        ACCESSIONS_KEY_TIMESTAMP: now,
        ACCESSIONS_KEY_BACTERIAL: bacterial,
        ACCESSIONS_KEY_VIRAL: viral,
        ACCESSIONS_KEY_ARCHAEA: archaea,
        ACCESSIONS_KEY_PLASMID: plasmid,
    }
    with path.open("w") as f:
        json.dump(data, f, indent=2)
    logger.info("Saved accessions (timestamp=%s) to %s", now, path)


def fetch_sequences(ids: list[str], max_retries: int = 3) -> list:
    """Fetch FASTA sequences for given IDs from NCBI. Retries on transient failures."""
    if not ids:
        return []
    id_str = ",".join(ids)
    last_error = None
    for attempt in range(max_retries):
        try:
            handle = Entrez.efetch(
                db="nucleotide",
                id=id_str,
                rettype="fasta",
                retmode="text",
            )
            records = list(SeqIO.parse(handle, "fasta"))
            handle.close()
            return records
        except Exception as e:
            last_error = e
            if attempt < max_retries - 1:
                time.sleep(1.0 * (attempt + 1))
            continue
    raise last_error


def download_genomes(
    num_organisms: int,
    output_dir: Path,
    *,
    num_archaea: int = 0,
    num_plasmid: int = 0,
    accessions_file: Path | None = None,
    save_accessions_to: Path | None = None,
) -> None:
    """Download bacterial and viral genomes (and optionally archaea, plasmid) into `output_dir`.

    Uses category subfolders: bacteria/, virus/, archaea/, plasmid/. Negative samples
    (bacteria, archaea, plasmid) support training classifiers to distinguish viruses.

    If accessions_file is set, load accession IDs from that JSON (and optional timestamp)
    instead of searching NCBI, so the run is reproducible. If save_accessions_to is set,
    write the accession list and current UTC timestamp to that file after searching.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    bacteria_dir = output_dir / BACTERIA_DIR
    virus_dir = output_dir / VIRUS_DIR
    bacteria_dir.mkdir(parents=True, exist_ok=True)
    virus_dir.mkdir(parents=True, exist_ok=True)

    if accessions_file is not None:
        data = load_accessions(accessions_file)
        bacterial_ids = data.get(ACCESSIONS_KEY_BACTERIAL, [])
        viral_ids = data.get(ACCESSIONS_KEY_VIRAL, [])
        archaea_ids = data.get(ACCESSIONS_KEY_ARCHAEA, [])
        plasmid_ids = data.get(ACCESSIONS_KEY_PLASMID, [])
        ts = data.get(ACCESSIONS_KEY_TIMESTAMP, "unknown")
        print(f"Using accessions from {accessions_file} (timestamp: {ts})")
        print(f"  Bacterial: {len(bacterial_ids)}, viral: {len(viral_ids)}, archaea: {len(archaea_ids)}, plasmid: {len(plasmid_ids)}")
    else:
        print("Searching for bacterial genomes...")
        bacterial_ids = search_genomes(DEFAULT_QUERIES["bacterial"], num_organisms)
        print(f"  Found {len(bacterial_ids)} IDs")

        print("Searching for viral genomes...")
        viral_ids = search_genomes(DEFAULT_QUERIES["viral"], num_organisms)
        print(f"  Found {len(viral_ids)} IDs")

        archaea_ids = []
        if num_archaea > 0:
            print("Searching for archaeal genomes...")
            archaea_ids = search_genomes(DEFAULT_QUERIES["archaea"], num_archaea)
            print(f"  Found {len(archaea_ids)} IDs")

        plasmid_ids = []
        if num_plasmid > 0:
            print("Searching for plasmids...")
            plasmid_ids = search_genomes(DEFAULT_QUERIES["plasmid"], num_plasmid)
            print(f"  Found {len(plasmid_ids)} IDs")

        if save_accessions_to is not None:
            save_accessions(
                save_accessions_to,
                bacterial_ids,
                viral_ids,
                archaea_ids,
                plasmid_ids,
            )
            print(f"Saved accession list to {save_accessions_to}")

    for i, gid in enumerate(bacterial_ids):
        path = bacteria_dir / f"bacterial_{i + 1}.fasta"
        print(f"Fetching bacterial {i + 1}/{len(bacterial_ids)}: {gid} -> {path}")
        logger.info("Download: accession=%s origin=bacterial path=%s", gid, path)
        time.sleep(0.4)
        try:
            records = fetch_sequences([gid])
            SeqIO.write(records, path, "fasta")
        except Exception as e:
            print(f"  Warning: failed to fetch {gid}: {e}")
            logger.warning("Download failed: accession=%s origin=bacterial error=%s", gid, e)

    for i, gid in enumerate(viral_ids):
        path = virus_dir / f"viral_{i + 1}.fasta"
        print(f"Fetching viral {i + 1}/{len(viral_ids)}: {gid} -> {path}")
        logger.info("Download: accession=%s origin=viral path=%s", gid, path)
        time.sleep(0.4)
        try:
            records = fetch_sequences([gid])
            SeqIO.write(records, path, "fasta")
        except Exception as e:
            print(f"  Warning: failed to fetch {gid}: {e}")
            logger.warning("Download failed: accession=%s origin=viral error=%s", gid, e)

    if archaea_ids:
        archaea_dir = output_dir / ARCHAEA_DIR
        archaea_dir.mkdir(parents=True, exist_ok=True)
        for i, gid in enumerate(archaea_ids):
            path = archaea_dir / f"archaea_{i + 1}.fasta"
            print(f"Fetching archaea {i + 1}/{len(archaea_ids)}: {gid} -> {path}")
            logger.info("Download: accession=%s origin=archaea path=%s", gid, path)
            time.sleep(0.4)
            try:
                records = fetch_sequences([gid])
                SeqIO.write(records, path, "fasta")
            except Exception as e:
                print(f"  Warning: failed to fetch {gid}: {e}")
                logger.warning("Download failed: accession=%s origin=archaea error=%s", gid, e)

    if plasmid_ids:
        plasmid_dir = output_dir / PLASMID_DIR
        plasmid_dir.mkdir(parents=True, exist_ok=True)
        for i, gid in enumerate(plasmid_ids):
            path = plasmid_dir / f"plasmid_{i + 1}.fasta"
            print(f"Fetching plasmid {i + 1}/{len(plasmid_ids)}: {gid} -> {path}")
            logger.info("Download: accession=%s origin=plasmid path=%s", gid, path)
            time.sleep(0.4)
            try:
                records = fetch_sequences([gid])
                SeqIO.write(records, path, "fasta")
            except Exception as e:
                print(f"  Warning: failed to fetch {gid}: {e}")
                logger.warning("Download failed: accession=%s origin=plasmid error=%s", gid, e)

    logger.info("Download step complete: genomes saved in %s", output_dir.resolve())
    print("Done. Genome FASTAs saved in", output_dir.resolve())


def _cli(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="Download bacterial and viral genomes (and optionally archaea, plasmid) from NCBI."
    )
    parser.add_argument(
        "--num-organisms",
        type=int,
        default=10,
        help="Number of bacterial and viral genomes each. Default: 10",
    )
    parser.add_argument(
        "--num-archaea",
        type=int,
        default=0,
        help="Number of archaeal genomes to download (negative samples). Default: 0",
    )
    parser.add_argument(
        "--num-plasmid",
        type=int,
        default=0,
        help="Number of plasmids to download (negative samples). Default: 0",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("output"),
        help="Directory for downloaded genome FASTAs. Default: output",
    )
    parser.add_argument(
        "--accessions-file",
        type=Path,
        default=None,
        metavar="PATH",
        help="Load accession IDs from JSON (skip NCBI search). File must contain bacterial, viral (and optionally archaea, plasmid) lists and may include a timestamp.",
    )
    parser.add_argument(
        "--save-accessions",
        type=Path,
        default=None,
        metavar="PATH",
        help="After searching NCBI, save accession list and UTC timestamp to this JSON for reproducible runs (use with --accessions-file later).",
    )
    args = parser.parse_args(argv)

    if args.accessions_file is not None and not args.accessions_file.exists():
        parser.error(f"--accessions-file not found: {args.accessions_file}")

    download_genomes(
        args.num_organisms,
        args.output_dir,
        num_archaea=args.num_archaea,
        num_plasmid=args.num_plasmid,
        accessions_file=args.accessions_file,
        save_accessions_to=args.save_accessions,
    )


def main() -> None:
    _cli()


if __name__ == "__main__":
    main()
