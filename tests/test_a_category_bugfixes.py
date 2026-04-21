"""Regression tests for the A-category bugfixes.

Each test targets a specific bug that the audit surfaced; they would have
failed on the pre-fix code and pass against the fix.
"""

from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path
from unittest.mock import patch

import pytest

from metagenome_generator.chunk_genomes import _collect_chunks_for_file
from metagenome_generator.genome_layout import VIRUS_DIR
from metagenome_generator.genome_pool import POOL_MANIFEST_NAME, materialize_from_pool
from metagenome_generator.similarity_filter import (
    filter_by_similarity,
    filter_test_against_train,
)
from metagenome_generator.temporal_split import fetch_accession_dates


BLAST_AVAILABLE = shutil.which("blastn") is not None and shutil.which("makeblastdb") is not None
requires_blast = pytest.mark.skipif(not BLAST_AVAILABLE, reason="BLAST+ not installed")


# ---------------------------------------------------------------------------
# A1: _collect_chunks_for_file must process every record in a multi-record FASTA
# ---------------------------------------------------------------------------


def _write_multirecord_virus_fasta(path: Path, n_records: int, seq_len: int = 400) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    seq = "ACGT" * (seq_len // 4)
    lines = []
    for i in range(n_records):
        lines.append(f">SEG{i}")
        lines.append(seq)
    path.write_text("\n".join(lines) + "\n")


def test_a1_fixed_length_consumes_all_records(tmp_path: Path) -> None:
    """Multi-segment virus FASTA must produce chunks from every segment, not just the first."""
    fp = tmp_path / VIRUS_DIR / "multi.fasta"
    _write_multirecord_virus_fasta(fp, n_records=3, seq_len=400)

    chunks = _collect_chunks_for_file(
        fp,
        sequence_length=100,
        reads_per_organism=None,
        prefix="multi",
    )

    # 3 records * 4 chunks per record = 12 chunks when the cap is removed. Before
    # the fix only record 0 was processed, yielding 4 chunks.
    assert len(chunks) == 12, f"expected chunks from every record, got {len(chunks)}"

    # Chunk IDs must be unique (regression: with a single shared prefix and
    # per-record idx restart, multi-record files emitted duplicate headers).
    ids = [rec.id for rec in chunks]
    assert len(ids) == len(set(ids)), f"chunk IDs collide across records: {ids}"

    # Records beyond the first must use a ``_seg{N}`` suffix; single-record files
    # keep the historical ``{prefix}_read_{idx}`` naming.
    assert any(rec.id.startswith("virus_multi_read_") for rec in chunks)
    assert any("multi_seg1_read_" in rec.id for rec in chunks)
    assert any("multi_seg2_read_" in rec.id for rec in chunks)


def test_a1_reads_per_organism_is_per_file_not_per_record(tmp_path: Path) -> None:
    """reads_per_organism caps across the whole file, even with multiple records."""
    fp = tmp_path / VIRUS_DIR / "multi.fasta"
    _write_multirecord_virus_fasta(fp, n_records=3, seq_len=400)

    chunks = _collect_chunks_for_file(
        fp,
        sequence_length=100,
        reads_per_organism=5,
        prefix="multi",
    )

    # Before the fix the count would be min(12, 5) via a single record.
    # After the fix we share the budget of 5 across records.
    assert len(chunks) == 5


def test_a1_variable_length_consumes_all_records(tmp_path: Path) -> None:
    """Variable-length chunking also must consume every record."""
    fp = tmp_path / VIRUS_DIR / "multi.fasta"
    _write_multirecord_virus_fasta(fp, n_records=2, seq_len=1000)

    chunks = _collect_chunks_for_file(
        fp,
        sequence_length=100,
        reads_per_organism=None,
        prefix="multi",
        min_length=200,
        max_length=300,
        seed=42,
    )

    # Each 1000-bp record yields at least one variable chunk of 200-300 bp, so at
    # least two chunks total. Before the fix only record 0 contributed.
    ids = [rec.id for rec in chunks]
    assert len(ids) == len(set(ids)), f"chunk IDs collide across records: {ids}"
    assert any(rid.startswith("virus_multi_contig_") for rid in ids)
    assert any("multi_seg1_contig_" in rid for rid in ids)


# ---------------------------------------------------------------------------
# A2: filter_by_similarity must deduplicate within the first batch
# ---------------------------------------------------------------------------


def _make_seqrecord(rec_id: str, seq: str):
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    return SeqRecord(Seq(seq), id=rec_id, description="")


def _random_dna(length: int, seed: int) -> str:
    import random as _random

    rng = _random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(length))


@requires_blast
def test_a2_first_batch_deduplicates_near_duplicates(tmp_path: Path) -> None:
    """Two identical sequences in the very first batch must not both survive."""
    # Non-repetitive (high-complexity) sequence so BLAST's default DUST masker
    # does not suppress the self-hits.
    base = _random_dna(400, seed=123)
    different = _random_dna(400, seed=456)
    records = [
        _make_seqrecord("dup_A", base),
        _make_seqrecord("dup_B", base),  # identical to dup_A
        _make_seqrecord("distinct", different),
    ]
    kept, stats = filter_by_similarity(
        records,
        target_count=5,
        similarity_threshold=90.0,
        min_coverage=0.8,
        batch_size=10,
        work_dir=tmp_path / "simfilter",
        num_threads=1,
    )
    kept_ids = {rec.id for rec in kept}
    # At most one of the two identical records should survive.
    assert not {"dup_A", "dup_B"}.issubset(kept_ids), (
        f"both near-duplicates from the first batch were kept: {kept_ids}"
    )
    assert stats["removed"] >= 1
    assert "distinct" in kept_ids


# ---------------------------------------------------------------------------
# A3: benchmark_recipe must handle multi-record viral FASTAs
# ---------------------------------------------------------------------------


def test_a3_write_concat_fasta_handles_multirecord(tmp_path: Path) -> None:
    from metagenome_generator.benchmark_recipe import _write_concat_fasta_and_lengths

    cache = tmp_path / "cache"
    cache.mkdir()
    # Two-segment virus in a single accession file.
    (cache / "NC_MULTI.1.fasta").write_text(
        ">NC_MULTI.1 segment 1\nACGTACGTAC\n>NC_MULTI.2 segment 2\nTTTTAAAA\n"
    )
    out_fasta = tmp_path / "concat.fasta"
    qlens = _write_concat_fasta_and_lengths(["NC_MULTI.1"], cache, out_fasta)

    # Both records must be accounted for; pre-fix SeqIO.read would have raised.
    assert "NC_MULTI.1" in qlens
    assert "NC_MULTI.2" in qlens
    assert qlens["NC_MULTI.1"] == 10
    assert qlens["NC_MULTI.2"] == 8
    text = out_fasta.read_text()
    assert "segment 1" in text and "segment 2" in text


# ---------------------------------------------------------------------------
# A4: filter_test_against_train must round-trip FASTQ
# ---------------------------------------------------------------------------


@requires_blast
def test_a4_filter_test_against_train_roundtrips_fastq(tmp_path: Path) -> None:
    # Non-repetitive sequence so BLAST's DUST masker does not suppress hits.
    seq_train = _random_dna(400, seed=11)
    seq_test_dup = seq_train  # identical -> should be removed
    seq_test_distinct = _random_dna(400, seed=22)

    def fq_record(rid: str, seq: str) -> str:
        q = "I" * len(seq)
        return f"@{rid}\n{seq}\n+\n{q}\n"

    train = tmp_path / "train.fastq"
    test = tmp_path / "test.fastq"
    out = tmp_path / "filtered.fastq"
    train.write_text(fq_record("train_1", seq_train))
    test.write_text(fq_record("test_dup", seq_test_dup) + fq_record("test_distinct", seq_test_distinct))

    n_removed, n_kept = filter_test_against_train(
        train, test, out,
        similarity_threshold=90.0,
        min_coverage=0.8,
        num_threads=1,
        work_dir=tmp_path / "work",
    )
    assert n_removed == 1
    assert n_kept == 1

    # Output must be a valid FASTQ (preserves qualities), not silently FASTA.
    lines = [ln for ln in out.read_text().splitlines() if ln]
    assert lines[0].startswith("@test_distinct")
    assert lines[2] == "+"
    assert lines[3] == "I" * len(seq_test_distinct)


# ---------------------------------------------------------------------------
# A5: genome-pool materialize must refuse dest==pool and pool-inside-dest
# ---------------------------------------------------------------------------


def _seed_pool(pool: Path) -> None:
    (pool / "bacteria").mkdir(parents=True)
    (pool / "virus").mkdir(parents=True)
    (pool / "bacteria" / "A1.fasta").write_text(">A1\nACGT\n")
    (pool / "virus" / "V1.fasta").write_text(">V1\nACGT\n")
    manifest = {"bacterial": ["A1"], "viral": ["V1"], "archaea": [], "plasmid": []}
    (pool / POOL_MANIFEST_NAME).write_text(json.dumps(manifest))


def test_a5_materialize_refuses_equal_dirs(tmp_path: Path) -> None:
    pool = tmp_path / "pool"
    _seed_pool(pool)

    with pytest.raises(ValueError, match="differ from --pool-dir"):
        materialize_from_pool(
            pool,
            pool,
            sample_seed=0,
            max_bacteria=None,
            max_virus=None,
            max_archaea=None,
            max_plasmid=None,
            use_symlinks=True,
            clean_dest=True,
        )
    # Pool must still be intact.
    assert (pool / POOL_MANIFEST_NAME).exists()
    assert (pool / "bacteria" / "A1.fasta").exists()


def test_a5_materialize_refuses_pool_inside_dest(tmp_path: Path) -> None:
    parent = tmp_path / "parent"
    pool = parent / "pool"
    _seed_pool(pool)

    with pytest.raises(ValueError, match="contains the pool directory"):
        materialize_from_pool(
            pool,
            parent,  # pool is a child of parent; clean_dest=True would wipe the pool
            sample_seed=0,
            max_bacteria=None,
            max_virus=None,
            max_archaea=None,
            max_plasmid=None,
            use_symlinks=True,
            clean_dest=True,
        )
    assert (pool / POOL_MANIFEST_NAME).exists()


# ---------------------------------------------------------------------------
# A6: fetch_accession_dates must retry on transient NCBI failures
# ---------------------------------------------------------------------------


def test_a6_fetch_accession_dates_retries_transient_failure() -> None:
    """A single transient esummary failure must not drop an entire batch."""
    from metagenome_generator import temporal_split

    call_count = {"n": 0}

    class _DummyHandle:
        def close(self) -> None:
            return None

    def flaky_esummary(**_kwargs):
        call_count["n"] += 1
        if call_count["n"] == 1:
            raise RuntimeError("simulated transient NCBI failure")
        return _DummyHandle()

    def fake_read(_handle):
        return [
            {"AccessionVersion": "NC_000001.1", "CreateDate": "2020/01/02", "Title": ""},
            {"AccessionVersion": "NC_000002.1", "CreateDate": "2021/06/07", "Title": ""},
        ]

    # Avoid sleeping in tests.
    with patch.object(temporal_split.Entrez, "esummary", side_effect=flaky_esummary), \
         patch.object(temporal_split.Entrez, "read", side_effect=fake_read), \
         patch.object(temporal_split.time, "sleep", lambda *_args, **_kwargs: None):
        dates = fetch_accession_dates(["NC_000001.1", "NC_000002.1"], max_retries=3)

    assert call_count["n"] >= 2, "retry loop did not trigger"
    assert dates == {
        "NC_000001.1": "2020/01/02",
        "NC_000002.1": "2021/06/07",
    }
