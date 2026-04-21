"""Regression tests for the C-category design/quirk fixes.

Each test targets a specific issue from the audit's C-category; they would
have failed on the pre-fix code and pass against the fix.
"""

from __future__ import annotations

import logging
from pathlib import Path
from unittest.mock import patch

import pytest

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from metagenome_generator.chunk_genomes import (
    _apply_viral_taxonomy_balance,
    chunk_sequence,
    normalize_train_split_percent,
)
from metagenome_generator.genome_layout import VIRUS_DIR
from metagenome_generator.similarity_filter import _build_blast_db
from metagenome_generator.temporal_split import _validate_split_date


# ---------------------------------------------------------------------------
# C1: filter_by_similarity — max_refill_rounds is a no-op kwarg, not dead code
# ---------------------------------------------------------------------------


def test_c1_max_refill_rounds_kwarg_still_accepted() -> None:
    """C1 removed the dead outer loop but must keep the kwarg for back-compat."""
    from metagenome_generator.similarity_filter import filter_by_similarity

    # Empty records: no BLAST needed, proves the signature still accepts the kwarg.
    result, stats = filter_by_similarity([], target_count=5, max_refill_rounds=7)
    assert result == []
    assert stats["kept"] == 0
    # The stats dict no longer has a "max rounds reached" warning path (empty
    # input short-circuits before the warning block).
    assert stats["warning"] is None


# ---------------------------------------------------------------------------
# C2: normalize_train_split_percent warns on the ambiguous 1.0 value
# ---------------------------------------------------------------------------


def test_c2_normalize_one_float_warns_and_returns_100(caplog) -> None:
    """Exactly 1.0 (float) is ambiguous between "1%" and "100%"; we warn and
    follow the fraction branch (→ 100.0)."""
    with caplog.at_level(logging.WARNING, logger="metagenome_generator.chunk_genomes"):
        out = normalize_train_split_percent(1.0)
    assert out == 100.0
    assert any("ambiguous" in rec.message.lower() for rec in caplog.records), (
        "expected a warning about the ambiguous 1.0 value"
    )


def test_c2_normalize_one_int_does_not_warn(caplog) -> None:
    """1 (int) is unambiguous ("1%" by convention: int is percentage). No warning."""
    with caplog.at_level(logging.WARNING, logger="metagenome_generator.chunk_genomes"):
        out = normalize_train_split_percent(1)
    # int 1 hits the 0<x<=1 branch → 100.0 as before, but must NOT warn (int
    # has no ambiguous float interpretation; users pass ints only as percent).
    # The branch still returns 100.0; the goal of C2 is the warning on float,
    # not to change numeric behavior.
    assert out == 100.0
    assert not any("ambiguous" in rec.message.lower() for rec in caplog.records)


def test_c2_normalize_percentages_unchanged() -> None:
    assert normalize_train_split_percent(80) == 80
    assert normalize_train_split_percent(80.0) == 80.0
    assert normalize_train_split_percent(0.8) == pytest.approx(80.0)
    assert normalize_train_split_percent(0.01) == pytest.approx(1.0)


# ---------------------------------------------------------------------------
# C3: _build_blast_db uses a per-fasta prefix so sibling DBs don't collide
# ---------------------------------------------------------------------------


def test_c3_build_blast_db_derives_prefix_from_fasta_stem(tmp_path: Path) -> None:
    """Before C3 the DB prefix was hard-coded to 'simfilter_db' so two builds
    in the same db_dir overwrote each other. Now the prefix defaults to
    '<fasta stem>_db', keeping distinct DBs distinct."""
    calls: list[list[str]] = []

    def fake_run(cmd, **_kwargs):
        calls.append(list(cmd))

        class _R:
            returncode = 0

        return _R()

    with patch("metagenome_generator.similarity_filter.subprocess.run", side_effect=fake_run):
        out_kept = _build_blast_db(tmp_path / "kept.fasta", tmp_path / "db")
        out_batch = _build_blast_db(tmp_path / "batch.fasta", tmp_path / "db")
        out_override = _build_blast_db(
            tmp_path / "ignored.fasta", tmp_path / "db", prefix="myprefix"
        )

    assert out_kept.name == "kept_db"
    assert out_batch.name == "batch_db"
    assert out_kept != out_batch
    assert out_override.name == "myprefix"
    assert len(calls) == 3


# ---------------------------------------------------------------------------
# C4: chunk_sequence internally respects reads_per_organism
# ---------------------------------------------------------------------------


def test_c4_chunk_sequence_respects_reads_per_organism() -> None:
    """chunk_sequence previously had no internal cap; only the caller enforced
    reads_per_organism. Direct callers that passed a cap got every possible
    read instead, risking silent over-generation."""
    rec = SeqRecord(Seq("ACGT" * 500), id="seq1", description="")
    # 2000 bp / 100 bp = 20 possible full-length reads.
    all_reads = list(chunk_sequence(rec, prefix="p", chunk_size=100))
    assert len(all_reads) == 20

    capped = list(chunk_sequence(rec, prefix="p", chunk_size=100, reads_per_organism=5))
    assert len(capped) == 5
    # Reads remain non-overlapping and correctly indexed.
    assert [r.id for r in capped] == [f"p_read_{i}" for i in range(5)]


# ---------------------------------------------------------------------------
# C5: _lineage_from_efetch uses explicit `is None` chaining so leaf elements
# don't get skipped by ElementTree's falsy-on-empty behavior.
# ---------------------------------------------------------------------------


def test_c5_lineage_parser_handles_childless_leaves() -> None:
    """Leaf Elements (TaxId, Rank, ScientificName) have no children and are
    therefore *falsy* under ElementTree's legacy __bool__. The old
    ``taxon.find("TaxId") or taxon.find("{*}TaxId")`` would incorrectly fall
    through to the namespaced variant even when the plain one matched.
    After C5 we use `is None` chaining and parsing works on a plain-XML
    (un-namespaced) TaxaSet document.
    """
    from metagenome_generator.viral_taxonomy import _lineage_from_efetch

    xml_bytes = b"""<TaxaSet>
  <Taxon>
    <TaxId>10310</TaxId>
    <ScientificName>Human alphaherpesvirus 2</ScientificName>
    <LineageEx>
      <Taxon>
        <TaxId>10239</TaxId>
        <ScientificName>Viruses</ScientificName>
        <Rank>superkingdom</Rank>
      </Taxon>
      <Taxon>
        <TaxId>10292</TaxId>
        <ScientificName>Herpesviridae</ScientificName>
        <Rank>family</Rank>
      </Taxon>
    </LineageEx>
  </Taxon>
</TaxaSet>
"""

    class _FakeHandle:
        def __init__(self, data: bytes):
            self._data = data
            self._pos = 0

        def read(self, *_a, **_k):
            d = self._data[self._pos :]
            self._pos = len(self._data)
            return d

        def close(self):  # pragma: no cover - not needed by ElementTree.parse
            pass

    captured = {}

    def fake_efetch(**kwargs):
        captured["ids"] = kwargs.get("id")
        import io

        return io.BytesIO(xml_bytes)

    with patch("metagenome_generator.viral_taxonomy.Entrez.efetch", side_effect=fake_efetch):
        mapping = _lineage_from_efetch(["10310"], level="family")

    assert mapping == {"10310": "Herpesviridae"}


# ---------------------------------------------------------------------------
# C6: _apply_viral_taxonomy_balance warns when most viral prefixes are unmapped
# ---------------------------------------------------------------------------


def test_c6_taxonomy_balance_warns_on_high_unknown_fraction(
    tmp_path: Path, caplog
) -> None:
    virus_dir = tmp_path / VIRUS_DIR
    virus_dir.mkdir(parents=True)
    prefixes_and_files: list[tuple[str, Path]] = []
    for i in range(4):
        fp = virus_dir / f"NC_{i:07d}.1.fasta"
        fp.write_text(f">NC_{i:07d}.1\n" + "ACGT" * 100 + "\n")
        prefixes_and_files.append((fp.stem, fp))

    # Taxonomy is keyed by *unversioned* accessions — classic format mismatch.
    viral_taxonomy = {
        "NC_0000000": "Herpesviridae",
        "NC_0000001": "Herpesviridae",
        "NC_0000002": "Poxviridae",
        "NC_0000003": "Poxviridae",
    }
    prefix_to_max_reads = {p: 10 for p, _ in prefixes_and_files}

    with caplog.at_level(logging.WARNING, logger="metagenome_generator.chunk_genomes"):
        _apply_viral_taxonomy_balance(
            prefixes_and_files,
            read_limits=[None] * 4,
            viral_taxonomy=viral_taxonomy,
            prefix_to_max_reads=prefix_to_max_reads,
        )

    assert any(
        "unmapped" in rec.message and "unknown" in rec.message
        for rec in caplog.records
    ), "expected a warning about unmapped viral prefixes"


# ---------------------------------------------------------------------------
# C7: _validate_split_date rejects impossible calendar dates
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "bad",
    [
        "2025-02-31",  # Feb 31st
        "2025-13-01",  # month 13
        "2025-00-10",  # month 0
        "2025-04-31",  # April has 30 days
        "2023-02-29",  # non-leap year
    ],
)
def test_c7_validate_split_date_rejects_invalid_calendar(bad: str) -> None:
    with pytest.raises(ValueError):
        _validate_split_date(bad)


@pytest.mark.parametrize(
    "good",
    [
        "2025-02-28",
        "2024-02-29",  # leap year
        "2020-12-31",
        "2000-01-01",
    ],
)
def test_c7_validate_split_date_accepts_valid_calendar(good: str) -> None:
    _validate_split_date(good)  # must not raise
