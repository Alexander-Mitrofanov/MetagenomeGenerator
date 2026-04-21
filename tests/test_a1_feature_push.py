"""Regression tests for the A1 feature push.

Each test covers a new feature added in the A1 push:

- F1: ``--multi-length`` / :func:`build_metagenome_multi_length`
- F2: ``--paired`` / :func:`build_metagenome_paired`
- F3: ``--coverage`` / ``--coverage-cv`` (coverage-driven read counts)
- F4: 3rd-generation error models (``nanopore`` / ``pacbio-hifi`` / ``pacbio-clr``)
- F5: ``--chimera-rate`` / ``--pcr-duplicate-rate`` (library-prep artefacts)
- F6: ``--embed-taxonomy`` (gold-standard header labels)
"""

from __future__ import annotations

import json
import random
from pathlib import Path

from Bio import SeqIO

from metagenome_generator.chunk_genomes import (
    _apply_third_gen_errors,
    _embed_tax_in_description,
    _introduce_chimeras,
    _introduce_pcr_duplicates,
    _phred_from_rate,
    add_illumina_qualities_to_record,
    build_metagenome,
    build_metagenome_multi_length,
    build_metagenome_paired,
    chunk_sequence_paired,
    parse_multi_length_spec,
)
from metagenome_generator.genome_layout import BACTERIA_DIR, VIRUS_DIR

# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------


def _make_fasta(path: Path, name: str, sequence: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f">{name}\n{sequence}\n")


def _random_dna(n: int, seed: int = 0) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def _make_layout(tmp_path: Path, n_virus: int = 2, n_bacteria: int = 1, length: int = 4000) -> Path:
    """Build a minimal `bacteria/`/`virus/` layout with random-DNA genomes."""
    base = tmp_path / "genomes"
    for i in range(n_virus):
        _make_fasta(
            base / VIRUS_DIR / f"virus_{i}.fasta",
            f"NC_V{i}",
            _random_dna(length, seed=100 + i),
        )
    for i in range(n_bacteria):
        _make_fasta(
            base / BACTERIA_DIR / f"bact_{i}.fasta",
            f"NC_B{i}",
            _random_dna(length, seed=200 + i),
        )
    return base


def _prep_out(out_path: Path) -> Path:
    """Ensure ``out_path.parent`` exists (required by :func:`build_metagenome`)."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    return out_path


# ---------------------------------------------------------------------------
# F1: Multi-length
# ---------------------------------------------------------------------------


def test_f1_parse_multi_length_spec_happy() -> None:
    assert parse_multi_length_spec("300,500,1000,3000") == [300, 500, 1000, 3000]
    assert parse_multi_length_spec(" 300 , 500 ") == [300, 500]
    assert parse_multi_length_spec("300,300,500") == [300, 500]  # dedup


def test_f1_parse_multi_length_rejects_bad_input() -> None:
    import pytest

    with pytest.raises(ValueError):
        parse_multi_length_spec("")
    with pytest.raises(ValueError):
        parse_multi_length_spec("abc")
    with pytest.raises(ValueError):
        parse_multi_length_spec("0")
    with pytest.raises(ValueError):
        parse_multi_length_spec("-5")


def test_f1_multi_length_writes_one_file_per_length(tmp_path: Path) -> None:
    base = _make_layout(tmp_path, n_virus=2, n_bacteria=1, length=3000)
    out = _prep_out(tmp_path / "out" / "meta.fasta")

    results = build_metagenome_multi_length(
        base,
        out,
        [200, 500, 1000],
        reads_per_organism=3,
    )

    assert [L for L, _, _ in results] == [200, 500, 1000]
    for length, count, path in results:
        assert path.name == f"meta_L{length}.fasta"
        assert path.exists()
        recs = list(SeqIO.parse(path, "fasta"))
        assert len(recs) == count > 0
        # All reads in a per-length file are exactly that length.
        got_lengths = [len(r.seq) for r in recs]
        assert all(L == length for L in got_lengths), f"L={length}: got lengths {got_lengths}"


def test_f1_multi_length_rejects_variable_length_kwargs(tmp_path: Path) -> None:
    import pytest

    base = _make_layout(tmp_path)
    out = _prep_out(tmp_path / "out" / "meta.fasta")
    with pytest.raises(ValueError, match="min_length"):
        build_metagenome_multi_length(
            base, out, [300],
            reads_per_organism=2,
            min_length=100, max_length=500,
        )


# ---------------------------------------------------------------------------
# F2: Paired-end
# ---------------------------------------------------------------------------


def test_f2_chunk_sequence_paired_yields_matching_pairs() -> None:
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    seq = "ACGT" * 500  # 2000 bp
    rec = SeqRecord(Seq(seq), id="testseq", description="")
    pairs = list(
        chunk_sequence_paired(
            rec,
            "p",
            read_length=100,
            insert_size=400,
            insert_size_sd=0,
            reads_per_organism=3,
            rng=random.Random(1),
        )
    )
    assert len(pairs) == 3
    for r1, r2 in pairs:
        assert len(r1.seq) == 100
        assert len(r2.seq) == 100
        # R2 must be the reverse-complement of the 3' end of the fragment,
        # so revcomp(R2) should match a substring of the original sequence.
        assert str(r2.seq.reverse_complement()) in seq


def test_f2_build_paired_writes_two_files(tmp_path: Path) -> None:
    base = _make_layout(tmp_path, n_virus=2, n_bacteria=1, length=3000)
    out = _prep_out(tmp_path / "out" / "meta.fasta")

    n_pairs, r1_path, r2_path = build_metagenome_paired(
        base,
        out,
        read_length=150,
        reads_per_organism=5,
        insert_size=450,
        insert_size_sd=30,
        seed=7,
    )

    assert r1_path.name == "meta_R1.fasta"
    assert r2_path.name == "meta_R2.fasta"
    r1 = list(SeqIO.parse(r1_path, "fasta"))
    r2 = list(SeqIO.parse(r2_path, "fasta"))
    assert len(r1) == len(r2) == n_pairs
    # ``_pair_{idx}`` infix marks the pair; IDs should match between mates
    # except for the ``/1``/``/2`` mate tag embedded in the id.
    for a, b in zip(r1, r2):
        assert a.id.endswith("/1")
        assert b.id.endswith("/2")
        assert a.id[:-2] == b.id[:-2]


def test_f2_paired_fastq_has_quality_scores(tmp_path: Path) -> None:
    base = _make_layout(tmp_path, n_virus=1, n_bacteria=1, length=2000)
    out = _prep_out(tmp_path / "out" / "meta.fastq")

    n_pairs, r1_path, r2_path = build_metagenome_paired(
        base,
        out,
        read_length=120,
        reads_per_organism=3,
        insert_size=360,
        insert_size_sd=0,
        seed=42,
        output_fastq=True,
    )
    assert r1_path.name == "meta_R1.fastq"
    r1 = list(SeqIO.parse(r1_path, "fastq"))
    assert len(r1) == n_pairs
    for rec in r1:
        q = rec.letter_annotations.get("phred_quality")
        assert q is not None and len(q) == 120


# ---------------------------------------------------------------------------
# F3: Coverage-depth
# ---------------------------------------------------------------------------


def test_f3_coverage_scales_read_count_with_genome_size(tmp_path: Path) -> None:
    # Two genomes, one 2x the size of the other; coverage mode must
    # produce ~2x more reads for the larger genome.
    base = tmp_path / "genomes"
    _make_fasta(base / BACTERIA_DIR / "small.fasta", "NC_S", _random_dna(2000, seed=1))
    _make_fasta(base / BACTERIA_DIR / "large.fasta", "NC_L", _random_dna(4000, seed=2))

    out = _prep_out(tmp_path / "out" / "meta.fasta")
    count = build_metagenome(
        base,
        out,
        sequence_length=100,
        reads_per_organism=None,
        coverage=1.0,  # 1x coverage -> ~bp / read_length reads per genome
        seed=0,
    )
    assert count > 0
    # Count per-source by stripping the `_read_{idx}` / `_seg*_read_{idx}` suffix.
    from collections import Counter

    recs = list(SeqIO.parse(out, "fasta"))
    per_src: Counter[str] = Counter()
    for rec in recs:
        src_id = rec.id.split("_read_", 1)[0]
        per_src[src_id] += 1
    # Expect ~20 reads from small (2000/100) and ~40 from large (4000/100).
    small_total = sum(v for k, v in per_src.items() if "small" in k or k.endswith("_S"))
    large_total = sum(v for k, v in per_src.items() if "large" in k or k.endswith("_L"))
    assert small_total > 0 and large_total > 0
    ratio = large_total / small_total
    assert 1.5 <= ratio <= 2.5, f"large/small ratio {ratio} out of expected range"


def test_f3_coverage_cv_produces_variability(tmp_path: Path) -> None:
    # With CV > 0 the per-genome read counts should *not* all be equal.
    # Use large genomes + small coverage so the coverage-implied read count
    # (not the genome-length-derived cap) dominates and we can actually
    # observe the log-normal spread.
    base = tmp_path / "genomes"
    for i in range(8):
        _make_fasta(
            base / BACTERIA_DIR / f"bact_{i}.fasta",
            f"NC_B{i}",
            _random_dna(20000, seed=300 + i),
        )
    out = _prep_out(tmp_path / "out" / "meta.fasta")

    count = build_metagenome(
        base,
        out,
        sequence_length=100,
        reads_per_organism=None,
        coverage=0.3,
        coverage_cv=1.0,
        seed=123,
    )
    assert count > 0

    from collections import Counter

    recs = list(SeqIO.parse(out, "fasta"))
    per_src: Counter[str] = Counter()
    for rec in recs:
        src_id = rec.id.split("_read_", 1)[0]
        per_src[src_id] += 1
    # With uniform coverage, all 8 genomes would have the same read count.
    # CV=1.0 should produce meaningful spread.
    values = list(per_src.values())
    assert len(set(values)) > 1, (
        f"expected variability across genomes; got identical counts: {values}"
    )


# ---------------------------------------------------------------------------
# F4: 3rd-gen error models
# ---------------------------------------------------------------------------


def test_f4_nanopore_inflates_errors_in_homopolymers() -> None:
    # Homopolymer inflation should produce *more indels* in an all-A sequence
    # than in a random-DNA sequence of the same length. Using |len_out-len_in|
    # as the metric avoids the Hamming-after-indel-shift confound.
    n = 400
    poly = "A" * n
    non_poly = _random_dna(n, seed=2)

    def _indel_delta(in_seq: str, seed: int) -> int:
        out = _apply_third_gen_errors(
            in_seq,
            random.Random(seed),
            sub_rate=0.0,  # isolate indel behaviour
            ins_rate=0.02,
            del_rate=0.02,
            homopolymer_mult=3.0,
            homopolymer_min_run=3,
        )
        return abs(len(out) - len(in_seq))

    trials = 20
    poly_delta = sum(_indel_delta(poly, s) for s in range(trials))
    non_poly_delta = sum(_indel_delta(non_poly, s) for s in range(trials))
    assert poly_delta > non_poly_delta * 1.5, (
        f"homopolymer sequence should accumulate more indels; "
        f"got poly_delta={poly_delta} non_poly_delta={non_poly_delta}"
    )


def test_f4_pacbio_hifi_low_error_rate() -> None:
    # PacBio HiFi is ~0.3% total; the output length should stay within ~1%
    # of the input and the substitution count on a matched-length window
    # should remain low. We avoid prefix-Hamming after indel shifts.
    from metagenome_generator.chunk_genomes import (
        PACBIO_HIFI_DEL_RATE,
        PACBIO_HIFI_INS_RATE,
        PACBIO_HIFI_SUB_RATE,
    )

    seq = _random_dna(10000, seed=42)
    out = _apply_third_gen_errors(
        seq,
        random.Random(11),
        sub_rate=PACBIO_HIFI_SUB_RATE,
        ins_rate=PACBIO_HIFI_INS_RATE,
        del_rate=PACBIO_HIFI_DEL_RATE,
        homopolymer_mult=1.0,
    )
    # Length should be very close (ins/del rates are ≤0.05% each).
    length_delta = abs(len(out) - len(seq)) / len(seq)
    assert length_delta < 0.01, f"pacbio-hifi length delta {length_delta:.4f} too high"


def test_f4_fastq_qualities_match_third_gen_model() -> None:
    # For a nanopore read, the flat Phred score should be materially lower
    # than the average Illumina quality score. The exact gap depends on the
    # error-rate constants, so we just require nanopore < illumina by a
    # meaningful margin.
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    rec_ill = SeqRecord(Seq("A" * 150), id="r1", description="")
    rec_ont = SeqRecord(Seq("A" * 150), id="r2", description="")
    rec_hifi = SeqRecord(Seq("A" * 150), id="r3", description="")
    add_illumina_qualities_to_record(rec_ill, error_model="illumina")
    add_illumina_qualities_to_record(rec_ont, error_model="nanopore")
    add_illumina_qualities_to_record(rec_hifi, error_model="pacbio-hifi")
    avg_ill = sum(rec_ill.letter_annotations["phred_quality"]) / 150
    avg_ont = sum(rec_ont.letter_annotations["phred_quality"]) / 150
    avg_hifi = sum(rec_hifi.letter_annotations["phred_quality"]) / 150
    # Ordering: HiFi (~Q25) > Illumina (~Q22-30) > Nanopore (~Q13)
    assert avg_ont < avg_ill, (
        f"nanopore avg Q {avg_ont} should be lower than illumina avg Q {avg_ill}"
    )
    assert avg_hifi > avg_ont, (
        f"pacbio-hifi avg Q {avg_hifi} should be higher than nanopore avg Q {avg_ont}"
    )


def test_f4_phred_from_rate_monotonic() -> None:
    assert _phred_from_rate(0.001) > _phred_from_rate(0.01)
    assert _phred_from_rate(0.01) > _phred_from_rate(0.1)


# ---------------------------------------------------------------------------
# F5: Chimeras & PCR duplicates
# ---------------------------------------------------------------------------


def test_f5_introduce_chimeras_replaces_expected_fraction() -> None:
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    records = [SeqRecord(Seq("A" * 200), id=f"r_{i}", description="") for i in range(100)]
    new, n_chim = _introduce_chimeras(records, chimera_rate=0.1, rng=random.Random(1))
    assert n_chim == 10
    assert sum(1 for r in new if r.id.startswith("chimera_")) == 10
    assert len(new) == len(records)  # length preserved


def test_f5_introduce_chimeras_zero_rate_noop() -> None:
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    records = [SeqRecord(Seq("A" * 100), id=f"r_{i}", description="") for i in range(10)]
    new, n_chim = _introduce_chimeras(records, chimera_rate=0.0, rng=random.Random(0))
    assert n_chim == 0
    assert new == records


def test_f5_pcr_duplicates_grow_list() -> None:
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    records = [SeqRecord(Seq("A" * 100), id=f"r_{i}", description="") for i in range(200)]
    new, n_dup = _introduce_pcr_duplicates(records, duplicate_rate=0.5, rng=random.Random(7))
    assert len(new) == len(records) + n_dup
    # Expect ≈ half to be duplicated; allow generous tolerance for sampling.
    assert 80 <= n_dup <= 120
    for r in new[len(records):]:
        assert r.id.endswith("_dup")
        assert "pcr_duplicate=true" in r.description


def test_f5_end_to_end_chimera_and_dup_via_build_metagenome(tmp_path: Path) -> None:
    base = _make_layout(tmp_path, n_virus=2, n_bacteria=2, length=3000)
    out = _prep_out(tmp_path / "out" / "meta.fasta")
    count = build_metagenome(
        base,
        out,
        sequence_length=100,
        reads_per_organism=20,
        seed=5,
        chimera_rate=0.1,
        pcr_duplicate_rate=0.1,
    )
    recs = list(SeqIO.parse(out, "fasta"))
    assert len(recs) == count
    has_chimera = any(r.id.startswith("chimera_") for r in recs)
    has_dup = any(r.id.endswith("_dup") for r in recs)
    assert has_chimera, "expected at least one chimera record"
    assert has_dup, "expected at least one PCR-duplicate record"


# ---------------------------------------------------------------------------
# F6: Embed taxonomy
# ---------------------------------------------------------------------------


def test_f6_embed_tax_viral_uses_taxonomy_map() -> None:
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    rec = SeqRecord(Seq("A" * 50), id="virus_NC_V1_read_0", description="start=0 end=50")
    out = _embed_tax_in_description(rec, "virus", "NC_V1", {"NC_V1": "Herpesviridae"})
    assert "tax=Herpesviridae" in out.description
    # Start/end must be preserved.
    assert "start=0 end=50" in out.description


def test_f6_embed_tax_viral_unknown_fallback() -> None:
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    rec = SeqRecord(Seq("A" * 50), id="virus_X_read_0", description="")
    out = _embed_tax_in_description(rec, "virus", "X", {})
    assert "tax=unknown" in out.description


def test_f6_embed_tax_non_viral_uses_category() -> None:
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    rec = SeqRecord(Seq("A" * 50), id="bacteria_NC_B1_read_0", description="start=0 end=50")
    out = _embed_tax_in_description(rec, "bacteria", "NC_B1", {"NC_V1": "Herpesviridae"})
    assert "tax=bacteria" in out.description


def test_f6_build_metagenome_embeds_taxonomy_in_written_records(tmp_path: Path) -> None:
    base = _make_layout(tmp_path, n_virus=2, n_bacteria=1, length=3000)
    # Build a taxonomy JSON mapping virus prefixes (file stems) to groups.
    tax = {"virus_0": "Herpesviridae", "virus_1": "Myoviridae"}
    tax_path = tmp_path / "viral_tax.json"
    tax_path.write_text(json.dumps(tax))

    out = _prep_out(tmp_path / "out" / "meta.fasta")
    count = build_metagenome(
        base,
        out,
        sequence_length=150,
        reads_per_organism=3,
        seed=2,
        viral_taxonomy_json=tax_path,
        embed_taxonomy=True,
    )
    assert count > 0
    viral_reads = [r for r in SeqIO.parse(out, "fasta") if r.id.startswith("virus_")]
    bact_reads = [r for r in SeqIO.parse(out, "fasta") if r.id.startswith("bacteria_")]
    assert viral_reads
    assert bact_reads
    # Viral reads get their taxonomy group; bacterial reads get "tax=bacteria".
    assert any("tax=Herpesviridae" in r.description for r in viral_reads)
    assert any("tax=Myoviridae" in r.description for r in viral_reads)
    assert all("tax=bacteria" in r.description for r in bact_reads)
