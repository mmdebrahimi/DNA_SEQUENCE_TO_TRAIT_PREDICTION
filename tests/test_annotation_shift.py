"""Tests for the circular annotation shift -- the structurally-matched null for the frame arm.

The load-bearing test here is `test_valid_offsets_match_an_independent_interval_derivation`.
`valid_shift_offsets` reads the CDS coverage mask, which is the same object `frame_features` uses
for `inside_cds`; checking it against that mask again would be checking it against itself. So the
expected set is re-derived by INTERVAL ARITHMETIC instead -- a row [s, e) straddles under delta
exactly when delta is in the open range (L - e, L - s) -- and the two derivations must agree
exactly. If they ever diverge, one of them is wrong and the test says which offsets differ.
"""

from __future__ import annotations

import pytest

from dna_decode.glm.annotation_shift import (
    AnnotationShiftError,
    circular_shift_cds,
    sample_shift_offset,
    shift_match_report,
    valid_shift_offsets,
)
from dna_decode.glm.corpus_window import CdsRow, corpus_windows, window_coverage

# A small synthetic replicon. Deliberately NOT uniform: two genes nearly adjacent, one long gene,
# both strands represented, and a gap at the origin so that some offsets are valid.
L = 1000
ROWS = [
    CdsRow(50, 200, "+", "g1"),
    CdsRow(205, 260, "-", "g2"),
    CdsRow(400, 700, "+", "g3"),
    CdsRow(750, 820, "-", "g4"),
]


def _straddling_by_interval_arithmetic(rows, genome_length: int) -> set[int]:
    """Offsets that make some row cross the origin, derived WITHOUT the coverage mask.

    Row [s, e) shifted by d occupies [s+d, e+d). It crosses a multiple of L exactly when
    s + d < L < e + d, i.e. L - e < d < L - s.
    """
    bad: set[int] = set()
    for r in rows:
        lo = genome_length - r.end0
        hi = genome_length - r.start0
        for d in range(max(lo + 1, 0), min(hi, genome_length)):
            bad.add(d)
    return bad


def test_valid_offsets_match_an_independent_interval_derivation():
    got = set(valid_shift_offsets(ROWS, L))
    bad = _straddling_by_interval_arithmetic(ROWS, L)
    expected = set(range(L)) - bad
    assert got == expected, (
        "mask-derived and interval-derived valid sets disagree on %r"
        % sorted(got.symmetric_difference(expected))[:20]
    )
    # Non-vacuity: the two sets must not be trivially everything.
    assert 0 < len(got) < L


def test_every_valid_offset_actually_shifts_without_raising():
    for d in valid_shift_offsets(ROWS, L):
        circular_shift_cds(ROWS, L, d)


def test_every_invalid_offset_is_refused():
    bad = _straddling_by_interval_arithmetic(ROWS, L)
    assert bad, "the fixture must have at least one straddling offset or this test is vacuous"
    for d in sorted(bad):
        with pytest.raises(AnnotationShiftError, match="straddle the origin"):
            circular_shift_cds(ROWS, L, d)


def test_shift_preserves_count_lengths_and_strands_exactly():
    for d in (0, 1, 300, 999):
        if d not in set(valid_shift_offsets(ROWS, L)):
            continue
        sh = circular_shift_cds(ROWS, L, d)
        rep = shift_match_report(ROWS, sh, L)
        assert rep["gene_count_match"]
        assert rep["gene_length_multiset_match"]
        assert rep["strand_composition_match"]
        assert rep["cds_covered_bases_match"]


def test_shift_by_zero_and_by_genome_length_are_the_identity():
    assert circular_shift_cds(ROWS, L, 0) == ROWS
    assert circular_shift_cds(ROWS, L, L) == ROWS


def test_shift_composes_back_to_the_original():
    d = 300
    assert d in set(valid_shift_offsets(ROWS, L))
    once = circular_shift_cds(ROWS, L, d)
    back = circular_shift_cds(once, L, L - d)
    assert back == ROWS


def test_shift_preserves_row_order_and_locus_tags():
    sh = circular_shift_cds(ROWS, L, 300)
    assert [r.locus_tag for r in sh] == [r.locus_tag for r in ROWS]
    assert [r.strand for r in sh] == [r.strand for r in ROWS]


def test_sample_offset_is_deterministic_and_respects_the_circular_floor():
    a = sample_shift_offset(ROWS, L, seed=7, min_abs_shift=100)
    b = sample_shift_offset(ROWS, L, seed=7, min_abs_shift=100)
    assert a == b, "same seed must give the same offset or a run is not reproducible"
    assert min(a, L - a) >= 100, "the floor is CIRCULAR: L-20 is as close to identity as 20"
    assert a in set(valid_shift_offsets(ROWS, L))


def test_sample_offset_refuses_rather_than_returning_a_near_identity():
    with pytest.raises(AnnotationShiftError, match="too dense"):
        sample_shift_offset(ROWS, L, seed=0, min_abs_shift=L)


def test_the_shift_is_not_a_no_op_on_the_window_set():
    """A control that produced the same windows would control nothing."""
    d = sample_shift_offset(ROWS, L, seed=1, min_abs_shift=100)
    w_true = corpus_windows(ROWS, L)
    w_shift = corpus_windows(circular_shift_cds(ROWS, L, d), L)
    m_true = window_coverage(w_true, L)
    m_shift = window_coverage(w_shift, L)
    assert m_true != m_shift, "the shifted mask is identical to the true one; this is not a null"


def test_match_report_detects_a_genuine_mismatch():
    """Non-vacuity: the report must be capable of returning False."""
    other = [CdsRow(10, 20, "+", "x")]
    rep = shift_match_report(ROWS, other, L)
    assert not rep["gene_count_match"]
    assert not rep["gene_length_multiset_match"]
    assert not rep["cds_covered_bases_match"]


def test_refuses_rows_that_are_not_half_open_intervals_in_range():
    with pytest.raises(AnnotationShiftError, match="half-open interval"):
        circular_shift_cds([CdsRow(200, 100, "+", "backwards")], L, 10)
    with pytest.raises(AnnotationShiftError, match="half-open interval"):
        circular_shift_cds([CdsRow(10, L + 5, "+", "past_end")], L, 10)


def test_refuses_an_empty_annotation_and_a_nonpositive_genome():
    with pytest.raises(AnnotationShiftError, match="no CDS rows"):
        circular_shift_cds([], L, 10)
    with pytest.raises(AnnotationShiftError, match="must be positive"):
        circular_shift_cds(ROWS, 0, 10)
    with pytest.raises(AnnotationShiftError, match="must be positive"):
        valid_shift_offsets(ROWS, 0)


def test_a_fully_coding_replicon_leaves_only_the_degenerate_identity_offset():
    """The density limit is real and must surface as a refusal, not as a bad shift.

    A gene spanning the whole replicon leaves exactly ONE valid offset: zero. That is not a
    loophole in the rule -- shifting by nothing cannot make anything straddle, and interval
    arithmetic agrees (for [0, L) the forbidden range is the OPEN interval (0, L), which excludes
    0). What matters is that the sampler will not hand that offset back as a control, because the
    identity is not a null.
    """
    full = [CdsRow(0, L, "+", "whole")]
    assert valid_shift_offsets(full, L) == [0]
    with pytest.raises(AnnotationShiftError, match="too dense"):
        sample_shift_offset(full, L, seed=0, min_abs_shift=1)
