"""Tests for the inverse-design falsifier.

Pins the DESIGN (non-circularity, the real baseline, the conformal honesty flag) rather than the headline
number — the design is what makes the verdict mean anything. Real-data tests skip without the D: cache.
"""
from __future__ import annotations

import statistics

import pytest

from scripts.forward_inverse_roundtrip import (
    PG,
    Candidate,
    empirical_null,
    fit_isotonic,
    run_inverse,
    single_nt_accessible,
)


# ---- the isotonic calibrator (the score -> effect map) ---------------------------------------------------

def test_isotonic_is_monotone_non_decreasing():
    f = fit_isotonic([1, 2, 3, 4, 5], [0.0, 1.0, 0.5, 3.0, 4.0])   # note the violation at x=3
    ys = [f(x) for x in (1, 2, 3, 4, 5)]
    assert all(a <= b + 1e-9 for a, b in zip(ys, ys[1:])), ys


def test_isotonic_pools_adjacent_violators_to_their_mean():
    f = fit_isotonic([1, 2], [1.0, 0.0])     # a pure violation -> both become the mean 0.5
    assert f(1) == pytest.approx(0.5)
    assert f(2) == pytest.approx(0.5)


def test_isotonic_recovers_an_already_monotone_signal():
    f = fit_isotonic([1, 2, 3], [0.0, 1.0, 2.0])
    assert f(1) == pytest.approx(0.0)
    assert f(3) == pytest.approx(2.0)


def test_isotonic_clamps_outside_the_fitted_range():
    f = fit_isotonic([1, 2, 3], [0.0, 1.0, 2.0])
    assert f(-99) == pytest.approx(0.0)      # never extrapolates wildly
    assert f(99) == pytest.approx(2.0)


# ---- non-circularity: the calibrator must never see the graded labels ------------------------------------

def _cands(n=40, seed=0):
    # a clean monotone score->effect relationship, one variant per position
    out = []
    for i in range(n):
        out.append(Candidate(f"A{i+1}V", i + 1, "A", "V", measured=-i / 10.0))
    return out


def test_run_inverse_grades_on_measured_labels_not_on_its_own_prediction():
    """The circularity the handoff's SME pass killed: grading with the generating model measures
    self-consistency. The graded field must be the MEASURED value of the proposed variant."""
    cands = _cands(40)
    cal, test = cands[::2], cands[1::2]
    score = lambda c: -c.pos            # noqa: E731  perfectly rank-correlated with `measured`
    r = run_inverse("x", cal, test, score, [-1.0], top_k=1)
    row = r["per_target"][0]
    proposed = next(c for c in test if c.mutant == row["proposed"])
    assert row["measured_effect"] == pytest.approx(proposed.measured)   # the WET-LAB value, not the score
    assert row["abs_err_top1"] == pytest.approx(abs(proposed.measured - (-1.0)))


def test_a_useless_oracle_does_not_beat_the_null():
    """Guard against a design that flatters any scorer: a constant score carries no information, so the
    inverse must not out-hit the empirical null with it."""
    cands = _cands(60)
    cal, test = cands[::2], cands[1::2]
    targets = list(statistics.quantiles([c.measured for c in cands], n=5))
    useless = run_inverse("useless", cal, test, lambda c: 0.0, targets, top_k=5)
    null = empirical_null(test, targets, top_k=5)
    assert useless["mean_abs_err_best_of_k"] >= null["mean_abs_err_best_of_k"] * 0.95


def test_a_perfect_oracle_beats_the_null():
    cands = _cands(60)
    cal, test = cands[::2], cands[1::2]
    targets = list(statistics.quantiles([c.measured for c in cands], n=5))
    perfect = run_inverse("perfect", cal, test, lambda c: -c.pos, targets, top_k=5)
    null = empirical_null(test, targets, top_k=5)
    assert perfect["mean_abs_err_best_of_k"] < null["mean_abs_err_best_of_k"]


# ---- the conformal honesty flag (coverage != informativeness) --------------------------------------------

def test_interval_informativeness_is_reported_separately_from_coverage():
    """J2's rail, restated: split-conformal coverage holds even for a USELESS model. So `brackets` is not
    evidence — the halfwidth relative to the effect span is. Both must be present in the output."""
    cands = _cands(60)
    cal, test = cands[::2], cands[1::2]
    r = run_inverse("x", cal, test, lambda c: 0.0, [-1.0], top_k=1)
    assert "interval_halfwidth_over_effect_span" in r
    assert "interval_is_informative" in r
    # a constant score => the calibrator is flat => wide interval relative to the span => NOT informative
    assert r["interval_is_informative"] is False


def test_a_perfect_oracle_yields_an_informative_interval():
    cands = _cands(60)
    cal, test = cands[::2], cands[1::2]
    r = run_inverse("x", cal, test, lambda c: -c.pos, [-1.0], top_k=1)
    assert r["interval_is_informative"] is True


# ---- the empirical null ----------------------------------------------------------------------------------

def test_empirical_null_best_of_k_is_better_than_a_single_draw():
    cands = _cands(50)
    n = empirical_null(cands, [-2.0], top_k=5)["per_target"][0]
    assert n["expected_abs_err_best_of_k"] < n["expected_abs_err_random_pick"]


def test_empirical_null_is_exact_not_sampled():
    """No RNG: the same input must give byte-identical output (a sampled null would drift run to run)."""
    cands = _cands(50)
    a = empirical_null(cands, [-1.0, -2.0], top_k=3)
    b = empirical_null(cands, [-1.0, -2.0], top_k=3)
    assert a == b


# ---- the candidate space ---------------------------------------------------------------------------------

def test_single_nt_accessible_excludes_unreachable_substitutions():
    # ATG(M) -> one nt change can reach L/V/I/T/K/R, but NOT e.g. W (needs 2 changes)
    acc = single_nt_accessible("M", "ATG")
    alts = {m[-1] for m in acc}
    assert "L" in alts and "V" in alts and "I" in alts
    assert "W" not in alts
    assert "M" not in alts                     # never the wild-type itself


def test_single_nt_accessible_excludes_nonsense():
    """A stop codon is not a fitness-tuning edit; the inverse must not propose one."""
    acc = single_nt_accessible("Y", "TAT")     # TAT(Y) -> TAA/TAG are stops, one nt away
    assert all(m[-1] != "*" for m in acc)


def test_single_nt_accessible_skips_a_codon_that_disagrees_with_the_protein():
    """Coordinate-frame guard: if the CDS codon does not translate to the stated residue, skip it rather
    than emit a wrong candidate."""
    assert single_nt_accessible("M", "GGG") == set()    # GGG is Gly, not Met -> refuse


# ---- real data -------------------------------------------------------------------------------------------

@pytest.mark.skipif(not PG.exists(), reason="ProteinGym cache (D:) not mounted")
def test_real_blatem_accessible_set_matches_the_committed_demo_scale():
    """The committed blatem_genome_demo reported 1,715 single-nt-accessible DMS variants; this falsifier
    scores the ~1.5k of them that are single-mutant + in-frame. Pin the scale so a parser regression that
    silently halves the candidate space is caught."""
    from scripts.forward_inverse_roundtrip import load_substrate

    target, dms, esm, cds = load_substrate()
    acc = single_nt_accessible(target, cds)
    scored = [m for m in dms if m in acc]
    assert 1400 <= len(scored) <= 1800, len(scored)
    assert len(esm) == len(target)             # one ESM column per residue
