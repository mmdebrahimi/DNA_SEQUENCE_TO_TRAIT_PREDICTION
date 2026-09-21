"""Guards for the CRyPTIC implied-boundary measurement.

The load-bearing property is that NO breakpoint or ECOFF is recalled from memory anywhere in this arm:
the clinical cut-off is RECOVERED from the shipped BINARY_PHENOTYPE as a pure step on the MIC ladder.
That recovery is only legitimate while the step really is pure, so `recover_cutoff` must REFUSE on an
impure one rather than pick the nearest rung -- most of these tests exist for that refusal.
"""
from __future__ import annotations

import importlib.util
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

ARTIFACT = ROOT / "wiki" / "tb_implied_boundary_2026-09-21.json"
BAR = ROOT / "wiki" / "tb_implied_boundary_acceptance_bar.json"


def _mod():
    spec = importlib.util.spec_from_file_location(
        "tb_implied_boundary", ROOT / "scripts" / "tb_implied_boundary.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


# --- the cut-off is recovered, never recalled -----------------------------------------------------

def test_a_clean_step_recovers_the_cutoff():
    m = _mod()
    rungs = np.array([-2.0, -2.0, -1.0, -1.0, 0.0, 0.0, 1.0])
    pheno = np.array(["S", "S", "S", "S", "R", "R", "R"])
    assert m.recover_cutoff(rungs, pheno) == -1.0


def test_an_impure_step_REFUSES_rather_than_guessing():
    """THE guard. A rung that is part S and part R means BINARY_PHENOTYPE is not a MIC threshold, so
    there is no cut-off to recover. Picking the nearest rung would silently manufacture the very
    number this arm exists to avoid recalling."""
    m = _mod()
    rungs = np.array([-2.0, -2.0, -1.0, -1.0, 0.0, 0.0])
    pheno = np.array(["S", "S", "S", "R", "R", "R"])      # -1.0 is mixed
    assert m.recover_cutoff(rungs, pheno) is None


def test_a_resistant_rung_below_a_susceptible_one_refuses():
    """Non-monotone phenotype: R at a LOW rung and S above it is not a threshold at all."""
    m = _mod()
    rungs = np.array([-2.0, -1.0, 0.0])
    pheno = np.array(["R", "S", "R"])
    assert m.recover_cutoff(rungs, pheno) is None


def test_an_all_resistant_ladder_has_no_recoverable_cutoff():
    m = _mod()
    assert m.recover_cutoff(np.array([0.0, 1.0]), np.array(["R", "R"])) is None


# --- censoring is parsed, never coerced to a point value ------------------------------------------

def test_every_censoring_form_parses_with_its_direction():
    """Operator-aware censoring: a bound is a real measurement of a bound. The direction must survive
    so the boundary rung can be ASSERTED uncensored rather than assumed."""
    m = _mod()
    assert m.parse_mic("0.5") == (pytest.approx(-1.0), "exact")
    assert m.parse_mic("<=0.06") == (pytest.approx(np.log2(0.06)), "left")
    assert m.parse_mic(">4") == (pytest.approx(2.0), "right")
    assert m.parse_mic("<0.03")[1] == "left"
    assert m.parse_mic(">=8")[1] == "right"


def test_the_ladder_is_log2_so_rungs_are_one_apart():
    m = _mod()
    a, _ = m.parse_mic("0.25")
    b, _ = m.parse_mic("0.5")
    assert b - a == pytest.approx(1.0)


# --- the permutation null actually controls -------------------------------------------------------

def test_an_unenriched_boundary_does_not_clear_its_own_null():
    """NON-VACUITY for the control: if carriage is spread evenly, the boundary must NOT be significant."""
    m = _mod()
    rng = np.random.default_rng(0)
    carries = rng.integers(0, 2, size=2000)
    boundary = np.zeros(2000, bool)
    boundary[:200] = True
    assert not m.permutation_gap(boundary, carries, np.random.default_rng(1), 200)["exceeds_perm_max"]


def test_a_strongly_enriched_boundary_does_clear_it():
    m = _mod()
    carries = np.zeros(2000, int)
    boundary = np.zeros(2000, bool)
    boundary[:200] = True
    carries[:150] = 1                                     # concentrated in the boundary rung
    assert m.permutation_gap(boundary, carries, np.random.default_rng(1), 200)["exceeds_perm_max"]


# --- the frozen verdict rule ----------------------------------------------------------------------

def _drug(bc, lc, clears, n=200, cutoff=-1.0):
    return {"cutoff_log2": cutoff, "boundary_n": n, "boundary_carriage": bc,
            "lower_carriage": lc, "permutation": {"exceeds_perm_max": clears}}


def test_all_four_frozen_branches():
    m = _mod()
    assert m.verdict_from_bar({"a": _drug(.4, .03, True), "b": _drug(.2, .01, True)}) == m.SUPPORTED
    assert m.verdict_from_bar({"a": _drug(.4, .03, True), "b": _drug(.2, .01, False)}) == m.WEAK
    assert m.verdict_from_bar({"a": _drug(.4, .03, True), "b": _drug(.01, .05, True)}) == m.FALSIFIED
    assert m.verdict_from_bar({"a": _drug(.4, .03, True, n=99)}) == m.INDETERMINATE


def test_an_unrecoverable_cutoff_is_indeterminate_not_falsified():
    """A refused cut-off means the question could not be ASKED on that drug -- which is a different
    fact from the claim failing, and must not be reported as a refutation."""
    m = _mod()
    assert m.verdict_from_bar({"a": _drug(.4, .03, True, cutoff=None)}) == m.INDETERMINATE


def test_one_drug_failing_falsifies_the_whole_claim():
    """The claim is stated over BOTH drugs, so a single failure is a failure -- not an average."""
    m = _mod()
    assert m.verdict_from_bar({"a": _drug(.4, .03, True), "b": _drug(.02, .02, True)}) == m.FALSIFIED


# --- the two composition statistics are COMPLEMENTARY ---------------------------------------------

def test_gene_share_is_vacuous_when_every_determinant_shares_one_gene():
    """Pinned so nobody reads rifampicin's `rpoB 1.0 vs 1.0` as 'no compositional shift'. Every RIF
    determinant is in rpoB, so the per-gene statistic CANNOT move there; its shift is at the allele
    level, which the top-1 comparison catches instead."""
    m = _mod()
    g = m._gene_share(Counter({"rpoB_p.Asp435Tyr": 20}), Counter({"rpoB_p.Ser450Leu": 2898}))
    assert g["boundary"] == {"rpoB": 1.0} and g["high_mic"] == {"rpoB": 1.0}
    assert g["largest_shift"] == 0.0 and not g["largest_shift_is_tied"]


def test_gene_share_catches_a_shift_the_top1_comparison_misses():
    """The isoniazid case: the commonest determinant is the SAME at both ends (so top-1 reports no
    shift) while the per-gene mix moves by ~0.29. Neither statistic alone detects both drugs."""
    m = _mod()
    boundary = Counter({"katG_p.Ser315Thr": 22, "inhA_c.-777C>T": 16, "inhA_c.-154G>A": 7})
    high = Counter({"katG_p.Ser315Thr": 4488, "inhA_c.-777C>T": 915, "inhA_c.-154G>A": 138})
    g = m._gene_share(boundary, high)
    assert g["boundary"]["inhA"] > g["high_mic"]["inhA"]
    # With exactly two genes the shares sum to 1, so BOTH shift by the same amount. The magnitude is
    # the real quantity; asserting a single winner would have pinned an arbitrary sort order.
    assert g["largest_shift"] > 0.2
    assert g["largest_shift_is_tied"] and g["largest_shift_genes"] == ["inhA", "katG"]


def test_gene_is_taken_lexically_so_no_biology_is_asserted():
    m = _mod()
    g = m._gene_share(Counter({"someGene_p.Xyz1Abc": 3}), Counter({"other_c.-1A>T": 1}))
    assert set(g["boundary"]) == {"someGene"} and set(g["high_mic"]) == {"other"}


# --- the committed run ----------------------------------------------------------------------------

@pytest.mark.skipif(not ARTIFACT.exists(), reason="run artifact not generated")
def test_the_committed_verdict_is_re_derivable_from_its_own_numbers():
    m = _mod()
    d = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    assert d["verdict"] == m.verdict_from_bar(d["results"]) == m.SUPPORTED


@pytest.mark.skipif(not ARTIFACT.exists(), reason="run artifact not generated")
def test_the_boundary_rung_is_uncensored_on_every_drug():
    """If a boundary rung were a censored bound, its carriage would be a statement about a bound
    rather than about that rung, and the headline would not mean what it says."""
    d = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    for drug, r in d["results"].items():
        assert r["boundary_is_uncensored"] is True, drug


@pytest.mark.skipif(not ARTIFACT.exists(), reason="run artifact not generated")
def test_each_gap_beats_its_permutation_MAXIMUM_not_merely_p95():
    d = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    for drug, r in d["results"].items():
        p = r["permutation"]
        assert p["observed_gap"] > p["perm_max"], drug


@pytest.mark.skipif(not ARTIFACT.exists(), reason="run artifact not generated")
def test_the_recovered_cutoffs_are_reported_with_their_provenance():
    """The cut-off provenance sentence is the reason this arm needs no external value; if it ever
    stops being true the artifact must not keep claiming it."""
    d = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    assert "RECOVERED" in d["cutoff_provenance"] and "memory" in d["cutoff_provenance"]
    for r in d["results"].values():
        assert r["cutoff_log2"] is not None and r["cutoff_mg_L"] > 0


@pytest.mark.skipif(not BAR.exists(), reason="bar not present")
def test_the_bar_was_frozen_after_the_mechanism_and_names_every_outcome():
    """Two bars earlier in this project were mis-specified by freezing a threshold BEFORE the analysis
    that says what a correct result can look like. This one records that it was frozen after."""
    b = json.loads(BAR.read_text(encoding="utf-8"))
    assert "FROZEN AFTER THE MECHANISM ANALYSIS" in b["status"]
    assert set(b["verdict_rule"]) == {"SUPPORTED", "WEAK_DIRECTIONAL", "FALSIFIED", "INDETERMINATE"}
    assert any("IN-DISTRIBUTION" in s for s in b["declared_in_advance"])


@pytest.mark.skipif(not ARTIFACT.exists(), reason="run artifact not generated")
def test_the_in_distribution_limit_survives_into_the_result():
    """The WHO catalogue was built partly from CRyPTIC. A reader must not take this for independent
    validation of the catalogue, so the limit ships with the numbers rather than only in the bar."""
    d = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    assert any("IN-DISTRIBUTION" in s for s in d["honest_limits"])
