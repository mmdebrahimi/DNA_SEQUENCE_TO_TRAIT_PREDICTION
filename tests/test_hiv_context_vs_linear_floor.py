"""Guards for the context-vs-linear-floor comparison on the HIV blind spot (2026-09-30)."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
_spec = importlib.util.spec_from_file_location(
    "hiv_context_vs_linear_floor", REPO / "scripts" / "hiv_context_vs_linear_floor.py")
C = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(C)

ART = REPO / "wiki" / "hiv_context_vs_linear_floor_2026-09-30.json"


def _art():
    if not ART.exists():
        pytest.skip("artifact absent (gitignored Stanford data missing at build time)")
    return json.loads(ART.read_text(encoding="utf-8"))


def test_the_bar_was_frozen_and_is_READ_not_restated():
    """A threshold restated in the measuring script can drift from the thing pre-registered."""
    bar = C.load_bar()
    assert bar["frozen_before_any_number_was_seen"] is True
    assert bar["frozen"] == "2026-09-30"
    # The CODE's thresholds must still agree with the frozen bar's prose. Asserted as a consistency
    # check, NOT by counting "0.02" in the source -- that substring also matches the bootstrap percentile
    # 0.025, which is the number-inside-a-number trap recorded in this repo's own lessons.
    assert C.MIN_GAIN == 0.02 and C.ANCHOR_TOL == 0.02
    assert f"+{C.MIN_GAIN}" in bar["verdict_rule"]["CONTEXT_BEATS_FLOOR"], \
        "the code's MIN_GAIN no longer matches the threshold the frozen bar states"
    assert str(C.ANCHOR_TOL) in bar["anchor"]["rule"], \
        "the code's ANCHOR_TOL no longer matches the frozen bar"


def test_the_prior_is_stated_and_matches_its_own_artifact():
    """The epistasis negative is prior evidence; citing it from memory would be the drift this repo hits."""
    bar = C.load_bar()
    prior = bar["prior_evidence_stated_up_front"]
    p = REPO / prior["artifact"]
    assert p.exists()
    assert json.loads(p.read_text(encoding="utf-8"))["verdict"] == prior["verdict"]
    assert "why_this_run_is_not_a_duplicate" in prior


def test_the_floor_is_READ_from_the_deployability_artifact():
    f = C.floor_targets()
    assert set(f) == {"NNRTI", "INSTI"}, "PI must be excluded -- its blind spot is unscoreable"
    assert f["NNRTI"]["gene"] == "RT" and f["INSTI"]["gene"] == "IN"


def test_the_linear_arm_REPRODUCES_the_measured_floor_exactly():
    """A context-vs-floor delta measured against a floor that does not reproduce is not a comparison."""
    d = _art()
    assert d["anchor"], "no anchor recorded"
    for cls, an in d["anchor"].items():
        assert an["reproduces"] is True, f"{cls}: {an}"
        assert abs(an["delta"]) <= 0.02, f"{cls}: floor drifted by {an['delta']}"


def test_every_arm_saw_the_same_rows_folds_and_mask():
    """Fairness is the whole comparison: only the estimator may differ."""
    d = _art()
    for cls, c in d["cells"].items():
        if c.get("status") != "SCORED":
            continue
        # one blind-spot mask per gene, shared by every arm -> one n and one R, not per-arm
        assert isinstance(c["blind_spot_n"], int) and isinstance(c["blind_spot_R"], int)
        arms = set(c["blind_spot_auroc"])
        assert {"linear_floor", "nonlinear"} <= arms
        assert set(c["delta_vs_floor"]) == arms - {"linear_floor"}


def test_the_verdict_is_the_frozen_rule_applied_mechanically():
    d = _art()
    bar = C.load_bar()
    recomputed = C.verdict(d["cells"], bar)
    assert recomputed["verdict"] == d["verdict"], "the artifact's verdict must be reproducible"
    assert d["verdict"] in bar["verdict_rule"], "the verdict must be one the frozen bar defines"


def test_context_did_not_beat_the_floor_and_the_deltas_are_reported_either_way():
    """The measured outcome: no context arm clears the bar. Deltas ship regardless of sign, so a future
    re-run that DOES win is visible rather than hidden behind a verdict string."""
    d = _art()
    assert d["verdict"] == "LINEAR_FLOOR_IS_THE_CEILING"
    for cls, c in d["cells"].items():
        if c.get("status") != "SCORED":
            continue
        for arm, dd in c["delta_vs_floor"].items():
            assert "delta" in dd and "ci_lo" in dd and "ci_hi" in dd
            assert dd["ci_positive"] is False, f"{cls}/{arm} is CI-positive -- the verdict must change"


def test_the_negative_is_not_merely_a_null_two_arms_are_significantly_WORSE():
    """Distinguishes 'added capacity does not pay' from 'no signal either way' -- a CI entirely below zero
    is a stronger statement than a delta near zero, and the memo must not round it to 'no difference'."""
    d = _art()
    worse = [(cls, arm) for cls, c in d["cells"].items() if c.get("status") == "SCORED"
             for arm, dd in c["delta_vs_floor"].items() if dd["ci_hi"] is not None and dd["ci_hi"] < 0]
    assert len(worse) >= 2, f"expected >=2 significantly-worse arms, got {worse}"


def test_the_pairwise_arm_selected_REAL_biology_which_is_what_makes_the_negative_strong():
    """`build_pairs` ranks by co-occurrence, and on RT it recovered the canonical TAM cluster (41L/210W/
    215Y) plus M184V. The negative is 'real interactions did not help', not 'we found only noise'."""
    d = _art()
    rt = d["cells"]["NNRTI"]
    tops = " ".join(rt.get("top_pairs_by_cooccurrence") or [])
    assert rt.get("n_pairs_kept", 0) > 0
    hits = sum(1 for m in ("41", "210", "215", "184") if m in tops)
    assert hits >= 3, f"expected canonical TAM/184 positions among the top pairs; got {tops}"


def test_the_prevalence_selection_caveat_is_recorded_not_buried():
    """NAMED LIMITATION: the pairs are chosen among the most PREVALENT tokens, which in an NNRTI cohort of
    treatment-experienced patients are NRTI TAMs -- not EFV drivers. So the pairwise arm is handicapped and
    the general claim rests on the NON-LINEAR arm, which sees every feature. If this caveat is ever dropped
    the negative would read stronger than it is."""
    d = _art()
    text = json.dumps(d.get("honest_scope", [])) + json.dumps(d.get("named_limitations", []))
    assert "prevalen" in text.lower(), "the pair-selection caveat must ship in the artifact"
    assert "nonlinear" in text.lower() or "non-linear" in text.lower()


def test_the_claim_is_scoped_to_this_feature_space_not_to_hiv_biology():
    """'Additive suffices for HIV' would overstate it: with ~1063 features and 844-1883 blind-spot
    isolates, capacity may simply be unaffordable out-of-distribution."""
    d = _art()
    scope = " ".join(d.get("honest_scope", [])).lower()
    assert "feature space" in scope or "sample size" in scope or "capacity" in scope
    assert "raw" in scope and "sequence" in scope, \
        "must state that this does not test a sequence model over raw residues"
