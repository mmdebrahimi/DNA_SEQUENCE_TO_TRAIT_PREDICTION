"""Tests for the as-given arm's bar.

The refusals matter for the same reason they did on the positional bar: a bar frozen after a run is
not a bar, and a bar anchored to retyped numbers is anchored to nothing.

**But this bar carries one thing its sibling did not, and it is the most important thing to test:
the number it grades is ALREADY VISIBLE.** 0.5086 was published as the positional run's ineligible
secondary before this arm existed. A bar that failed to say so would read as a blind
pre-registration, which would be a false claim about the evidence's strength -- so `verify()`
refuses to freeze unless the positional artifact actually carries that rung as an INELIGIBLE one,
and tests assert the hazard is present with its figures rather than as an adjective.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts.freeze_asgiven_arm_bar import (
    ARM,
    CHECK_ORDER,
    OUT,
    OUTCOMES,
    PRIMARY,
    build_bar,
    classify,
    reachable_outcomes,
    verify,
)
from scripts.freeze_positional_arm_bar import MARGIN, BarRefused, ceiling_arms_and_criterion

REPO = Path(__file__).resolve().parents[1]
RECONCILE = REPO / "data/promoteratlas/local_reconcile_2026-11-01.json"
LEDGER = REPO / "soraya_runs/arms-ledger.json"
POS = REPO / "wiki/positional_arm_result_2026-11-01.json"


def _copies(tmp_path: Path):
    rec = json.loads(RECONCILE.read_text(encoding="utf-8"))
    led = json.loads(LEDGER.read_text(encoding="utf-8"))
    pos = json.loads(POS.read_text(encoding="utf-8"))
    rp, lp, pp = tmp_path / "rec.json", tmp_path / "led.json", tmp_path / "pos.json"
    for p, o in ((rp, rec), (lp, led), (pp, pos)):
        p.write_text(json.dumps(o), encoding="utf-8")
    return rec, led, pos, rp, lp, pp


def _w(p: Path, o) -> Path:
    p.write_text(json.dumps(o), encoding="utf-8")
    return p


# --------------------------------------------------------------------------------------------------
# the closed outcome set
# --------------------------------------------------------------------------------------------------
def test_the_reachable_outcome_set_equals_the_closed_set():
    assert reachable_outcomes() == set(OUTCOMES)


def test_equality_not_containment_is_the_check():
    padded = list(CHECK_ORDER) + [("a_gate_nobody_calls", "AN_OUTCOME_NOBODY_REACHES")]
    assert reachable_outcomes(padded) > set(OUTCOMES)
    assert reachable_outcomes(padded) != set(OUTCOMES)


def test_validity_gates_precede_the_decision_and_split_equality_is_first():
    gates = [g for g, _ in CHECK_ORDER]
    assert gates[0] == "split_equality"
    assert gates[-1] == "primary_decision"


# --------------------------------------------------------------------------------------------------
# classify -- exhaustive, ordered, and the (0, margin] absorption is pre-committed
# --------------------------------------------------------------------------------------------------
def test_an_interval_above_the_margin_is_an_ADD():
    assert classify(0.03, 0.09)[0] == "ASGIVEN_ADDS"


def test_an_interval_below_zero_is_a_NO_GAIN():
    assert classify(-0.05, -0.01)[0] == "ASGIVEN_NO_GAIN"


def test_an_interval_containing_zero_is_a_TIE():
    assert classify(-0.0079, 0.0017)[0] == "ASGIVEN_TIES"


def test_a_SUB_MARGIN_GAIN_is_a_TIE_and_that_is_decided_in_advance():
    """The case that would otherwise invite a post-hoc re-read: a real but sub-margin gain."""
    v, why = classify(0.005, 0.015)
    assert v == "ASGIVEN_TIES"
    assert "neither clears" in why
    # and it must NOT be reported as an add even though the whole interval is positive
    assert v != "ASGIVEN_ADDS"


def test_the_margin_is_a_strict_boundary():
    assert classify(MARGIN, 0.5)[0] == "ASGIVEN_TIES"
    assert classify(MARGIN + 1e-9, 0.5)[0] == "ASGIVEN_ADDS"
    assert classify(-0.5, -1e-9)[0] == "ASGIVEN_NO_GAIN"
    assert classify(-0.5, 0.0)[0] == "ASGIVEN_TIES"


@pytest.mark.parametrize("lo,hi", [(-0.3, -0.2), (-0.1, 0.1), (0.005, 0.015), (0.03, 0.2)])
def test_classify_is_total_and_every_outcome_is_in_the_closed_set(lo, hi):
    v, why = classify(lo, hi)
    assert v in OUTCOMES
    assert why.strip()


# --------------------------------------------------------------------------------------------------
# verify -- the refusals
# --------------------------------------------------------------------------------------------------
def test_verify_passes_on_the_committed_anchors():
    v = verify()
    assert v["reconcile"]["verdict"] == "RECONCILED"
    assert v["arms"][ARM]["status"] in ("planned", "measured")


def test_verify_refuses_a_mutated_comparator(tmp_path):
    rec, _l, _p, rp, lp, pp = _copies(tmp_path)
    rec["atlas_median"] = 0.5100
    with pytest.raises(BarRefused, match="atlas median drifted"):
        verify(reconcile_path=_w(rp, rec), ledger_path=lp, positional_result=pp)


def test_verify_refuses_an_unreconciled_comparator(tmp_path):
    rec, _l, _p, rp, lp, pp = _copies(tmp_path)
    rec["verdict"] = "MISMATCHED"
    with pytest.raises(BarRefused, match="the comparator is"):
        verify(reconcile_path=_w(rp, rec), ledger_path=lp, positional_result=pp)


def test_verify_refuses_a_mutated_arm_status(tmp_path):
    _r, led, _p, rp, lp, pp = _copies(tmp_path)
    arms, _c = ceiling_arms_and_criterion(led)
    arms[ARM]["status"] = "measured"
    with pytest.raises(BarRefused, match="not 'planned'"):
        verify(reconcile_path=rp, ledger_path=_w(lp, led), positional_result=pp)


def test_verify_refuses_an_arm_that_already_carries_a_metric(tmp_path):
    _r, led, _p, rp, lp, pp = _copies(tmp_path)
    arms, _c = ceiling_arms_and_criterion(led)
    arms[ARM]["metric"] = 0.5086
    with pytest.raises(BarRefused, match="already carries metric"):
        verify(reconcile_path=rp, ledger_path=_w(lp, led), positional_result=pp)


def test_verify_refuses_a_missing_arm(tmp_path):
    _r, led, _p, rp, lp, pp = _copies(tmp_path)
    arms, _c = ceiling_arms_and_criterion(led)
    del arms[ARM]
    with pytest.raises(BarRefused, match="is NOT declared"):
        verify(reconcile_path=rp, ledger_path=_w(lp, led), positional_result=pp)


def test_verify_refuses_when_the_positional_artifact_lacks_the_rung(tmp_path):
    """Without that rung the hazard cannot be stated WITH A FIGURE, and the bar must not ship."""
    _r, _l, pos, rp, lp, pp = _copies(tmp_path)
    del pos["rungs"][PRIMARY]
    with pytest.raises(BarRefused, match="no scored"):
        verify(reconcile_path=rp, ledger_path=lp, positional_result=_w(pp, pos))


def test_verify_refuses_if_the_rung_was_not_INELIGIBLE_in_the_positional_arm(tmp_path):
    """The arm's whole premise is secondary-there, primary-here."""
    _r, _l, pos, rp, lp, pp = _copies(tmp_path)
    pos["rungs"][PRIMARY]["eligible_as_the_arms_result"] = True
    with pytest.raises(BarRefused, match="not recorded as INELIGIBLE"):
        verify(reconcile_path=rp, ledger_path=lp, positional_result=_w(pp, pos))


def test_verify_refuses_an_absent_anchor(tmp_path):
    with pytest.raises(BarRefused, match="anchor source absent"):
        verify(reconcile_path=tmp_path / "nope.json")


# --------------------------------------------------------------------------------------------------
# the written bar -- the hazard must be a FIGURE, not an adjective
# --------------------------------------------------------------------------------------------------
@pytest.mark.skipif(not OUT.is_file(), reason="bar not frozen")
def test_the_bar_states_the_observed_before_declared_hazard_with_its_numbers():
    bar = json.loads(OUT.read_text(encoding="utf-8"))
    h = bar["NOT_A_BLIND_PRE_REGISTRATION"]
    pos = json.loads(POS.read_text(encoding="utf-8"))["rungs"][PRIMARY]

    assert "not a blind pre-registration" in h["statement"].lower()
    assert h["already_observed_rung"] == PRIMARY
    assert h["already_observed_median"] == pos["median_spearman"]
    assert h["already_observed_gain"] == pos["gain_bootstrap"]["median"]
    assert h["already_observed_ci"] == [pos["gain_bootstrap"]["ci_lo"],
                                        pos["gain_bootstrap"]["ci_hi"]]
    assert "LOWER" in h["what_it_costs"], "the cost of the hazard must be stated, not implied"
    assert h["what_survives"].strip()


@pytest.mark.skipif(not OUT.is_file(), reason="bar not frozen")
def test_the_bar_declares_one_primary_and_no_correction_and_no_secondaries():
    m = json.loads(OUT.read_text(encoding="utf-8"))["multiplicity"]
    assert m["n_primaries"] == 1
    assert m["primary"] == PRIMARY
    assert m["correction"].startswith("NONE")
    assert m["ci_percentiles"] == [2.5, 97.5]
    assert m["no_secondaries_declared"] is True


@pytest.mark.skipif(not OUT.is_file(), reason="bar not frozen")
def test_the_bars_margin_is_the_inherited_value():
    bar = json.loads(OUT.read_text(encoding="utf-8"))
    assert bar["inherited_not_chosen"]["margin"]["value"] == MARGIN == 0.02
    assert "INHERITED" in bar["inherited_not_chosen"]["margin"]["derivation"]


@pytest.mark.skipif(not OUT.is_file(), reason="bar not frozen")
def test_the_bar_pins_the_committed_comparator_verbatim():
    bar = json.loads(OUT.read_text(encoding="utf-8"))
    rec = json.loads(RECONCILE.read_text(encoding="utf-8"))
    comp = bar["inherited_not_chosen"]["comparator"]
    assert comp["atlas_per_seed"] == rec["atlas_per_seed"]
    assert comp["must_be_recomputed_in_run"] is False
    assert bar["inherited_not_chosen"]["split"]["n_test_per_seed"] == rec["n_test_per_seed"]


@pytest.mark.skipif(not OUT.is_file(), reason="bar not frozen")
def test_the_bar_predicts_a_TIE_and_pre_commits_the_sub_margin_rule():
    bar = json.loads(OUT.read_text(encoding="utf-8"))
    f = bar["forecast"]
    assert f["predicted_outcome"] == "ASGIVEN_TIES"
    assert sum((f["p_asgiven_ties"], f["p_asgiven_adds"], f["p_asgiven_no_gain"])) == pytest.approx(1.0)
    assert "nearly free" in f["my_forecast_record_in_this_project"], (
        "a forecast of an already-visible number must say so")
    assert any("sub-margin gain returns ASGIVEN_TIES" in p for p in bar["pre_commitments"])
    assert any("NOT to re-report" in p for p in bar["pre_commitments"])


@pytest.mark.skipif(not OUT.is_file(), reason="bar not frozen")
def test_the_bar_names_the_1_bit_orientation_inference_as_a_limit():
    limits = " ".join(json.loads(OUT.read_text(encoding="utf-8"))["honest_limits_fixed_in_advance"])
    assert "1-BIT INFERENCE" in limits.upper()
    assert "blind pre-registration" in limits
    assert "0.941" in limits


def test_build_bar_is_a_pure_function_of_the_verified_anchors():
    a, b = build_bar(verify()), build_bar(verify())
    a.pop("frozen"), b.pop("frozen")
    assert a == b
