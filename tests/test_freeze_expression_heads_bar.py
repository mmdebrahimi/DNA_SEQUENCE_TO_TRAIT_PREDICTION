"""Tests for the expression-heads acceptance bar.

The refusals matter for the usual reason: a bar frozen after a run is not a bar, and a bar anchored
to retyped numbers is anchored to nothing.

Two tests here carry more weight than the rest.

**The route rule must be able to choose EITHER way.** On the real probe it rejects the cheap route,
because one head's rank agreement (0.9888) falls below the margin-derived tolerance (0.999). A rule
that could only ever reject would be a disguised constant, so it is also driven with a synthetic
probe whose heads all clear the tolerance, and it must then accept.

**The secondary must be unable to reach the verdict.** `classify` takes an interval and nothing
else, so the only way a max-over-windows rung could enter a verdict is if a caller passed it in.
The bar declares it ineligible and the test drives a secondary that WOULD win to show the verdict
function has no door for it.
"""

from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from scripts.freeze_expression_heads_bar import (
    ARM,
    CHECK_ORDER,
    NULL_SIGNAL_RATIO,
    OUT,
    OUTCOMES,
    PRIMARY_REDUCTION,
    SCORING_DIRECT,
    SCORING_RIDGE_INELIGIBLE,
    ROUTE_CHEAP,
    ROUTE_RANK_TOLERANCE,
    ROUTE_TRUE,
    SECONDARY_REDUCTION,
    build_bar,
    choose_route,
    classify,
    reachable_outcomes,
    verify,
)
from scripts.freeze_positional_arm_bar import MARGIN, BarRefused, ceiling_arms_and_criterion

REPO = Path(__file__).resolve().parents[1]
RECONCILE = REPO / "data/promoteratlas/local_reconcile_2026-11-01.json"
LEDGER = REPO / "soraya_runs/arms-ledger.json"
PROBE = REPO / "wiki/expression_heads_feasibility_2026-11-02.json"


def _as_at_freeze_time(led: dict) -> dict:
    """The ledger as it stood when this bar was frozen: the arm declared but NOT yet measured.

    WHY THIS IS NEEDED AND WHY IT IS NOT A RELAXATION. `verify()` refuses to freeze a bar for an
    arm that already carries a number -- correctly, because a bar written after the result is not a
    pre-registration. The arm was recorded MEASURED on 2026-11-04 at 0.1178, so every test calling
    `verify()` against the LIVE ledger began failing with "the arm's status is 'measured', not
    'planned'", and waiting would never fix it: the question `verify()` asks was answered at freeze
    time and cannot be re-asked of a later ledger.

    The guard is untouched and still fires on the real file -- see
    `test_verify_REFUSES_against_the_LIVE_ledger_now_that_the_arm_is_measured`, which is what makes
    this a scoping decision rather than a way to dodge a red test.
    """
    led = json.loads(json.dumps(led))
    arms, _crit = ceiling_arms_and_criterion(led)
    if ARM in arms:
        arms[ARM].update(status="planned", metric=None, evidence="", resolved_at="")
    return led


def _copies(tmp_path: Path, *, at_freeze_time: bool = True):
    rec = json.loads(RECONCILE.read_text(encoding="utf-8"))
    led = json.loads(LEDGER.read_text(encoding="utf-8"))
    if at_freeze_time:
        led = _as_at_freeze_time(led)
    pr = json.loads(PROBE.read_text(encoding="utf-8"))
    rp, lp, pp = tmp_path / "rec.json", tmp_path / "led.json", tmp_path / "probe.json"
    for p, o in ((rp, rec), (lp, led), (pp, pr)):
        p.write_text(json.dumps(o), encoding="utf-8")
    return rec, led, pr, rp, lp, pp


def _verified(tmp_path: Path):
    """`verify()` against the freeze-time snapshot, for the tests that used a bare call."""
    _r, _l, _p, rp, lp, pp = _copies(tmp_path)
    return verify(reconcile_path=rp, ledger_path=lp, probe_path=pp)


def _w(p: Path, o) -> Path:
    p.write_text(json.dumps(o), encoding="utf-8")
    return p


# --------------------------------------------------------------------------------------------------
# the closed outcome set
# --------------------------------------------------------------------------------------------------
def test_the_reachable_outcome_set_EQUALS_the_closed_set():
    assert reachable_outcomes() == set(OUTCOMES)


def test_equality_and_not_containment_is_the_check():
    """The non-vacuity proof: padding the check order must break the equality."""
    padded = list(CHECK_ORDER) + [("a_gate_nobody_runs", "AN_OUTCOME_NOBODY_DECLARED")]
    assert reachable_outcomes(padded) > set(OUTCOMES)
    assert reachable_outcomes(padded) != set(OUTCOMES)


def test_dropping_a_gate_also_breaks_the_equality_in_the_other_direction():
    trimmed = [c for c in CHECK_ORDER if c[0] != "head_validity"]
    assert reachable_outcomes(trimmed) < set(OUTCOMES)


def test_validity_gates_run_before_the_decision_and_split_equality_is_first():
    gates = [g for g, _ in CHECK_ORDER]
    assert gates[0] == "split_equality"
    assert gates[-1] == "primary_decision"
    assert gates.index("head_validity") < gates.index("primary_decision")


# --------------------------------------------------------------------------------------------------
# classify -- total, ordered, and the sub-margin band is pre-committed
# --------------------------------------------------------------------------------------------------
def test_an_interval_above_the_margin_ADDS():
    assert classify(0.03, 0.09)[0] == "HEADS_ADDS"


def test_an_interval_below_zero_is_NO_GAIN():
    assert classify(-0.30, -0.05)[0] == "HEADS_NO_GAIN"


def test_an_interval_containing_zero_TIES():
    assert classify(-0.02, 0.01)[0] == "HEADS_TIES"


def test_a_real_but_SUB_MARGIN_gain_TIES_and_that_is_decided_in_advance():
    v, why = classify(0.005, 0.015)
    assert v == "HEADS_TIES"
    assert "neither clears" in why
    assert v != "HEADS_ADDS"


def test_the_margin_is_a_strict_boundary():
    assert classify(MARGIN, 0.5)[0] == "HEADS_TIES"
    assert classify(MARGIN + 1e-9, 0.5)[0] == "HEADS_ADDS"
    assert classify(-0.5, -1e-9)[0] == "HEADS_NO_GAIN"
    assert classify(-0.5, 0.0)[0] == "HEADS_TIES"


@pytest.mark.parametrize("lo,hi", [(-0.4, -0.3), (-0.1, 0.1), (0.001, 0.019), (0.05, 0.2)])
def test_classify_is_total_and_lands_in_the_closed_set(lo, hi):
    v, why = classify(lo, hi)
    assert v in OUTCOMES and why.strip()


def test_the_verdict_function_has_NO_DOOR_for_a_secondary():
    """A max-over-windows rung can only enter a verdict if a caller hands it in.

    `classify` takes an interval and a margin; it has no rung argument, so the ineligibility of the
    secondary is structural rather than a convention someone has to remember.
    """
    import inspect

    params = set(inspect.signature(classify).parameters)
    assert params == {"ci_lo", "ci_hi", "margin"}
    assert SECONDARY_REDUCTION not in params
    # and a secondary interval that WOULD win cannot be distinguished from any other interval,
    # which is exactly why eligibility has to be enforced upstream by the scorer
    assert classify(0.10, 0.20)[0] == "HEADS_ADDS"


# --------------------------------------------------------------------------------------------------
# the route rule -- BOTH directions
# --------------------------------------------------------------------------------------------------
def _probe_with_rhos(rhos: dict, *, cheap=0.24, true_=4.64) -> dict:
    return {
        "routes": {name: {"per_reduction": {PRIMARY_REDUCTION:
                                            {"spearman_cheap_vs_true": r}}}
                   for name, r in rhos.items()},
        "cost_projection": {"cheap_route": {"hours_all_heads": cheap},
                            "true_route": {"hours_all_heads": true_}},
        "geometry": {"n_windows": 65, "arithmetic": "11008 = 128 x 86"},
    }


def test_the_rule_ACCEPTS_the_cheap_route_when_every_head_clears_the_tolerance():
    """Without this the rule would be a disguised constant."""
    out = choose_route(_probe_with_rhos({"a": 0.9995, "b": 0.9999, "c": 0.99991}))
    assert out["chosen"] == ROUTE_CHEAP
    assert out["rejected_alternative"] == ROUTE_TRUE
    assert "reaches the margin-derived tolerance" in out["why"]


def test_the_rule_REJECTS_the_cheap_route_when_ONE_head_falls_short():
    out = choose_route(_probe_with_rhos({"a": 0.9999, "b": 0.9888, "c": 0.9999}))
    assert out["chosen"] == ROUTE_TRUE
    assert out["worst_head"] == "b"
    assert "'probably fine' is not a frozen bar" in out["why"]


def test_one_head_short_is_enough_because_the_metric_is_a_MAXIMUM_over_heads():
    """Four of five clearing it is still a rejection."""
    out = choose_route(_probe_with_rhos({"a": 1.0, "b": 1.0, "c": 1.0, "d": 1.0, "e": 0.95}))
    assert out["chosen"] == ROUTE_TRUE


def test_the_tolerance_is_a_strict_boundary():
    assert choose_route(_probe_with_rhos({"a": ROUTE_RANK_TOLERANCE}))["chosen"] == ROUTE_CHEAP
    just_under = ROUTE_RANK_TOLERANCE - 1e-6
    assert choose_route(_probe_with_rhos({"a": just_under}))["chosen"] == ROUTE_TRUE


def test_the_tolerance_derivation_cites_the_MARGIN_and_not_the_observed_numbers():
    out = choose_route(_probe_with_rhos({"a": 0.9999}))
    d = out["tolerance_derivation"]
    assert "DERIVED FROM THE MARGIN" in d
    assert "not from these numbers" in d
    assert "%.2f" % MARGIN in d


def test_the_route_record_names_the_faithful_route_and_why_it_is_faithful():
    out = choose_route(_probe_with_rhos({"a": 0.9999}))
    assert out["faithful_route_is"] == ROUTE_TRUE
    assert "what APPLYING the head means" in out["why_the_true_route_is_the_faithful_one"]


def test_choose_route_refuses_a_probe_with_no_routes_or_no_agreement_figure():
    with pytest.raises(BarRefused, match="no route measurements"):
        choose_route({"routes": {}})
    with pytest.raises(BarRefused, match="no .* rank agreement"):
        choose_route({"routes": {"a": {"per_reduction": {}}},
                      "cost_projection": {}, "geometry": {}})


# --------------------------------------------------------------------------------------------------
# verify -- the refusals
# --------------------------------------------------------------------------------------------------
def test_verify_passes_on_the_committed_anchors(tmp_path):
    """Against the ledger AS AT FREEZE TIME -- see `_as_at_freeze_time`. The arm has since been
    recorded measured at 0.1178, and `verify()` refuses that by design."""
    v = _verified(tmp_path)
    assert v["reconcile"]["verdict"] == "RECONCILED"
    assert v["n_heads"] == 5
    assert v["arms"][ARM]["status"] == "planned"


def test_verify_REFUSES_against_the_LIVE_ledger_now_that_the_arm_is_measured():
    """The other half, asserted rather than left implicit: the guard is still live on the real
    file. If the live ledger ever read `planned` for this arm again, this fails and says so."""
    with pytest.raises(BarRefused, match="not 'planned'"):
        verify()


def test_verify_refuses_a_mutated_comparator(tmp_path):
    rec, _l, _p, rp, lp, pp = _copies(tmp_path)
    rec["atlas_median"] = 0.5200
    with pytest.raises(BarRefused, match="atlas median drifted"):
        verify(reconcile_path=_w(rp, rec), ledger_path=lp, probe_path=pp)


def test_verify_refuses_an_unreconciled_comparator(tmp_path):
    rec, _l, _p, rp, lp, pp = _copies(tmp_path)
    rec["verdict"] = "MISMATCHED"
    with pytest.raises(BarRefused, match="not RECONCILED"):
        verify(reconcile_path=_w(rp, rec), ledger_path=lp, probe_path=pp)


def test_verify_refuses_a_comparator_recording_mismatches(tmp_path):
    rec, _l, _p, rp, lp, pp = _copies(tmp_path)
    rec["mismatched"] = ["seed3"]
    with pytest.raises(BarRefused, match="records mismatches"):
        verify(reconcile_path=_w(rp, rec), ledger_path=lp, probe_path=pp)


def test_verify_refuses_an_arm_already_measured(tmp_path):
    _r, led, _p, rp, lp, pp = _copies(tmp_path)
    arms, _c = ceiling_arms_and_criterion(led)
    arms[ARM]["status"] = "measured"
    with pytest.raises(BarRefused, match="not 'planned'"):
        verify(reconcile_path=rp, ledger_path=_w(lp, led), probe_path=pp)


def test_verify_refuses_an_arm_that_already_carries_a_metric(tmp_path):
    _r, led, _p, rp, lp, pp = _copies(tmp_path)
    arms, _c = ceiling_arms_and_criterion(led)
    arms[ARM]["metric"] = 0.42
    with pytest.raises(BarRefused, match="already carries metric"):
        verify(reconcile_path=rp, ledger_path=_w(lp, led), probe_path=pp)


def test_verify_refuses_a_missing_arm(tmp_path):
    _r, led, _p, rp, lp, pp = _copies(tmp_path)
    arms, _c = ceiling_arms_and_criterion(led)
    del arms[ARM]
    with pytest.raises(BarRefused, match="is NOT declared"):
        verify(reconcile_path=rp, ledger_path=_w(lp, led), probe_path=pp)


def test_verify_refuses_a_probe_that_did_not_establish_feasibility(tmp_path):
    _r, _l, pr, rp, lp, pp = _copies(tmp_path)
    pr["verdict"] = "INDETERMINATE_COLLAPSED_HEAD"
    with pytest.raises(BarRefused, match="did not establish feasibility"):
        verify(reconcile_path=rp, ledger_path=lp, probe_path=_w(pp, pr))


def test_verify_refuses_a_probe_of_the_wrong_schema(tmp_path):
    _r, _l, pr, rp, lp, pp = _copies(tmp_path)
    pr["schema"] = "something-else-v9"
    with pytest.raises(BarRefused, match="schema"):
        verify(reconcile_path=rp, ledger_path=lp, probe_path=_w(pp, pr))


def test_verify_refuses_a_probe_whose_head_set_SHRANK_between_loading_and_measurement(tmp_path):
    """A head must not be able to go missing quietly -- the divisor is the head count."""
    _r, _l, pr, rp, lp, pp = _copies(tmp_path)
    pr["routes"].pop(sorted(pr["routes"])[0])
    with pytest.raises(BarRefused, match="go missing"):
        verify(reconcile_path=rp, ledger_path=lp, probe_path=_w(pp, pr))


def test_verify_refuses_a_probe_recording_DRIFTED_base_weights(tmp_path):
    _r, _l, pr, rp, lp, pp = _copies(tmp_path)
    name = sorted(pr["heads"])[0]
    pr["heads"][name]["base_tensors_identical"] = 123
    with pytest.raises(BarRefused, match="FINE-TUNED-encoder"):
        verify(reconcile_path=rp, ledger_path=lp, probe_path=_w(pp, pr))


def test_verify_refuses_a_probe_recording_a_wrong_fc1_shape(tmp_path):
    _r, _l, pr, rp, lp, pp = _copies(tmp_path)
    name = sorted(pr["heads"])[0]
    pr["heads"][name]["fc1_shape"] = [64, 19200]
    with pytest.raises(BarRefused, match=r"not \[64, 11008\]"):
        verify(reconcile_path=rp, ledger_path=lp, probe_path=_w(pp, pr))


def test_verify_refuses_an_absent_anchor(tmp_path):
    with pytest.raises(BarRefused, match="anchor source absent"):
        verify(reconcile_path=tmp_path / "nope.json")


# --------------------------------------------------------------------------------------------------
# purity + the divisor
# --------------------------------------------------------------------------------------------------
def test_build_bar_is_a_pure_function_of_the_verified_anchors(tmp_path):
    a, b = build_bar(_verified(tmp_path)), build_bar(_verified(tmp_path))
    a.pop("frozen"), b.pop("frozen")
    assert a == b


def test_the_bonferroni_divisor_EQUALS_the_number_of_declared_heads(tmp_path):
    v = _verified(tmp_path)
    bar = build_bar(v)
    m = bar["multiplicity"]
    assert m["n_primaries"] == v["n_heads"] == 5
    assert "alpha/%d" % v["n_heads"] in m["correction"]
    assert sorted(m["primaries"]) == v["head_names"]
    assert m["divisor_is_read_from_the_probe"] is True


def test_the_divisor_FOLLOWS_the_head_count_so_it_cannot_go_stale(tmp_path):
    """Drop a head from BOTH probe sections and the divisor must move with it."""
    _r, _l, pr, rp, lp, pp = _copies(tmp_path)
    name = sorted(pr["heads"])[0]
    pr["heads"].pop(name)
    pr["routes"].pop(name)
    v = verify(reconcile_path=rp, ledger_path=lp, probe_path=_w(pp, pr))
    bar = build_bar(v)
    assert v["n_heads"] == 4
    assert bar["multiplicity"]["n_primaries"] == 4
    assert "alpha/4" in bar["multiplicity"]["correction"]


def test_the_correction_is_over_heads_and_NOT_over_heads_times_reductions(tmp_path):
    m = build_bar(_verified(tmp_path))["multiplicity"]
    assert m["n_primaries"] == 5, "ten would be heads x reductions"
    assert "costume" in m["why_not_heads_times_reductions"]


# --------------------------------------------------------------------------------------------------
# the written bar
# --------------------------------------------------------------------------------------------------
@pytest.mark.skipif(not OUT.is_file(), reason="bar not frozen")
def test_the_frozen_bar_declares_one_primary_reduction_and_an_INELIGIBLE_secondary():
    r = json.loads(OUT.read_text(encoding="utf-8"))["reductions"]
    assert r["primary"] == PRIMARY_REDUCTION == "mean"
    assert r["secondary"] == SECONDARY_REDUCTION == "max"
    assert r["secondary_is_eligible_as_the_arms_result"] is False
    assert "promotable" in r["why"]


@pytest.mark.skipif(not OUT.is_file(), reason="bar not frozen")
def test_the_frozen_bar_chose_the_FAITHFUL_route_and_recorded_what_that_costs():
    w = json.loads(OUT.read_text(encoding="utf-8"))["windowing"]
    assert w["chosen"] == ROUTE_TRUE
    assert w["worst_rank_agreement"] < w["tolerance"]
    assert w["true_route_hours"] > w["cheap_route_hours"]
    assert "NOT allowed to define the arm" in w["why"]


@pytest.mark.skipif(not OUT.is_file(), reason="bar not frozen")
def test_the_frozen_bar_inherits_the_margin_and_the_comparator_verbatim():
    bar = json.loads(OUT.read_text(encoding="utf-8"))
    rec = json.loads(RECONCILE.read_text(encoding="utf-8"))
    inh = bar["inherited_not_chosen"]
    assert inh["margin"]["value"] == MARGIN == 0.02
    assert "INHERITED" in inh["margin"]["derivation"]
    assert inh["comparator"]["atlas_median"] == 0.5076
    assert inh["comparator"]["atlas_per_seed"] == rec["atlas_per_seed"]
    assert inh["comparator"]["must_be_recomputed_in_run"] is False
    assert inh["split"]["n_test_per_seed"] == rec["n_test_per_seed"]


@pytest.mark.skipif(not OUT.is_file(), reason="bar not frozen")
def test_the_frozen_bar_records_a_forecast_that_sums_to_one_and_declines_its_own_prior():
    f = json.loads(OUT.read_text(encoding="utf-8"))["forecast"]
    assert f["predicted_outcome"] in OUTCOMES
    assert f["p_heads_no_gain"] + f["p_heads_ties"] + f["p_heads_adds"] == pytest.approx(1.0)
    assert "5 of 9" in f["my_forecast_record_in_this_project"]
    assert "Do not defer" in f["my_forecast_record_in_this_project"]


@pytest.mark.skipif(not OUT.is_file(), reason="bar not frozen")
def test_the_frozen_bar_names_every_limit_the_arm_must_be_read_under():
    limits = " ".join(json.loads(OUT.read_text(encoding="utf-8"))["honest_limits_fixed_in_advance"])
    assert "ZERO-SHOT TRANSFER" in limits
    assert "GEOMETRY MISMATCH" in limits
    assert "TRANSLATION" in limits, "the tslpred endpoint mismatch must be a stated limit"
    assert "CIRCULARITY IS SCOPED" in limits and "GSE108535" in limits
    assert "0.941" in limits, "the replicate ceiling must be present so a tie is not read as near it"
    assert "asymmetric AGAINST this arm" in limits, (
        "the incumbent is a FITTED ridge and this arm is not fitted at all; the direction of that "
        "asymmetry must be stated"
    )


@pytest.mark.skipif(not OUT.is_file(), reason="bar not frozen")
def test_the_frozen_bar_pre_commits_the_sub_margin_rule_and_the_route():
    pc = json.loads(OUT.read_text(encoding="utf-8"))["pre_commitments"]
    joined = " ".join(pc)
    assert "sub-margin gain returns HEADS_TIES" in joined
    assert "NOT eligible as the arm's result" in joined
    assert "NOT to be switched after seeing a score" in joined
    assert "bounds TRANSFER" in joined
    assert "never a scored null" in joined


@pytest.mark.skipif(not OUT.is_file(), reason="bar not frozen")
def test_the_frozen_bar_states_the_metric_rule_as_a_max_over_heads():
    bar = json.loads(OUT.read_text(encoding="utf-8"))
    assert "max over the 5 declared heads" in bar["metric_rule"]
    assert "never across arms" in bar["metric_rule"]
    assert bar["reachable_equals_declared"] is True


@pytest.mark.skipif(not OUT.is_file(), reason="bar not frozen")
def test_the_frozen_bar_requires_a_NON_CONSTANT_head_output_as_a_validity_control():
    vc = json.loads(OUT.read_text(encoding="utf-8"))["validity_controls"]
    assert "NON-CONSTANT" in vc["head_validity"]
    assert "clean null" in vc["head_validity"]
    assert set(vc) == {"split_equality", "geometry", "head_validity", "permutation_null"}


# --------------------------------------------------------------------------------------------------
# the amended null -- the defect that would have cost a 4.64 h run
# --------------------------------------------------------------------------------------------------
@pytest.mark.skipif(not OUT.is_file(), reason="bar not frozen")
def test_the_bar_does_NOT_use_the_shuffled_TRAINING_label_null_on_a_one_feature_arm():
    """The inherited null is informative at 1,280 features and VACUOUS at one.

    With a single feature, StandardScaler + Ridge can only scale and sign-flip, so shuffling the
    training labels cannot destroy the ranking -- only its sign. The gate then either fails for a
    reason unrelated to the heads or passes by sign cancellation.
    """
    bar = json.loads(OUT.read_text(encoding="utf-8"))
    gates = [g for g, _ in bar["check_order"]]
    assert "permutation_null" in gates
    assert "shuffled_label_null" not in gates
    assert "shuffled_label_null" not in bar["validity_controls"]


@pytest.mark.skipif(not OUT.is_file(), reason="bar not frozen")
def test_the_null_is_PREDICTOR_side_valid_and_carries_BOTH_its_conditions():
    n = json.loads(OUT.read_text(encoding="utf-8"))["validity_controls"]["permutation_null"]
    assert "FIXED predictor" in n["form"] and "PERMUTED TEST labels" in n["form"]
    assert "NO refit" in n["form"]
    assert "|median null| < 0.02" == n["condition_median"]
    assert "x |real median Spearman|" in n["condition_magnitude"]
    # the magnitude condition is what catches a null that cannot destroy its own signal
    assert "cannot destroy the" in n["condition_magnitude"]


@pytest.mark.skipif(not OUT.is_file(), reason="bar not frozen")
def test_the_bar_records_WHY_the_inherited_null_was_replaced_with_its_measured_numbers():
    n = json.loads(OUT.read_text(encoding="utf-8"))["validity_controls"]["permutation_null"]
    w = n["why_not_the_inherited_shuffled_TRAINING_label_null"]
    assert "ONE feature" in w and "SIGN-FLIP" in w
    assert "0.8814" in w, "the measured pathology must be a figure, not an adjective"
    assert "0.0088" in w, "the 128-feature control case must be recorded so the siblings are scoped"
    assert "0.0371" in n["measured_replacement"]
    assert "1/sqrt(n_test)" in n["measured_replacement"], (
        "the replacement's magnitudes must be anchored to the sampling scale they sit at"
    )


@pytest.mark.skipif(not OUT.is_file(), reason="bar not frozen")
def test_the_magnitude_ceiling_separates_the_broken_null_from_the_fixed_one():
    """Scale-free by construction: ratio 1.00 for the rejected null, 0.04 for the kept one."""
    real, broken, fixed = 0.8814, 0.8814, 0.0371
    assert broken / real > NULL_SIGNAL_RATIO, "the rejected null must trip the ceiling"
    assert fixed / real < NULL_SIGNAL_RATIO, "the kept null must clear it comfortably"


# --------------------------------------------------------------------------------------------------
# the declared scoring definition
# --------------------------------------------------------------------------------------------------
@pytest.mark.skipif(not OUT.is_file(), reason="bar not frozen")
def test_the_bar_DECLARES_how_a_head_scalar_becomes_a_spearman_and_how_sign_is_handled():
    s = json.loads(OUT.read_text(encoding="utf-8"))["scoring"]
    assert s["definition"] == SCORING_DIRECT
    assert "sign AS GIVEN" in s["formula"]
    assert s["an_anticorrelated_head_reports_a_NEGATIVE_spearman"] is True


@pytest.mark.skipif(not OUT.is_file(), reason="bar not frozen")
def test_the_ridge_scoring_form_is_INELIGIBLE_because_it_would_end_the_zero_shot_claim():
    s = json.loads(OUT.read_text(encoding="utf-8"))["scoring"]
    assert s["ineligible_alternative"] == SCORING_RIDGE_INELIGIBLE
    assert "fits the SIGN from the training labels" in s["why_the_ridge_form_is_ineligible"]
    assert "stop being zero-shot" in s["why_the_ridge_form_is_ineligible"]


@pytest.mark.skipif(not OUT.is_file(), reason="bar not frozen")
def test_the_bar_keeps_the_scoring_half_of_the_incumbents_pipeline_and_says_which_half():
    c = json.loads(OUT.read_text(encoding="utf-8"))["scoring"]["comparability_note"]
    assert "log1p(expression)" in c and "RAW expression" in c
    assert "fits nothing" in c, "a zero-shot arm inherits the scoring half only"


# --------------------------------------------------------------------------------------------------
# the amendment record -- a bar may be corrected before first use, not after
# --------------------------------------------------------------------------------------------------
@pytest.mark.skipif(not OUT.is_file(), reason="bar not frozen")
def test_the_amendment_records_its_own_outcome_neutrality_and_a_revert_pre_commitment():
    a = json.loads(OUT.read_text(encoding="utf-8"))["amendments_before_first_score"]
    assert len(a) >= 1
    e = a[0]
    assert "TRIVIALLY SATISFIED" in e["outcome_neutrality"]
    assert "no head has been scored" in e["outcome_neutrality"]
    assert "after it, this same edit would be" in e["outcome_neutrality"].replace("\n", " ")
    assert "YES" in e["new_threshold_introduced"], (
        "a new threshold WAS added; declaring it openly is the point, since the usual rule forbids "
        "one when amending after results exist"
    )
    assert "reverted" in e["revert_pre_commitment"]


@pytest.mark.skipif(not OUT.is_file(), reason="bar not frozen")
def test_the_amended_bar_pre_commits_the_scoring_and_the_null():
    joined = " ".join(json.loads(OUT.read_text(encoding="utf-8"))["pre_commitments"])
    assert "DIRECT Spearman with the sign as given" in joined
    assert "PREDICTOR-side-valid" in joined
    assert "BOTH its median and its magnitude condition" in joined
