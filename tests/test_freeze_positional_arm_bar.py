"""Tests for the positional arm's bar-freezing script.

The point of the freeze step is that a bar frozen AFTER a run is not a bar, and a bar anchored to
RETYPED numbers is anchored to nothing. So the tests that matter most are the refusals: `verify()`
must reject a mutated comparator, a mutated arm status and a missing arm, and each refusal is
exercised against a real mutated copy of the committed artifact rather than a stub.

The second-most-important test is the reachable-vs-closed outcome EQUALITY. Containment would pass
a bar carrying an unreachable outcome -- a dead validity branch nobody would ever notice -- which is
exactly the defect the plan asked this check to catch.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts.freeze_positional_arm_bar import (
    ARM,
    CHECK_ORDER,
    CO_PRIMARIES,
    FAMILY_ALPHA,
    INCUMBENT_ARM,
    LEDGER,
    MARGIN,
    OUT,
    OUTCOMES,
    RECONCILE,
    SECONDARIES,
    BarRefused,
    bonferroni_percentiles,
    build_bar,
    classify,
    derive_parity_tolerance,
    reachable_outcomes,
    verify,
)

REPO = Path(__file__).resolve().parents[1]
FEASIBILITY = REPO / "wiki/positional_arm_feasibility_2026-11-01.json"


# --------------------------------------------------------------------------------------------------
# the closed outcome set -- EQUALITY, not containment
# --------------------------------------------------------------------------------------------------
def test_the_reachable_outcome_set_equals_the_closed_set_exactly():
    assert reachable_outcomes() == set(OUTCOMES)


def test_equality_is_the_check_not_containment():
    """A bar with an unreachable outcome is a dead validity branch; containment would pass it."""
    padded = list(CHECK_ORDER) + [("a_gate_nobody_calls", "OUTCOME_NOBODY_CAN_REACH")]
    assert reachable_outcomes(padded) > set(OUTCOMES), "strictly larger, so containment holds"
    assert reachable_outcomes(padded) != set(OUTCOMES), "but equality fails, which is the check"


def test_every_check_order_gate_maps_to_a_listed_outcome():
    for gate, result in CHECK_ORDER:
        for name in ((result,) if isinstance(result, str) else result):
            assert name in OUTCOMES, "%s returns unlisted outcome %s" % (gate, name)


def test_validity_gates_come_before_the_decision_and_split_equality_is_first():
    gates = [g for g, _ in CHECK_ORDER]
    assert gates[0] == "split_equality", "every other number is meaningless if the splits differ"
    assert gates.index("co_primary_decision") == len(gates) - 1, "validity first"
    for v in ("split_equality", "parity", "shuffled_label_null"):
        assert gates.index(v) < gates.index("co_primary_decision")


# --------------------------------------------------------------------------------------------------
# Bonferroni arithmetic
# --------------------------------------------------------------------------------------------------
def test_bonferroni_percentiles_for_two_co_primaries():
    lo, hi = bonferroni_percentiles(2)
    assert (lo, hi) == pytest.approx((1.25, 98.75))


def test_bonferroni_percentiles_for_the_full_ladder():
    lo, hi = bonferroni_percentiles(6)
    assert (lo, hi) == pytest.approx((0.41666667, 99.58333333))


def test_bonferroni_at_one_comparison_is_the_plain_95_percent_interval():
    assert bonferroni_percentiles(1) == pytest.approx((2.5, 97.5))


def test_more_comparisons_give_a_strictly_wider_interval():
    widths = [bonferroni_percentiles(n)[1] - bonferroni_percentiles(n)[0] for n in (1, 2, 4, 6, 12)]
    assert widths == sorted(widths), "correcting for more comparisons must never narrow the interval"


def test_bonferroni_percentiles_are_symmetric_about_fifty():
    for n in (1, 2, 3, 6):
        lo, hi = bonferroni_percentiles(n)
        assert lo + hi == pytest.approx(100.0)


def test_bonferroni_percentiles_refuses_degenerate_inputs():
    with pytest.raises(ValueError, match="n_comparisons"):
        bonferroni_percentiles(0)
    for bad in (0.0, 1.0, -0.1, 1.5):
        with pytest.raises(ValueError, match="family_alpha"):
            bonferroni_percentiles(2, family_alpha=bad)


# --------------------------------------------------------------------------------------------------
# the parity tolerance is DERIVED from the measurement, not preferred
# --------------------------------------------------------------------------------------------------
def test_the_parity_tolerance_derives_1e_5_from_the_achieved_deviation():
    """4.352e-07 measured -> 10x is 4.352e-06 -> smallest power of ten at or above it is 1e-05."""
    assert derive_parity_tolerance(4.352e-07) == pytest.approx(1e-5)


@pytest.mark.parametrize("measured,want", [
    (1e-9, 1e-8), (5e-9, 1e-7), (1e-7, 1e-6), (4.352e-7, 1e-5), (9e-7, 1e-5), (1e-6, 1e-5),
])
def test_the_tolerance_rule_is_monotone_and_always_above_the_measurement(measured, want):
    got = derive_parity_tolerance(measured)
    assert got == pytest.approx(want)
    assert got >= 10 * measured, "the tolerance must leave at least a 10x margin"


def test_the_tolerance_refuses_a_non_finite_or_non_positive_measurement():
    for bad in (0.0, -1e-7, float("nan"), float("inf")):
        with pytest.raises(ValueError, match="measured_deviation"):
            derive_parity_tolerance(bad)


# --------------------------------------------------------------------------------------------------
# classify -- co-primaries only, ordered, first match wins
# --------------------------------------------------------------------------------------------------
def test_classify_returns_adds_when_a_co_primary_interval_clears_the_margin():
    v, why = classify({"bins10_mean_oriented": (0.03, 0.09), "autocorr": (-0.01, 0.01)})
    assert v == "POSITIONAL_ADDS"
    assert "bins10_mean_oriented" in why


def test_classify_returns_adds_when_only_the_autocorr_co_primary_clears():
    v, why = classify({"bins10_mean_oriented": (-0.05, -0.01), "autocorr": (0.025, 0.07)})
    assert v == "POSITIONAL_ADDS"
    assert "autocorr" in why


def test_classify_returns_no_gain_only_when_both_intervals_sit_below_the_margin():
    v, _ = classify({"bins10_mean_oriented": (-0.3, -0.1), "autocorr": (-0.2, 0.005)})
    assert v == "POSITIONAL_NO_GAIN"


def test_classify_returns_ambiguous_when_an_interval_straddles_the_margin():
    v, _ = classify({"bins10_mean_oriented": (0.01, 0.05), "autocorr": (-0.1, -0.05)})
    assert v == "POSITIONAL_AMBIGUOUS"


def test_the_margin_is_a_strict_boundary():
    """An interval whose lower bound merely TOUCHES the margin does not clear it."""
    v, _ = classify({"bins10_mean_oriented": (MARGIN, 0.5), "autocorr": (-1.0, -0.5)})
    assert v == "POSITIONAL_AMBIGUOUS"
    v, _ = classify({"bins10_mean_oriented": (MARGIN + 1e-9, 0.5), "autocorr": (-1.0, -0.5)})
    assert v == "POSITIONAL_ADDS"


def test_classify_cannot_see_the_secondaries_at_all():
    """Ineligibility is enforced by construction, not by convention."""
    losing = {"bins10_mean_oriented": (-0.3, -0.1), "autocorr": (-0.3, -0.1)}
    base, _ = classify(losing)
    with_winners = dict(losing, **{s: (0.9, 0.99) for s in SECONDARIES})
    assert classify(with_winners)[0] == base == "POSITIONAL_NO_GAIN"


def test_classify_refuses_a_missing_co_primary_rather_than_grading_a_subset():
    with pytest.raises(ValueError, match="missing"):
        classify({"bins10_mean_oriented": (0.1, 0.2)})
    with pytest.raises(ValueError, match="missing"):
        classify({})


def test_every_classify_outcome_is_in_the_closed_set():
    cases = [
        {"bins10_mean_oriented": (0.1, 0.2), "autocorr": (0.1, 0.2)},
        {"bins10_mean_oriented": (-0.2, -0.1), "autocorr": (-0.2, -0.1)},
        {"bins10_mean_oriented": (-0.1, 0.1), "autocorr": (-0.1, 0.1)},
    ]
    got = {classify(c)[0] for c in cases}
    assert got == {"POSITIONAL_ADDS", "POSITIONAL_NO_GAIN", "POSITIONAL_AMBIGUOUS"}
    assert got <= set(OUTCOMES)


# --------------------------------------------------------------------------------------------------
# verify -- the refusals, against real mutated copies
# --------------------------------------------------------------------------------------------------
def _copies(tmp_path: Path):
    """Anchor copies with the arm reset to PLANNED -- the only state a freeze is legal in.

    **These tests used to read the LIVE ledger as their fixture, and that was my defect.** Once the
    arm was recorded as `measured` (Step 6 of the plan), `verify()` began refusing -- correctly,
    because a bar may only be frozen while no number exists -- and five tests went red *because the
    code works*. A unit test of `verify()`'s logic must not be coupled to a transient live state
    that the project legitimately and permanently moves past. The reset is what makes these tests
    about the FUNCTION rather than about today's ledger; the live post-measurement behaviour is
    pinned separately below.
    """
    rec = json.loads(RECONCILE.read_text(encoding="utf-8"))
    led = json.loads(LEDGER.read_text(encoding="utf-8"))
    led["criteria"][_crit(led)][ARM].update(
        status="planned", metric=None, evidence="", resolved_at="")
    rp, lp = tmp_path / "reconcile.json", tmp_path / "ledger.json"
    rp.write_text(json.dumps(rec), encoding="utf-8")
    lp.write_text(json.dumps(led), encoding="utf-8")
    return rec, led, rp, lp


def _write(path: Path, obj) -> Path:
    path.write_text(json.dumps(obj), encoding="utf-8")
    return path


def _crit(led: dict) -> str:
    return [k for k in led["criteria"] if k.startswith("ceiling:")][0]


@pytest.mark.skipif(not FEASIBILITY.is_file(), reason="feasibility artifact absent; run the probe")
def test_verify_passes_on_the_committed_anchors_in_their_pre_freeze_state(tmp_path):
    _rec, _led, rp, lp = _copies(tmp_path)
    v = verify(reconcile_path=rp, ledger_path=lp)
    assert v["reconcile"]["verdict"] == "RECONCILED"
    assert v["arms"][ARM]["status"] == "planned"


@pytest.mark.skipif(not FEASIBILITY.is_file(), reason="feasibility artifact absent; run the probe")
def test_verify_now_REFUSES_on_the_LIVE_ledger_because_the_arm_is_measured():
    """The post-run reality, pinned deliberately: a bar cannot be re-frozen after the fact.

    This is the mirror of the test above. Together they say: `verify()` accepts the anchors as they
    stood when the bar was frozen, and refuses them now that a number exists -- which is the whole
    point of the status check, not an accident of ordering.
    """
    live = json.loads(LEDGER.read_text(encoding="utf-8"))
    status = live["criteria"][_crit(live)][ARM]["status"]
    if status != "measured":
        pytest.skip("the arm is %r, not yet measured" % status)
    with pytest.raises(BarRefused, match="not 'planned'|already carries metric"):
        verify()


@pytest.mark.skipif(not FEASIBILITY.is_file(), reason="feasibility artifact absent; run the probe")
def test_verify_refuses_a_mutated_atlas_median(tmp_path):
    rec, _led, rp, lp = _copies(tmp_path)
    rec["atlas_median"] = 0.5100
    with pytest.raises(BarRefused, match="atlas median drifted"):
        verify(reconcile_path=_write(rp, rec), ledger_path=lp)


@pytest.mark.skipif(not FEASIBILITY.is_file(), reason="feasibility artifact absent; run the probe")
def test_verify_refuses_a_mutated_gc_median(tmp_path):
    rec, _led, rp, lp = _copies(tmp_path)
    rec["gc_median"] = 0.3200
    with pytest.raises(BarRefused, match="GC median drifted"):
        verify(reconcile_path=_write(rp, rec), ledger_path=lp)


@pytest.mark.skipif(not FEASIBILITY.is_file(), reason="feasibility artifact absent; run the probe")
def test_verify_refuses_an_unreconciled_or_mismatched_comparator(tmp_path):
    rec, _led, rp, lp = _copies(tmp_path)
    rec["verdict"] = "MISMATCHED"
    with pytest.raises(BarRefused, match="not RECONCILED"):
        verify(reconcile_path=_write(rp, rec), ledger_path=lp)

    rec, _led, rp, lp = _copies(tmp_path)
    rec["mismatched"] = [3]
    with pytest.raises(BarRefused, match="mismatched seeds"):
        verify(reconcile_path=_write(rp, rec), ledger_path=lp)


@pytest.mark.skipif(not FEASIBILITY.is_file(), reason="feasibility artifact absent; run the probe")
def test_verify_refuses_a_comparator_with_the_wrong_seed_or_tile_count(tmp_path):
    rec, _led, rp, lp = _copies(tmp_path)
    rec["n_seeds"] = 5
    with pytest.raises(BarRefused, match="seeds, expected"):
        verify(reconcile_path=_write(rp, rec), ledger_path=lp)

    rec, _led, rp, lp = _copies(tmp_path)
    rec["n_tiles"] = 44105
    with pytest.raises(BarRefused, match="tiles, expected"):
        verify(reconcile_path=_write(rp, rec), ledger_path=lp)


@pytest.mark.skipif(not FEASIBILITY.is_file(), reason="feasibility artifact absent; run the probe")
def test_verify_refuses_a_comparator_missing_a_per_seed_entry(tmp_path):
    rec, _led, rp, lp = _copies(tmp_path)
    rec["atlas_per_seed"].pop("7")
    with pytest.raises(BarRefused, match="atlas_per_seed has 9 entries"):
        verify(reconcile_path=_write(rp, rec), ledger_path=lp)


@pytest.mark.skipif(not FEASIBILITY.is_file(), reason="feasibility artifact absent; run the probe")
def test_verify_refuses_a_mutated_arm_status(tmp_path):
    """A bar may only be frozen while NO number exists."""
    _rec, led, rp, lp = _copies(tmp_path)
    led["criteria"][_crit(led)][ARM]["status"] = "measured"
    with pytest.raises(BarRefused, match="not 'planned'"):
        verify(reconcile_path=rp, ledger_path=_write(lp, led))


@pytest.mark.skipif(not FEASIBILITY.is_file(), reason="feasibility artifact absent; run the probe")
def test_verify_refuses_an_arm_that_already_carries_a_metric(tmp_path):
    _rec, led, rp, lp = _copies(tmp_path)
    led["criteria"][_crit(led)][ARM]["metric"] = 0.52
    with pytest.raises(BarRefused, match="already carries metric"):
        verify(reconcile_path=rp, ledger_path=_write(lp, led))


@pytest.mark.skipif(not FEASIBILITY.is_file(), reason="feasibility artifact absent; run the probe")
def test_verify_refuses_a_missing_arm(tmp_path):
    _rec, led, rp, lp = _copies(tmp_path)
    del led["criteria"][_crit(led)][ARM]
    with pytest.raises(BarRefused, match="is NOT declared"):
        verify(reconcile_path=rp, ledger_path=_write(lp, led))


@pytest.mark.skipif(not FEASIBILITY.is_file(), reason="feasibility artifact absent; run the probe")
def test_verify_refuses_a_drifted_incumbent_or_anchor_arm(tmp_path):
    for arm, pat in ((INCUMBENT_ARM, "in the ledger, expected"),
                     ("conv-3width-13k-params-150bp", "conv anchor reads"),
                     ("gc-single-feature", "gc anchor reads")):
        _rec, led, rp, lp = _copies(tmp_path)
        led["criteria"][_crit(led)][arm]["metric"] = 0.99
        with pytest.raises(BarRefused, match=pat):
            verify(reconcile_path=rp, ledger_path=_write(lp, led))


@pytest.mark.skipif(not FEASIBILITY.is_file(), reason="feasibility artifact absent; run the probe")
def test_verify_refuses_a_failed_feasibility_probe(tmp_path):
    _rec, _led, rp, lp = _copies(tmp_path)
    feas = json.loads(FEASIBILITY.read_text(encoding="utf-8"))
    feas["parity"]["passed"] = False
    fp = _write(tmp_path / "feas.json", feas)
    with pytest.raises(BarRefused, match="feasibility verdict"):
        verify(reconcile_path=rp, ledger_path=lp, feasibility_path=fp)


def test_verify_refuses_an_absent_anchor_source(tmp_path):
    with pytest.raises(BarRefused, match="anchor source absent"):
        verify(reconcile_path=tmp_path / "nope.json")


@pytest.mark.skipif(not FEASIBILITY.is_file(), reason="feasibility artifact absent; run the probe")
def test_verify_refuses_a_ledger_without_exactly_one_ceiling_criterion(tmp_path):
    _rec, led, rp, lp = _copies(tmp_path)
    led["criteria"]["ceiling: a second one that should not exist"] = {}
    with pytest.raises(BarRefused, match="expected exactly one"):
        verify(reconcile_path=rp, ledger_path=_write(lp, led))


# --------------------------------------------------------------------------------------------------
# the written bar
# --------------------------------------------------------------------------------------------------
@pytest.mark.skipif(not OUT.is_file(), reason="bar not frozen yet")
def test_the_frozen_bars_margin_equals_the_ledger_derived_inherited_value():
    """An inherited number must not be able to become a re-chosen one."""
    bar = json.loads(OUT.read_text(encoding="utf-8"))
    assert bar["inherited_not_chosen"]["margin"]["value"] == MARGIN == 0.02
    assert "INHERITED" in bar["inherited_not_chosen"]["margin"]["derivation"]
    # and the margin must be the same number the outcome definitions are written against
    for key in ("POSITIONAL_ADDS", "POSITIONAL_NO_GAIN", "POSITIONAL_AMBIGUOUS"):
        assert "0.02" in bar["outcomes_defined_before_running"][key]


@pytest.mark.skipif(not OUT.is_file(), reason="bar not frozen yet")
def test_the_frozen_bar_pins_the_committed_comparator_values_verbatim():
    bar = json.loads(OUT.read_text(encoding="utf-8"))
    rec = json.loads(RECONCILE.read_text(encoding="utf-8"))
    comp = bar["inherited_not_chosen"]["comparator"]
    assert comp["atlas_per_seed"] == rec["atlas_per_seed"]
    assert comp["gc_per_seed"] == rec["gc_per_seed"]
    assert comp["atlas_median"] == 0.5076
    assert comp["must_be_recomputed_in_run"] is False
    assert bar["inherited_not_chosen"]["split"]["n_test_per_seed"] == rec["n_test_per_seed"]


@pytest.mark.skipif(not OUT.is_file(), reason="bar not frozen yet")
def test_the_frozen_bar_declares_the_closed_set_and_their_equality():
    bar = json.loads(OUT.read_text(encoding="utf-8"))
    assert sorted(bar["closed_outcome_set"]) == sorted(OUTCOMES)
    assert bar["reachable_equals_closed"] is True
    assert sorted(bar["reachable_outcome_set"]) == sorted(OUTCOMES)
    assert bar["frozen_before_any_score_exists"] is True


@pytest.mark.skipif(not OUT.is_file(), reason="bar not frozen yet")
def test_the_frozen_bar_records_the_multiplicity_correction_and_the_divisor_discrepancy():
    m = json.loads(OUT.read_text(encoding="utf-8"))["multiplicity"]
    assert m["co_primaries"] == list(CO_PRIMARIES)
    assert m["secondaries"] == list(SECONDARIES)
    assert m["co_primary_alpha_each"] == pytest.approx(FAMILY_ALPHA / 2)
    assert m["secondary_alpha_each"] == pytest.approx(FAMILY_ALPHA / 6)
    assert m["secondaries_are_ineligible_as_the_arms_result"] is True
    # the discrepancy between the plan's prose and the count is RECORDED, not silently applied
    assert "5-comparison" in m["divisor_discrepancy_recorded"]
    assert "outcome-neutral" in m["divisor_discrepancy_recorded"]


@pytest.mark.skipif(not OUT.is_file(), reason="bar not frozen yet")
def test_the_frozen_bar_carries_both_scoping_corrections_and_the_null_more_seeds_rule():
    bar = json.loads(OUT.read_text(encoding="utf-8"))
    assert len(bar["scoping_corrections_fixed_in_advance"]) == 2
    joined = " ".join(bar["scoping_corrections_fixed_in_advance"])
    assert "-35/-10" in joined and "CANNOT" in joined
    assert "MORE SEEDS" in bar["validity_controls"]["shuffled_label_null"]["on_failure"]
    assert any("MORE SEEDS" in p for p in bar["pre_commitments"])
    assert any("0.5076 is NEVER retro-edited" in p for p in bar["pre_commitments"])


@pytest.mark.skipif(not OUT.is_file(), reason="bar not frozen yet")
def test_the_frozen_bar_names_the_no_gain_over_read_and_the_1_bit_orientation_inference():
    limits = " ".join(json.loads(OUT.read_text(encoding="utf-8"))["honest_limits_fixed_in_advance"])
    assert "NOT that positional structure is absent" in limits
    assert "1-bit inference" in limits
    assert "689 duplicate sequences" in limits
    assert "0.941" in limits


@pytest.mark.skipif(not OUT.is_file(), reason="bar not frozen yet")
def test_the_frozen_bar_records_a_forecast_with_its_record_and_the_do_not_defer_note():
    f = json.loads(OUT.read_text(encoding="utf-8"))["forecast"]
    assert f["recorded_before_any_number_exists"] is True
    assert f["predicted_outcome"] in OUTCOMES
    ps = [f["p_positional_adds"], f["p_positional_no_gain"], f["p_positional_ambiguous"]]
    assert sum(ps) == pytest.approx(1.0), "the three verdict probabilities must be a distribution"
    assert "3 of 7" in f["my_forecast_record_in_this_project"]
    assert "DO NOT defer" in f["my_forecast_record_in_this_project"]


@pytest.mark.skipif(not FEASIBILITY.is_file(), reason="feasibility artifact absent; run the probe")
def test_build_bar_is_a_pure_function_of_the_verified_anchors(tmp_path):
    """Re-deriving the bar from the same anchors must reproduce it except for the freeze stamp."""
    _rec, _led, rp, lp = _copies(tmp_path)
    v = verify(reconcile_path=rp, ledger_path=lp)
    a, b = build_bar(v), build_bar(v)
    a.pop("frozen"), b.pop("frozen")
    assert a == b


@pytest.mark.skipif(not OUT.is_file(), reason="bar not frozen yet")
def test_the_bar_on_disk_still_matches_a_rebuild_from_its_pre_freeze_anchors(tmp_path):
    """The frozen bar must not have drifted from the anchors it was built on.

    Compared with `frozen` and the probe-derived throughput/free-RAM stripped: those are timing and
    host measurements that legitimately differ between the freeze run and a rebuild.
    """
    _rec, _led, rp, lp = _copies(tmp_path)
    rebuilt = build_bar(verify(reconcile_path=rp, ledger_path=lp))
    on_disk = json.loads(OUT.read_text(encoding="utf-8"))
    volatile = ("throughput_seq_per_second", "free_ram_bytes_at_freeze")
    for d in (rebuilt, on_disk):
        d.pop("frozen", None)
        for k in volatile:
            d["measured_and_inherited_from_the_probe"].pop(k, None)
    assert rebuilt == on_disk
