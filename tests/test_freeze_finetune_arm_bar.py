"""Tests for the fine-tune arm's frozen acceptance bar.

The bar's whole value is that it cannot be bent after a number exists, so most of these tests are
REFUSALS: each one mutates a single anchor and asserts the freeze declines. A bar that builds
against a drifted comparator or an already-scored arm is not a bar.

The ceiling tests are the load-bearing ones. The plan replaces a fixed 6-hour constant with a
window-derived ceiling, which is exactly the move that looks like motivated reasoning -- so the
rule is driven in BOTH directions and must route this arm to Kaggle on a short window just as the
old constant did.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts.freeze_positional_arm_bar import MARGIN, BarRefused
from scripts.freeze_expression_heads_bar import NULL_SIGNAL_RATIO
import scripts.freeze_finetune_arm_bar as fb

REPO = Path(__file__).resolve().parents[1]


# --------------------------------------------------------------------- anchor fixtures


@pytest.fixture
def anchors(tmp_path):
    """Real anchors copied to tmp, so a test can mutate one without touching the repo.

    The arm is PINNED back to `planned` in the temp copy. 2026-11-17: the fine-tuned arm was
    recorded MEASURED, and the freeze script correctly refuses to build a bar over a scored arm --
    `a bar frozen after a run is not a bar`, which is the whole anti-bar-shopping guard. That
    refusal is right and is NOT relaxed anywhere: `test_refuses_an_arm_that_is_no_longer_planned`
    still sets `measured` explicitly and still asserts the refusal, and it overrides this pin.

    What was wrong was the FIXTURE, not the guard. These tests ask what the freeze does given a
    planned arm, so `planned` is their PRECONDITION and must be established here rather than
    inherited from live state -- inheriting it meant 44 tests passed only until the arm was run,
    i.e. only until the project reached the end state the whole arc was aiming at.
    """
    out = {}
    for key, src in (("reconcile", fb.RECONCILE), ("ledger", fb.LEDGER), ("probe", fb.PROBE)):
        dst = tmp_path / src.name
        dst.write_text(src.read_text(encoding="utf-8"), encoding="utf-8")
        out[key] = dst

    led = json.loads(out["ledger"].read_text(encoding="utf-8"))
    crit = next(k for k in led["criteria"] if fb.ARM in led["criteria"][k])
    arm = led["criteria"][crit][fb.ARM]
    arm.update(status="planned", metric=None, evidence="")
    out["ledger"].write_text(json.dumps(led, indent=2), encoding="utf-8")
    return out


def _edit(path: Path, fn):
    d = json.loads(path.read_text(encoding="utf-8"))
    fn(d)
    path.write_text(json.dumps(d), encoding="utf-8")


def _verify(anchors):
    return fb.verify(reconcile_path=anchors["reconcile"], ledger_path=anchors["ledger"],
                     probe_path=anchors["probe"])


def _arm(ledger: dict) -> dict:
    crit = [k for k in ledger["criteria"] if k.startswith("ceiling:")][0]
    return ledger["criteria"][crit]


# --------------------------------------------------------------------- the unmutated path


def test_the_real_anchors_verify(anchors):
    """The baseline. Without this the refusal tests could all pass on a broken verify()."""
    v = _verify(anchors)
    assert v["frozen_epochs"] > 0
    assert v["probe"]["schema"] == fb.PROBE_SCHEMA


def test_the_bar_builds_from_the_real_anchors(anchors):
    bar = fb.build_bar(_verify(anchors))
    assert bar["schema"] == "finetune-arm-acceptance-bar-v1"
    assert bar["arm"] == fb.ARM


def test_the_bar_is_a_pure_function_of_its_anchors(anchors):
    v = _verify(anchors)
    a = json.dumps(fb.build_bar(v), sort_keys=True)
    b = json.dumps(fb.build_bar(v), sort_keys=True)
    assert a == b


# --------------------------------------------------------------------- comparator refusals


def test_refuses_an_unreconciled_comparator(anchors):
    _edit(anchors["reconcile"], lambda d: d.update(verdict="MISMATCH"))
    with pytest.raises(BarRefused, match="not RECONCILED"):
        _verify(anchors)


def test_refuses_a_drifted_comparator_median(anchors):
    """0.5076 is the number the whole arm is graded against. If it moves, the bar is void."""
    _edit(anchors["reconcile"], lambda d: d.update(atlas_median=0.52))
    with pytest.raises(BarRefused, match="drifted"):
        _verify(anchors)


def test_refuses_a_comparator_with_mismatched_seeds(anchors):
    _edit(anchors["reconcile"], lambda d: d.update(mismatched=["seed3"]))
    with pytest.raises(BarRefused, match="mismatches"):
        _verify(anchors)


def test_refuses_a_missing_comparator_file(anchors):
    anchors["reconcile"].unlink()
    with pytest.raises(BarRefused, match="anchor source absent"):
        _verify(anchors)


# --------------------------------------------------------------------- arm refusals


def test_refuses_an_undeclared_arm(anchors):
    _edit(anchors["ledger"], lambda d: _arm(d).pop(fb.ARM))
    with pytest.raises(BarRefused, match="NOT declared"):
        _verify(anchors)


def test_refuses_an_arm_that_already_carries_a_metric(anchors):
    """The decisive refusal. Freezing a bar around a number that already exists is choosing the
    bar that the number passes."""
    _edit(anchors["ledger"], lambda d: _arm(d)[fb.ARM].update(metric=0.61))
    with pytest.raises(BarRefused, match="already carries metric"):
        _verify(anchors)


def test_refuses_an_arm_that_is_no_longer_planned(anchors):
    _edit(anchors["ledger"], lambda d: _arm(d)[fb.ARM].update(status="measured"))
    with pytest.raises(BarRefused, match="not 'planned'"):
        _verify(anchors)


def test_refuses_a_missing_incumbent(anchors):
    """With no incumbent in the ledger there is nothing to compare against."""
    _edit(anchors["ledger"], lambda d: _arm(d).pop(fb.INCUMBENT_ARM))
    with pytest.raises(BarRefused, match="incumbent"):
        _verify(anchors)


# --------------------------------------------------------------------- probe refusals


def test_refuses_a_probe_with_the_wrong_schema(anchors):
    _edit(anchors["probe"], lambda d: d.update(schema="something-else-v1"))
    with pytest.raises(BarRefused, match="schema"):
        _verify(anchors)


def test_refuses_an_INDETERMINATE_probe(anchors):
    _edit(anchors["probe"], lambda d: d.update(verdict="INDETERMINATE_NO_CURVE"))
    with pytest.raises(BarRefused, match="did not establish feasibility"):
        _verify(anchors)


def test_refuses_a_probe_whose_curve_never_turned(anchors):
    """A monotonically rising curve means the budget was binding, so a null from this arm could
    not be read at all -- which is the premise the whole budget-versus-representation
    pre-commitment rests on."""
    _edit(anchors["probe"], lambda d: d["convergence"]["early_stop"].update(turned=False))
    with pytest.raises(BarRefused, match="did NOT turn"):
        _verify(anchors)


@pytest.mark.parametrize("key", ["learning_rate", "batch_size", "patience"])
def test_refuses_a_probe_missing_a_hyperparameter_to_freeze(anchors, key):
    _edit(anchors["probe"], lambda d: d["hyperparameters_used"].update({key: None}))
    with pytest.raises(BarRefused, match=key):
        _verify(anchors)


def test_refuses_a_probe_trained_on_the_wrong_transform(anchors):
    """ridge_spearman fits on log1p and scores against raw; deviating makes the number
    incomparable with the committed incumbent tables, by that function's own docstring."""
    _edit(anchors["probe"], lambda d: d["target_transform"].update(trained_on="raw expression"))
    with pytest.raises(BarRefused, match="log1p"):
        _verify(anchors)


def test_refuses_a_probe_with_no_deviation_list(anchors):
    _edit(anchors["probe"], lambda d: d.pop("deviations_from_pretraining_recipe"))
    with pytest.raises(BarRefused, match="deviation"):
        _verify(anchors)


def test_refuses_a_probe_with_no_epoch_budget(anchors):
    _edit(anchors["probe"], lambda d: d["local_projection"].pop("n_epochs"))
    with pytest.raises(BarRefused, match="epoch count"):
        _verify(anchors)


def test_refuses_a_probe_with_no_local_projection_to_route_against(anchors):
    v = _verify(anchors)
    probe = dict(v["probe"])
    probe.pop("local_projection")
    with pytest.raises(BarRefused, match="no local projection"):
        fb.choose_route(probe, fb.derive_ceiling(216.0))


# --------------------------------------------------------------------- the derived ceiling


def test_a_ten_day_window_routes_the_arm_LOCAL(anchors):
    """The declared window's direction."""
    bar = fb.build_bar(_verify(anchors), window=fb.WINDOW)
    assert bar["route"]["chosen"] == fb.ROUTE_LOCAL
    assert bar["route"]["fits"] is True


def test_a_SHORT_window_routes_the_arm_to_KAGGLE(anchors):
    """The other direction, and the one that proves the ceiling is not a disguised constant.

    Tonight's ~11 h window yields a ~6.6 h ceiling, which rejects this arm's 23.82 h local
    projection exactly as the fixed 6-hour constant did. If this ever passes LOCAL, the
    derivation has stopped being a derivation.
    """
    short = dict(fb.WINDOW, start="2026-11-03T20:00:00", end="2026-11-04T07:00:00")
    bar = fb.build_bar(_verify(anchors), window=short)
    assert bar["route"]["chosen"] == fb.ROUTE_KAGGLE
    assert bar["route"]["fits"] is False


def test_the_KAGGLE_route_declares_the_horizon_inapplicable_rather_than_refusing(anchors):
    """A short window must still FREEZE a bar -- with a Kaggle route -- not raise.

    The first version of `build_bar` sized a local horizon unconditionally and raised here,
    conflating "this window cannot host the LOCAL run" (true, and the reason the route is Kaggle)
    with "the bar cannot be frozen" (false).
    """
    short = dict(fb.WINDOW, start="2026-11-03T20:00:00", end="2026-11-04T07:00:00")
    h = fb.build_bar(_verify(anchors), window=short)["horizon"]
    assert h["applicable"] is False
    assert h["diagnostic_horizon_epochs"] is None
    assert h["frozen_budget_epochs"] == _verify(anchors)["frozen_epochs"]


def test_the_KAGGLE_route_records_that_losing_the_diagnostics_WEAKENS_the_arm(anchors):
    """The honesty rail on the rejected route. Without the long horizon the budget-limited
    reading stays closed by inference, which is strictly weaker -- and must be visible in the
    bar rather than discovered at reporting time."""
    short = dict(fb.WINDOW, start="2026-11-03T20:00:00", end="2026-11-04T07:00:00")
    bar = fb.build_bar(_verify(anchors), window=short)
    cost = bar["horizon"]["what_this_costs"]
    assert "UNAVAILABLE" in cost
    assert "INFERENCE" in cost
    assert "weakening" in cost
    assert "UNAVAILABLE" in \
        bar["budget_versus_representation"]["ineligible_diagnostics"]["what"]


def test_the_ceiling_is_monotone_in_the_window():
    """More window must never buy a smaller ceiling."""
    prev = 0.0
    for hours in (11.0, 24.0, 72.0, 216.0, 480.0):
        c = fb.derive_ceiling(hours)["ceiling_hours"]
        assert c > prev
        prev = c


def test_the_ceiling_is_strictly_below_the_window():
    """Committing the whole window leaves nothing for the measured 34% throughput swing or for a
    reboot, which is the failure mode the fraction exists to prevent."""
    for hours in (11.0, 216.0):
        c = fb.derive_ceiling(hours)
        assert c["ceiling_hours"] < c["hours_available"]


def test_the_commit_fraction_survives_the_measured_throughput_swing():
    """The derivation's own arithmetic, checked rather than asserted in prose: a run sized at the
    ceiling must still fit the window if the host runs at the slower of the two measured rates."""
    hours = fb.window_hours(fb.WINDOW)
    c = fb.derive_ceiling(hours)
    worst_case = c["ceiling_hours"] * (1 + fb.MEASURED_THROUGHPUT_SWING)
    assert worst_case <= hours


def test_a_nonpositive_window_is_refused():
    with pytest.raises(BarRefused):
        fb.derive_ceiling(0.0)
    with pytest.raises(BarRefused):
        fb.derive_ceiling(-5.0)


def test_an_inverted_window_is_refused():
    with pytest.raises(BarRefused, match="ends at or before"):
        fb.window_hours({"name": "bad", "start": "2026-11-13T00:00:00",
                         "end": "2026-11-04T00:00:00"})


@pytest.mark.parametrize("frac", [0.0, -0.1, 1.5])
def test_an_out_of_range_commit_fraction_is_refused(frac):
    with pytest.raises(BarRefused, match="commit fraction"):
        fb.derive_ceiling(216.0, fraction=frac)


def test_window_hours_is_computed_not_restated():
    """2026-11-04 to 2026-11-13 is nine days. The plan's round 'about 240 hours' is not the
    number this bar uses, and a bar that restated it would be unverifiable."""
    assert fb.window_hours(fb.WINDOW) == 216.0


def test_the_ceiling_records_what_it_replaces_and_why(anchors):
    bar = fb.build_bar(_verify(anchors))
    rep = bar["ceiling"]["replaces"]
    assert "LOCAL_WALL_CEILING_SECONDS" in rep["constant"]
    assert "operational" in rep["its_own_stated_justification"].lower()
    assert "attendance" in rep["why_it_no_longer_binds"].lower()


def test_the_route_records_that_only_the_ceiling_moved(anchors):
    """The honesty rail on the route inversion: the probe chose Kaggle, this bar chooses local,
    and the reader must be able to see that no MEASUREMENT changed between them."""
    bar = fb.build_bar(_verify(anchors))
    assert bar["route"]["probe_had_chosen"] == fb.ROUTE_KAGGLE
    assert bar["route"]["chosen"] == fb.ROUTE_LOCAL
    why = bar["route"]["why_the_probe_chose_differently"]
    assert "Nothing measured changed" in why
    assert "%.2f" % bar["route"]["local_projection_hours"] in why


# --------------------------------------------------------------------- horizon sizing


def test_the_horizon_is_derived_from_the_ceiling_and_fits_it(anchors):
    bar = fb.build_bar(_verify(anchors))
    h = bar["horizon"]
    assert h["projected_hours_at_horizon"] <= bar["ceiling"]["ceiling_hours"]
    assert h["diagnostic_horizon_epochs"] > h["frozen_budget_epochs"]


def test_a_smaller_window_buys_a_shorter_horizon(anchors):
    """The horizon must follow the ceiling, or an overrun is possible on a slower host."""
    v = _verify(anchors)
    big = fb.size_horizon(fb.derive_ceiling(216.0), v["probe"], v["frozen_epochs"])
    small = fb.size_horizon(fb.derive_ceiling(80.0), v["probe"], v["frozen_epochs"])
    assert small["diagnostic_horizon_epochs"] < big["diagnostic_horizon_epochs"]


def test_a_window_too_small_for_the_frozen_BUDGET_is_refused(anchors):
    """If the arm's own verdict does not fit the window, the bar must not be frozen at all --
    rather than silently freezing a budget the run cannot reach."""
    v = _verify(anchors)
    with pytest.raises(BarRefused, match="fewer than"):
        fb.size_horizon(fb.derive_ceiling(20.0), v["probe"], v["frozen_epochs"])


def test_the_horizon_refuses_a_probe_with_no_per_epoch_cost(anchors):
    v = _verify(anchors)
    probe = json.loads(json.dumps(v["probe"]))
    probe["local_projection"].pop("seconds_per_epoch")
    with pytest.raises(BarRefused, match="per-epoch cost"):
        fb.size_horizon(fb.derive_ceiling(216.0), probe, v["frozen_epochs"])


# --------------------------------------------------------------------- outcome closure


def test_reachable_equals_declared():
    assert fb.reachable_outcomes() == set(fb.OUTCOMES)


def test_reachability_is_non_vacuous_under_PADDING():
    """Adding a declared outcome no gate can produce must break the equality.

    Containment would not catch this, which is why the check is equality: a bar promising a
    verdict it cannot deliver is worse than one promising fewer.
    """
    padded = set(fb.OUTCOMES) | {"FINETUNE_SPECTACULAR"}
    assert fb.reachable_outcomes() != padded


def test_reachability_is_non_vacuous_under_TRIMMING():
    """Removing a gate must break the equality, proving the check reads the order it claims to."""
    trimmed = tuple(g for g in fb.CHECK_ORDER if g[0] != "shuffled_label_null")
    assert fb.reachable_outcomes(trimmed) != set(fb.OUTCOMES)
    assert "INDETERMINATE_NULL_FAIL" not in fb.reachable_outcomes(trimmed)


def test_build_bar_refuses_if_the_outcome_set_drifts(anchors, monkeypatch):
    """The equality is ENFORCED at freeze time, not merely tested here."""
    monkeypatch.setattr(fb, "OUTCOMES", fb.OUTCOMES + ("FINETUNE_UNREACHABLE",))
    with pytest.raises(BarRefused, match="reachable outcomes"):
        fb.build_bar(_verify(anchors))


def test_liveness_is_the_FIRST_gate(anchors):
    """The ordering inversion against the heads arm, pinned. A dead run has no numbers, so asking
    about split equality first would classify a silent death as a data problem."""
    assert fb.CHECK_ORDER[0][0] == "run_completed"
    assert fb.CHECK_ORDER[0][1] == "INDETERMINATE_RUN_DIED"
    bar = fb.build_bar(_verify(anchors))
    assert bar["check_order"][0] == ["run_completed", "INDETERMINATE_RUN_DIED"]
    assert "DEAD_NO_SENTINEL" in bar["why_liveness_is_checked_first"]


def test_the_budget_exhausted_outcome_is_distinct_from_a_verdict(anchors):
    bar = fb.build_bar(_verify(anchors))
    assert "INDETERMINATE_BUDGET_EXHAUSTED" in bar["outcomes"]
    assert "Neither is a verdict about the arm" in bar["why_liveness_is_checked_first"]


# --------------------------------------------------------------------- the decision band


def test_classify_is_total_on_the_interval():
    for lo, hi in [(-1.0, -0.5), (-0.01, 0.01), (0.0, 0.03), (0.03, 0.09), (0.021, 0.5)]:
        o, _ = fb.classify(lo, hi)
        assert o in fb.OUTCOMES


def test_a_gain_smaller_than_the_margin_is_a_TIE():
    """Decided in advance, so it cannot be re-read as a win later."""
    o, _ = fb.classify(0.005, 0.019)
    assert o == "FINETUNE_TIES"


def test_a_gain_clearing_the_margin_ADDS():
    o, _ = fb.classify(MARGIN + 1e-6, 0.2)
    assert o == "FINETUNE_ADDS"


def test_an_interval_entirely_below_zero_is_NO_GAIN():
    o, why = fb.classify(-0.3, -0.05)
    assert o == "FINETUNE_NO_GAIN"
    assert "LOSS" in why


def test_the_margin_is_exactly_at_the_boundary_a_TIE():
    """Inclusive at the margin: `ci_lo > margin` is strict, so landing exactly on it ties."""
    assert fb.classify(MARGIN, 0.2)[0] == "FINETUNE_TIES"


# --------------------------------------------------------------------- budget vs representation


def test_the_bar_names_the_probe_curve_FIELD_not_just_the_claim(anchors):
    """A pre-commitment that points at a field a reader can open is checkable; one that
    paraphrases a curve is not."""
    ev = fb.build_bar(_verify(anchors))["budget_versus_representation"]["evidence_today"]
    assert ev["curve_field"] == "convergence.val_spearman_by_epoch"
    assert ev["probe_field"] == "convergence.early_stop"
    probe = json.loads(fb.PROBE.read_text(encoding="utf-8"))
    cur = probe
    for part in ev["curve_field"].split("."):
        cur = cur[part]
    assert isinstance(cur, list) and len(cur) > 1, "the named field must actually resolve"


def test_the_bar_records_the_subsample_scope_limit(anchors):
    ev = fb.build_bar(_verify(anchors))["budget_versus_representation"]["evidence_today"]
    assert "SUBSAMPLE" in ev["why_this_is_not_enough"]
    assert "2,160" in ev["why_this_is_not_enough"] or "2160" in ev["why_this_is_not_enough"]


def test_the_longer_horizons_are_pre_declared_INELIGIBLE_by_name(anchors):
    d = fb.build_bar(_verify(anchors))["budget_versus_representation"]["ineligible_diagnostics"]
    assert d["status"].startswith("INELIGIBLE")
    assert "max-over-windows" in d["why_pre_declared"]
    assert "pre-declared" in d["status"] or "pre-declared" in d["why_pre_declared"]


def test_the_two_readings_are_both_named(anchors):
    bar = fb.build_bar(_verify(anchors))
    text = bar["budget_versus_representation"]["the_distinction"]
    assert "REPRESENTATION" in text
    assert "BUDGET" in text


# --------------------------------------------------------------------- the null


def test_the_null_is_declared_NON_vacuous_here(anchors):
    n = fb.build_bar(_verify(anchors))["null"]
    assert n["kind"] == "shuffled_training_labels"
    assert n["is_it_vacuous_here"] is False


def test_the_bar_says_WHY_the_two_arms_nulls_DIFFER(anchors):
    """The plan's requirement, and the reason it matters: an uncorrected difference between two
    arms' nulls reads as one arm having been held to a weaker standard."""
    n = fb.build_bar(_verify(anchors))["null"]
    why = n["why_the_heads_arm_needed_a_DIFFERENT_null"]
    assert "SIGN-FLIP" in why
    assert "0.8814" in why, "the measured vacuity must be quoted, not alluded to"
    assert "0.0371" in why, "the replacement null's measured magnitude must be quoted too"
    assert "weaker standard" in n["why_this_is_recorded_here"]


def test_the_null_ratio_is_IMPORTED_not_restated(anchors):
    """A restated copy drifts. Importing means the two arms cannot silently diverge."""
    nr = fb.build_bar(_verify(anchors))["inherited_not_chosen"]["null_signal_ratio"]
    assert nr["value"] == NULL_SIGNAL_RATIO
    assert "freeze_expression_heads_bar" in nr["imported_from"]


def test_the_tighter_sampling_figure_is_a_DIAGNOSTIC_not_the_gate(anchors):
    """Re-deriving a gate against the number it is about to judge is the logged failure."""
    n = fb.build_bar(_verify(anchors))["null"]
    assert "not the gate" in n["sampling_diagnostic_not_a_gate"]


# --------------------------------------------------------------------- inherited anchors


def test_the_margin_is_inherited_verbatim(anchors):
    m = fb.build_bar(_verify(anchors))["inherited_not_chosen"]["margin"]
    assert m["value"] == MARGIN == 0.02
    assert "INHERITED" in m["derivation"]


def test_the_comparator_median_is_read_from_the_anchor_not_retyped(anchors):
    bar = fb.build_bar(_verify(anchors))
    rec = json.loads(fb.RECONCILE.read_text(encoding="utf-8"))
    assert bar["inherited_not_chosen"]["comparator"]["median_spearman"] == rec["atlas_median"]


def test_the_target_transform_is_inherited_and_says_why_it_matters(anchors):
    tt = fb.build_bar(_verify(anchors))["inherited_not_chosen"]["target_transform"]
    assert tt["trained_on"] == "log1p(expression)"
    assert "raw" in str(tt["scored_against"]).lower()
    assert "ridge_spearman" in tt["inherited_from"]
    assert "confound" in tt["why_it_is_not_a_free_parameter"]


def test_the_frozen_hyperparameters_come_from_the_probe(anchors):
    bar = fb.build_bar(_verify(anchors))
    probe = json.loads(fb.PROBE.read_text(encoding="utf-8"))
    hp = bar["frozen_hyperparameters"]
    for k in ("learning_rate", "batch_size", "patience"):
        assert hp[k] == probe["hyperparameters_used"][k]
    assert hp["n_epochs_budget"] == probe["local_projection"]["n_epochs"]


def test_the_frozen_budget_is_a_MAXIMUM_not_a_requirement(anchors):
    """Early stopping may end a seed sooner, and that is not a budget failure -- otherwise every
    well-converged seed would classify INDETERMINATE_BUDGET_EXHAUSTED."""
    es = fb.build_bar(_verify(anchors))["frozen_hyperparameters"]["early_stopping"]
    assert "MAXIMUM" in es["rule"]


# --------------------------------------------------------------------- multiplicity


def test_no_multiplicity_correction_and_the_bar_explains_the_contrast(anchors):
    m = fb.build_bar(_verify(anchors))["multiplicity"]
    assert m["correction"] == "none"
    assert m["n_declared_comparisons"] == 1
    assert "alpha/5" in m["why"], "the heads arm's divisor must be named for the contrast to read"
    assert "smaller hypothesis space" in m["contrast_is_deliberate"]


# --------------------------------------------------------------------- scope honesty


def test_the_bar_bounds_what_a_null_would_mean(anchors):
    bar = fb.build_bar(_verify(anchors))
    text = bar["what_a_null_here_bounds"]
    assert "THIS rung" in text
    assert "does not bound fine-tuning in general" in text


def test_the_window_records_that_attendance_is_the_binding_constraint(anchors):
    w = fb.build_bar(_verify(anchors))["window"]
    assert "NONE" in w["attendance"]
    assert "Kaggle" in w["attendance"], "the reason Kaggle is unavailable must be the attendance"
