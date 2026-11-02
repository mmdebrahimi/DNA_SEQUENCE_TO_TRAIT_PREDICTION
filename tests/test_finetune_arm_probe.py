"""Tests for the fine-tune feasibility and compute-route probe.

Three rails matter more than the arithmetic here, and each exists because the alternative is a
number nobody can interpret:

* **The route rule is exercised in BOTH directions.** A selector that can only ever return one
  answer is not a rule. This host projects over the ceiling and therefore routes to Kaggle, so the
  local branch would never be reached by a real run and is driven synthetically.

* **An unrecorded hyperparameter deviation RAISES.** The fine-tune is anchored to the checkpoint's
  own pretraining recipe; a silent substitution would make the arm a weaker test of the same
  question while looking identical in the artifact.

* **A curve that never turned must say so.** If the validation metric was still improving at the
  last epoch, the epoch budget it suggests is a LOWER BOUND, and a null under that budget cannot be
  read as representation-limited. That distinction is the whole reason this probe gates the arm.
"""

from __future__ import annotations

import json
import pathlib

import pytest

from scripts import finetune_arm_probe as F
from scripts.expression_heads_probe import ProbeRefused

ARTIFACT = pathlib.Path(__file__).resolve().parents[1] / "wiki/finetune_arm_feasibility_2026-11-02.json"


# --------------------------------------------------------------------------------------------------
# cost arithmetic -- known answers
# --------------------------------------------------------------------------------------------------
def test_seconds_per_epoch_is_sequences_over_rate():
    assert F.seconds_per_epoch(1000, 50.0) == pytest.approx(20.0)


@pytest.mark.parametrize("n,rate", [(0, 50.0), (-1, 50.0)])
def test_seconds_per_epoch_refuses_a_nonsense_cohort(n, rate):
    with pytest.raises(ValueError, match="n_train must be"):
        F.seconds_per_epoch(n, rate)


@pytest.mark.parametrize("rate", [0.0, -3.0])
def test_seconds_per_epoch_refuses_a_nonsense_rate(rate):
    with pytest.raises(ValueError, match="throughput must be positive"):
        F.seconds_per_epoch(1000, rate)


def test_project_arm_cost_multiplies_epochs_and_seeds():
    """1000 sequences at 50/s = 20 s/epoch; 4 epochs x 3 seeds = 240 s."""
    p = F.project_arm_cost(n_train=1000, train_seq_per_second=50.0, n_epochs=4, n_seeds=3)
    assert p["seconds_per_epoch"] == pytest.approx(20.0)
    assert p["seconds_per_seed"] == pytest.approx(80.0)
    assert p["seconds_total"] == pytest.approx(240.0)
    assert p["hours_total"] == round(240.0 / 3600.0, 2)


@pytest.mark.parametrize("epochs,seeds", [(0, 10), (5, 0), (-1, 10)])
def test_project_arm_cost_refuses_nonsense_budgets(epochs, seeds):
    with pytest.raises(ValueError, match="must both be >= 1"):
        F.project_arm_cost(n_train=100, train_seq_per_second=10.0,
                           n_epochs=epochs, n_seeds=seeds)


# --------------------------------------------------------------------------------------------------
# the route rule -- BOTH directions, so it cannot be vacuous
# --------------------------------------------------------------------------------------------------
def test_a_projection_UNDER_the_ceiling_routes_LOCAL():
    route, why = F.choose_route(3600.0, ceiling_seconds=F.LOCAL_WALL_CEILING_SECONDS)
    assert route == F.ROUTE_LOCAL
    assert "within the named ceiling" in why


def test_a_projection_OVER_the_ceiling_routes_KAGGLE():
    route, why = F.choose_route(10 * 3600.0, ceiling_seconds=F.LOCAL_WALL_CEILING_SECONDS)
    assert route == F.ROUTE_KAGGLE
    assert "EXCEEDS the named ceiling" in why
    assert "abandoned part-way" in why, "the reason for the ceiling must travel with the decision"


def test_the_route_rule_is_not_vacuous_the_two_branches_are_reachable_and_distinct():
    lo, _ = F.choose_route(1.0)
    hi, _ = F.choose_route(1e9)
    assert {lo, hi} == {F.ROUTE_LOCAL, F.ROUTE_KAGGLE}


def test_the_ceiling_is_a_strict_boundary_and_exactly_at_it_stays_local():
    at, _ = F.choose_route(float(F.LOCAL_WALL_CEILING_SECONDS))
    over, _ = F.choose_route(float(F.LOCAL_WALL_CEILING_SECONDS) + 1.0)
    assert at == F.ROUTE_LOCAL and over == F.ROUTE_KAGGLE


def test_choose_route_refuses_a_nonpositive_projection():
    with pytest.raises(ValueError, match="local_seconds must be positive"):
        F.choose_route(0.0)


# --------------------------------------------------------------------------------------------------
# convergence sizing -- the epoch axis is the target
# --------------------------------------------------------------------------------------------------
def test_sizing_aims_for_the_TARGET_EPOCHS_and_derives_the_subsample_to_fit():
    s = F.size_convergence_run(train_seq_per_second=40.0, budget_seconds=900.0,
                               n_available=44106)
    assert s["feasible"] is True
    assert s["max_epochs"] == F.TARGET_CURVE_EPOCHS
    assert s["hit_target_epochs"] is True
    assert s["projected_seconds"] <= 900.0 * 1.05
    assert "epoch axis" in s["why"]


def test_a_BIGGER_budget_buys_more_DATA_at_the_same_epoch_count_not_more_epochs():
    """The direction that was corrected: shape lives on the epoch axis, level on the data axis."""
    small = F.size_convergence_run(train_seq_per_second=40.0, budget_seconds=450.0,
                                   n_available=44106)
    big = F.size_convergence_run(train_seq_per_second=40.0, budget_seconds=900.0,
                                 n_available=44106)
    assert small["max_epochs"] == big["max_epochs"] == F.TARGET_CURVE_EPOCHS
    assert big["n_subsample"] > small["n_subsample"]


def test_sizing_is_capped_by_what_the_cohort_actually_holds():
    s = F.size_convergence_run(train_seq_per_second=1e6, budget_seconds=900.0, n_available=500)
    assert s["n_subsample"] <= 500


def test_sizing_falls_back_toward_the_FLOOR_before_declaring_infeasible():
    """A budget that cannot buy the target but can buy the floor still yields a usable curve."""
    s = F.size_convergence_run(train_seq_per_second=10.0, budget_seconds=60.0,
                               n_available=44106, target_epochs=200, min_epochs=6)
    assert s["feasible"] is True
    assert s["max_epochs"] == 6 and s["hit_target_epochs"] is False


def test_sizing_declares_INFEASIBLE_rather_than_returning_a_two_epoch_curve():
    s = F.size_convergence_run(train_seq_per_second=0.5, budget_seconds=5.0, n_available=44106)
    assert s["feasible"] is False
    assert "too few to say anything about convergence" in s["why"]


def test_sizing_refuses_a_target_below_its_own_floor():
    with pytest.raises(ValueError, match="target_epochs must be >= min_epochs"):
        F.size_convergence_run(train_seq_per_second=40.0, budget_seconds=900.0,
                               n_available=100, target_epochs=2, min_epochs=6)


@pytest.mark.parametrize("rate,budget", [(0.0, 900.0), (40.0, 0.0), (-1.0, 900.0)])
def test_sizing_refuses_nonsense_inputs(rate, budget):
    with pytest.raises(ValueError, match="must both be positive"):
        F.size_convergence_run(train_seq_per_second=rate, budget_seconds=budget, n_available=100)


# --------------------------------------------------------------------------------------------------
# the pretraining recipe and the deviation rail
# --------------------------------------------------------------------------------------------------
_RECIPE = {"learning_rate": 0.001, "patience": 25, "batch_size": 1024, "n_epochs": 2000}


def test_the_recipe_is_read_verbatim_from_the_checkpoints_own_config():
    cfg = dict(_RECIPE, device="cuda", architecture="DNATransformer", subset_fraction=1)
    assert F.read_pretraining_recipe(cfg) == _RECIPE


@pytest.mark.parametrize("cfg", [None, {}])
def test_a_missing_config_is_REFUSED_rather_than_invented(cfg):
    with pytest.raises(ProbeRefused, match="no `config`"):
        F.read_pretraining_recipe(cfg)


def test_a_config_missing_a_needed_key_is_refused_and_names_the_key():
    with pytest.raises(ProbeRefused, match="learning_rate"):
        F.read_pretraining_recipe({"patience": 25, "batch_size": 8, "n_epochs": 10})


def test_an_inherited_value_is_NOT_reported_as_a_deviation():
    used = {"learning_rate": {"value": 0.001, "why": ""}}
    assert F.deviations_from_recipe(_RECIPE, used) == []


def test_a_JUSTIFIED_deviation_is_recorded_with_both_values_and_its_reason():
    used = {"batch_size": {"value": 32, "why": "measured fastest on 4 CPU threads"}}
    out = F.deviations_from_recipe(_RECIPE, used)
    assert out == [{"parameter": "batch_size", "pretraining_value": 1024,
                    "probe_value": 32, "why": "measured fastest on 4 CPU threads"}]


def test_an_UNJUSTIFIED_deviation_RAISES_so_a_silent_substitution_cannot_ship():
    """The honesty rail: deviating is allowed, deviating quietly is not."""
    used = {"batch_size": {"value": 32, "why": ""}}
    with pytest.raises(ProbeRefused, match="no justification recorded"):
        F.deviations_from_recipe(_RECIPE, used)


def test_the_deviation_rail_is_non_vacuous_adding_a_reason_makes_the_same_call_pass():
    used = {"batch_size": {"value": 32, "why": "because measured"}}
    assert len(F.deviations_from_recipe(_RECIPE, used)) == 1


# --------------------------------------------------------------------------------------------------
# curve_turned -- the distinction the arm's interpretability rests on
# --------------------------------------------------------------------------------------------------
def test_a_curve_that_peaked_and_then_stalled_HAS_turned():
    c = F.curve_turned([0.1, 0.3, 0.5, 0.4, 0.45], patience=2)
    assert c["turned"] is True
    assert c["best_epoch"] == 2 and c["epochs_after_best"] == 2
    assert "bounded by observation" in c["why"]


def test_a_curve_still_RISING_at_its_last_epoch_has_NOT_turned_and_says_what_that_costs():
    c = F.curve_turned([0.1, 0.2, 0.3, 0.4, 0.5], patience=2)
    assert c["turned"] is False
    assert c["best_epoch"] == 4 and c["epochs_after_best"] == 0
    assert "LOWER bound" in c["why"]
    assert "could not be distinguished from a representation-limited one" in c["why"], (
        "the consequence -- that a null under this budget cannot be told apart from a "
        "representation-limited one -- must be stated, not left to the reader"
    )


def test_a_curve_whose_peak_is_too_recent_to_judge_has_not_turned_either():
    c = F.curve_turned([0.1, 0.5, 0.4], patience=2)
    assert c["epochs_after_best"] == 1 and c["turned"] is False


def test_an_empty_curve_reports_None_rather_than_False():
    """None means UNMEASURED; False would claim the curve was observed and did not turn."""
    c = F.curve_turned([], patience=2)
    assert c["turned"] is None


# --------------------------------------------------------------------------------------------------
# the verdict set
# --------------------------------------------------------------------------------------------------
def _conv(curve, *, patience=2, shared=0):
    return {"peaks_shared_between_splits": shared,
            "early_stop": F.curve_turned(curve, patience=patience)}


def test_an_unaffordable_budget_is_INDETERMINATE():
    v, _ = F.classify({"feasible": False, "why": "too small"}, None, F.ROUTE_LOCAL)
    assert v == "INDETERMINATE_BUDGET_TOO_SMALL_TO_OBSERVE_CONVERGENCE"


def test_skipping_the_curve_is_INDETERMINATE_and_forbids_the_arm():
    v, why = F.classify({"feasible": True}, None, F.ROUTE_KAGGLE)
    assert v == "INDETERMINATE_CONVERGENCE_NOT_MEASURED"
    assert "must not run" in why


def test_a_LEAKED_inner_split_is_INDETERMINATE_not_a_curve():
    v, why = F.classify({"feasible": True}, _conv([0.1, 0.5, 0.4, 0.3], shared=7),
                        F.ROUTE_KAGGLE)
    assert v == "INDETERMINATE_SPLIT_LEAKED"
    assert "7 peaks" in why


def test_a_curve_that_did_not_turn_is_FEASIBLE_but_flagged():
    v, why = F.classify({"feasible": True}, _conv([0.1, 0.2, 0.3]), F.ROUTE_KAGGLE)
    assert v == "FEASIBLE_BUT_CURVE_DID_NOT_TURN"
    assert "LOWER BOUND" in why


def test_a_turned_curve_is_the_clean_FEASIBLE_verdict():
    v, why = F.classify({"feasible": True}, _conv([0.1, 0.5, 0.4, 0.3]), F.ROUTE_KAGGLE)
    assert v == "FEASIBLE_CURVE_TURNED"
    assert F.ROUTE_KAGGLE in why


def test_the_five_verdict_branches_are_all_distinguishable():
    seen = {
        F.classify({"feasible": False, "why": ""}, None, F.ROUTE_LOCAL)[0],
        F.classify({"feasible": True}, None, F.ROUTE_LOCAL)[0],
        F.classify({"feasible": True}, _conv([0.1, 0.5, 0.4, 0.3], shared=1), F.ROUTE_LOCAL)[0],
        F.classify({"feasible": True}, _conv([0.1, 0.2, 0.3]), F.ROUTE_LOCAL)[0],
        F.classify({"feasible": True}, _conv([0.1, 0.5, 0.4, 0.3]), F.ROUTE_LOCAL)[0],
    }
    assert len(seen) == 5
    assert all(v.startswith("INDETERMINATE_") or v.startswith("FEASIBLE_") for v in seen)


# --------------------------------------------------------------------------------------------------
# refusals through main()
# --------------------------------------------------------------------------------------------------
def test_main_returns_nonzero_and_writes_nothing_when_the_clone_is_absent(tmp_path):
    out = tmp_path / "x.json"
    assert F.main(["--pa-repo", str(tmp_path), "--out", str(out)]) == 2
    assert not out.exists()


@pytest.mark.skipif(not ARTIFACT.is_file(), reason="probe has not been run")
def test_main_returns_nonzero_on_a_moved_checkpoint_hash(monkeypatch, tmp_path):
    monkeypatch.setattr(F, "EXPECTED_SHA256", "0" * 64)
    out = tmp_path / "x.json"
    assert F.main(["--out", str(out), "--skip-convergence"]) == 2
    assert not out.exists()


# --------------------------------------------------------------------------------------------------
# the committed artifact
# --------------------------------------------------------------------------------------------------
@pytest.mark.skipif(not ARTIFACT.is_file(), reason="probe has not been run")
def test_the_artifact_anchors_every_hyperparameter_to_the_checkpoints_own_recipe():
    rec = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    assert rec["pretraining_recipe"] == {"learning_rate": 0.001, "patience": 25,
                                         "batch_size": 1024, "n_epochs": 2000}
    used, dev = rec["hyperparameters_used"], rec["deviations_from_pretraining_recipe"]
    assert used["learning_rate"] == rec["pretraining_recipe"]["learning_rate"], (
        "the learning rate is inherited, so it must NOT appear as a deviation"
    )
    assert "learning_rate" not in {d["parameter"] for d in dev}
    for d in dev:
        assert d["why"].strip(), "%s deviates with no reason" % d["parameter"]
        assert d["pretraining_value"] != d["probe_value"]


@pytest.mark.skipif(not ARTIFACT.is_file(), reason="probe has not been run")
def test_the_artifact_projects_cost_at_the_SAME_batch_size_it_measured():
    """Catches the inconsistency this probe had on its first run: fastest rate, different batch."""
    rec = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    thr = rec["throughput"]
    best_bs = str(thr["best_batch_size"])
    assert thr["per_batch_size"][best_bs]["train_seq_per_second"] == thr["best_train_seq_per_second"]
    assert rec["hyperparameters_used"]["batch_size"] == thr["best_batch_size"], (
        "the batch the arm trains at must be the batch the projected rate was measured at"
    )
    assert "SAME_row" in "".join(thr.keys()) or thr.get("rate_and_batch_come_from_the_SAME_row")


@pytest.mark.skipif(not ARTIFACT.is_file(), reason="probe has not been run")
def test_the_artifact_records_the_tractability_reason_and_the_real_parameter_count():
    ck = json.loads(ARTIFACT.read_text(encoding="utf-8"))["checkpoint"]
    assert ck["n_parameters_via_model_parameters"] == 1_847_428
    assert "1.2B" in ck["tractability_note"]
    assert "NOT comparable" in ck["tractability_note"]


@pytest.mark.skipif(not ARTIFACT.is_file(), reason="probe has not been run")
def test_the_artifact_chose_its_route_from_the_measurement_with_the_ceiling_named():
    r = json.loads(ARTIFACT.read_text(encoding="utf-8"))["route"]
    assert r["chosen"] in (F.ROUTE_LOCAL, F.ROUTE_KAGGLE)
    assert r["ceiling_seconds"] == F.LOCAL_WALL_CEILING_SECONDS
    assert set(r["alternatives_considered"]) == {F.ROUTE_LOCAL, F.ROUTE_KAGGLE}
    # the key NAMES the ceiling as operational; the value must give the operational REASON and
    # say the number is a constant a reader may disagree with, not a per-run judgement
    why = r["ceiling_is_operational_not_statistical"]
    assert "overnight run is recoverable" in why
    assert "abandoned part-way" in why
    assert "constant, not a judgement made per run" in why


@pytest.mark.skipif(not ARTIFACT.is_file(), reason="probe has not been run")
def test_the_artifact_states_the_two_readings_and_ties_them_to_the_curve_FIELD():
    rec = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    t = rec["the_two_readings_that_must_not_be_conflated"]
    assert "REPRESENTATION" in t and "TRAINING BUDGET" in t
    assert "convergence.early_stop.turned" in t, (
        "the distinction must point at the FIELD that settles it, not be left as prose"
    )
    assert rec["comparator"]["frozen_incumbent_median_spearman"] == 0.5076


@pytest.mark.skipif(not ARTIFACT.is_file(), reason="probe has not been run")
def test_the_artifacts_inner_split_is_peak_blocked_with_no_shared_peaks():
    c = json.loads(ARTIFACT.read_text(encoding="utf-8"))["convergence"]
    if c is None:
        pytest.skip("convergence was skipped in the committed run")
    assert c["peaks_shared_between_splits"] == 0, (
        "a peak in both halves makes the curve something other than a held-out measurement"
    )
    assert c["n_train"] > 0 and c["n_val"] > 0
    assert len(c["val_spearman_by_epoch"]) >= F.MIN_CURVE_EPOCHS
    assert "SUBSAMPLE curve" in c["scope"]
    assert "NOT the full-cohort epoch count" in c["scope"]
