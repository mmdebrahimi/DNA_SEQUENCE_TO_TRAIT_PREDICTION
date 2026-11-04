"""Tests for the fine-tune arm's scorer.

SEVEN OUTCOMES, ALL DRIVEN. The bar declares seven and `reachable_outcomes()` asserts the check
order can produce exactly those seven; these tests additionally produce each one from a synthetic
run artifact, because a reachability proof over a static table says nothing about whether the
scorer's own branches can get there.

TWO DOCUMENTED DEPARTURES FROM THE PLAN'S STEP 2, PINNED HERE RATHER THAN QUIETLY DROPPED. The plan
said two things the implementation deliberately does not do, and both changes are argued in the
scorer's module docstring:

  1. The plan said `run_completed` must refuse a `BOUNDED_BY_WALL_CLOCK` status outright, and said
     the tests should "assert a BOUNDED run cannot produce a decision outcome". The horizon is far
     longer than the frozen budget, so a run bounded during the INELIGIBLE DIAGNOSTIC TAIL has
     every seed past its budget and a fully readable verdict -- voiding it would discard a result
     over the loss of epochs the bar already declared can never be the arm's result. So the test
     here asserts the OPPOSITE of the plan's wording, and a separate test pins that a bound which
     actually COST a seed or an epoch is still caught -- by `budget_reached`, the gate that owns it.
  2. The plan said `split_equality` reads the HEADS bar's `n_test_per_seed`. The implementation
     reads `data/promoteratlas/local_reconcile_2026-11-01.json`, because the heads bar's own
     `comparator.source` field names that artifact and its copies are byte-identical to it -- so
     the heads bar is a DERIVED copy and this is the PRIMARY. The plan's "refuse when the two bars'
     frac/n_seeds disagree" test maps onto `verify()` refusing when the RECONCILE artifact's
     `frac`/`n_seeds` disagree with the fine-tune bar's comparator block, which is what is pinned.

Both are recorded so a reader who diffs the plan against the tests finds the reason here instead of
concluding the tests were bent to fit the code.

THE FIXTURES ARE BUILT FROM THE REAL COMMITTED COMPARATOR, not from invented per-seed numbers. A
hand-written `atlas_per_seed` would let a split-equality or pairing bug pass, because the test's
reference and the scorer's reference would both be the test author's invention. Gains are applied
ON TOP of the committed atlas values, so the pairing the bar names is exercised for real.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

import scripts.finetune_arm_score as S
from scripts.finetune_arm_score import GATE_FUNCS, GATE_SEQUENCE, ci_percentiles
from scripts.freeze_finetune_arm_bar import CHECK_ORDER, OUTCOMES, reachable_outcomes
from scripts.positional_arm_score import ScoreError

REPO = Path(__file__).resolve().parents[1]
BAR_PATH = REPO / "wiki/finetune_arm_acceptance_bar.json"
RECONCILE_PATH = REPO / "data/promoteratlas/local_reconcile_2026-11-01.json"

#: `data/` is gitignored, so the committed comparator is absent on a fresh checkout.
pytestmark = pytest.mark.skipif(
    not (BAR_PATH.is_file() and RECONCILE_PATH.is_file()),
    reason="the frozen bar or the committed reconcile comparator is absent here")


def _load(p: Path) -> dict:
    return json.loads(Path(p).read_text(encoding="utf-8"))


# Tight spreads so the bootstrap is genuinely exercised (not a degenerate constant) while the
# resulting interval still lands unambiguously in one band.
ADDS_GAINS = [0.050, 0.052, 0.048, 0.051, 0.049, 0.050, 0.053, 0.047, 0.050, 0.051]
TIES_GAINS = [0.000, 0.010, -0.010, 0.005, -0.005, 0.000, 0.015, -0.015, 0.000, 0.010]
NO_GAIN_GAINS = [-0.100, -0.095, -0.105, -0.098, -0.102, -0.100, -0.097, -0.103, -0.100, -0.099]

PASSING_NULL = [0.01, -0.02, 0.015]
#: max|null| 0.40 against an arm median near 0.56 is a ratio of ~0.72, over the bar's 0.5 ceiling.
FAILING_NULL = [0.40, 0.10, 0.05]

#: A null small in absolute terms, for the NO_GAIN fixtures whose signal is lower.
SMALL_NULL = [0.004, -0.003, 0.002]


def make_run(*, gains=None, status="COMPLETE", n_epochs_read=7, stopped_early=False,
             null_values=PASSING_NULL, horizon=27, schema="finetune-arm-run-v1",
             drop_seeds=(), n_val_delta=None, shared_seeds=(), value_none_seeds=(),
             omit_null=False, per_seed_epochs=None):
    """A synthetic `finetune-arm-run-v1` artifact, keyed off the REAL committed comparator."""
    rec = _load(RECONCILE_PATH)
    gains = list(ADDS_GAINS if gains is None else gains)
    ordered = sorted(rec["atlas_per_seed"].items(), key=lambda kv: int(kv[0]))
    per_seed = {}
    for i, (key, atlas) in enumerate(ordered):
        if key in drop_seeds:
            continue
        read = (per_seed_epochs or {}).get(key, n_epochs_read)
        value = None if key in value_none_seeds else round(float(atlas) + gains[i], 6)
        n_val = int(rec["n_test_per_seed"][key])
        if n_val_delta and key == "0":
            n_val += n_val_delta
        per_seed[key] = {
            "frozen_budget": {
                "value": value,
                "n_epochs_read": read,
                "stopped_early": bool(stopped_early),
                "best_epoch": min(5, read),
            },
            "n_val": n_val,
            "peaks_shared_between_splits": key in shared_seeds,
            "diagnostic_tail": {"epochs": list(range(8, horizon + 1)), "eligible": False},
        }
    run = {
        "schema": schema,
        "status": status,
        "horizon_epochs": horizon,
        "per_seed": per_seed,
        "ineligible_diagnostics_note": "epochs past the frozen budget are pre-declared INELIGIBLE",
    }
    if not omit_null:
        run["null"] = {
            "kind": "shuffled_label",
            "values_at_frozen_budget": list(null_values),
            "untrained_baseline_per_replicate": [0.02, 0.01, 0.03],
            "how_to_read_a_HIGH_null": "compare against the untrained baseline first",
            "scope_limit": "bounds label-shuffling only",
        }
    return run


def score(*, job_state=None, **kw) -> dict:
    """Run the scorer over a synthetic run against the REAL bar and comparator."""
    v = {"bar": _load(BAR_PATH), "run": make_run(**kw), "reconcile": _load(RECONCILE_PATH)}
    return S.run_scoring(v, job_state=job_state)


DECISIONS = ("FINETUNE_ADDS", "FINETUNE_TIES", "FINETUNE_NO_GAIN")


# --------------------------------------------------------------- structural coupling to the bar


def test_reachable_outcomes_equals_the_bars_declared_outcomes():
    """Equality, never containment: containment would let a declared outcome be unreachable."""
    assert reachable_outcomes() == set(OUTCOMES)


def test_reachable_equality_is_non_vacuous_under_PADDING():
    """Add an outcome the order cannot produce -- the equality must break."""
    padded = CHECK_ORDER + (("extra_gate", "INDETERMINATE_NOT_DECLARED"),)
    assert reachable_outcomes(check_order=padded) != set(OUTCOMES)


def test_reachable_equality_is_non_vacuous_under_TRIMMING():
    """Drop a gate -- its outcome becomes unreachable and the equality must break."""
    trimmed = tuple(e for e in CHECK_ORDER if e[0] != "run_completed")
    reach = reachable_outcomes(check_order=trimmed)
    assert reach != set(OUTCOMES)
    assert "INDETERMINATE_RUN_DIED" not in reach


def test_every_gate_in_the_bars_order_has_an_implementation():
    """A gate named in the frozen order with no implementation would raise only after the run
    artifact had been read -- i.e. at the worst possible moment."""
    assert set(GATE_FUNCS) == {g for g, _ in GATE_SEQUENCE}


def test_the_scorer_runs_the_gates_in_the_bars_declared_ORDER():
    expected = tuple(g for g, o in CHECK_ORDER if isinstance(o, str))
    assert tuple(g for g, _ in GATE_SEQUENCE) == expected


def test_the_primary_decision_row_is_not_treated_as_a_gate():
    """`primary_decision` carries a 3-tuple of outcomes, not a single one; it is the classifier,
    not a gate, and giving it a gate function would make the decision refusable."""
    assert "primary_decision" not in GATE_FUNCS
    assert "primary_decision" not in {g for g, _ in GATE_SEQUENCE}


# --------------------------------------------------------------- the derived interval width


def test_ci_percentiles_are_derived_from_the_real_frozen_bar():
    assert ci_percentiles(_load(BAR_PATH)) == [2.5, 97.5]


def test_ci_percentiles_REFUSE_a_corrected_bar():
    bar = _load(BAR_PATH)
    bar["multiplicity"] = dict(bar["multiplicity"], correction="bonferroni")
    with pytest.raises(ScoreError, match="single UNCORRECTED comparison"):
        ci_percentiles(bar)


def test_ci_percentiles_REFUSE_more_than_one_declared_comparison():
    bar = _load(BAR_PATH)
    bar["multiplicity"] = dict(bar["multiplicity"], n_declared_comparisons=5)
    with pytest.raises(ScoreError):
        ci_percentiles(bar)


# --------------------------------------------------------------- the three decision outcomes


def test_FINETUNE_ADDS_when_the_whole_interval_clears_the_margin():
    out = score(gains=ADDS_GAINS)
    assert out["verdict"] == "FINETUNE_ADDS"
    assert out["scored"]["gain_bootstrap"]["ci_lo"] > out["margin"]


def test_FINETUNE_TIES_when_the_interval_straddles():
    out = score(gains=TIES_GAINS)
    assert out["verdict"] == "FINETUNE_TIES"
    b = out["scored"]["gain_bootstrap"]
    assert b["ci_lo"] <= out["margin"] and b["ci_hi"] >= 0.0


def test_FINETUNE_NO_GAIN_when_the_whole_interval_lies_below_zero():
    out = score(gains=NO_GAIN_GAINS, null_values=SMALL_NULL)
    assert out["verdict"] == "FINETUNE_NO_GAIN"
    assert out["scored"]["gain_bootstrap"]["ci_hi"] < 0.0


def test_a_real_gain_smaller_than_the_margin_is_a_TIE_not_a_win():
    """The margin band is absorbed into TIES in advance; this is the whole point of fixing it."""
    small = [0.012, 0.011, 0.013, 0.012, 0.012, 0.011, 0.013, 0.012, 0.012, 0.013]
    out = score(gains=small)
    assert out["scored"]["gain_bootstrap"]["ci_lo"] > 0.0
    assert out["verdict"] == "FINETUNE_TIES"


# --------------------------------------------------------------- INDETERMINATE_RUN_DIED


def test_RUN_DIED_on_an_unrecognised_run_schema():
    assert score(schema="something-else-v9")["verdict"] == "INDETERMINATE_RUN_DIED"


def test_RUN_DIED_when_the_run_records_no_per_seed_results():
    out = score(drop_seeds=tuple(str(i) for i in range(10)))
    assert out["verdict"] == "INDETERMINATE_RUN_DIED"


def test_RUN_DIED_when_the_launcher_reports_DEAD_NO_SENTINEL():
    """Whatever is on disk is then whatever survived -- not a result."""
    out = score(job_state="DEAD_NO_SENTINEL")
    assert out["verdict"] == "INDETERMINATE_RUN_DIED"
    assert "DEAD_NO_SENTINEL" in out["gates"]["run_completed"]["detail"]


def test_RUN_DIED_when_the_launcher_refused_to_start():
    """Exit 9 means the payload never ran, which is a different fact from running and failing."""
    out = score(job_state="REFUSED_LOCK_HELD")
    assert out["verdict"] == "INDETERMINATE_RUN_DIED"
    assert "never ran" in out["gates"]["run_completed"]["detail"]


def test_a_COMPLETE_job_state_does_not_trip_run_completed():
    assert score(job_state="COMPLETE")["gates"]["run_completed"]["passed"] is True


# --------------------------------------------------------------- INDETERMINATE_BUDGET_EXHAUSTED


def test_BUDGET_EXHAUSTED_when_a_seed_is_missing_and_it_does_NOT_score_the_rest():
    """Scoring the seeds it happens to have would silently change the declared design."""
    out = score(drop_seeds=("7",))
    assert out["verdict"] == "INDETERMINATE_BUDGET_EXHAUSTED"
    assert "7" in out["gates"]["budget_reached"]["missing_seeds"]
    assert "scored" not in out
    assert "arm_metric" not in out


def test_BUDGET_EXHAUSTED_when_a_seed_fell_short_without_stopping_early():
    out = score(per_seed_epochs={"3": 4})
    assert out["verdict"] == "INDETERMINATE_BUDGET_EXHAUSTED"
    assert any(s.startswith("3:") for s in out["gates"]["budget_reached"]["short_seeds"])


def test_BUDGET_EXHAUSTED_when_a_seed_has_no_frozen_budget_value():
    """A null value is a MISSING measurement, not a measurement of zero."""
    assert score(value_none_seeds=("2",))["verdict"] == "INDETERMINATE_BUDGET_EXHAUSTED"


def test_a_seed_that_stopped_early_INSIDE_the_budget_PASSES():
    """The budget is a MAXIMUM. The opposite reading fails every well-converged run."""
    out = score(n_epochs_read=5, stopped_early=True)
    assert out["gates"]["budget_reached"]["passed"] is True
    assert out["verdict"] in DECISIONS


def test_the_budget_rule_is_recorded_in_the_artifact_as_a_MAXIMUM():
    assert "MAXIMUM" in score()["gates"]["budget_reached"]["rule"]


# --------------------------------------------------------------- INDETERMINATE_SPLIT_MISMATCH


def test_SPLIT_MISMATCH_when_a_seeds_validation_size_differs_from_the_committed_reference():
    assert score(n_val_delta=-1)["verdict"] == "INDETERMINATE_SPLIT_MISMATCH"


def test_SPLIT_MISMATCH_when_peaks_are_shared_between_the_two_splits():
    """A shared peak is leakage; equal sizes do not rule it out, so it is checked separately."""
    out = score(shared_seeds=("4",))
    assert out["verdict"] == "INDETERMINATE_SPLIT_MISMATCH"
    assert "4" in str(out["gates"]["split_equality"]["detail"])


def test_split_equality_names_the_PRIMARY_reference_it_read():
    g = score()["gates"]["split_equality"]
    assert g["committed_reference"] == "data/promoteratlas/local_reconcile_2026-11-01.json"


# --------------------------------------------------------------- INDETERMINATE_NULL_FAIL


def test_NULL_FAIL_when_the_null_is_large_relative_to_the_arms_own_signal():
    out = score(null_values=FAILING_NULL)
    assert out["verdict"] == "INDETERMINATE_NULL_FAIL"
    assert out["gates"]["shuffled_label_null"]["max_abs_over_signal"] >= 0.5


def test_a_MISSING_null_is_not_a_passing_null():
    """The dangerous failure: an absent null reads as 'nothing went wrong'."""
    out = score(omit_null=True)
    assert out["verdict"] == "INDETERMINATE_NULL_FAIL"
    assert out["gates"]["shuffled_label_null"]["n_replicates"] == 0


def test_the_null_gate_reports_the_untrained_baseline_beside_it():
    """A high null has two causes that look identical in the number; only the baseline separates
    them, so it must travel with the null rather than be looked up later."""
    g = score(null_values=FAILING_NULL)["gates"]["shuffled_label_null"]
    assert g["untrained_baseline_per_replicate"] == [0.02, 0.01, 0.03]
    assert "untrained baseline" in g["detail"]


def test_the_null_is_scaled_against_the_signal_not_an_absolute_floor():
    g = score()["gates"]["shuffled_label_null"]
    assert g["ratio_ceiling"] == 0.5
    assert g["arm_median_signal"] > 0.5


# --------------------------------------------------------------- the corrected BOUNDED semantics


def test_a_run_bounded_in_the_ineligible_tail_still_decides():
    """DEPARTS FROM THE PLAN'S WORDING, DELIBERATELY -- see the module docstring. Every seed is
    past the frozen budget, so the verdict is readable, and voiding it would throw away a result
    over epochs the bar already declared INELIGIBLE."""
    out = score(status="BOUNDED_BY_WALL_CLOCK")
    assert out["gates"]["run_completed"]["passed"] is True
    assert out["verdict"] in DECISIONS
    assert out["run_status"] == "BOUNDED_BY_WALL_CLOCK"


def test_a_bound_that_actually_COST_a_seed_is_caught_by_budget_reached():
    """The case the plan was protecting against, caught by the gate that owns it."""
    out = score(status="BOUNDED_BY_WALL_CLOCK", drop_seeds=("9",))
    assert out["verdict"] == "INDETERMINATE_BUDGET_EXHAUSTED"


def test_a_bound_that_cost_EPOCHS_is_caught_by_budget_reached():
    out = score(status="BOUNDED_BY_WALL_CLOCK", n_epochs_read=3)
    assert out["verdict"] == "INDETERMINATE_BUDGET_EXHAUSTED"


def test_run_completed_records_why_bounded_is_not_fatal_there():
    assert "INELIGIBLE" in score()["gates"]["run_completed"]["bounded_is_not_fatal_here"]


# --------------------------------------------------------------- gate order and short-circuiting


def test_an_earlier_gate_failure_short_circuits_the_later_gates():
    """Two things broken at once must report the FIRST in the bar's order, and must not publish
    partial findings from gates that never ran."""
    out = score(schema="wrong-v0", n_val_delta=-5, null_values=FAILING_NULL)
    assert out["verdict"] == "INDETERMINATE_RUN_DIED"
    assert set(out["gates"]) == {"run_completed"}


def test_no_arm_metric_is_published_when_any_gate_fails():
    for kw in ({"schema": "x"}, {"drop_seeds": ("1",)}, {"n_val_delta": 3},
               {"null_values": FAILING_NULL}):
        out = score(**kw)
        assert "arm_metric" not in out and "scored" not in out, kw


def test_a_passing_run_publishes_the_metric_and_all_four_gates():
    out = score()
    assert set(out["gates"]) == {g for g, _ in GATE_SEQUENCE}
    assert out["arm_metric"] == out["scored"]["median_spearman"]


# --------------------------------------------------------------- re-derivability


@pytest.mark.parametrize("kw", [
    {}, {"gains": TIES_GAINS}, {"gains": NO_GAIN_GAINS, "null_values": SMALL_NULL},
    {"schema": "x"}, {"drop_seeds": ("5",)}, {"n_val_delta": 2}, {"null_values": FAILING_NULL},
])
def test_the_verdict_is_rederivable_from_the_artifact_plus_the_bar_alone(kw):
    """Nothing outside the emitted artifact and the frozen bar may be needed to reach the verdict,
    or a reader cannot check it."""
    out = score(**kw)
    verdict, _detail = S.verdict_from_artifact(out, _load(BAR_PATH))
    assert verdict == out["verdict"]
    assert out["verdict_is_rederivable_from_this_artifact"] is True


def test_a_verdict_is_always_one_of_the_declared_outcomes():
    for kw in ({}, {"gains": TIES_GAINS}, {"schema": "x"}, {"drop_seeds": ("0",)},
               {"n_val_delta": 1}, {"omit_null": True}):
        assert score(**kw)["verdict"] in OUTCOMES, kw


def test_rederivation_refuses_a_gate_that_was_never_evaluated():
    out = score()
    out["gates"].pop("split_equality")
    verdict, detail = S.verdict_from_artifact(out, _load(BAR_PATH))
    assert verdict == "INDETERMINATE_SPLIT_MISMATCH"
    assert "never evaluated" in detail


# --------------------------------------------------------------- verify() refusals


def _anchors(tmp_path, *, bar_edit=None, rec_edit=None, run_edit=None):
    bar, rec, run = _load(BAR_PATH), _load(RECONCILE_PATH), make_run()
    for obj, edit in ((bar, bar_edit), (rec, rec_edit), (run, run_edit)):
        if edit:
            obj.update(edit)
    paths = {}
    for name, obj in (("bar", bar), ("rec", rec), ("run", run)):
        p = tmp_path / ("%s.json" % name)
        p.write_text(json.dumps(obj), encoding="utf-8")
        paths[name] = p
    return paths


def test_verify_accepts_the_real_anchors(tmp_path):
    p = _anchors(tmp_path)
    v = S.verify(bar_path=p["bar"], run_path=p["run"], reconcile_path=p["rec"])
    assert v["reconcile"]["verdict"] == "RECONCILED"


def test_verify_REFUSES_a_frac_disagreement_rather_than_comparing_two_designs(tmp_path):
    """The plan's 'two references disagree' refusal, at the surface that holds the reference."""
    p = _anchors(tmp_path, rec_edit={"frac": 0.5})
    with pytest.raises(ScoreError, match="different designs"):
        S.verify(bar_path=p["bar"], run_path=p["run"], reconcile_path=p["rec"])


def test_verify_REFUSES_an_n_seeds_disagreement(tmp_path):
    p = _anchors(tmp_path, rec_edit={"n_seeds": 8})
    with pytest.raises(ScoreError, match="n_seeds"):
        S.verify(bar_path=p["bar"], run_path=p["run"], reconcile_path=p["rec"])


def test_verify_REFUSES_a_comparator_that_did_not_reconcile(tmp_path):
    p = _anchors(tmp_path, rec_edit={"verdict": "MISMATCHED"})
    with pytest.raises(ScoreError, match="not RECONCILED"):
        S.verify(bar_path=p["bar"], run_path=p["run"], reconcile_path=p["rec"])


def test_verify_REFUSES_a_comparator_recording_mismatches(tmp_path):
    p = _anchors(tmp_path, rec_edit={"mismatched": ["3"]})
    with pytest.raises(ScoreError, match="mismatches"):
        S.verify(bar_path=p["bar"], run_path=p["run"], reconcile_path=p["rec"])


def test_verify_REFUSES_comparator_median_drift_from_the_bars_inherited_value(tmp_path):
    """The bar inherited a number; if the comparator now reads differently, one of them moved."""
    p = _anchors(tmp_path, rec_edit={"atlas_median": 0.6})
    with pytest.raises(ScoreError, match="inherits"):
        S.verify(bar_path=p["bar"], run_path=p["run"], reconcile_path=p["rec"])


def test_verify_REFUSES_an_unrecognised_bar_schema(tmp_path):
    p = _anchors(tmp_path, bar_edit={"schema": "not-the-bar-v0"})
    with pytest.raises(ScoreError, match="bar schema"):
        S.verify(bar_path=p["bar"], run_path=p["run"], reconcile_path=p["rec"])


def test_verify_REFUSES_an_absent_anchor(tmp_path):
    p = _anchors(tmp_path)
    with pytest.raises(ScoreError, match="anchor source absent"):
        S.verify(bar_path=p["bar"], run_path=tmp_path / "nope.json", reconcile_path=p["rec"])


# --------------------------------------------------------------- what the artifact must carry


def test_the_artifact_records_the_interval_derivation_not_just_the_number():
    out = score()
    assert out["ci_percentiles"] == [2.5, 97.5]
    assert "DERIVED" in out["ci_percentiles_derivation"]


def test_the_artifact_names_the_comparator_source_and_why_it_is_the_primary():
    c = score()["comparator"]
    assert c["source"] == "data/promoteratlas/local_reconcile_2026-11-01.json"
    assert "byte-identical" in c["why_the_primary_not_the_heads_bar"]


def test_the_artifact_marks_the_diagnostic_tail_ineligible():
    out = score()
    assert out["ineligible_diagnostics_note"]
    tails = out["diagnostic_tail_per_seed"]
    assert tails and all(t["eligible"] is False for t in tails.values())


def test_the_scored_block_records_that_the_frozen_rule_was_READ_not_recomputed():
    """Two implementations of one early-stopping rule disagree silently on exactly the seeds that
    stopped early, which is the hardest case to notice."""
    assert "READ from the run artifact" in score()["scored"]["frozen_budget_rule_source"]


def test_the_per_seed_pairing_uses_the_committed_atlas_values():
    out = score(gains=ADDS_GAINS)
    rec = _load(RECONCILE_PATH)
    for key, row in out["scored"]["per_seed"].items():
        assert row["atlas_committed"] == rec["atlas_per_seed"][key]
        assert row["gain"] == pytest.approx(
            row["frozen_budget_spearman"] - row["atlas_committed"], abs=1e-4)


def test_seed_wins_are_counted_against_the_committed_comparator():
    assert score(gains=ADDS_GAINS)["scored"]["seed_wins_over_atlas"] == 10
    assert score(gains=NO_GAIN_GAINS,
                 null_values=SMALL_NULL)["scored"]["seed_wins_over_atlas"] == 0
