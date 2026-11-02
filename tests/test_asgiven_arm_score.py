"""Tests for the as-given arm's scorer.

Every gate branch is driven with a stub scorer and proven REACHABLE, as on the positional arm, and
the parity gate is proven NON-VACUOUS against a real perturbation. The synthetic fixtures are
IMPORTED from the positional scorer's tests rather than copied, because two sets of drifting
fixtures would be two chances to test different things while looking the same.

One test here has no sibling on the positional arm and is the most interesting: the as-given rung
was already scored once, as that arm's INELIGIBLE secondary. **The two runs must agree on the same
number.** If they do not, one of them is wrong, and a disagreement would be far more informative
than either value.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from scripts.freeze_asgiven_arm_bar import OUTCOMES, PRIMARY
from scripts.positional_arm_probe import RAM_HEADROOM
from scripts.positional_arm_score import ScoreError, median_of
from tests.test_positional_arm_score import (
    CHANNELS,
    N_SEEDS,
    FRAC,
    make_block,
    make_cohort,
    stub_scorer,
)
import scripts.asgiven_arm_score as sc
from scripts.positional_arm_score import leave_peak_out_indices

REPO = Path(__file__).resolve().parents[1]
POS_RESULT = REPO / "wiki/positional_arm_result_2026-11-01.json"
ASG_RESULT = REPO / "wiki/asgiven_arm_result_2026-11-02.json"


def make_asgiven_bar(tiles, *, atlas=0.50, n_test=None, parity_tol=1e-5):
    """A bar in the AS-GIVEN shape: one primary, plain 95%, tolerance under validity_controls."""
    if n_test is None:
        n_test = {}
        for s in range(N_SEEDS):
            _tr, te, _h = leave_peak_out_indices(tiles, frac=FRAC, seed=s)
            n_test[str(s)] = len(te)
    per_seed = {str(s): round(atlas + 0.001 * s, 4) for s in range(N_SEEDS)}
    return {
        "arm": "asgiven-position-resolved-promoteratlas-150bp",
        "criterion": "ceiling: synthetic",
        "inherited_not_chosen": {
            "margin": {"value": 0.02},
            "comparator": {"atlas_per_seed": per_seed,
                           "atlas_median": median_of(per_seed.values()),
                           "source": "synthetic"},
            "split": {"frac": FRAC, "n_seeds": N_SEEDS, "n_test_per_seed": n_test},
        },
        "validity_controls": {"parity": {"tolerance": parity_tol}},
        "multiplicity": {"ci_percentiles": [2.5, 97.5], "n_primaries": 1, "primary": PRIMARY},
        "NOT_A_BLIND_PRE_REGISTRATION": {"already_observed_median": 0.5086},
        "honest_limits_fixed_in_advance": ["synthetic"],
        "forecast": {"predicted_outcome": "ASGIVEN_TIES"},
    }


def _passing_parity(bar):
    return [bar["inherited_not_chosen"]["comparator"]["atlas_per_seed"][str(s)]
            for s in range(N_SEEDS)]


def go(tiles, st, block, bar, *, scorer, null_extension=0, emb_cache=None):
    return sc.run(block_dir=Path("unused"), bar=bar, null_extension=null_extension, n_boot=200,
                  scorer=scorer, cohort=(tiles, st), block=block, emb_cache=emb_cache,
                  verbose=False)


# --------------------------------------------------------------------------------------------------
# every gate branch is reachable
# --------------------------------------------------------------------------------------------------
def test_split_mismatch_is_reachable_and_stops_the_run():
    tiles, st = make_cohort()
    bar = make_asgiven_bar(tiles, n_test={str(s): 1 for s in range(N_SEEDS)})
    rec = go(tiles, st, make_block(), bar, scorer=stub_scorer(parity_by_seed=[0.5] * N_SEEDS))
    assert rec["verdict"] == "INDETERMINATE_SPLIT_MISMATCH"
    assert "parity" not in rec["gates"]
    assert "rung" not in rec


def test_parity_failure_is_reachable_and_stops_the_run():
    tiles, st = make_cohort()
    bar = make_asgiven_bar(tiles)
    rec = go(tiles, st, make_block(), bar, scorer=stub_scorer(parity_by_seed=[0.90] * N_SEEDS))
    assert rec["verdict"] == "INDETERMINATE_PARITY_FAIL"
    assert "shuffled_label_null" not in rec["gates"]
    assert "rung" not in rec


def test_null_failure_is_reachable_and_stops_the_run():
    tiles, st = make_cohort()
    bar = make_asgiven_bar(tiles)
    rec = go(tiles, st, make_block(), bar,
             scorer=stub_scorer(parity_by_seed=_passing_parity(bar), null_value=0.25),
             null_extension=0)
    assert rec["verdict"] == "INDETERMINATE_NULL_FAIL"
    assert "rung" not in rec
    assert "MORE SEEDS" in rec["gates"]["shuffled_label_null"]["detail"]


def test_a_perturbed_block_trips_the_parity_gate(tmp_path):
    """NON-VACUITY. An always-passing parity gate would look identical from the artifact."""
    from scripts.positional_arm_score import build_rung

    tiles, st = make_cohort()
    bar = make_asgiven_bar(tiles)
    committed = build_rung(make_block()[0], make_block()[2], "global_mean",
                           [t.strand for t in tiles])
    cache = tmp_path / "emb.npy"
    np.save(cache, committed.astype(np.float32))
    atlas = _passing_parity(bar)

    clean = go(tiles, st, make_block(), bar, scorer=stub_scorer(parity_by_seed=atlas),
               emb_cache=cache)
    assert clean["gates"]["parity"]["passed"] is True

    dirty = go(tiles, st, make_block(perturb=(3, 7, 0.5)), bar,
               scorer=stub_scorer(parity_by_seed=atlas), emb_cache=cache)
    assert dirty["verdict"] == "INDETERMINATE_PARITY_FAIL"
    assert dirty["gates"]["parity"]["embedding_level"]["passed"] is False


def test_the_null_runs_on_THIS_arms_primary_not_the_positional_one():
    """`ScoreCtx.primary_rung` is the one thing that differs; if it were ignored the null would be
    computed on a rung this arm does not declare."""
    tiles, st = make_cohort()
    bar = make_asgiven_bar(tiles)
    rec = go(tiles, st, make_block(), bar, scorer=stub_scorer(parity_by_seed=_passing_parity(bar)))
    assert rec["gates"]["shuffled_label_null"]["rung"] == PRIMARY
    assert rec["primary_rung"] == PRIMARY


# --------------------------------------------------------------------------------------------------
# the three terminal verdicts
# --------------------------------------------------------------------------------------------------
@pytest.mark.parametrize("value,want", [
    (0.60, "ASGIVEN_ADDS"),       # gain ~0.099, interval well above the margin
    (0.30, "ASGIVEN_NO_GAIN"),    # gain ~-0.20, interval below zero
    (0.5015, "ASGIVEN_TIES"),     # gain ~0.0005, straddles zero
    (0.51, "ASGIVEN_TIES"),       # gain ~0.010: a REAL but SUB-MARGIN gain -> still a tie
])
def test_all_three_terminal_verdicts_are_reachable(value, want):
    tiles, st = make_cohort()
    bar = make_asgiven_bar(tiles)
    rec = go(tiles, st, make_block(), bar,
             scorer=stub_scorer(parity_by_seed=_passing_parity(bar),
                                by_width={CHANNELS * 10: value}))
    assert rec["verdict"] == want, "gain %+.4f" % rec["rung"]["gain_bootstrap"]["median"]
    assert rec["verdict"] in OUTCOMES


def test_the_union_of_reachable_verdicts_equals_the_bars_closed_set():
    tiles, st = make_cohort()
    bar = make_asgiven_bar(tiles)
    seen = {
        go(tiles, st, make_block(), make_asgiven_bar(tiles, n_test={str(s): 1 for s in range(N_SEEDS)}),
           scorer=stub_scorer(parity_by_seed=[0.5] * N_SEEDS))["verdict"],
        go(tiles, st, make_block(), bar,
           scorer=stub_scorer(parity_by_seed=[0.9] * N_SEEDS))["verdict"],
        go(tiles, st, make_block(), bar,
           scorer=stub_scorer(parity_by_seed=_passing_parity(bar), null_value=0.25),
           null_extension=0)["verdict"],
    }
    for v in (0.60, 0.30, 0.5015):
        seen.add(go(tiles, st, make_block(), bar,
                    scorer=stub_scorer(parity_by_seed=_passing_parity(bar),
                                       by_width={CHANNELS * 10: v}))["verdict"])
    assert seen == set(OUTCOMES)


def test_a_rung_over_the_ram_budget_yields_TIES_not_a_no_gain_claim():
    """A rung nobody could afford to score is not evidence of no gain."""
    tiles, st = make_cohort()
    bar = make_asgiven_bar(tiles)
    import scripts.asgiven_arm_score as mod
    old = mod.RAM_HEADROOM
    try:
        mod.RAM_HEADROOM = 1e12
        rec = go(tiles, st, make_block(), bar,
                 scorer=stub_scorer(parity_by_seed=_passing_parity(bar)))
    finally:
        mod.RAM_HEADROOM = old
    assert rec["rung"]["skipped"] is True
    assert rec["verdict"] == "ASGIVEN_TIES"
    assert "not evidence of no gain" in rec["rung"]["why"]


# --------------------------------------------------------------------------------------------------
# re-derivability
# --------------------------------------------------------------------------------------------------
@pytest.mark.parametrize("value", [0.60, 0.30, 0.5015, 0.51])
def test_the_verdict_is_rederivable_from_the_emitted_artifact_alone(value):
    tiles, st = make_cohort()
    bar = make_asgiven_bar(tiles)
    rec = go(tiles, st, make_block(), bar,
             scorer=stub_scorer(parity_by_seed=_passing_parity(bar),
                                by_width={CHANNELS * 10: value}))
    again, _ = sc.verdict_from_artifact(json.loads(json.dumps(rec)), bar)
    assert again == rec["verdict"]
    assert rec["verdict_is_rederivable_from_this_artifact"] is True


def test_the_rederivation_treats_a_never_evaluated_gate_as_that_gates_failure():
    bar = make_asgiven_bar(make_cohort()[0])
    assert sc.verdict_from_artifact({"gates": {}}, bar)[0] == "INDETERMINATE_SPLIT_MISMATCH"
    assert sc.verdict_from_artifact(
        {"gates": {"split_equality": {"passed": True}}}, bar)[0] == "INDETERMINATE_PARITY_FAIL"
    assert sc.verdict_from_artifact(
        {"gates": {"split_equality": {"passed": True}, "parity": {"passed": True}}},
        bar)[0] == "INDETERMINATE_NULL_FAIL"


def test_run_refuses_a_block_whose_row_count_differs_from_the_cohort():
    tiles, st = make_cohort()
    bar = make_asgiven_bar(tiles)
    with pytest.raises(ScoreError, match="pairing is impossible"):
        go(tiles, st, make_block(n=199), bar,
           scorer=stub_scorer(parity_by_seed=_passing_parity(bar)))


# --------------------------------------------------------------------------------------------------
# the two runs must agree -- this rung was already scored once, as an ineligible secondary
# --------------------------------------------------------------------------------------------------
@pytest.mark.skipif(not (ASG_RESULT.is_file() and POS_RESULT.is_file()),
                    reason="both artifacts needed; run the as-given scorer")
def test_the_as_given_median_REPRODUCES_the_positional_runs_secondary():
    """If these disagree one of them is wrong, and the disagreement would matter more than either."""
    asg = json.loads(ASG_RESULT.read_text(encoding="utf-8"))
    pos = json.loads(POS_RESULT.read_text(encoding="utf-8"))
    if asg["rung"].get("skipped"):
        pytest.skip("the primary was skipped on the RAM pre-check")

    here = asg["rung"]
    there = pos["rungs"][PRIMARY]
    assert here["median_spearman"] == there["median_spearman"], (
        "the same rung on the same splits scored %.4f as a primary and %.4f as a secondary"
        % (here["median_spearman"], there["median_spearman"]))
    assert here["per_seed"] == there["per_seed"], "per-seed values differ between the two runs"
    assert here["seed_wins_over_atlas"] == there["seed_wins_over_atlas"]
    assert here["n_features"] == there["n_features"]


@pytest.mark.skipif(not ASG_RESULT.is_file(), reason="run the as-given scorer")
def test_the_committed_result_carries_the_hazard_and_a_closed_set_verdict():
    asg = json.loads(ASG_RESULT.read_text(encoding="utf-8"))
    assert asg["verdict"] in OUTCOMES
    assert asg["primary_rung"] == PRIMARY
    assert asg["comparator"]["recomputed_in_run"] is False
    h = asg["not_a_blind_pre_registration"]
    assert "not a blind pre-registration" in h["statement"].lower()
    assert h["already_observed_median"] == 0.5086
    assert asg["multiplicity"]["n_primaries"] == 1


@pytest.mark.skipif(not ASG_RESULT.is_file(), reason="run the as-given scorer")
def test_the_committed_result_passed_every_gate_or_says_which_one_did_not():
    asg = json.loads(ASG_RESULT.read_text(encoding="utf-8"))
    if asg["verdict"].startswith("INDETERMINATE"):
        failed = [g for g, d in asg["gates"].items() if not d["passed"]]
        assert failed, "an INDETERMINATE verdict with every gate passing is incoherent"
        return
    for name, g in asg["gates"].items():
        assert g["passed"] is True, "%s did not pass but the verdict is terminal" % name
    assert asg["gates"]["parity"]["derived_median"] == asg["comparator"]["atlas_median"], (
        "parity passed but the derived median does not equal the committed one")
