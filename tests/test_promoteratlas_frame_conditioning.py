"""Tests for the frame-conditioning experiment.

The strongest test here is `test_the_recorded_verdict_is_reproducible_from_the_artifacts_own_numbers`:
it rebuilds the grader's input out of the committed artifact and re-grades it. That makes the
artifact self-checking -- a figure edited by hand, or a verdict copied from an earlier run, stops
agreeing with the numbers printed beside it.

The two unit tests cover the claims the experiment's design rests on: that the strengthening
diagnostic really is monotone in added columns (so the plumbing gate cannot fire on a rounding
artifact), and that the stability summary actually detects a statistic straddling its own bar (which
is what this run turned out to need).
"""

from __future__ import annotations

import json
import pathlib

import numpy as np
import pytest

from dna_decode.glm.frame_verdict import CHECK_ORDER, FrameInput, grade
from scripts.promoteratlas_frame_conditioning import _ols_train_r2, _stability

REPO = pathlib.Path(__file__).resolve().parents[1]
ART = REPO / "wiki/promoteratlas_frame_conditioning_2026-11-01.json"
BAR = REPO / "wiki/promoteratlas_frame_acceptance_bar.json"
ARM = REPO / "wiki/promoteratlas_tile_arm_2026-11-01.json"


# --------------------------------------------------------------------------------- unit tests


def test_adding_a_column_cannot_lower_the_ols_training_fit():
    """Why the plumbing gate reads OLS and not ridge: this monotonicity is exact for OLS.

    A ridge penalty can move a training fit a hair in either direction, so a strict inequality on
    the ridge figure could report a plumbing fault that is really a rounding artifact.
    """
    rng = np.random.default_rng(0)
    X = rng.normal(size=(300, 3))
    y = np.abs(rng.normal(size=300)) * 5.0
    base = _ols_train_r2(X, y)
    for _ in range(5):
        X2 = np.hstack([X, rng.normal(size=(300, 1))])
        assert _ols_train_r2(X2, y) >= base - 1e-9
        X = X2
        base = _ols_train_r2(X, y)


def test_a_pure_noise_column_still_cannot_lower_it():
    rng = np.random.default_rng(7)
    X = rng.normal(size=(200, 2))
    y = np.abs(X[:, 0] * 3.0 + 10.0)
    with_noise = np.hstack([X, rng.normal(size=(200, 1))])
    assert _ols_train_r2(with_noise, y) >= _ols_train_r2(X, y) - 1e-9


def test_stability_detects_a_statistic_straddling_its_bar():
    st = _stability([0.05, -0.02, 0.01, 0.09, -0.04], bar=0.02)
    assert st["straddles_bar"] is True
    assert st["straddles_zero"] is True
    assert st["n_seeds_below_bar"] == 3
    assert st["n_seeds_negative"] == 2


def test_stability_reports_a_clean_statistic_as_not_straddling():
    st = _stability([0.25, 0.27, 0.26, 0.30], bar=0.02)
    assert st["straddles_bar"] is False
    assert st["straddles_zero"] is False
    assert st["n_seeds_below_bar"] == 0


def test_stability_survives_nan_seeds_without_claiming_them():
    st = _stability([0.1, float("nan"), 0.3], bar=0.02)
    assert st["n_seeds"] == 2


# ------------------------------------------------------------------------------- the artifact

pytestmark = pytest.mark.skipif(not ART.is_file(), reason="experiment has not been run")


@pytest.fixture(scope="module")
def art() -> dict:
    return json.loads(ART.read_text(encoding="utf-8"))


def test_the_anchor_reproduced_the_committed_arm_before_anything_else_was_read(art):
    a = art["anchor"]
    committed = float(json.loads(ARM.read_text(encoding="utf-8"))["arms"]["native150"]["median"])
    assert a["committed_atlas_native150_median"] == committed
    assert a["abs_delta"] <= a["tolerance"]
    assert a["reconciled"] is True


def test_the_tile_accounting_reconciles_exactly(art):
    acc = art["accounting"]
    assert acc["n_tiles_with_a_frame"] + acc["n_tiles_excluded"] == acc["n_tiles_scored_by_the_arm"]
    assert acc["reconciles"] is True
    assert acc["n_tiles_scored_by_the_arm"] == 44106
    assert (
        acc["distinct_unique"] + acc["distinct_ambiguous"] + acc["distinct_not_found"]
        == acc["distinct_sequences"]
    )


def test_the_verdict_is_one_the_grader_can_actually_reach(art):
    assert art["verdict"] in CHECK_ORDER


def test_the_recorded_verdict_is_reproducible_from_the_artifacts_own_numbers(art):
    """Rebuild the grader's input from the artifact and re-grade. Makes the artifact self-checking."""
    c2, c3 = art["c2_matched_strength_series"], art["c3_stratified"]
    fi = FrameInput(
        harness_reconciled=art["anchor"]["reconciled"],
        unique_locate_rate=art["accounting"]["unique_locate_rate"],
        frame_alone_spearman=c2["frame_alone_median"],
        train_fit_gc=c2["train_r2_gc_median"],
        train_fit_gc_frame=c2["train_r2_gc_frame_median"],
        atlas=c2["atlas_median"],
        gc_baseline=c2["gc_median"],
        gc_plus_frame=c2["gc_plus_frame_median"],
        gc_plus_matched=max(c2["gc_plus_dinucleotide_median"], c2["gc_plus_kmer4_median"]),
        atlas_in_window=c3["atlas_in_window_median"],
        gc_in_window=c3["gc_in_window_median"],
        atlas_out_window=c3["atlas_out_window_median"],
        gc_out_window=c3["gc_out_window_median"],
        random_stratum_min_margin=c3["random_stratum_min_margin_median"],
        margin=art["inherited_margin"],
        # READ FROM THE ARTIFACT, never a literal: this test's whole claim is reproducibility
        # from the artifact's own numbers, so a hand-typed value here would make it vacuous.
        best_matched_alone=art["c2_matched_strength_series"]["matched_control_premise"][
            "best_matched_alone"],
    )
    # THE RECORDED VERDICT NO LONGER REPRODUCES, AND THAT IS THE FINDING -- not a test to relax.
    #
    # `art["verdict"]` is FRAME_PARTIAL, which is what the EIGHT-outcome check order returned on
    # 2026-11-05 and is preserved as the historical record. The current NINE-outcome order adds
    # INDETERMINATE_CONTROL_PREMISE_UNMET, which reads a boolean the runner had computed and only
    # PRINTED since the beginning: control_is_at_least_as_strong_as_the_frame. On this run it is
    # False -- best_matched_alone 0.2895 against frame_alone 0.3170, a shortfall of 0.0275 or 8.7%
    # of the frame's own signal -- so the matched control was never matched, and the substantive
    # verdict rested on an unmet premise.
    #
    # The predicate is the RUNNER'S OWN (`max(dinuc, kmer4) >= frame_alone`); no tolerance is
    # invented here, which is why an 8.7% shortfall demotes rather than being waved through. If a
    # tolerance is wanted, that is a user call and the shortfall is stated above so it can be made.
    assert art["verdict"] == "FRAME_PARTIAL"
    assert grade(fi).verdict == "INDETERMINATE_CONTROL_PREMISE_UNMET"
    assert not grade(fi).is_valid


def test_the_margin_matches_the_frozen_bar(art):
    bar = json.loads(BAR.read_text(encoding="utf-8"))
    assert art["inherited_margin"] == bar["inherited_margin"]
    assert art["primary_frame_granularity"] == bar["primary_frame_granularity"]


# ------------------------------------------------------- the two controls the review demanded


def test_the_matched_strength_control_reports_its_OWN_premise(art):
    """A control whose premise is unverified cannot support an attribution claim."""
    mc = art["c2_matched_strength_series"]["matched_control_premise"]
    for key in ("frame_alone_median", "dinucleotide_alone_median", "kmer4_alone_median",
                "best_matched_alone", "control_is_at_least_as_strong_as_the_frame"):
        assert key in mc, key
    assert isinstance(mc["control_is_at_least_as_strong_as_the_frame"], bool)


def test_the_random_stratification_null_is_reported_beside_the_real_strata(art):
    c3 = art["c3_stratified"]
    assert "random_stratum_min_margin_median" in c3
    assert "stability_random_stratum_min_margin" in c3


def test_each_stratum_reports_its_expression_spread_beside_its_margin(art):
    """Without the spread a reader cannot tell range restriction from absorption."""
    c3 = art["c3_stratified"]
    for key in ("expression_spread_in_window", "expression_spread_out_window"):
        s = c3[key]
        assert s["n"] > 0
        for f in ("median", "q1", "q3", "iqr", "sd"):
            assert f in s, (key, f)


def test_each_stratums_baseline_was_re_measured_on_that_stratum(art):
    """Pre-commitment (c): never compared against the full-cohort figure."""
    c3 = art["c3_stratified"]
    for key in ("gc_in_window_median", "gc_out_window_median"):
        assert key in c3
    assert c3["gc_in_window_median"] != c3["gc_out_window_median"], (
        "one GC number for both strata would mean the baseline was not re-measured per stratum"
    )


# ------------------------------------------------------------------------- stability honesty


def test_the_artifact_says_whether_this_run_RESOLVED_its_own_verdict(art):
    vs = art["verdict_stability"]
    assert isinstance(vs["resolved"], bool)
    assert "coin flip" in vs["reading"]
    assert "user decision" in vs["reading"]


def test_the_stability_block_is_present_for_every_statistic_the_verdict_turns_on(art):
    c3 = art["c3_stratified"]
    for key in ("stability_margin_out_window", "stability_margin_in_window"):
        st = c3[key]
        for f in ("median", "min", "max", "n_seeds_below_bar", "straddles_bar"):
            assert f in st, (key, f)


def test_both_readings_of_the_stratum_margin_are_reported(art):
    """A gap between them is itself evidence the quantity is unstable."""
    blk = art["c3_stratified"]["median_of_differences_vs_difference_of_medians"]
    assert "median_of_per_seed_out_window_margins" in blk
    assert "difference_of_out_window_medians" in blk


# ----------------------------------------------------------------------- no float dump in wiki


def test_the_wiki_artifact_carries_no_float_dump(art):
    """A 102,755-decimal dump broke two contract-number guards on 2026-11-01."""
    longest = 0

    def walk(o):
        nonlocal longest
        if isinstance(o, list):
            if o and all(isinstance(x, (int, float)) for x in o):
                longest = max(longest, len(o))
            for x in o:
                walk(x)
        elif isinstance(o, dict):
            for v in o.values():
                walk(v)

    walk(art)
    assert longest <= 20, "a %d-element numeric array reached the wiki artifact" % longest


def test_the_per_seed_detail_stayed_in_gitignored_data(art):
    detail = REPO / "data/promoteratlas/frame_conditioning_detail_2026-11-01.json"
    if detail.is_file():
        d = json.loads(detail.read_text(encoding="utf-8"))
        assert len(d["per_seed"]["atlas"]) == art["n_seeds"]


def test_the_honest_limits_came_from_the_frozen_bar_not_from_prose(art):
    bar = json.loads(BAR.read_text(encoding="utf-8"))
    assert art["honest_limits"] == bar["honest_limits_fixed_in_advance"]
