"""Tests for the positional arm's scoring runner.

Every test drives the REAL `run()` control flow with a synthetic cohort, a synthetic stored block
and a stub scorer, so each `INDETERMINATE_*` branch is proven REACHABLE rather than assumed to be.
A gate that cannot be reached is a gate that cannot protect anything, and this project has a
recorded instance of a validity branch that could never fire.

The non-vacuity test that matters most perturbs the stored block and asserts the parity gate trips:
without it, a parity gate that always passed would look identical from the artifact.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from dna_decode.glm.genomewide import Tile
from scripts.freeze_positional_arm_bar import CO_PRIMARIES, MARGIN, OUTCOMES, SECONDARIES
from scripts.positional_arm_extract import as_bins, feature_layout, rebin
from scripts.positional_arm_score import (
    SPEARMAN_PARITY_TOL,
    ScoreError,
    build_rung,
    leave_peak_out_indices,
    median_of,
    order_statistic_median_ci,
    paired_bootstrap,
    run,
    verdict_from_artifact,
)

CHANNELS = 4
N_BINS = 10
N_LAGS = 3
N_TILES = 200
N_SEEDS = 3
FRAC = 0.25
LAYOUT = feature_layout(N_LAGS, channels=CHANNELS, n_bins=N_BINS)
TOTAL = LAYOUT["_total"]          # 4*10 + 4*10 + 4*3 = 92


# --------------------------------------------------------------------------------------------------
# synthetic fixtures
# --------------------------------------------------------------------------------------------------
def make_cohort(n: int = N_TILES, n_peaks: int = 20):
    rng = np.random.default_rng(5)
    tiles = []
    for i in range(n):
        tiles.append(Tile(
            seq="ACGT" * 10, expression=float(abs(rng.normal()) + 0.1),
            peak="peak_%02d" % (i % n_peaks), strand="+" if i % 2 else "-",
            category="tile", active="yes",
        ))

    class St:
        path = "synthetic"
    return tiles, St()


def make_block(n: int = N_TILES, *, perturb: tuple | None = None):
    rng = np.random.default_rng(11)
    arr = rng.normal(size=(n, TOTAL)).astype(np.float32)
    if perturb is not None:
        row, col, delta = perturb
        arr[row, col] += delta
    side = {
        "block": "synthetic", "status": "COMPLETE", "checkpoint_sha256": "stub",
        "geometry": {"n_rows": n, "n_features": TOTAL, "dtype": "float32",
                     "n_bins": N_BINS, "lags": list(range(1, N_LAGS + 1)),
                     "channels": CHANNELS, "sequence_length": 150},
    }
    return arr, side, LAYOUT


def make_bar(tiles, *, atlas=0.50, n_test=None, parity_tol=1e-5, headroom=1.5):
    if n_test is None:
        n_test = {}
        for s in range(N_SEEDS):
            _tr, te, _h = leave_peak_out_indices(tiles, frac=FRAC, seed=s)
            n_test[str(s)] = len(te)
    per_seed = {str(s): round(atlas + 0.001 * s, 4) for s in range(N_SEEDS)}
    return {
        "arm": "positional-non-mean-pooled-promoteratlas-150bp",
        "criterion": "ceiling: synthetic",
        "inherited_not_chosen": {
            "margin": {"value": MARGIN},
            "comparator": {"atlas_per_seed": per_seed, "atlas_median": median_of(per_seed.values()),
                           "source": "synthetic"},
            "split": {"frac": FRAC, "n_seeds": N_SEEDS, "n_test_per_seed": n_test},
        },
        "measured_and_inherited_from_the_probe": {
            "parity_tolerance": parity_tol, "ram_headroom_factor": headroom,
        },
        "multiplicity": {"co_primary_ci_percentiles": [1.25, 98.75],
                         "secondary_ci_percentiles": [0.41666667, 99.58333333]},
        "honest_limits_fixed_in_advance": ["synthetic"],
        "scoping_corrections_fixed_in_advance": ["synthetic"],
        "forecast": {"predicted_outcome": "POSITIONAL_NO_GAIN"},
    }


def stub_scorer(*, parity_by_seed, null_value=0.0, by_width=None, calls=None):
    """Keyed on what `run` can actually vary: the shuffle flag and the matrix width.

    width 4 -> the parity rung (global mean), 40 -> the three bins10 rungs, 12 -> autocorr,
    8 -> bins2, 20 -> bins5.
    """
    widths = {4: None, 40: 0.60, 12: 0.60, 8: 0.55, 20: 0.57}
    widths.update(by_width or {})
    seen = {"i": 0}

    def score(Xtr, ytr, Xte, yte, shuffle_seed=None):
        w = int(np.asarray(Xtr).shape[1])
        preds = [0.0] * len(Xte)
        if calls is not None:
            calls.append({"width": w, "shuffled": shuffle_seed is not None})
        if shuffle_seed is not None:
            return float(null_value), preds
        if w == 4:
            v = parity_by_seed[seen["i"] % len(parity_by_seed)]
            seen["i"] += 1
            return float(v), preds
        if w not in widths or widths[w] is None:
            raise AssertionError("stub scorer got an unexpected width %d" % w)
        return float(widths[w]), preds

    return score


def go(tiles, st, block, bar, *, scorer, null_extension=0, emb_cache=None):
    return run(block_dir=Path("unused"), bar=bar, null_extension=null_extension, n_boot=200,
               scorer=scorer, cohort=(tiles, st), block=block, emb_cache=emb_cache, verbose=False)


# --------------------------------------------------------------------------------------------------
# rung construction
# --------------------------------------------------------------------------------------------------
def test_every_declared_rung_builds_with_the_declared_width():
    tiles, _st = make_cohort()
    arr, _side, lay = make_block()
    strands = [t.strand for t in tiles]
    want = {
        "global_mean": CHANNELS,
        "bins10_mean_oriented": CHANNELS * 10, "bins10_mean_asgiven": CHANNELS * 10,
        "bins10_max_oriented": CHANNELS * 10,
        "bins2_mean_oriented": CHANNELS * 2, "bins5_mean_oriented": CHANNELS * 5,
        "autocorr": CHANNELS * N_LAGS,
    }
    for name, w in want.items():
        X = build_rung(arr, lay, name, strands)
        assert X.shape == (N_TILES, w), name
        assert np.all(np.isfinite(X)), name


def test_build_rung_refuses_an_unknown_rung():
    tiles, _ = make_cohort()
    arr, _s, lay = make_block()
    with pytest.raises(ScoreError, match="unknown rung"):
        build_rung(arr, lay, "bins7_mean_oriented", [t.strand for t in tiles])


def test_oriented_and_asgiven_rungs_differ_on_a_mixed_strand_cohort():
    """If they were equal the orientation diagnostic would be measuring nothing."""
    tiles, _ = make_cohort()
    arr, _s, lay = make_block()
    strands = [t.strand for t in tiles]
    a = build_rung(arr, lay, "bins10_mean_oriented", strands)
    b = build_rung(arr, lay, "bins10_mean_asgiven", strands)
    assert not np.allclose(a, b)
    # and on an all-plus cohort they must coincide exactly
    allplus = build_rung(arr, lay, "bins10_mean_oriented", ["+"] * len(tiles))
    assert np.allclose(allplus, b)


def test_rebin_and_orient_commute_which_is_why_the_cheap_order_is_used():
    """The mean within a group is order-invariant, so grouping then reversing equals the reverse."""
    from dna_decode.glm.positional_readout import flatten_bins, orient

    arr, _s, lay = make_block()
    strands = ["+" if i % 3 else "-" for i in range(N_TILES)]
    bm = as_bins(arr[:, 0:CHANNELS * 10], lay)
    rebin_then_orient = flatten_bins(orient(rebin(bm, 5), strands))
    orient_then_rebin = flatten_bins(rebin(orient(bm, strands), 5))
    assert np.allclose(rebin_then_orient, orient_then_rebin, atol=1e-12)


# --------------------------------------------------------------------------------------------------
# statistics helpers
# --------------------------------------------------------------------------------------------------
def test_median_of_handles_even_and_odd_lengths():
    assert median_of([1, 2, 3]) == 2.0
    assert median_of([1, 2, 3, 4]) == 2.5
    with pytest.raises(ScoreError, match="empty"):
        median_of([])


def test_paired_bootstrap_brackets_a_clearly_positive_gain():
    b = paired_bootstrap([0.08] * 10, percentiles=[1.25, 98.75], n_boot=500)
    assert b["ci_lo"] == pytest.approx(0.08, abs=1e-9)
    assert b["mean"] == pytest.approx(0.08)
    assert b["n_seeds"] == 10


def test_paired_bootstrap_is_deterministic_and_wider_at_a_stricter_level():
    g = list(np.random.default_rng(3).normal(0.05, 0.03, size=10))
    a = paired_bootstrap(g, percentiles=[1.25, 98.75], n_boot=2000)
    b = paired_bootstrap(g, percentiles=[1.25, 98.75], n_boot=2000)
    assert a == b
    strict = paired_bootstrap(g, percentiles=[0.41666667, 99.58333333], n_boot=2000)
    assert (strict["ci_hi"] - strict["ci_lo"]) >= (a["ci_hi"] - a["ci_lo"])


def test_paired_bootstrap_refuses_empty_or_non_finite_gains():
    with pytest.raises(ScoreError, match="no gains"):
        paired_bootstrap([], percentiles=[2.5, 97.5])
    with pytest.raises(ScoreError, match="non-finite"):
        paired_bootstrap([0.1, float("nan")], percentiles=[2.5, 97.5])


def test_order_statistic_median_ci_reports_its_achieved_coverage():
    lo, hi, cov = order_statistic_median_ci(list(range(21)))
    assert lo is not None and hi is not None
    assert cov >= 0.95
    assert lo <= 10 <= hi
    assert order_statistic_median_ci([1.0]) == (None, None, None)


# --------------------------------------------------------------------------------------------------
# GATE 1: split equality -- reachable, and checked first
# --------------------------------------------------------------------------------------------------
def test_split_mismatch_is_reachable_and_stops_the_run():
    tiles, st = make_cohort()
    bar = make_bar(tiles, n_test={str(s): 1 for s in range(N_SEEDS)})
    rec = go(tiles, st, make_block(), bar,
             scorer=stub_scorer(parity_by_seed=[0.5] * N_SEEDS))
    assert rec["verdict"] == "INDETERMINATE_SPLIT_MISMATCH"
    assert rec["gates"]["split_equality"]["passed"] is False
    assert rec["gates"]["split_equality"]["mismatched_seeds"] == list(range(N_SEEDS))
    assert "parity" not in rec["gates"], "the run must stop at the first validity failure"
    assert rec["rungs"] == {}


def test_split_equality_passes_on_the_recomputed_committed_sizes():
    tiles, st = make_cohort()
    bar = make_bar(tiles)
    rec = go(tiles, st, make_block(), bar,
             scorer=stub_scorer(parity_by_seed=[bar["inherited_not_chosen"]["comparator"]
                                                ["atlas_per_seed"][str(s)]
                                                for s in range(N_SEEDS)]))
    assert rec["gates"]["split_equality"]["passed"] is True
    assert all(v["match"] for v in rec["gates"]["split_equality"]["per_seed"].values())


def test_splits_are_reused_byte_identically_across_rungs():
    """A rung scored on a different split is not comparable to the incumbent."""
    tiles, _ = make_cohort()
    a = [leave_peak_out_indices(tiles, frac=FRAC, seed=s) for s in range(N_SEEDS)]
    b = [leave_peak_out_indices(tiles, frac=FRAC, seed=s) for s in range(N_SEEDS)]
    assert a == b
    for (tr, te, _h) in a:
        assert set(tr).isdisjoint(te)
        assert len(tr) + len(te) == N_TILES


def test_the_split_holds_out_whole_peaks():
    tiles, _ = make_cohort()
    tr, te, _ = leave_peak_out_indices(tiles, frac=FRAC, seed=0)
    train_peaks = {tiles[i].peak for i in tr}
    test_peaks = {tiles[i].peak for i in te}
    assert train_peaks.isdisjoint(test_peaks), "peak disjointness is the leakage unit"


# --------------------------------------------------------------------------------------------------
# GATE 2: parity -- reachable, and NON-VACUOUS
# --------------------------------------------------------------------------------------------------
def test_parity_failure_is_reachable_and_stops_the_run():
    tiles, st = make_cohort()
    bar = make_bar(tiles)
    rec = go(tiles, st, make_block(), bar,
             scorer=stub_scorer(parity_by_seed=[0.90] * N_SEEDS))  # nowhere near committed 0.50
    assert rec["verdict"] == "INDETERMINATE_PARITY_FAIL"
    assert rec["gates"]["parity"]["passed"] is False
    assert rec["gates"]["parity"]["spearman_level"]["deviant_seeds"] == list(range(N_SEEDS))
    assert "shuffled_label_null" not in rec["gates"]
    assert rec["rungs"] == {}


def test_a_perturbed_stored_block_trips_the_embedding_level_parity_gate(tmp_path):
    """NON-VACUITY. Without this, an always-passing parity gate looks identical from the artifact."""
    tiles, st = make_cohort()
    bar = make_bar(tiles)
    committed = build_rung(make_block()[0], LAYOUT, "global_mean", [t.strand for t in tiles])
    cache = tmp_path / "emb.npy"
    np.save(cache, committed.astype(np.float32))

    atlas = [bar["inherited_not_chosen"]["comparator"]["atlas_per_seed"][str(s)]
             for s in range(N_SEEDS)]

    clean = go(tiles, st, make_block(), bar,
               scorer=stub_scorer(parity_by_seed=atlas), emb_cache=cache)
    assert clean["gates"]["parity"]["passed"] is True
    assert clean["gates"]["parity"]["embedding_level"]["passed"] is True

    dirty = go(tiles, st, make_block(perturb=(3, 7, 0.5)), bar,
               scorer=stub_scorer(parity_by_seed=atlas), emb_cache=cache)
    assert dirty["verdict"] == "INDETERMINATE_PARITY_FAIL"
    assert dirty["gates"]["parity"]["embedding_level"]["passed"] is False
    assert dirty["gates"]["parity"]["embedding_level"]["max_abs_deviation"] > 0


def test_the_spearman_parity_tolerance_is_a_boundary():
    tiles, st = make_cohort()
    bar = make_bar(tiles)
    atlas = [bar["inherited_not_chosen"]["comparator"]["atlas_per_seed"][str(s)]
             for s in range(N_SEEDS)]
    inside = [round(v + SPEARMAN_PARITY_TOL * 0.5, 4) for v in atlas]
    outside = [round(v + SPEARMAN_PARITY_TOL * 5, 4) for v in atlas]
    assert go(tiles, st, make_block(), bar,
              scorer=stub_scorer(parity_by_seed=inside))["gates"]["parity"]["passed"] is True
    assert go(tiles, st, make_block(), bar,
              scorer=stub_scorer(parity_by_seed=outside))["verdict"] == "INDETERMINATE_PARITY_FAIL"


# --------------------------------------------------------------------------------------------------
# GATE 3: the shuffled-label null -- reachable, and it extends rather than weakens
# --------------------------------------------------------------------------------------------------
def _passing_parity(bar):
    return [bar["inherited_not_chosen"]["comparator"]["atlas_per_seed"][str(s)]
            for s in range(N_SEEDS)]


def test_null_failure_is_reachable_and_stops_the_run():
    tiles, st = make_cohort()
    bar = make_bar(tiles)
    rec = go(tiles, st, make_block(), bar,
             scorer=stub_scorer(parity_by_seed=_passing_parity(bar), null_value=0.25),
             null_extension=0)
    assert rec["verdict"] == "INDETERMINATE_NULL_FAIL"
    assert rec["gates"]["shuffled_label_null"]["passed"] is False
    assert rec["gates"]["shuffled_label_null"]["pooled_median"] == pytest.approx(0.25)
    assert rec["rungs"] == {}
    assert "MORE SEEDS" in rec["gates"]["shuffled_label_null"]["detail"]


def test_a_failing_null_spends_more_seeds_rather_than_moving_the_bar():
    """The bar's pre-commitment, exercised: the extension changes the SEED COUNT, not the bar."""
    tiles, st = make_cohort()
    bar = make_bar(tiles)
    calls: list[dict] = []
    rec = go(tiles, st, make_block(), bar,
             scorer=stub_scorer(parity_by_seed=_passing_parity(bar), null_value=0.25, calls=calls),
             null_extension=40)
    g = rec["gates"]["shuffled_label_null"]
    assert g["extended"] is True
    assert g["n_extension_seeds"] == 40
    assert g["n_pooled"] == 40
    assert g["bar"] == "|median| < %.2f" % MARGIN, "the bar text must be unchanged by extending"
    assert sum(1 for c in calls if c["shuffled"]) == N_SEEDS + 40
    assert g["median_ci_achieved_coverage"] >= 0.95


def test_the_null_does_not_extend_when_it_already_passes():
    tiles, st = make_cohort()
    bar = make_bar(tiles)
    rec = go(tiles, st, make_block(), bar,
             scorer=stub_scorer(parity_by_seed=_passing_parity(bar), null_value=0.001),
             null_extension=40)
    g = rec["gates"]["shuffled_label_null"]
    assert g["passed"] is True
    assert g["extended"] is False
    assert g["n_extension_seeds"] == 0


def test_the_null_bar_is_two_sided():
    """A large NEGATIVE shuffled-label Spearman is just as much a broken pipeline."""
    tiles, st = make_cohort()
    bar = make_bar(tiles)
    rec = go(tiles, st, make_block(), bar,
             scorer=stub_scorer(parity_by_seed=_passing_parity(bar), null_value=-0.25),
             null_extension=0)
    assert rec["verdict"] == "INDETERMINATE_NULL_FAIL"


def test_the_null_statistic_uses_the_primary_rung_and_shuffles_only_train_labels():
    tiles, st = make_cohort()
    bar = make_bar(tiles)
    calls: list[dict] = []
    go(tiles, st, make_block(), bar,
       scorer=stub_scorer(parity_by_seed=_passing_parity(bar), null_value=0.0, calls=calls))
    shuffled = [c for c in calls if c["shuffled"]]
    assert shuffled, "the null must actually be computed"
    assert all(c["width"] == CHANNELS * 10 for c in shuffled), "null runs on the primary rung"


def test_the_artifact_records_the_null_prose_versus_implementation_discrepancy():
    tiles, st = make_cohort()
    bar = make_bar(tiles)
    rec = go(tiles, st, make_block(), bar,
             scorer=stub_scorer(parity_by_seed=_passing_parity(bar)))
    prov = rec["gates"]["shuffled_label_null"]["statistic_provenance"]
    assert "gain over GC" in prov and "no GC term" in prov
    assert "0.0203" in prov and "+0.0517" in prov


# --------------------------------------------------------------------------------------------------
# the ladder and the verdict
# --------------------------------------------------------------------------------------------------
# The committed atlas per-seed values here are 0.5000 / 0.5010 / 0.5020, median 0.5010.
# AMBIGUOUS needs the corrected interval to CONTAIN the margin, not merely to fall short of it:
# a rung at 0.51 has a gain of ~0.010, which is entirely BELOW +0.02, so a tight interval there is
# correctly NO_GAIN. A rung at 0.521 has gains 0.021/0.020/0.019, so its interval straddles +0.02.
@pytest.mark.parametrize("bins10,autoc,want", [
    (0.60, 0.40, "POSITIONAL_ADDS"),        # bins10 clears the margin outright
    (0.40, 0.60, "POSITIONAL_ADDS"),        # autocorr clears on its own
    (0.30, 0.30, "POSITIONAL_NO_GAIN"),     # both well below
    (0.51, 0.30, "POSITIONAL_NO_GAIN"),     # a gain of 0.010 is short of the margin, not ambiguous
    (0.521, 0.30, "POSITIONAL_AMBIGUOUS"),  # gains 0.019-0.021 straddle the margin
])
def test_all_three_terminal_verdicts_are_reachable(bins10, autoc, want):
    tiles, st = make_cohort()
    bar = make_bar(tiles)
    rec = go(tiles, st, make_block(), bar,
             scorer=stub_scorer(parity_by_seed=_passing_parity(bar),
                                by_width={40: bins10, 12: autoc}))
    assert rec["verdict"] == want
    assert rec["verdict"] in OUTCOMES


def test_the_union_of_reachable_verdicts_equals_the_bars_closed_set():
    """Together with the three gate tests above, every one of the six outcomes is reachable."""
    tiles, st = make_cohort()
    bar = make_bar(tiles)
    seen = set()
    seen.add(go(tiles, st, make_block(), make_bar(tiles, n_test={str(s): 1 for s in range(N_SEEDS)}),
                scorer=stub_scorer(parity_by_seed=[0.5] * N_SEEDS))["verdict"])
    seen.add(go(tiles, st, make_block(), bar,
                scorer=stub_scorer(parity_by_seed=[0.9] * N_SEEDS))["verdict"])
    seen.add(go(tiles, st, make_block(), bar,
                scorer=stub_scorer(parity_by_seed=_passing_parity(bar), null_value=0.25),
                null_extension=0)["verdict"])
    for b10, ac in ((0.60, 0.40), (0.30, 0.30), (0.521, 0.30)):
        seen.add(go(tiles, st, make_block(), bar,
                    scorer=stub_scorer(parity_by_seed=_passing_parity(bar),
                                       by_width={40: b10, 12: ac}))["verdict"])
    assert seen == set(OUTCOMES)


def test_secondaries_are_scored_but_cannot_change_the_verdict():
    tiles, st = make_cohort()
    bar = make_bar(tiles)
    rec = go(tiles, st, make_block(), bar,
             scorer=stub_scorer(parity_by_seed=_passing_parity(bar),
                                by_width={40: 0.30, 12: 0.30, 8: 0.99, 20: 0.99}))
    assert rec["verdict"] == "POSITIONAL_NO_GAIN"
    for s in ("bins2_mean_oriented", "bins5_mean_oriented"):
        assert rec["rungs"][s]["median_spearman"] == pytest.approx(0.99)
        assert rec["rungs"][s]["eligible_as_the_arms_result"] is False
    for c in CO_PRIMARIES:
        assert rec["rungs"][c]["eligible_as_the_arms_result"] is True


def test_every_rung_reports_per_seed_values_wins_and_a_paired_ci():
    tiles, st = make_cohort()
    bar = make_bar(tiles)
    rec = go(tiles, st, make_block(), bar,
             scorer=stub_scorer(parity_by_seed=_passing_parity(bar)))
    for name in list(CO_PRIMARIES) + list(SECONDARIES):
        r = rec["rungs"][name]
        assert len(r["per_seed"]) == N_SEEDS
        assert 0 <= r["seed_wins_over_atlas"] <= N_SEEDS
        assert set(r["gain_bootstrap"]) >= {"ci_lo", "ci_hi", "median", "mean", "ci_percentiles"}
        assert r["n_features"] > 0


def test_co_primaries_and_secondaries_get_their_declared_interval_levels():
    tiles, st = make_cohort()
    bar = make_bar(tiles)
    rec = go(tiles, st, make_block(), bar,
             scorer=stub_scorer(parity_by_seed=_passing_parity(bar)))
    for c in CO_PRIMARIES:
        assert rec["rungs"][c]["gain_bootstrap"]["ci_percentiles"] == [1.25, 98.75]
    for s in SECONDARIES:
        assert rec["rungs"][s]["gain_bootstrap"]["ci_percentiles"] == [0.41666667, 99.58333333]


# --------------------------------------------------------------------------------------------------
# the RAM pre-check: SKIPPED and recorded, never silently dropped
# --------------------------------------------------------------------------------------------------
def test_a_rung_over_the_ram_budget_is_recorded_as_skipped_rather_than_attempted():
    tiles, st = make_cohort()
    bar = make_bar(tiles, headroom=1e12)       # make every rung unaffordable
    calls: list[dict] = []
    rec = go(tiles, st, make_block(), bar,
             scorer=stub_scorer(parity_by_seed=_passing_parity(bar), calls=calls))
    for name in list(CO_PRIMARIES) + list(SECONDARIES):
        r = rec["rungs"][name]
        assert r["skipped"] is True
        assert r["peak_design_bytes"] > 0
        assert r["fits_in_ram"] is not True
        assert "silently dropped" in r["why"]
    assert not any(c["width"] in (40, 12, 8, 20) and not c["shuffled"] for c in calls), \
        "a skipped rung must not have been scored"


def test_a_skipped_co_primary_yields_ambiguous_not_a_gain_claim():
    """A rung nobody could afford to score is not evidence of no gain."""
    tiles, st = make_cohort()
    bar = make_bar(tiles, headroom=1e12)
    rec = go(tiles, st, make_block(), bar,
             scorer=stub_scorer(parity_by_seed=_passing_parity(bar)))
    assert rec["verdict"] == "POSITIONAL_AMBIGUOUS"
    assert "not scored" in rec["verdict_detail"]


# --------------------------------------------------------------------------------------------------
# the verdict is a pure function of the emitted artifact plus the frozen anchors
# --------------------------------------------------------------------------------------------------
@pytest.mark.parametrize("bins10,autoc", [(0.60, 0.40), (0.30, 0.30), (0.521, 0.30)])
def test_the_verdict_is_rederivable_from_the_emitted_artifact_alone(bins10, autoc):
    tiles, st = make_cohort()
    bar = make_bar(tiles)
    rec = go(tiles, st, make_block(), bar,
             scorer=stub_scorer(parity_by_seed=_passing_parity(bar),
                                by_width={40: bins10, 12: autoc}))
    # round-trip through JSON so only what was WRITTEN is available to the re-derivation
    again, _why = verdict_from_artifact(json.loads(json.dumps(rec)), bar)
    assert again == rec["verdict"]
    assert rec["verdict_is_rederivable_from_this_artifact"] is True


def test_the_rederivation_reproduces_every_indeterminate_branch_too():
    tiles, st = make_cohort()
    bar = make_bar(tiles)
    cases = [
        go(tiles, st, make_block(), make_bar(tiles, n_test={str(s): 1 for s in range(N_SEEDS)}),
           scorer=stub_scorer(parity_by_seed=[0.5] * N_SEEDS)),
        go(tiles, st, make_block(), bar, scorer=stub_scorer(parity_by_seed=[0.9] * N_SEEDS)),
        go(tiles, st, make_block(), bar,
           scorer=stub_scorer(parity_by_seed=_passing_parity(bar), null_value=0.25),
           null_extension=0),
    ]
    for rec in cases:
        again, _ = verdict_from_artifact(json.loads(json.dumps(rec)), bar)
        assert again == rec["verdict"]


def test_the_rederivation_treats_a_never_evaluated_gate_as_that_gates_failure():
    """An artifact truncated before a gate must not re-derive as a terminal verdict."""
    bar = make_bar(make_cohort()[0])
    assert verdict_from_artifact({"gates": {}}, bar)[0] == "INDETERMINATE_SPLIT_MISMATCH"
    assert verdict_from_artifact(
        {"gates": {"split_equality": {"passed": True}}}, bar)[0] == "INDETERMINATE_PARITY_FAIL"
    assert verdict_from_artifact(
        {"gates": {"split_equality": {"passed": True}, "parity": {"passed": True}}},
        bar)[0] == "INDETERMINATE_NULL_FAIL"


# --------------------------------------------------------------------------------------------------
# bookkeeping the artifact must carry
# --------------------------------------------------------------------------------------------------
def test_the_artifact_records_the_comparator_as_not_recomputed():
    tiles, st = make_cohort()
    bar = make_bar(tiles)
    rec = go(tiles, st, make_block(), bar,
             scorer=stub_scorer(parity_by_seed=_passing_parity(bar)))
    assert rec["comparator"]["recomputed_in_run"] is False
    assert rec["comparator"]["atlas_per_seed"] == bar["inherited_not_chosen"]["comparator"]["atlas_per_seed"]
    assert rec["margin"] == MARGIN


def test_the_artifact_carries_the_scoping_corrections_and_honest_limits_forward():
    tiles, st = make_cohort()
    bar = make_bar(tiles)
    rec = go(tiles, st, make_block(), bar,
             scorer=stub_scorer(parity_by_seed=_passing_parity(bar)))
    assert rec["scoping_corrections"] == bar["scoping_corrections_fixed_in_advance"]
    assert rec["honest_limits"] == bar["honest_limits_fixed_in_advance"]
    assert rec["forecast_recorded_in_bar"] == bar["forecast"]


def test_run_refuses_a_block_whose_row_count_differs_from_the_cohort():
    tiles, st = make_cohort()
    bar = make_bar(tiles)
    with pytest.raises(ScoreError, match="pairing is impossible"):
        go(tiles, st, make_block(n=N_TILES - 1), bar,
           scorer=stub_scorer(parity_by_seed=_passing_parity(bar)))


# --------------------------------------------------------------------------------------------------
# main exit codes
# --------------------------------------------------------------------------------------------------
def test_main_returns_2_when_the_bar_is_absent(monkeypatch, tmp_path):
    import scripts.positional_arm_score as sc

    monkeypatch.setattr(sc, "BAR_PATH", tmp_path / "no_bar.json")
    assert sc.main([]) == 2


def test_main_returns_1_when_scoring_raises(monkeypatch):
    import scripts.positional_arm_score as sc

    monkeypatch.setattr(sc, "run", lambda **_: (_ for _ in ()).throw(ScoreError("nope")))
    assert sc.main([]) == 1


# --------------------------------------------------------------------------------------------------
# the decomposition: the gate ORDER is data, and it is SOURCED from the frozen bar
# --------------------------------------------------------------------------------------------------
def test_gate_sequence_is_derived_from_the_frozen_bars_check_order():
    """The order must not be restated anywhere. It is the bar's, and this proves the derivation."""
    from scripts.freeze_positional_arm_bar import CHECK_ORDER
    from scripts.positional_arm_score import GATE_SEQUENCE

    want = tuple((g, o) for g, o in CHECK_ORDER if isinstance(o, str))
    assert GATE_SEQUENCE == want
    assert [g for g, _ in GATE_SEQUENCE] == ["split_equality", "parity", "shuffled_label_null"]


def test_gate_sequence_matches_the_frozen_bar_artifact_on_disk():
    """A reordered refactor must fail rather than silently score before a gate."""
    from scripts.freeze_positional_arm_bar import OUT as BAR_OUT
    from scripts.positional_arm_score import GATE_SEQUENCE

    if not BAR_OUT.is_file():
        pytest.skip("bar not frozen")
    declared = json.loads(BAR_OUT.read_text(encoding="utf-8"))["verdict_check_order"]
    gates = [g for g, _ in GATE_SEQUENCE]
    assert declared[:len(gates)] == gates, (
        "the runner's gate order and the frozen bar's verdict_check_order disagree")
    assert declared[-1] == "co_primary_decision", "the decision must come last"


def test_every_gate_in_the_sequence_has_an_implementation_and_no_orphans():
    from scripts.positional_arm_score import GATE_FUNCS, GATE_SEQUENCE

    assert set(GATE_FUNCS) == {g for g, _ in GATE_SEQUENCE}
    assert all(callable(f) for f in GATE_FUNCS.values())


def test_the_gates_are_invoked_in_the_declared_order_and_stop_at_the_first_failure(monkeypatch):
    """Observes the real dispatch rather than trusting the tuple."""
    import scripts.positional_arm_score as sc

    tiles, st = make_cohort()
    bar = make_bar(tiles)
    called: list[str] = []
    wrapped = {}
    for name, fn in sc.GATE_FUNCS.items():
        def mk(n, f):
            def inner(ctx):
                called.append(n)
                return f(ctx)
            return inner
        wrapped[name] = mk(name, fn)
    monkeypatch.setattr(sc, "GATE_FUNCS", wrapped)

    go(tiles, st, make_block(), bar, scorer=stub_scorer(parity_by_seed=_passing_parity(bar)))
    assert called == ["split_equality", "parity", "shuffled_label_null"]

    # a split mismatch must stop the chain at the FIRST gate
    called.clear()
    go(tiles, st, make_block(), make_bar(tiles, n_test={str(s): 1 for s in range(N_SEEDS)}),
       scorer=stub_scorer(parity_by_seed=[0.5] * N_SEEDS))
    assert called == ["split_equality"], "a failed gate must not let later gates run"


def test_each_extracted_gate_returns_a_dict_with_passed_and_detail():
    """The contract `run`'s loop depends on: it reads only `passed` and `detail`."""
    import scripts.positional_arm_score as sc

    tiles, st = make_cohort()
    bar = make_bar(tiles)
    rec = go(tiles, st, make_block(), bar,
             scorer=stub_scorer(parity_by_seed=_passing_parity(bar)))
    for name in sc.GATE_FUNCS:
        g = rec["gates"][name]
        assert isinstance(g["passed"], bool)
        assert isinstance(g["detail"], str) and g["detail"].strip()


def test_the_primary_rung_is_built_exactly_once_across_the_null_gate_and_the_ladder():
    """Before the extraction the null gate and the ladder shared one local. `ScoreCtx.primary_matrix`
    preserves that, and this test is what stops the decomposition quietly doubling the work."""
    import scripts.positional_arm_score as sc

    tiles, st = make_cohort()
    bar = make_bar(tiles)
    built: list[str] = []
    real = sc.build_rung

    def counting(block, layout, name, strands):
        built.append(name)
        return real(block, layout, name, strands)

    import unittest.mock as mock
    with mock.patch.object(sc, "build_rung", counting):
        go(tiles, st, make_block(), bar,
           scorer=stub_scorer(parity_by_seed=_passing_parity(bar)))
    primary = CO_PRIMARIES[0]
    assert built.count(primary) == 1, "the primary rung was built %d times" % built.count(primary)
    assert built.count("global_mean") == 1


def test_run_stays_decomposed_and_below_the_complexity_ceiling():
    """Guards the refactor itself. `run` measured cyclomatic ~58 over 304 lines before Step 5."""
    import ast

    src = (Path(__file__).resolve().parents[1] / "scripts/positional_arm_score.py").read_text(
        encoding="utf-8")
    tree = ast.parse(src)
    found = {}
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef):
            c = 1 + sum(
                isinstance(k, (ast.If, ast.For, ast.While, ast.ExceptHandler, ast.BoolOp,
                               ast.IfExp, ast.Assert))
                for k in ast.walk(node))
            found[node.name] = c
    assert found["run"] < 25, "run is cyclomatic ~%d; the gate chain has been re-inlined" % found["run"]
    for gate in ("gate_split_equality", "gate_parity", "gate_shuffled_label_null", "score_ladder"):
        assert gate in found, "%s was removed; the decomposition is gone" % gate
        assert found[gate] < 25, "%s is cyclomatic ~%d" % (gate, found[gate])
