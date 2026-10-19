"""Pure-part tests for the mechanism-decomposition arm.

The fit itself needs a ~3M-row sweep output on D:, so what is pinned here is everything that can be
wrong WITHOUT the data: the feature vector's shape and ordering, the paired statistics, and — most
importantly — that each verdict is REACHABLE. A verdict ladder nothing can reach is theatre, and this
repo has shipped one before (a filter that removed nothing, reported as an active control).
"""
from __future__ import annotations

import importlib.util
import math
from pathlib import Path

import pytest

SPEC = importlib.util.spec_from_file_location(
    "promcalc_decomposition_arm",
    Path(__file__).resolve().parent.parent / "scripts" / "promcalc_decomposition_arm.py",
)
M = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)


def _site(tx: float = 100.0, **kw) -> dict:
    d = {k: 0.0 for k in M.DG_FIELDS}
    d.update({"Tx_rate": tx, "spacer_len": 15, "hex35_start": 25, "hex10_start": 46})
    d.update(kw)
    return d


# --------------------------------------------------------------------------------------------------
# features
# --------------------------------------------------------------------------------------------------
def test_feature_vector_length_matches_its_names():
    """The names are the only documentation of a positional vector; a drift makes every column lie."""
    assert len(M.tile_features([_site()])) == len(M.feature_names())


def test_an_empty_tile_yields_the_same_length_and_flags_itself():
    """A tile with no candidate promoter is a real observation, not an error. It must stay
    distinguishable from a promoter whose energies happen to be zero -- hence the flag."""
    empty = M.tile_features([])
    assert len(empty) == len(M.feature_names())
    names = M.feature_names()
    assert empty[names.index("no_sites_flag")] == 1.0
    assert M.tile_features([_site()])[names.index("no_sites_flag")] == 0.0


def test_the_winning_site_is_the_one_with_the_HIGHEST_rate():
    """`best_*` must track max Tx_rate. Picking the first or last site would still produce a
    plausible number on real data while measuring something else entirely."""
    names = M.feature_names()
    sites = [_site(tx=1.0, dG_10=-9.9), _site(tx=500.0, dG_10=-1.1), _site(tx=50.0, dG_10=-5.5)]
    f = M.tile_features(sites)
    assert f[names.index("best_dG_10")] == pytest.approx(-1.1)
    assert f[names.index("best_log1p_tx")] == pytest.approx(math.log1p(500.0))


def test_n_sites_is_recorded_and_distribution_stats_span_all_sites():
    names = M.feature_names()
    f = M.tile_features([_site(tx=1.0, dG_total=-3.0), _site(tx=2.0, dG_total=1.0)])
    assert f[names.index("n_sites")] == 2.0
    assert f[names.index("dG_total_min")] == pytest.approx(-3.0)
    assert f[names.index("dG_total_max")] == pytest.approx(1.0)
    assert f[names.index("dG_total_mean")] == pytest.approx(-1.0)


def test_top_k_rates_are_sorted_descending_and_pad_when_short():
    names = M.feature_names()
    f = M.tile_features([_site(tx=5.0), _site(tx=50.0)])
    assert f[names.index("log1p_tx1")] == pytest.approx(math.log1p(50.0))
    assert f[names.index("log1p_tx2")] == pytest.approx(math.log1p(5.0))
    assert f[names.index("log1p_tx3")] == 0.0


# --------------------------------------------------------------------------------------------------
# paired statistics
# --------------------------------------------------------------------------------------------------
def test_paired_reports_a_ci_that_can_include_zero():
    """The CI is the point of this arm's reporting: two learned arms on this substrate had positive
    medians whose CIs included zero, and that is why they did not count."""
    noisy = M.paired([0.31, 0.29, 0.33, 0.28, 0.34], [0.30, 0.30, 0.30, 0.30, 0.30])
    assert noisy["includes_zero"] is True
    clean = M.paired([0.40] * 5, [0.30] * 5)
    assert clean["includes_zero"] is False
    assert clean["median"] == pytest.approx(0.10)
    assert clean["wins"] == 5


# --------------------------------------------------------------------------------------------------
# the verdict ladder -- every rung reachable
# --------------------------------------------------------------------------------------------------
def _arms(decomp: float, max_tx: float = 0.3968, gc0: float = 0.3000, ctrl: float = 0.0) -> dict:
    n = 10
    a = {
        "gc": {"median": 0.3112, "per_seed": [gc0] + [0.31] * (n - 1)},
        "max_tx": {"median": max_tx, "per_seed": [max_tx] * n},
        "decomposition": {"median": decomp, "per_seed": [decomp] * n},
        "control_feature_permuted": {"median": ctrl, "per_seed": [ctrl] * n},
    }
    a["_paired_vs_max_tx"] = M.paired(a["decomposition"]["per_seed"], a["max_tx"]["per_seed"])
    return a


def test_verdict_adds_when_the_margin_and_seed_count_are_met():
    v, why = M.verdict(_arms(0.4268))
    assert v == "DECOMPOSITION_ADDS", why


def test_verdict_weak_below_the_derived_margin():
    v, why = M.verdict(_arms(0.4068))  # +0.010 -- inside the band that swallowed the learned arms
    assert v == "WEAK", why


def test_verdict_no_gain_when_the_decomposition_does_not_beat_the_scalar():
    v, why = M.verdict(_arms(0.3900))
    assert v == "NO_GAIN", why


def test_reconcile_failure_on_the_gc_anchor_blocks_any_grade():
    v, why = M.verdict(_arms(0.4268, gc0=0.2000))
    assert v == "INDETERMINATE_RECONCILE_FAILED"
    assert "not reading the same data" in why


def test_reconcile_failure_when_the_SCALAR_does_not_reproduce():
    """Without this gate a feature-set comparison would silently confound features with pipeline."""
    v, why = M.verdict(_arms(0.4268, max_tx=0.3000))
    assert v == "INDETERMINATE_RECONCILE_FAILED"
    assert "confound" in why


def test_a_scoring_control_blocks_any_grade():
    v, why = M.verdict(_arms(0.4268, ctrl=0.30))
    assert v == "INDETERMINATE_CONTROL_NOT_CLEAN"
    assert "fabricates signal" in why


# --------------------------------------------------------------------------------------------------
# the pre-registration itself
# --------------------------------------------------------------------------------------------------
def test_the_control_is_feature_permutation_NOT_shuffled_labels():
    """Load-bearing. A shuffled-LABEL control cannot collapse here: it randomises the fitted
    coefficient's sign but not its ranking power, so it returns +-0.4 by construction. That error was
    made and caught on this exact substrate on 2026-10-12; pinning the choice stops its return."""
    assert M.PREREGISTERED["control"] == "feature_permutation"
    assert "shuffled" not in M.PREREGISTERED["control"]


def test_the_margin_is_documented_as_derived_from_measured_effects():
    """A pre-registered threshold nobody analysed is an assertion with a timestamp."""
    assert M.PREREGISTERED["margin"] == 0.02
    d = M.PREREGISTERED["margin_is"]
    assert "DERIVED" in d and "0.0686" in d


def test_both_reconcile_anchors_name_a_source():
    for k, v in M.PREREGISTERED["reconcile"].items():
        assert v["source"], f"{k} has no provenance"
        assert v["expected"] and v["tolerance"]
