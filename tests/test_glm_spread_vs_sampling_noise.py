"""The across-seed spread is sampling noise, and this pins BOTH the benchmark and the retraction.

This probe exists because a seed-role conclusion was drawn twice from the same data and was wrong twice.
A held-out AUROC between `n1` and `n2` items has a known null sampling sd -- the Mann-Whitney variance,
`sqrt((n1 + n2 + 1) / (12 * n1 * n2))`. At the committed geometry that is 0.0409. The committed grid's ten
nulls have an observed across-seed sd of 0.0410, and at `n_generated=400` the formula gives 0.0264 against
an observed 0.0266. Two agreements to three decimals across a geometry change that moved the value by a
third is what a score dominated by test-set sampling noise looks like, and it is not something a seed-role
mechanism predicts.

**The retraction this pins.** A previous run reported the stream-decoupled arm A2 at an sd ratio of 1.515
with a CI of [1.191, 1.912] in `both` mode and read it as stream sharing suppressing the spread. Two things
undo it: arm A was measured at 60 replicates and came in low (0.0352 against 0.0448 at 240 replicates), and
A2's inflation appears at only ONE of five offset pairs -- the pair the earlier run happened to use. With
arm A measured properly it is consistent with the pure sampling benchmark in BOTH modes, and the arms
bracket it.

So the honest conclusion is not that a seed-role split would widen the error bar, nor that it would narrow
it, but that **there is nothing for it to remove**: the only lever on the error bar is test-set size.
"""

from __future__ import annotations

import importlib.util
import json
import statistics
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "glm_spread_vs_sampling_noise.py"
OUT = ROOT / "wiki" / "glm_spread_vs_sampling_noise_2026-11-21.json"


@pytest.fixture(scope="module")
def mod():
    spec = importlib.util.spec_from_file_location("_spread", SCRIPT)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


# ------------------------------------------------------------------ the benchmark

def test_the_benchmark_reproduces_the_committed_geometry_value(mod):
    n1, n2 = mod.expected_test_split(600, 120)
    assert (n1, n2) == (300, 60)
    assert mod.mann_whitney_null_sd(n1, n2) == pytest.approx(0.0409, abs=5e-4)


def test_the_benchmark_reproduces_the_off_geometry_n400_value(mod):
    """The second agreement is what makes the first one more than a coincidence."""
    n1, n2 = mod.expected_test_split(600, 400)
    assert mod.mann_whitney_null_sd(n1, n2) == pytest.approx(0.0264, abs=5e-4)


def test_the_benchmark_falls_as_either_class_grows(mod):
    """Monotonicity in both arguments; a formula that did not would not explain the sweep."""
    base = mod.mann_whitney_null_sd(300, 60)
    assert mod.mann_whitney_null_sd(600, 60) < base
    assert mod.mann_whitney_null_sd(300, 200) < base
    assert mod.mann_whitney_null_sd(900, 300) < mod.mann_whitney_null_sd(300, 200)


def test_the_benchmark_is_symmetric_in_its_two_classes(mod):
    assert mod.mann_whitney_null_sd(300, 60) == pytest.approx(mod.mann_whitney_null_sd(60, 300))


def test_the_expected_split_keeps_the_class_proportion(mod):
    """`heldout_auroc` trains on half the shuffled union, so the test half keeps the overall ratio."""
    n1, n2 = mod.expected_test_split(1800, 600)
    assert n1 + n2 == 1200
    assert n1 / (n1 + n2) == pytest.approx(0.75, abs=0.01)


# ------------------------------------------------------------------ the sd interval

def test_an_sd_gets_its_own_interval(mod):
    """Both earlier runs compared sd values with no interval on the sd ITSELF; at ten replicates the
    sd of an sd estimate is roughly 24% of its value, wide enough to swallow most of these effects."""
    vals = [0.50 + 0.01 * ((i % 7) - 3) for i in range(60)]
    out = mod.sd_ci(vals)
    assert out["ci95"][0] < out["sd"] < out["ci95"][1]
    assert out["n"] == 60


def test_more_replicates_give_a_tighter_sd_interval(mod):
    """The reason arm A was re-run at 240: its own sd was the loose term in every ratio."""
    base = [0.50 + 0.01 * ((i % 7) - 3) for i in range(10)]
    wide = mod.sd_ci(base)
    tight = mod.sd_ci(base * 24)
    w = wide["ci95"][1] - wide["ci95"][0]
    t = tight["ci95"][1] - tight["ci95"][0]
    assert t < w / 2, (w, t)


# ------------------------------------------------------------------ the design

def test_the_original_offset_pair_is_included_in_the_sweep(mod):
    """A replication attempt that omitted the original choice could not test it."""
    assert (10_000, 20_000) in mod.OFFSET_PAIRS
    assert len(mod.OFFSET_PAIRS) >= 4


def test_arm_A_is_run_at_many_more_replicates_than_the_offset_arms(mod):
    assert mod.N_MAIN >= 4 * mod.N_OFFSET // 2
    assert mod.N_MAIN >= 200


def test_both_outcomes_are_named_before_the_run(mod):
    p = mod.PREDICTION
    assert p["registered_before_any_number"] is True
    assert "if_A2_replicates_across_offsets" in p
    assert "if_A2_does_not_replicate" in p
    src = SCRIPT.read_text(encoding="utf-8")
    assert "STREAM_SHARING_SUPPRESSES_THE_SPREAD__ERROR_BARS_BIASED_LOW" in src
    assert "NO_STREAM_EFFECT__V2_RESULT_WAS_ONE_OFFSET_CHOICE" in src


# ------------------------------------------------------------------ all grid pairs, not a chosen three

def test_every_grid_pair_is_computed_not_a_hand_picked_subset(mod):
    """Stating a pattern from three chosen pairs is how a sign becomes a mechanism.

    Five series (four checkpoints plus the null) give C(5,2) = 10 pairs: 4 against the null and 6
    between checkpoints.
    """
    grid = {"grid": []}
    for model, base in (("zeroshot", 0.70), ("epoch1", 0.72), ("epoch2", 0.74), ("epoch3", 0.76)):
        for s in range(10):
            grid["grid"].append({"model": model, "epoch": None, "seed": s,
                                 "kmer": base + 0.01 * s, "both": base + 0.01 * s,
                                 "null": {"kmer": 0.50 + 0.01 * s, "both": 0.60 + 0.01 * s}})
    out = mod.all_grid_pairs(grid)
    for mode in ("kmer", "both"):
        assert len(out[mode]["pairs"]) == 10, sorted(out[mode]["pairs"])
        assert out[mode]["grouped"]["vs_null"]["n_pairs"] == 4
        assert out[mode]["grouped"]["vs_checkpoint"]["n_pairs"] == 6


# ------------------------------------------------------------------ the committed artifact

@pytest.mark.skipif(not OUT.exists(), reason="probe has not been run on this host")
def test_arm_A_is_consistent_with_the_pure_sampling_benchmark_in_both_modes():
    """THE result. If the coupled spread matches what sampling noise alone predicts, no seed role
    carries the spread and a role split has nothing to remove."""
    rec = json.loads(OUT.read_text(encoding="utf-8"))
    for mode in ("kmer", "both"):
        r = rec["replication"][mode]
        assert r["A_consistent_with_benchmark"] is True, (mode, r["sd_A_ci95"], r["benchmark_sd"])
        assert r["n_A"] >= 200


@pytest.mark.skipif(not OUT.exists(), reason="probe has not been run on this host")
def test_the_stream_effect_does_not_replicate_across_offsets():
    """The retraction, pinned as a retraction. One of five offset pairs is not a mechanism."""
    rec = json.loads(OUT.read_text(encoding="utf-8"))
    for mode in ("kmer", "both"):
        r = rec["replication"][mode]
        assert r["replicates"] is False, (mode, r["n_offsets_resolvably_above_A"])
        assert r["n_offsets_resolvably_above_A"] <= 1
        assert not rec["verdict_by_mode"][mode].startswith("STREAM_SHARING_SUPPRESSES")


@pytest.mark.skipif(not OUT.exists(), reason="probe has not been run on this host")
def test_the_floor_spreads_bracket_the_benchmark_across_the_sweep():
    """The clean test: the floor's true AUROC is 0.5, exactly the null the formula assumes."""
    rec = json.loads(OUT.read_text(encoding="utf-8"))
    b = rec["benchmark_vs_sweep_floors"]
    if "by_config" not in b:
        pytest.skip("sweep artifact was absent when this ran")
    assert 0.8 < b["floor_over_predicted_mean"] < 1.3, b["floor_over_predicted_mean"]
    assert "24%" in b["caveat"], "the width of the acceptable band must be stated, not implied"


@pytest.mark.skipif(not OUT.exists(), reason="probe has not been run on this host")
def test_the_refit_variance_caveat_is_declared():
    """The formula assumes a fixed scoring function; the classifier is refitted per split, which adds
    variance it does not model. That makes the benchmark a LOWER bound, and saying so is what stops
    'observed above predicted' being read as a seed-role effect."""
    rec = json.loads(OUT.read_text(encoding="utf-8"))
    joined = " ".join(rec["honest_limits"])
    assert "refitted" in joined and "LOWER bound" in joined


@pytest.mark.skipif(not OUT.exists(), reason="probe has not been run on this host")
def test_the_pairing_pattern_is_reported_over_all_pairs_with_counts():
    """Not 'pairing helps model-vs-null' -- the count of pairs where it is RESOLVABLE, out of how many."""
    rec = json.loads(OUT.read_text(encoding="utf-8"))
    for mode in ("kmer", "both"):
        g = rec["all_grid_pairs"][mode]["grouped"]
        for kind in ("vs_null", "vs_checkpoint"):
            assert g[kind]["n_pairs"] >= 4
            assert 0 <= g[kind]["n_with_resolvable_benefit"] <= g[kind]["n_pairs"]
            assert isinstance(g[kind]["median_ratio"], float)


@pytest.mark.skipif(not OUT.exists(), reason="probe has not been run on this host")
def test_the_arm_A2_values_are_shipped_so_an_outlier_is_visible():
    """A single extreme replicate can carry an sd at n=60; the sorted values make that checkable."""
    rec = json.loads(OUT.read_text(encoding="utf-8"))
    for key, v in rec["arm_A2_by_offset_pair"].items():
        vals = v["sorted_values_both"]
        assert len(vals) >= 50, key
        assert vals == sorted(vals)
        assert statistics.stdev(vals) > 0
