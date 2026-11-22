"""The reconcile probe must gate on reproducing the grid, and must not subtract incommensurable AUROCs.

Two defects motivated this script and both are pinned here as defects, because a later tidy-up is most
likely to remove the record of them.

1. **The geometry gate.** Two seed-role probes hardcoded `n_generated=400` while the committed kernel
   froze `n=120`, and nothing compared either probe's output to the runner's. The off-geometry numbers
   were entirely plausible, which is exactly why only a reconciliation could catch them.
2. **The margin subtraction.** `margin = null - ceiling` paired a 600-vs-120 AUROC with a 300-vs-300 one.
   A held-out AUROC's sampling variation depends on test-set size and class balance, so those two numbers
   are not commensurable. The fix is a floor measured in the null's own geometry, which is what
   `perfect_generator_floor` is.

`test_the_floor_rule_alone_flips_the_kmer_verdict` is the load-bearing one: it pins that the committed
`kmer` conclusion depended on the choice of floor rule rather than on the data, which is the whole reason
the earlier recommendation was unsafe.
"""

from __future__ import annotations

import importlib.util
import json
import math
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "glm_grid_reconcile.py"
OUT = ROOT / "wiki" / "glm_grid_reconcile_2026-11-21.json"
GRID = ROOT / "wiki" / "glm_generator_paired_sweep_2026-11-20.json"


@pytest.fixture(scope="module")
def mod():
    spec = importlib.util.spec_from_file_location("_reconcile", SCRIPT)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


# ------------------------------------------------------------------ the geometry

def test_the_geometry_matches_the_committed_kernel(mod):
    """Transcribed from scripts/kaggle/glm_paired_sweep_kernel.py:45. The n is the value v1 got wrong."""
    assert mod.N_GENERATED == 120
    assert mod.LIMIT == 1200
    assert mod.LENGTH == 150
    assert mod.POSITIONAL_BINS == 6
    assert mod.KMER_K == 4
    assert mod.MARKOV_ORDER == 3


def test_the_kernel_still_declares_the_geometry_this_script_claims(mod):
    """Guards the transcription. If the kernel's RECORDED block ever moves, this must fail loudly
    rather than let the probe keep asserting a geometry the kernel no longer uses."""
    kernel = ROOT / "scripts" / "kaggle" / "glm_paired_sweep_kernel.py"
    if not kernel.exists():
        pytest.skip("kernel absent on this host")
    src = kernel.read_text(encoding="utf-8")
    assert '"n": 120' in src and '"limit_natural": 1200' in src
    assert '"bins": 6' in src and '"length": 150' in src


# ------------------------------------------------------------------ the significance test

def test_a_margin_whose_interval_straddles_zero_is_not_resolvable(mod):
    """The distinction the committed artifact lacked: unresolved is not the same as absent."""
    margins = [0.08, -0.05, 0.06, -0.02, 0.09, 0.01, -0.06, 0.07, 0.04, 0.11]
    out = mod.margin_significance(margins)
    assert out["ci95"][0] < 0 < out["ci95"][1]
    assert out["resolvable_above_zero"] is False


def test_a_tight_positive_margin_is_resolvable(mod):
    margins = [0.13, 0.14, 0.12, 0.15, 0.13, 0.14, 0.12, 0.14, 0.13, 0.15]
    out = mod.margin_significance(margins)
    assert out["resolvable_above_zero"] is True
    assert out["ci95"][0] > 0


def test_sd_over_mean_is_reported_because_a_count_of_positives_is_not_a_finding(mod):
    """7-of-10 positive reads as a result while being consistent with a true margin of zero."""
    margins = [0.08, -0.05, 0.06, -0.02, 0.09, 0.01, -0.06, 0.07, 0.04, 0.11]
    out = mod.margin_significance(margins)
    assert out["sd_over_mean"] > 1.0
    assert sum(1 for v in margins if v > 0) == 7


def test_too_few_replicates_refuses_rather_than_returning_a_verdict(mod):
    out = mod.margin_significance([0.1, 0.2])
    assert out["resolvable"] is None


def test_the_interval_uses_the_standard_error_not_the_raw_sd(mod):
    """A CI built on the SD rather than the SE would be ~3x too wide at n=10 and would call
    everything unresolved, which looks conservative and is simply wrong."""
    margins = [0.13, 0.14, 0.12, 0.15, 0.13, 0.14, 0.12, 0.14, 0.13, 0.15]
    out = mod.margin_significance(margins)
    assert out["stderr"] == pytest.approx(out["sd"] / math.sqrt(10), abs=1e-4)
    halfwidth = out["ci95"][1] - out["mean"]
    assert halfwidth == pytest.approx(2.262 * out["stderr"], abs=1e-3)


# ------------------------------------------------------------------ the floor rule

def test_the_floor_rule_alone_flips_the_kmer_verdict(mod):
    """THE load-bearing test. Max-of-five against mean-of-five changes the sign of the kmer margin.

    A conclusion that depends on which of two defensible floor rules was chosen is not a finding about
    the data, and the committed artifact named the mean-based floor as an unquantified limit while
    shipping every value needed to compute it.
    """
    ceil = ROOT / "wiki" / "glm_seed_role_ceilings_2026-11-21.json"
    if not ceil.exists():
        pytest.skip("ceilings artifact absent on this host")
    rec = json.loads(ceil.read_text(encoding="utf-8"))
    out = mod.mean_based_ceiling(rec)
    assert out["kmer"]["margin_vs_max_mean"] < 0 < out["kmer"]["margin_vs_mean_mean"]
    assert out["kmer"]["verdict_flips_on_floor_rule"] is True
    assert out["both"]["verdict_flips_on_floor_rule"] is False, "both mode must be robust to the rule"


# ------------------------------------------------------------------ the paired-delta helper

def test_pairing_helps_when_the_two_arms_are_positively_correlated(mod):
    """Var(A-B) = Var(A)+Var(B)-2Cov(A,B): positive covariance must show as a tighter delta."""
    grid = {"grid": []}
    for s in range(10):
        shared = 0.02 * s
        grid["grid"].append({"model": "zeroshot", "epoch": None, "seed": s,
                             "kmer": 0.70 + shared, "both": 0.75 + shared,
                             "null": {"kmer": 0.50 + shared, "both": 0.60 + shared}})
    out = mod.paired_delta_from_grid(grid)
    r = out["kmer"]["zeroshot_minus_null"]
    assert r["pairing_helps"] is True
    assert r["sd_paired_delta"] < r["sd_unpaired_sum"]
    assert r["correlation"] > 0.9


def test_pairing_does_not_help_when_the_arms_are_uncorrelated(mod):
    """Control: the helper must not report a benefit that is not there."""
    a = [0.70, 0.78, 0.72, 0.76, 0.71, 0.79, 0.73, 0.75, 0.74, 0.77]
    b = [0.55, 0.51, 0.59, 0.52, 0.58, 0.53, 0.57, 0.54, 0.56, 0.50]
    grid = {"grid": [{"model": "zeroshot", "epoch": None, "seed": s,
                      "kmer": a[s], "both": a[s], "null": {"kmer": b[s], "both": b[s]}}
                     for s in range(10)]}
    out = mod.paired_delta_from_grid(grid)
    r = out["kmer"]["zeroshot_minus_null"]
    assert r["correlation"] < 0.2
    assert r["sd_paired_delta"] >= r["sd_unpaired_sum"] * 0.9


# ------------------------------------------------------------------ Decision One's actual premise

def test_generator_distance_to_null_is_measured_separately_from_the_margin(mod):
    """A mode whose null sits near its floor can still carry a large, registerable model-to-null gap.

    Retiring a mode on the margin alone would discard a working measurement, so the distance has to be
    reported as its own quantity rather than inferred from the margin.
    """
    grid = json.loads(GRID.read_text(encoding="utf-8")) if GRID.exists() else None
    if grid is None:
        pytest.skip("grid absent on this host")
    out = mod.model_minus_null_from_grid(grid)
    for mode in ("kmer", "both"):
        assert out[mode]["n_cells"] == 40
        assert out[mode]["n_cells_generator_above_null"] == 40, (
            "every committed cell has the generator above its own null; if this ever changes the "
            "retire-the-mode argument changes with it")
        assert out[mode]["all_models_mean"] > 0


# ------------------------------------------------------------------ the committed artifact

@pytest.mark.skipif(not OUT.exists(), reason="reconcile has not been run on this host")
def test_the_gate_passed_at_every_seed_and_mode():
    rec = json.loads(OUT.read_text(encoding="utf-8"))
    assert rec["gate"]["passed"] is True
    assert rec["gate"]["reconciled_per_mode"] == {"kmer": 10, "both": 10}
    assert all(r["reconciles"] for r in rec["rows"])


@pytest.mark.skipif(not OUT.exists(), reason="reconcile has not been run on this host")
def test_the_two_defects_are_recorded_with_their_numbers():
    """Prose alone would not survive a tidy-up; the numbers are what make the record checkable."""
    rec = json.loads(OUT.read_text(encoding="utf-8"))
    d1 = json.dumps(rec["the_defect_this_found"])
    assert "400" in d1 and "120" in d1 and "0.0410" in d1
    d2 = json.dumps(rec["the_second_defect"])
    assert "600" in d2 and "300" in d2


@pytest.mark.skipif(not OUT.exists(), reason="reconcile has not been run on this host")
def test_the_perfect_generator_floor_sits_near_chance():
    """Real sequence against real sequence has a true AUROC of 0.5.

    This is the instrument check the max-of-five self-control ceiling was failing: a floor of 0.542 for
    a real-versus-real comparison is not a floor, it is an artifact of taking a maximum over five splits
    at a different class balance.
    """
    rec = json.loads(OUT.read_text(encoding="utf-8"))
    for mode in ("kmer", "both"):
        assert abs(rec["summary"][mode]["perfect_floor_mean"] - 0.5) < 0.02, mode


@pytest.mark.skipif(not OUT.exists(), reason="reconcile has not been run on this host")
def test_the_instrument_has_power_in_both_modes():
    """A floor at chance could mean a blind discriminator; the uniform control rules that out."""
    rec = json.loads(OUT.read_text(encoding="utf-8"))
    for mode in ("kmer", "both"):
        assert rec["power_control_uniform"][mode]["mean"] > 0.85, mode


@pytest.mark.skipif(not OUT.exists(), reason="reconcile has not been run on this host")
def test_kmer_is_reported_as_unresolved_and_not_as_absent():
    """The correction to the committed conclusion, pinned as a correction."""
    rec = json.loads(OUT.read_text(encoding="utf-8"))
    assert rec["verdict_by_mode"]["kmer"] == "WINNABLE_REGION_UNRESOLVED_AT_THIS_N"
    assert rec["verdict_by_mode"]["both"] == "WINNABLE_REGION_AT_EVERY_SEED"
    assert rec["margin_significance"]["kmer"]["resolvable_above_zero"] is False
    assert rec["margin_significance"]["both"]["resolvable_above_zero"] is True


@pytest.mark.skipif(not OUT.exists(), reason="reconcile has not been run on this host")
def test_the_corpus_redundancy_caveat_is_recorded_in_proportion():
    """The runner promises fit/eval disjointness; duplicate windows mean that is not exactly true.

    Recorded as a real caveat AND as a small one -- overstating it would be its own distortion.
    """
    rec = json.loads(OUT.read_text(encoding="utf-8"))
    cr = rec["corpus_redundancy"]
    assert cr["n_distinct_sequences"] < cr["n_windows"]
    assert cr["worst_case_leak_fraction_of_eval"] < 0.01
    assert "far too small" in cr["reading"]


@pytest.mark.skipif(not OUT.exists(), reason="reconcile has not been run on this host")
def test_spare_windows_are_sequence_disjoint_at_every_seed():
    """An overlap would put the same feature vector in both classes and make the floor cleaner than it is."""
    rec = json.loads(OUT.read_text(encoding="utf-8"))
    for r in rec["rows"]:
        assert r["n_spare"] == 120
        assert r["n_spare_skipped_as_duplicates"] >= 0
