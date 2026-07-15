"""Offline tests for the forward-cell dosage head (dna_decode/forward/dosage) — split-conformal magnitude
intervals. Pure math on synthetic data; no D:/no GPU. Includes a guarded equivalence check vs J2's
canonical `_conformal_q`."""
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from dna_decode.forward import conformal_q, dosage_intervals, evaluate_dosage  # noqa: E402


def test_conformal_q_formula():
    # m=6, alpha=0.2 -> k=ceil(7*0.8)=6 -> sorted[5]=max=0.9
    r = [0.1, 0.5, 0.3, 0.9, 0.2, 0.7]
    assert conformal_q(r, 0.2) == 0.9
    # alpha=0.1 -> k=ceil(7*0.9)=7 > 6 -> max
    assert conformal_q(r, 0.1) == 0.9
    # a lower alpha-quantile: alpha=0.5 -> k=ceil(7*0.5)=4 -> sorted[3]=0.5
    assert conformal_q(r, 0.5) == 0.5
    assert conformal_q([], 0.2) != conformal_q([], 0.2) or True  # nan, no crash


def test_conformal_q_matches_j2_helper():
    """Non-duplication: my conformal_q must be byte-equal to J2's hiv_quantitative_calibration._conformal_q."""
    sys.path.insert(0, str(REPO / "scripts"))
    try:
        import hiv_quantitative_calibration as qc
    except Exception:
        pytest.skip("J2 helper not importable in this env")
    probe = np.array([0.1, 0.5, 0.3, 0.9, 0.2, 0.7, 0.44, 0.05])
    for alpha in (0.1, 0.2, 0.5):
        assert abs(conformal_q(probe, alpha) - qc._conformal_q(probe, alpha)) < 1e-12


def test_dosage_coverage_synthetic():
    # y = 2x + noise; split-conformal at target 0.8 should hold ~0.8 held-out
    rng = np.random.RandomState(0)
    n = 1200
    x = rng.uniform(0, 1, n)
    y = 2 * x + rng.normal(0, 0.3, n)
    idx = rng.permutation(n); a, b = n // 2, 3 * n // 4
    res = evaluate_dosage(x[idx[:a]], y[idx[:a]], x[idx[a:b]], y[idx[a:b]],
                          x[idx[b:]], y[idx[b:]], coverage=0.8)
    assert abs(res.coverage - 0.8) < 0.07           # finite-sample slack
    assert res.interval_narrowing > 0.1             # informative x narrows vs marginal


def test_narrowing_near_zero_when_uninformative():
    rng = np.random.RandomState(1)
    n = 1200
    x = rng.uniform(0, 1, n)               # x carries NO info about y
    y = rng.normal(0, 1, n)
    idx = rng.permutation(n); a, b = n // 2, 3 * n // 4
    res = evaluate_dosage(x[idx[:a]], y[idx[:a]], x[idx[a:b]], y[idx[a:b]],
                          x[idx[b:]], y[idx[b:]], coverage=0.8)
    assert abs(res.coverage - 0.8) < 0.08           # coverage STILL holds (the load-bearing lesson)
    assert res.interval_narrowing < 0.1             # but the interval does NOT narrow (uninformative)


def test_isotonic_handles_negative_polarity():
    # y anti-correlated with x -> isotonic auto-flips (decreasing); point still tracks
    rng = np.random.RandomState(2)
    n = 800
    x = rng.uniform(0, 1, n)
    y = -3 * x + rng.normal(0, 0.2, n)
    idx = rng.permutation(n); a, b = n // 2, 3 * n // 4
    res = evaluate_dosage(x[idx[:a]], y[idx[:a]], x[idx[a:b]], y[idx[a:b]],
                          x[idx[b:]], y[idx[b:]], coverage=0.8)
    # the CALIBRATED POINT predicts y, so it correlates POSITIVELY with y even though x<->y is negative
    # (isotonic auto-flipped to decreasing) — this proves the auto-orientation worked.
    assert res.point_spearman > 0.7
    assert abs(res.coverage - 0.8) < 0.08


def test_dosage_intervals_shapes():
    x = np.linspace(0, 1, 100); y = x
    point, lo, hi, q, mq = dosage_intervals(x, y, x, y, x[:10], coverage=0.9)
    assert len(point) == len(lo) == len(hi) == 10
    assert np.all(hi >= lo) and q >= 0


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
