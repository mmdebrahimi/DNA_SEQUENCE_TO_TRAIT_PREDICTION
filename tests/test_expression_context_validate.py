"""Unit tests for the pure endpoint of the expression_context independent-cohort validator.

Pins the promotion-gate verdict matrix (PROMOTE iff s_upgrades==0 AND r_rescues>=1 AND n_S>=15) including
the inert-detector guard, and the Wilson 95% upper bound. No network / BLAST / cohort needed.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.expression_context_validate import compute_rescue_endpoint, wilson_upper  # noqa: E402


def _cohort(n_R, n_S, r_signals, s_signals):
    """Build (signals, labels): r_signals true-R fire, s_signals true-S fire."""
    signals, labels = [], []
    for i in range(n_R):
        labels.append("R"); signals.append(i < r_signals)
    for i in range(n_S):
        labels.append("S"); signals.append(i < s_signals)
    return signals, labels


def test_promote_zero_s_upgrades_one_rescue_enough_s():
    sig, lab = _cohort(n_R=12, n_S=15, r_signals=3, s_signals=0)
    ep = compute_rescue_endpoint(sig, lab)
    assert ep["verdict"] == "PROMOTE"
    assert ep["r_rescues"] == 3 and ep["s_upgrades"] == 0 and ep["n_S"] == 15
    assert ep["abstain_rescue_rate"] == 0.25


def test_hold_on_any_s_upgrade():
    sig, lab = _cohort(n_R=12, n_S=15, r_signals=3, s_signals=1)
    ep = compute_rescue_endpoint(sig, lab)
    assert ep["verdict"] == "HOLD"
    assert any("s_upgrades=1" in r for r in ep["hold_reasons"])


def test_hold_on_inert_detector_zero_rescues():
    """Zero S-upgrades but ALSO zero R-rescues = inert detector -> HOLD (the guard)."""
    sig, lab = _cohort(n_R=12, n_S=15, r_signals=0, s_signals=0)
    ep = compute_rescue_endpoint(sig, lab)
    assert ep["verdict"] == "HOLD"
    assert any("inert" in r for r in ep["hold_reasons"])


def test_hold_on_insufficient_s():
    sig, lab = _cohort(n_R=12, n_S=10, r_signals=3, s_signals=0)
    ep = compute_rescue_endpoint(sig, lab)
    assert ep["verdict"] == "HOLD"
    assert any("n_S=10" in r for r in ep["hold_reasons"])


def test_wilson_upper_bounds():
    assert wilson_upper(0, 0) is None
    u15 = wilson_upper(0, 15)
    u10 = wilson_upper(0, 10)
    assert 0.15 < u15 < 0.25            # 0/15 -> ~0.20
    assert u10 > u15                     # smaller n -> wider bound (~0.28)
    assert wilson_upper(0, 50) < u15     # larger n -> narrower


if __name__ == "__main__":
    import pytest
    sys.exit(pytest.main([__file__, "-q"]))
