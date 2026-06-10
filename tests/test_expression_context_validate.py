"""Unit tests for the pure endpoint of the expression_context independent-cohort validator.

Pins the promotion-gate verdict matrix (PROMOTE iff s_upgrades==0 AND r_rescues>=1 AND n_S>=15) including
the inert-detector guard, and the Wilson 95% upper bound. No network / BLAST / cohort needed.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.expression_context_validate import compute_rescue_endpoint, wilson_upper  # noqa: E402


def _cohort(n_target_R, n_acquired_R, n_S, target_signals, acquired_signals, s_signals):
    """Build (signals, labels, intrinsic_only). target_R = intrinsic-only-R (the signal's target);
    acquired_R = R with a strong acquired carbapenemase (NON-target)."""
    signals, labels, io = [], [], []
    for i in range(n_target_R):
        labels.append("R"); io.append(True); signals.append(i < target_signals)
    for i in range(n_acquired_R):
        labels.append("R"); io.append(False); signals.append(i < acquired_signals)
    for i in range(n_S):
        labels.append("S"); io.append(False); signals.append(i < s_signals)
    return signals, labels, io


def test_promote_needs_target_subset_support():
    sig, lab, io = _cohort(n_target_R=12, n_acquired_R=0, n_S=15, target_signals=3, acquired_signals=0, s_signals=0)
    ep = compute_rescue_endpoint(sig, lab, io)
    assert ep["verdict"] == "PROMOTE"
    assert ep["n_target_R"] == 12 and ep["target_R_rescues"] == 3 and ep["s_upgrades"] == 0


def test_hold_on_any_s_upgrade():
    sig, lab, io = _cohort(n_target_R=12, n_acquired_R=0, n_S=15, target_signals=3, acquired_signals=0, s_signals=1)
    ep = compute_rescue_endpoint(sig, lab, io)
    assert ep["verdict"] == "HOLD"
    assert any("s_upgrades=1" in r for r in ep["hold_reasons"])


def test_hold_underpowered_target_is_not_falsification():
    """The REAL independent-cohort shape: n_target_R=1 (14/15 R acquired-carbapenemase) -> HOLD UNDERPOWERED,
    explicitly NOT a falsification. Even if every acquired-R fired, the verdict still HOLDs on n_target_R."""
    sig, lab, io = _cohort(n_target_R=1, n_acquired_R=14, n_S=15, target_signals=0, acquired_signals=0, s_signals=0)
    ep = compute_rescue_endpoint(sig, lab, io)
    assert ep["verdict"] == "HOLD"
    assert ep["n_target_R"] == 1 and ep["n_acquired_R"] == 14
    assert any("UNDERPOWERED" in r and "NOT a falsification" in r for r in ep["hold_reasons"])


def test_hold_on_inert_detector_when_adequately_powered():
    """Adequately-powered target subset but detector fires on none -> HOLD (inert, the real guard)."""
    sig, lab, io = _cohort(n_target_R=12, n_acquired_R=0, n_S=15, target_signals=0, acquired_signals=0, s_signals=0)
    ep = compute_rescue_endpoint(sig, lab, io)
    assert ep["verdict"] == "HOLD"
    assert any("inert" in r for r in ep["hold_reasons"])


def test_acquired_R_rescue_does_not_count_for_promotion():
    """A rescue on an acquired-carbapenemase R (non-target) must NOT satisfy the gate."""
    sig, lab, io = _cohort(n_target_R=1, n_acquired_R=11, n_S=15, target_signals=0, acquired_signals=11, s_signals=0)
    ep = compute_rescue_endpoint(sig, lab, io)
    assert ep["verdict"] == "HOLD"                         # n_target_R=1 underpowered despite 11 r_rescues
    assert ep["r_rescues"] == 11 and ep["target_R_rescues"] == 0


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
