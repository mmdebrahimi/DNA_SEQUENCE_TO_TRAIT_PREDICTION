"""Unit tests for the Oxford live-run driver's pure control logic (no subprocess)."""
from __future__ import annotations

import scripts.run_oxford_revalidation as drv


# --------------------------------------------------------------------------- #
# worst_exit (severity 3 > 1 > 2 > 0)
# --------------------------------------------------------------------------- #
def test_worst_exit_hard_fail_dominates():
    assert drv.worst_exit([0, 1, 2, 3]) == 3


def test_worst_exit_degraded_over_gate():
    assert drv.worst_exit([0, 2, 1]) == 1   # 1 ranks above 2


def test_worst_exit_gate_over_clean():
    assert drv.worst_exit([0, 2, 0]) == 2


def test_worst_exit_all_clean():
    assert drv.worst_exit([0, 0]) == 0
    assert drv.worst_exit([]) == 0


# --------------------------------------------------------------------------- #
# roll_up_allowed (driver gating — the primary invariant)
# --------------------------------------------------------------------------- #
def test_rollup_allowed_all_clean():
    assert drv.roll_up_allowed([0, 0, 0], allow_degraded=False) is True


def test_rollup_blocked_on_hard_fail():
    assert drv.roll_up_allowed([0, 3], allow_degraded=False) is False
    assert drv.roll_up_allowed([0, 3], allow_degraded=True) is False   # hard-fail blocks even degraded


def test_rollup_degraded_gated():
    assert drv.roll_up_allowed([0, 1], allow_degraded=False) is False
    assert drv.roll_up_allowed([0, 1], allow_degraded=True) is True


def test_rollup_blocked_on_gate_refusal():
    assert drv.roll_up_allowed([0, 2], allow_degraded=True) is False


def test_rollup_nothing_ran():
    assert drv.roll_up_allowed([], allow_degraded=True) is False


# --------------------------------------------------------------------------- #
# preflight_blocks + plan_steps
# --------------------------------------------------------------------------- #
def test_preflight_blocks():
    assert drv.preflight_blocks("FAIL", allow_degraded=False) is True
    assert drv.preflight_blocks("FAIL", allow_degraded=True) is False
    assert drv.preflight_blocks("PASS", allow_degraded=False) is False


def test_plan_steps_order_and_fanout():
    steps = drv.plan_steps(["ciprofloxacin", "ceftriaxone"])
    assert steps[0] == "w0_probe" and steps[1] == "build_labels" and steps[2] == "preflight"
    assert "revalidate:ciprofloxacin" in steps and "revalidate:ceftriaxone" in steps
    assert steps[-1] == "rollup"
