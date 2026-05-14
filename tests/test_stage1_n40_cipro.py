"""Unit tests for scripts/stage1_n40_cipro.py helpers.

Pins the gate-formula verdict bucketing + paired-bootstrap-CI behavior. The
end-to-end LOSO run is tested via dry-run on the 12-strain mini cohort
(separately, not in CI — see scripts/stage1_n40_cipro.py docstring).
"""
from __future__ import annotations

import numpy as np
import pytest

from scripts.stage1_n40_cipro import (
    GATE_THRESHOLD_PP,
    paired_bootstrap_ci,
    per_mlst_breakdown,
    verdict_label,
)


class TestVerdictLabel:
    def test_clean_pass_at_5pp(self):
        assert "CLEAN PASS" in verdict_label(5.0)
        assert "CLEAN PASS" in verdict_label(7.5)

    def test_clean_pass_below_5pp_not_clean(self):
        # 4.999 pp should be NOISY, not CLEAN
        v = verdict_label(4.999)
        assert "NOISY PASS" in v
        assert "CLEAN PASS" not in v

    def test_noisy_pass_in_3_to_5_window(self):
        assert "NOISY PASS" in verdict_label(3.0)
        assert "NOISY PASS" in verdict_label(4.0)
        assert "NOISY PASS" in verdict_label(4.99)

    def test_fail_below_threshold(self):
        assert "FAIL" in verdict_label(2.99)
        assert "FAIL" in verdict_label(0.0)
        assert "FAIL" in verdict_label(-5.0)

    def test_threshold_is_3pp(self):
        # Sanity-pin the locked threshold constant
        assert GATE_THRESHOLD_PP == 3.0


class TestPairedBootstrapCI:
    def test_identical_scores_give_zero_centered_ci(self):
        # If both classifiers produce identical scores, the gap is always 0
        rng = np.random.default_rng(0)
        n = 30
        y = rng.integers(0, 2, size=n)
        # Ensure ≥1 of each class for AUROC validity
        y[0] = 0
        y[1] = 1
        scores = rng.random(n).astype(np.float32)
        mean, lo, hi, _ = paired_bootstrap_ci(y, scores, scores, n_iterations=200, seed=42)
        assert mean == pytest.approx(0.0, abs=1e-6)
        assert lo == pytest.approx(0.0, abs=1e-6)
        assert hi == pytest.approx(0.0, abs=1e-6)

    def test_a_strictly_better_than_b_gives_positive_ci(self):
        # Build scores where A is perfectly correlated with y, B is anti-correlated
        n = 30
        y = np.array([0] * 15 + [1] * 15)
        scores_a = y.astype(np.float32) + np.random.default_rng(1).normal(0, 0.01, size=n).astype(np.float32)
        scores_b = (1 - y).astype(np.float32) + np.random.default_rng(2).normal(0, 0.01, size=n).astype(np.float32)
        mean, lo, hi, _ = paired_bootstrap_ci(y, scores_a, scores_b, n_iterations=500, seed=42)
        # A's AUROC ≈ 1.0, B's AUROC ≈ 0.0 → gap ≈ 1.0; CI must be tightly positive
        assert mean > 0.9
        assert lo > 0.5

    def test_returns_four_values_three_floats_and_int(self):
        # Step 2 of Stage1_Refactor_And_Test_Hardening_Plan: paired_bootstrap_ci now
        # returns 4-tuple (mean, lo, hi, n_effective) — n_effective is the count of
        # non-degenerate resamples that contributed to the CI.
        rng = np.random.default_rng(0)
        n = 20
        y = np.array([0] * 10 + [1] * 10)
        a = rng.random(n).astype(np.float32)
        b = rng.random(n).astype(np.float32)
        result = paired_bootstrap_ci(y, a, b, n_iterations=100, seed=1)
        assert len(result) == 4
        mean, lo, hi, n_eff = result
        assert isinstance(mean, float)
        assert isinstance(lo, float)
        assert isinstance(hi, float)
        assert isinstance(n_eff, int)

    def test_handles_degenerate_resample_silently(self):
        # All same class → AUROC undefined → bootstrap should still return a result (skipping bad folds)
        n = 10
        y = np.array([0] * 5 + [1] * 5)
        a = np.linspace(0, 1, n).astype(np.float32)
        b = np.linspace(1, 0, n).astype(np.float32)
        # Few iterations; if any single-class resample comes through it's skipped
        result = paired_bootstrap_ci(y, a, b, n_iterations=50, seed=3)
        # 4-tuple now; n_effective may be < 50 if some resamples were degenerate
        assert len(result) == 4


class TestPerMlstBreakdown:
    def test_groups_by_mlst_and_counts_rs(self):
        rows = per_mlst_breakdown(
            strain_ids=["s1", "s2", "s3", "s4"],
            mlsts=["ST10", "ST10", "ST131", "ST131"],
            y=np.array([1, 0, 1, 1]),
        )
        # Two MLSTs in sorted order
        assert len(rows) == 2
        st10 = next(r for r in rows if r["mlst"] == "ST10")
        st131 = next(r for r in rows if r["mlst"] == "ST131")
        assert st10 == {"mlst": "ST10", "n": 2, "r": 1, "s": 1}
        assert st131 == {"mlst": "ST131", "n": 2, "r": 2, "s": 0}

    def test_singleton_mlsts(self):
        rows = per_mlst_breakdown(
            strain_ids=["s1", "s2", "s3"],
            mlsts=["A", "B", "C"],
            y=np.array([0, 1, 0]),
        )
        assert len(rows) == 3
        assert all(r["n"] == 1 for r in rows)
