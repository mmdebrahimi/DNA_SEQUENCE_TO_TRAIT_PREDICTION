"""Tests for HIV quantitative calibration (`scripts/hiv_quantitative_calibration.py`, Family B).
Pure conformal-math tests run offline (numpy core -> CI-safe); the real-data smoke SKIPS when the
gitignored HIV DataSets are absent.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))

import hiv_quantitative_calibration as qc  # noqa: E402


def test_conformal_q_finite_sample_quantile():
    # m=9 residuals; alpha=0.1 -> k = ceil(10*0.9) = 9 -> the 9th (max) sorted value
    r = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])
    assert qc._conformal_q(r, alpha=0.1) == pytest.approx(0.9)
    # alpha=0.5 -> k = ceil(10*0.5) = 5 -> 5th sorted value
    assert qc._conformal_q(r, alpha=0.5) == pytest.approx(0.5)


def test_conformal_q_insufficient_calib_returns_max():
    r = np.array([0.2, 0.5])           # m=2; alpha=0.05 -> k=ceil(3*0.95)=3 > m -> max
    assert qc._conformal_q(r, alpha=0.05) == pytest.approx(0.5)


def test_calibrate_drug_achieves_nominal_coverage_on_gaussian():
    rng = np.random.default_rng(0)
    n = 2000
    y = rng.normal(size=n)
    oof = y + rng.normal(scale=0.5, size=n)     # unbiased noisy predictor
    out = qc.calibrate_drug(y, oof, repeats=20, seed=0)
    assert abs(out["cover_90"] - 0.90) <= 0.03   # split-conformal is coverage-valid
    assert abs(out["cover_80"] - 0.80) <= 0.03
    assert out["calibrated_90"] is True
    assert out["fold_factor_90"] > 1.0           # interval has positive width


def test_calibrate_drug_deterministic_with_seed():
    y = np.linspace(0, 2, 300); oof = y + 0.1
    assert qc.calibrate_drug(y, oof, repeats=5, seed=3) == qc.calibrate_drug(y, oof, repeats=5, seed=3)


def test_prereg_constants_frozen():
    assert qc.COVER_TOL == 0.05 and qc.PASS_FRACTION == 0.5 and qc.TARGETS == (0.90, 0.80)


def test_module_does_not_import_frozen_amr_surface():
    src = (REPO / "scripts" / "hiv_quantitative_calibration.py").read_text(encoding="utf-8")
    assert "amr_rules" not in src and "calibrated_amr_rules" not in src and "shipped_decoder_surface" not in src


_NNRTI = REPO / "data" / "raw" / "hiv" / "NNRTI_DataSet.txt"


@pytest.mark.skipif(not _NNRTI.exists(), reason="HIV NNRTI dataset absent (gitignored)")
def test_run_all_real_smoke_nnrti():
    res = qc.run_all(_NNRTI.parent, classes=["NNRTI"], seed=0)
    efv = res["per_class"]["NNRTI"]["per_drug"]["EFV"]
    assert efv["powered"] and 0.80 <= efv["cover_90"] <= 0.98   # roughly nominal on real data
    assert res["verdict"] in ("CALIBRATED_INTERVALS", "MISCALIBRATED")


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(pytest.main([__file__, "-q"]))
