"""Tests for CRyPTIC TB MIC calibration (`scripts/tb_mic_calibration.py`, Family B MIC extension).
Pure-function tests offline (numpy/sklearn/pyarrow are core -> CI-safe); the real-data test SKIPS when the
gitignored feature cache is absent.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))

import tb_mic_calibration as tm  # noqa: E402


def test_parse_mic_exact_and_censored():
    assert tm.parse_mic("1.0") == (0.0, "exact")           # log2(1)=0
    assert tm.parse_mic("4.0") == (2.0, "exact")           # log2(4)=2
    lo = tm.parse_mic("<=0.25"); assert lo[1] == "left" and lo[0] == pytest.approx(-2.0)
    hi = tm.parse_mic(">8"); assert hi[1] == "right" and hi[0] == pytest.approx(3.0)
    assert tm.parse_mic("NA") is None and tm.parse_mic("") is None


def test_calibrate_drug_achieves_nominal_coverage_on_resolved():
    rng = np.random.default_rng(0)
    feat, mic = {}, {}
    for i in range(240):
        uid = f"g{i}"
        has = i < 120
        feat[uid] = ["rpoB_S450L"] if has else []
        # high MIC (index ~3) if determinant present, low (~ -1) if not, + noise on the ladder
        idx = (3 if has else -1) + int(rng.integers(-1, 2))
        val = 2.0 ** idx
        mic[uid] = f"{val:g}"
    out = tm.calibrate_drug(feat, mic, seed=0)
    assert out["powered"] is True
    assert out["n_determinant_features"] >= 1
    assert abs(out["cover_resolved_90"] - 0.90) <= 0.06     # split-conformal coverage-valid
    assert out["interval_fold_factor"] >= 1.0


def test_calibrate_drug_unpowered_on_tiny_input():
    out = tm.calibrate_drug({"a": ["x"]}, {"a": "1.0"}, seed=0)
    assert out["powered"] is False


def test_prereg_constants_frozen():
    assert tm.COVER_TOL == 0.05 and tm.TARGET == 0.90 and tm.MIN_SUPPORT == 10


def test_module_does_not_modify_frozen_bacterial_surface():
    src = (REPO / "scripts" / "tb_mic_calibration.py").read_text(encoding="utf-8")
    # reuses tb_amr (TB rule) + hiv harness READ-only; must not touch the frozen bacterial surface
    assert "calibrated_amr_rules" not in src and "shipped_decoder_surface" not in src


_CACHE = REPO / "data" / "processed" / "tb_mic_features_cache.json"


@pytest.mark.skipif(not _CACHE.exists(), reason="TB MIC feature cache absent (gitignored; needs parquet stream)")
def test_run_from_cache_real():
    res = tm.run(tm.DEFAULT_DUMP, tm.DEFAULT_REUSE, _CACHE, force=False)
    assert res["verdict"] in ("CALIBRATED_MIC_INTERVALS", "MIC_CENSORING_OR_MODEL_DOMINATES", "NO_POWERED_DRUGS")
    assert "rifampicin" in res["per_drug"] and "isoniazid" in res["per_drug"]


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(pytest.main([__file__, "-q"]))
