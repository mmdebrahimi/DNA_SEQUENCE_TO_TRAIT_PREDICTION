"""Tests for the HIV epistasis-vs-additive world-model experiment (`scripts/hiv_epistasis.py`).

Pure-function tests run offline on synthetic data (numpy + sklearn are CORE deps -> CI-safe). The
end-to-end real-data test SKIPS when the gitignored Stanford HIVDB DataSets are absent (CI has no data).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

import hiv_epistasis as he  # noqa: E402


# --- metrics ---
def test_spearman_monotone_and_anti():
    x = np.arange(20, dtype=float)
    assert he._spearman(x, x * 3 + 1) == pytest.approx(1.0, abs=1e-9)      # monotone up
    assert he._spearman(x, -x) == pytest.approx(-1.0, abs=1e-9)            # monotone down
    assert abs(he._spearman(x, np.zeros(20))) < 1e-9 or np.isnan(he._spearman(x, np.zeros(20)))


def test_spearman_too_few_is_nan():
    assert np.isnan(he._spearman(np.arange(3.0), np.arange(3.0)))


def test_r2_perfect_and_mean():
    y = np.array([1.0, 2, 3, 4])
    assert he._r2(y, y) == pytest.approx(1.0)
    assert he._r2(y, np.full(4, y.mean())) == pytest.approx(0.0)  # predicting the mean => R2 0


# --- feature construction ---
def _rows(cells):
    """cells: list of dicts {colname: value}. Build rows with P10/P20/P30 position columns."""
    return [dict(c) for c in cells]


def test_build_presence_and_min_muts_filter():
    # 12 isolates: P10=A in 11 of them (>=MIN_MUTS=10 -> kept); P20=C in only 3 (< 10 -> dropped)
    rows = []
    for i in range(12):
        r = {"P10": "A" if i < 11 else "-", "P20": "C" if i < 3 else "-"}
        rows.append(r)
    X, feats, counts = he.build_presence(rows, ["P10", "P20"], min_muts=10)
    assert "10A" in feats and "20C" not in feats          # min-count filter works
    assert counts["10A"] == 11
    assert X.shape == (12, 1)
    assert X[:11, 0].sum() == 11 and X[11, 0] == 0


def test_build_pairs_cooccurrence_and_cap():
    # two features co-occurring in 10 isolates -> a valid pair; a third co-occurs only twice -> excluded
    rows = []
    for i in range(15):
        rows.append({"P10": "A" if i < 12 else "-", "P20": "C" if i < 12 else "-",
                     "P30": "D" if i < 2 else "-"})
    X, feats, counts = he.build_presence(rows, ["P10", "P20", "P30"], min_muts=1)
    pairs, names = he.build_pairs(X, feats, counts, top_k=10, min_cooc=10, max_pairs=50)
    assert "10A:20C" in names or "20C:10A" in names      # the co-occurring pair is kept
    assert all("30D" not in n for n in names)            # 30D co-occurs < 10 -> excluded


def test_augment_products():
    X = np.array([[1.0, 1.0], [1.0, 0.0], [0.0, 1.0]])
    Xi = he.augment(X, [(0, 1)])
    assert Xi.shape == (3, 3)
    assert list(Xi[:, 2]) == [1.0, 0.0, 0.0]             # product column = AND of the two


# --- paired bootstrap verdict ---
def test_paired_bootstrap_ci_positive_when_interaction_better():
    rng = np.random.default_rng(1)
    y = rng.normal(size=400)
    add = y + rng.normal(scale=1.2, size=400)            # weak predictor
    inter = y + rng.normal(scale=0.3, size=400)          # clearly better predictor
    out = he.paired_bootstrap(y, add, inter, boot=500, seed=0)
    assert out["int_rho"] > out["add_rho"]
    assert out["ci_positive"] is True and out["ci_lo"] > 0


def test_paired_bootstrap_not_positive_when_equal():
    rng = np.random.default_rng(2)
    y = rng.normal(size=400)
    same = y + rng.normal(scale=0.8, size=400)
    out = he.paired_bootstrap(y, same, same.copy(), boot=500, seed=0)
    assert out["ci_positive"] is False                   # identical predictors => no gain
    assert abs(out["delta_rho"]) < 1e-6


def test_paired_bootstrap_deterministic_with_seed():
    y = np.linspace(0, 1, 200)
    a = y + 0.1; b = y + 0.05
    assert he.paired_bootstrap(y, a, b, boot=200, seed=7) == he.paired_bootstrap(y, a, b, boot=200, seed=7)


# --- verdict thresholds (freeze the PASS rule) ---
def test_pass_fraction_frozen_at_half():
    assert he.PASS_FRACTION == 0.5
    assert he.N_MIN == 30 and he.MIN_MUTS == 10 and he.MIN_COOC == 10


# --- frozen-surface decoupling (the epistasis arm must not touch the deployed decoder) ---
def test_module_does_not_import_frozen_amr_surface():
    src = (REPO / "scripts" / "hiv_epistasis.py").read_text(encoding="utf-8")
    assert "amr_rules" not in src
    assert "calibrated_amr_rules" not in src
    assert "shipped_decoder_surface" not in src


# --- end-to-end on real data (skips when the gitignored HIV data is absent) ---
_NNRTI = REPO / "data" / "raw" / "hiv" / "NNRTI_DataSet.txt"


@pytest.mark.skipif(not _NNRTI.exists(), reason="Stanford HIVDB NNRTI dataset not present (gitignored)")
def test_run_class_real_data_smoke():
    res = he.run_class(_NNRTI, "NNRTI", max_drugs=1, seed=0)
    m = next(iter(res["per_drug"].values()))
    assert m["powered"] is True
    assert 0.0 <= m["add_rho"] <= 1.0                    # additive baseline orders isolates well
    assert "delta_rho" in m and "ci_lo" in m and "top_interactions" in m


@pytest.mark.skipif(not _NNRTI.exists(), reason="Stanford HIVDB NNRTI dataset not present (gitignored)")
def test_run_all_verdict_consistent_with_its_counts():
    res = he.run_all(_NNRTI.parent, classes=["NNRTI"], max_drugs=1, seed=0)
    frac = res["fraction_ci_positive"]
    expected = ("PASS_EPISTASIS_BEATS_ADDITIVE" if (res["n_powered_cells"] > 0 and frac >= he.PASS_FRACTION)
                else ("FAIL_ADDITIVE_SUFFICES" if res["n_powered_cells"] > 0 else "NO_POWERED_CELLS"))
    assert res["verdict"] == expected                    # emitted verdict matches its own counts


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(pytest.main([__file__, "-q"]))
