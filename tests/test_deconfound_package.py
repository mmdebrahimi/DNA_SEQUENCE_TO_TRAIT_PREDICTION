"""Prove the de-confounding methodology is importable from the INSTALLABLE package (promoted 2026-07-02)."""
import numpy as np

from dna_decode.deconfound import (
    Candidate,
    GATE_KEYS,
    cluster_from_distance,
    group_centered_association,
    group_centered_biomarker_t,
    group_centered_spearman,
    r2,
    score,
    within_group_r2,
)


def test_public_api_imports():
    # the promoted surface is importable straight from the package (not scripts/)
    assert callable(within_group_r2) and callable(group_centered_biomarker_t)
    assert callable(group_centered_spearman) and callable(group_centered_association)


def test_within_group_r2_deconfounds():
    rng = np.random.default_rng(0)
    groups = np.repeat(np.arange(4), 80)
    X = rng.integers(0, 2, (320, 3)).astype(float)
    base = np.repeat(rng.normal(0, 5, 4), 80)
    y = base + 1.5 * X[:, 0] + rng.normal(0, 0.3, 320)   # within-group signal on feature 0
    w, used = within_group_r2(X, y, groups)
    assert used == 4 and w > 0.3


def test_biomarker_survives_and_collapses():
    rng = np.random.default_rng(1)
    grp = np.repeat(["A", "B", "C"], 100)
    g = rng.integers(0, 2, 300)
    base = np.repeat(rng.normal(0, 3, 3), 100)
    y = base + 2.0 * g + rng.normal(0, 0.3, 300)         # within-group carrier effect
    assert group_centered_biomarker_t(y, g, grp)["within_lineage_t"] > 5
    g2 = (grp == "A").astype(int)                         # pure structure confound
    y2 = np.where(grp == "A", -5.0, 0.0) + rng.normal(0, 0.3, 300)
    assert abs(group_centered_biomarker_t(y2, g2, grp)["within_lineage_t"]) < 2


def test_scorecard_and_cluster():
    c = Candidate("X", "yeast", "growth", {k: "pass" for k in GATE_KEYS}, depth_estimate=1000)
    assert score(c)["verdict"] == "PASS"
    d = np.array([[0, 1, 9, 9], [1, 0, 9, 9], [9, 9, 0, 1], [9, 9, 1, 0]], float)
    assert len(set(cluster_from_distance(d, 2))) == 2
    assert abs(r2(np.array([1.0, 2, 3]), np.array([1.0, 2, 3])) - 1.0) < 1e-9
