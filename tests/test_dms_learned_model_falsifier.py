"""Offline tests for the DMS learned-model falsifier (pure functions + join logic; no D: dependency)."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
import dms_learned_model_falsifier as fz  # noqa: E402


def test_spearman_monotonic():
    x = [1, 2, 3, 4, 5, 6]
    assert abs(fz.spearman(x, [2, 4, 6, 8, 10, 12]) - 1.0) < 1e-9      # perfectly monotone
    assert abs(fz.spearman(x, [12, 10, 8, 6, 4, 2]) + 1.0) < 1e-9      # perfectly anti-monotone


def test_spearman_shuffle_is_near_zero():
    rng = np.random.default_rng(0)
    x = rng.normal(size=500)
    y = rng.normal(size=500)   # independent
    assert abs(fz.spearman(x, y)) < 0.15


def test_spearman_too_few_points_is_nan():
    assert np.isnan(fz.spearman([1, 2], [3, 4]))


def test_load_dms_filters_multi_mutants(tmp_path):
    p = tmp_path / "assay.csv"
    p.write_text(
        "mutant,mutated_sequence,DMS_score,DMS_score_bin\n"
        "A1C,MSEQ,0.5,1\n"
        "M2K,MSEQ,-1.2,0\n"
        "A1C:M2K,MSEQ,9.9,1\n"      # multi-mutant -> must be skipped
        "BAD,MSEQ,notanum,0\n",     # unparseable score -> skipped
        encoding="utf-8",
    )
    d = fz.load_dms(str(p))
    assert d == {"A1C": 0.5, "M2K": -1.2}


def test_join_correlation_recovers_signal():
    # synthetic: AlphaMissense score anti-correlated with DMS fitness (pathogenic -> low fitness)
    am_prot = {f"A{i}C": float(i) for i in range(30)}
    dms = {f"A{i}C": float(30 - i) for i in range(30)}
    xs = [am_prot[m] for m in dms]
    ys = [dms[m] for m in dms]
    assert fz.spearman(xs, ys) < -0.9   # strong negative, as biology predicts
