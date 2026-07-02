"""Pin the DepMap multimodal clade-centered continuous association on synthetic data (no network/big files)."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from scripts.depmap_multimodal import _clade_centered_rho  # noqa: E402


def test_recovers_within_lineage_continuous_signal():
    rng = np.random.default_rng(0)
    lin = np.repeat(["A", "B", "C"], 100)
    x = rng.normal(size=300)                                  # a continuous feature (e.g. expression)
    base = np.repeat(rng.normal(0, 3, 3), 100)                # between-lineage structure
    y = base - 0.8 * x + rng.normal(0, 0.3, 300)              # within-lineage: high x -> low y (sensitize)
    g, w, n = _clade_centered_rho(y, x, lin)
    assert n == 300
    assert w < -0.5                                           # de-confounded negative association recovered


def test_null_on_pure_structure():
    rng = np.random.default_rng(1)
    lin = np.repeat(["A", "B", "C"], 100)
    x = np.repeat(rng.normal(0, 3, 3), 100) + rng.normal(0, 0.1, 300)  # x tracks lineage only
    y = np.repeat(rng.normal(0, 3, 3), 100) + rng.normal(0, 0.3, 300)  # y tracks lineage only
    _, w, _ = _clade_centered_rho(y, x, lin)
    assert abs(w) < 0.3                                       # no within-lineage association


def test_nan_handling():
    lin = np.array(["A"] * 30)
    x = np.array([np.nan] * 15 + list(range(15)), dtype=float)
    y = np.array(list(range(15)) + [np.nan] * 15, dtype=float)
    g, w, n = _clade_centered_rho(y, x, lin)
    assert n == 0                                             # no overlapping non-NaN pairs
