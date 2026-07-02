"""Pin the DepMap decoder's attribution + de-confounding logic on synthetic data (no data/network)."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from scripts.depmap_decoder import univariate_top, within_lineage_biomarker  # noqa: E402


def test_univariate_top_finds_driver():
    rng = np.random.default_rng(0)
    n = 300
    g = rng.integers(0, 2, n).astype(float)
    y = -3.0 * g + rng.normal(0, 0.5, n)                 # gene 0 drives y (mutant -> lower)
    G = np.column_stack([g, rng.normal(size=n), rng.normal(size=n)])
    top = univariate_top(y, G, np.array(["DRIVER", "n1", "n2"]))
    assert top[0][0] == "DRIVER" and top[0][1] < 0        # driver is #1, negative t


def test_within_lineage_biomarker_survives_when_real():
    # gene separates y INSIDE every lineage (real mechanism) -> large within-lineage t
    rng = np.random.default_rng(1)
    lin = np.repeat(["A", "B", "C"], 100)
    g = rng.integers(0, 2, 300).astype(float)
    base = np.repeat(rng.normal(0, 3, 3), 100)            # big between-lineage offset
    y = base + 2.0 * g + rng.normal(0, 0.3, 300)          # within-lineage effect of g
    r = within_lineage_biomarker(y, g, lin)
    assert r["within_lineage_t"] > 5                      # survives de-confounding
    assert len(r["per_lineage_delta_lfc"]) == 3


def test_within_lineage_biomarker_collapses_on_pure_confound():
    # gene ONLY marks a sensitive lineage (no within-lineage effect) -> within-lineage t ~ 0
    rng = np.random.default_rng(2)
    lin = np.repeat(["A", "B", "C"], 100)
    g = (lin == "A").astype(float)                        # gene present iff lineage A
    y = np.where(lin == "A", -5.0, 0.0) + rng.normal(0, 0.3, 300)   # A is sensitive
    r = within_lineage_biomarker(y, g, lin)
    assert abs(r["within_lineage_t"]) < 2                 # no within-lineage signal (pure confound)
