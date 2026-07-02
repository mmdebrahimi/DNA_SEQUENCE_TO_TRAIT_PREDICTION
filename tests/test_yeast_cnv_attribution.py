"""Pin the yeast copy-number attribution logic on synthetic data (no data/network)."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from scripts.yeast_cnv_attribution import _perm_null, clade_centered_spearman  # noqa: E402


def test_clade_centered_recovers_within_clade_dosage():
    rng = np.random.default_rng(0)
    clades = np.repeat(np.arange(4), 80)
    c = rng.integers(1, 10, 320).astype(float)                      # copy number
    base = np.repeat(rng.normal(0, 5, 4), 80)                        # big between-clade offset (structure)
    y = base + 0.6 * c + rng.normal(0, 0.3, 320)                     # within-clade dose-response
    rho, _, cr = clade_centered_spearman(y, c, clades)
    assert rho > 0.6                                                 # de-confounded dose-response recovered
    perms = _perm_null(y, cr, clades, n=100)
    assert np.percentile(np.abs(perms), 95) < 0.4                    # permutation null is small


def test_clade_centered_null_on_pure_structure():
    rng = np.random.default_rng(1)
    clades = np.repeat(np.arange(4), 80)
    c = np.repeat(rng.normal(0, 3, 4), 80) + rng.normal(0, 0.1, 320)  # copy tracks clade only (confound)
    y = np.repeat(rng.normal(0, 3, 4), 80) + rng.normal(0, 0.3, 320)  # phenotype tracks clade only
    rho, _, _ = clade_centered_spearman(y, c, clades)
    assert abs(rho) < 0.3                                            # no within-clade dose-response


def test_direction_sign():
    rng = np.random.default_rng(2)
    clades = np.repeat(np.arange(3), 100)
    c = rng.integers(1, 8, 300).astype(float)
    y = -0.5 * c + rng.normal(0, 0.3, 300)                           # more copies -> LOWER phenotype
    rho, _, _ = clade_centered_spearman(y, c, clades)
    assert rho < 0                                                   # negative direction detected
