"""Tests for the GEUVADIS organism-multimodal de-confounding evaluator.

Pure-logic tests (no data/network): the pooled-vs-within-population inflation detector
must expose a population-structure confound (the R3 lesson). Real-data loader smoke
skips when the gitignored GEUVADIS cache is absent.
"""
import os

import numpy as np
import pytest

from dna_decode.organism_multimodal.deconfound_eval import evaluate


def test_pooled_vs_within_exposes_population_confound():
    rng = np.random.default_rng(1)
    pops = np.array(sum(([p] * 90 for p in "CEU FIN GBR TSI YRI".split()), []))
    means = {p: v for p, v in zip("CEU FIN GBR TSI YRI".split(), [-2, -1, 0, 1, 2])}
    base = np.array([means[p] for p in pops])
    y_true = base + rng.normal(0, 1.0, len(pops))
    y_pred = base + rng.normal(0, 1.0, len(pops))   # predictor knows only population mean
    r = evaluate(y_true, y_pred, pops, n_perm=100)
    assert r.pooled_rho > 0.5                 # confound inflates pooled
    assert abs(r.within_rho_mean) < 0.2       # within-population collapses to null
    assert r.inflation > 0.4


def test_within_group_true_signal_survives():
    # a genuine within-population signal is NOT killed by the de-confounding
    rng = np.random.default_rng(2)
    pops = np.array(sum(([p] * 90 for p in "CEU FIN GBR TSI YRI".split()), []))
    y_true = rng.normal(0, 1, len(pops))
    y_pred = y_true + rng.normal(0, 0.5, len(pops))   # real per-individual signal
    r = evaluate(y_true, y_pred, pops, n_perm=100)
    assert r.within_rho_mean > 0.6
    assert abs(r.inflation) < 0.2             # no population confound -> pooled ~ within


def test_small_group_is_nan_not_crash():
    y = np.arange(20.0)
    g = np.array(["A"] * 18 + ["B"] * 2)      # B has n=2 < 5
    r = evaluate(y, y + 0.1, g, n_perm=10)
    assert np.isnan(r.per_group_rho["B"]) and not np.isnan(r.per_group_rho["A"])


@pytest.mark.skipif(not os.path.exists("D:/dna_decode_cache/geuvadis/GD462.GeneQuantRPKM.txt.gz"),
                    reason="GEUVADIS cache absent")
def test_real_data_loads_and_maps():
    from dna_decode.organism_multimodal.geuvadis_data import (
        load_expr, parse_sample_population, canon_pop)
    E = load_expr("D:/dna_decode_cache/geuvadis/GD462.GeneQuantRPKM.txt.gz")
    S = parse_sample_population("D:/dna_decode_cache/geuvadis/E-GEUV-1.sdrf.txt")
    assert len(E.samples) == 462 and len(E.genes) > 20000
    assert all(s in S for s in E.samples)
    assert set(canon_pop(S[s]) for s in E.samples) == {"CEU", "FIN", "GBR", "TSI", "YRI"}
