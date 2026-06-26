"""Frozen kill-criterion (compute_probe_verdict) — bucket coverage on synthetic metric values."""
from __future__ import annotations

import importlib

import numpy as np

probe = importlib.import_module("scripts.functional_alphabet_probe")
compute_probe_verdict = probe.compute_probe_verdict


def test_beats_kmer():
    v = compute_probe_verdict(func_wl=0.80, kmer_wl=0.55, gap=0.25, gap_p=0.01, n_pairs=43)
    assert v["verdict"] == "BEATS_KMER"


def test_ties_positive_but_not_significant():
    v = compute_probe_verdict(func_wl=0.60, kmer_wl=0.55, gap=0.05, gap_p=0.30, n_pairs=43)
    assert v["verdict"] == "TIES"


def test_fails_when_not_better():
    v = compute_probe_verdict(func_wl=0.50, kmer_wl=0.62, gap=-0.12, gap_p=0.90, n_pairs=43)
    assert v["verdict"] == "FAILS"


def test_underpowered_below_min_pairs():
    v = compute_probe_verdict(func_wl=0.9, kmer_wl=0.4, gap=0.5, gap_p=0.01, n_pairs=3)
    assert v["verdict"] == "UNDERPOWERED"


def test_underpowered_on_nan():
    v = compute_probe_verdict(func_wl=float("nan"), kmer_wl=float("nan"), gap=float("nan"),
                              gap_p=float("nan"), n_pairs=0)
    assert v["verdict"] == "UNDERPOWERED"


def test_verdict_schema_keys():
    v = compute_probe_verdict(0.7, 0.5, 0.2, 0.02, 43)
    assert set(v) >= {"verdict", "func_wl", "kmer_wl", "gap", "gap_p", "n_pairs"}


def test_paired_within_lineage_degenerate_returns_nan_gap():
    """Zero shared lineages -> NaN concordance + gap (the N=40-smoke regime), not a crash."""
    y = np.array([1, 0, 1, 0])
    mlst = ["A", "B", "C", "D"]  # all unique -> no shared R/S lineage
    out = probe._paired_within_lineage([0.1, 0.2, 0.3, 0.4], [0.4, 0.3, 0.2, 0.1], y, mlst, n_perm=50)
    assert out["n_pairs"] == 0
    assert not np.isfinite(out["gap"])
