"""Tests for the supervised-vs-catalog blind-spot experiment (offline; real-data run skip-guarded).

The full experiment needs the gitignored Stanford dataset; here we unit-test the pure feature builder
(one-hot correctness + rare-singleton guard) on synthetic rows so CI stays green without the data.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

pytest.importorskip("sklearn")

spec = importlib.util.spec_from_file_location(
    "hiv_supervised_vs_catalog", REPO / "scripts" / "hiv_supervised_vs_catalog.py")
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def _rows(muts_per_row, pcols):
    """Build synthetic dataset rows: each row maps Pxxx -> residue letter ('-' = WT/consensus)."""
    rows = []
    for muts in muts_per_row:
        r = {c: "-" for c in pcols}
        for p, aa in muts:
            r[f"P{p}"] = aa
        rows.append(r)
    return rows


def test_build_onehot_binary_and_singleton_guard():
    pcols = [f"P{p}" for p in range(1, 6)]
    prot = "AAAAA"                        # WT = A at positions 1..5
    # pos2->N appears in 2 rows (kept), pos4->C appears in 1 row (singleton -> dropped)
    rows = _rows([[(2, "N")], [(2, "N")], [(4, "C")]], pcols)
    X, feats = mod.build_onehot(rows, pcols, drifted=set(), prot=prot)
    assert set(np.unique(X)) <= {0.0, 1.0}                 # strictly one-hot
    assert (2, "N") in feats                                # present in >=2 rows -> kept
    assert (4, "C") not in feats                            # singleton -> dropped
    j = feats.index((2, "N"))
    assert X[0, j] == 1.0 and X[1, j] == 1.0 and X[2, j] == 0.0


def test_build_onehot_excludes_wt_and_drifted():
    pcols = [f"P{p}" for p in range(1, 5)]
    prot = "AAAA"
    # pos1 stays WT ('A') in all rows -> never a feature; pos3->G twice (kept); drifted excludes pos2
    rows = _rows([[(1, "A"), (2, "N"), (3, "G")], [(2, "N"), (3, "G")]], pcols)
    X, feats = mod.build_onehot(rows, pcols, drifted={2}, prot=prot)
    assert (1, "A") not in feats                            # WT residue is never a feature
    assert all(p != 2 for p, _ in feats)                   # drifted position excluded
    assert (3, "G") in feats


def test_drug_and_cutoff_reused_from_prior():
    # keeps the comparison apples-to-apples with the ESM closed-negative
    assert mod.DRUG == "EFV" and mod.CUTOFF == 3.0
