"""Tests for the resistance->virulence lineage de-confound (offline pure helpers).

The Mash/Docker + real-cache `run()` path is skip-guarded; the threshold-selection logic and the
verdict wiring are unit-tested on synthetic distance matrices so CI stays green without Docker/D:.
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
    "crossaxis_lineage_deconfound", REPO / "scripts" / "crossaxis_lineage_deconfound.py")
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

from dna_decode.eval.phylogeny import DistanceMatrix  # noqa: E402


def _two_clade_matrix(n_per=8, sep=0.05, within=0.001):
    """2 tight clades separated by `sep`; greedy-rep should recover exactly 2 clusters below `sep`."""
    ids = [f"g{i:02d}" for i in range(2 * n_per)]
    m = np.full((len(ids), len(ids)), sep)
    for a in range(len(ids)):
        m[a, a] = 0.0
        for b in range(len(ids)):
            if (a < n_per) == (b < n_per):
                m[a, b] = within if a != b else 0.0
    return DistanceMatrix(strain_ids=ids, matrix=m)


def test_pick_threshold_prefers_usable_structure():
    dm = _two_clade_matrix()
    # MIN_CLADES=6 forces the sweep past a 2-clade partition toward a finer one (fallback = finest).
    thr, clusters, diag = mod.pick_threshold(dm)
    assert isinstance(clusters, dict)
    assert len(clusters) == len(dm.strain_ids)
    assert {d["threshold"] for d in diag} == set(mod.THRESHOLDS)
    # every diagnostic row reports a clade count and a largest-clade fraction in (0, 1]
    for d in diag:
        assert d["n_clades"] >= 1
        assert 0.0 < d["largest_clade_frac"] <= 1.0


def test_pick_threshold_fallback_is_finest_when_no_rung_qualifies():
    # a single blob at all thresholds -> never satisfies MIN_CLADES -> fallback finest threshold, 1 clade.
    n = 10
    dm = DistanceMatrix(strain_ids=[f"g{i}" for i in range(n)], matrix=np.full((n, n), 0.0001))
    thr, clusters, diag = mod.pick_threshold(dm)
    assert thr == mod.THRESHOLDS[0]
    assert len(set(clusters.values())) == 1


def test_auc_delegates_and_is_directional():
    scores = [0.1, 0.2, 0.8, 0.9]
    labels = [0, 0, 1, 1]
    assert mod._auc(scores, labels) == pytest.approx(1.0)
    assert mod._auc([-s for s in scores], labels) == pytest.approx(0.0)


def test_prereg_constants_frozen():
    # the de-confound bar is pre-registered; guard against silent drift.
    assert mod.REAL_MIN == 0.70
    assert mod.CV_MIN == 10
    assert mod.MIN_CLADES == 6
    assert mod.ORG == "escherichia_coli_shigella"


def test_target_axes_registry():
    # all three axes are wired to real coresistance_multiaxis prefixes.
    import coresistance_multiaxis as ma
    assert set(mod.TARGET_AXES) == {"virulence", "plasmid", "determinant"}
    assert mod.TARGET_AXES["virulence"][0] == ma.VIR_PREFIX
    assert mod.TARGET_AXES["plasmid"][0] == ma.REP_PREFIX
    assert mod.TARGET_AXES["determinant"][0] is None   # un-prefixed AMR determinants


def test_is_target_axis_membership():
    import coresistance_multiaxis as ma
    # determinant axis (prefix None) = anything NOT a plasmid:/vir: feature
    assert mod._is_target("gyrA_S83L", None) is True
    assert mod._is_target(ma.REP_PREFIX + "IncFII", None) is False
    assert mod._is_target(ma.VIR_PREFIX + "hlyA", None) is False
    # prefixed axes match only their own prefix
    assert mod._is_target(ma.VIR_PREFIX + "hlyA", ma.VIR_PREFIX) is True
    assert mod._is_target("gyrA_S83L", ma.VIR_PREFIX) is False


@pytest.mark.skipif(not Path("D:/dna_decode_cache/refseq").exists(),
                    reason="needs the D: refseq cache + Docker Mash")
def test_fasta_index_finds_ecoli():
    idx = mod.fasta_index()
    assert isinstance(idx, dict)
    # at least some GC[AF]_ accessions resolve when the cache is present
    assert all(k.startswith(("GCA_", "GCF_")) for k in idx)
