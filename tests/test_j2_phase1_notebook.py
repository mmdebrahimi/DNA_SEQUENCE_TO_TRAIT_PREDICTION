"""Offline tests for the J2 Phase 1 free-GPU runner (no torch, no network, no GPU).

Covers the pure logic (spearman / parse_variant / load_dms / shard / merge) + a DRIFT GUARD that pins the
self-contained notebook copy's scoring core to the canonical scripts/esm_zeroshot_dms.py.
"""
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent


def _load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


NB = _load(ROOT / "notebooks" / "j2_phase1_esm2_proteingym.py", "j2nb")
CANON = _load(ROOT / "scripts" / "esm_zeroshot_dms.py", "esm_zeroshot_dms_canon")


def test_spearman_monotonic():
    assert NB.spearman([1, 2, 3, 4, 5], [10, 20, 30, 40, 50]) == pytest.approx(1.0)
    assert NB.spearman([1, 2, 3, 4, 5], [50, 40, 30, 20, 10]) == pytest.approx(-1.0)


def test_spearman_shuffle_near_zero():
    rng = np.random.default_rng(0)
    x = rng.normal(size=400)
    y = rng.normal(size=400)  # independent
    assert abs(NB.spearman(x, y)) < 0.15


def test_parse_variant():
    assert NB.parse_variant("M1A") == ("M", 1, "A")
    assert NB.parse_variant("K103N") == ("K", 103, "N")


def test_load_dms_filters_multi_mutants(tmp_path):
    p = tmp_path / "assay.csv"
    p.write_text("mutant,DMS_score\nM1A,1.0\nK2R,2.0\nM1A:K2R,9.9\nA3G;D4E,8.8\nbad,x\n", encoding="utf-8")
    d = NB.load_dms(str(p))
    assert d == {"M1A": 1.0, "K2R": 2.0}          # multi-mutants + unparseable dropped


def test_shard_partition_disjoint_covers_all_deterministic():
    ids = [f"DMS_{i}" for i in range(23)]
    n = 3
    shards = [NB.shard_assays(ids, i, n) for i in range(n)]
    # disjoint
    for a in range(n):
        for b in range(a + 1, n):
            assert not (set(shards[a]) & set(shards[b]))
    # covers all
    assert set().union(*shards) == set(ids)
    # deterministic + strided on sorted order
    assert shards[0] == sorted(ids)[0::n]
    assert NB.shard_assays(ids, 1, n) == NB.shard_assays(ids, 1, n)


def test_shard_out_of_range():
    with pytest.raises(ValueError):
        NB.shard_assays(["a", "b"], 3, 2)


def test_merge_pools_and_dedups(tmp_path):
    s0 = tmp_path / "s0.json"
    s1 = tmp_path / "s1.json"
    json.dump({"results": [{"dms_id": "A", "rho": 0.6, "rho_shuf": 0.01},
                            {"dms_id": "B", "rho": -0.5, "rho_shuf": 0.00}]}, s0.open("w"))
    json.dump({"results": [{"dms_id": "C", "rho": 0.4, "rho_shuf": 0.02}]}, s1.open("w"))
    m = NB.merge([str(s0), str(s1)])
    assert m["n_assays"] == 3
    assert m["median_abs_rho"] == pytest.approx(0.5)          # median(|0.6|,|-0.5|,|0.4|)
    assert m["ok"] is True                                    # >=0.45 and shuffled <0.05


def test_merge_dedups_disjoint_collision(tmp_path):
    s0 = tmp_path / "s0.json"
    s1 = tmp_path / "s1.json"
    json.dump({"results": [{"dms_id": "A", "rho": 0.6, "rho_shuf": 0.0}]}, s0.open("w"))
    json.dump({"results": [{"dms_id": "A", "rho": 0.6, "rho_shuf": 0.0}]}, s1.open("w"))  # same id
    m = NB.merge([str(s0), str(s1)])
    assert m["n_assays"] == 1                                 # dedup by dms_id


def test_drift_guard_scoring_core_matches_canonical():
    """The self-contained notebook copy's pure scoring core must equal scripts/esm_zeroshot_dms.py."""
    rng = np.random.default_rng(1)
    x = rng.normal(size=60).tolist()
    y = rng.normal(size=60).tolist()
    assert NB.spearman(x, y) == pytest.approx(CANON.spearman(x, y), abs=1e-12)
    for v in ("M1A", "K103N", "W42F"):
        assert NB.parse_variant(v) == CANON.parse_variant(v)


def test_pass_constants_match_canonical():
    assert (NB.PASS_BAR, NB.STRETCH, NB.SHUFFLE_MAX) == (CANON.PASS_BAR, CANON.STRETCH, CANON.SHUFFLE_MAX)
