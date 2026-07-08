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


# ---- J2 Phase-2 (beat 0.48) additions: windowing + ensemble ----
def test_window_for_position_math():
    # short protein: whole-sequence path, token index == residue index
    assert NB.window_for_position(5, 500, 1022) == (0, 500, 5)
    # long protein, maxlen=1022: window is always exactly maxlen residues, loc in [1, maxlen]
    assert NB.window_for_position(1, 2000, 1022) == (0, 1022, 1)          # left edge
    assert NB.window_for_position(2000, 2000, 1022) == (978, 2000, 1022)  # right edge (clamped)
    assert NB.window_for_position(1000, 2000, 1022) == (488, 1510, 512)   # centered
    for p in (1, 250, 999, 1500, 2000):
        s, e, loc = NB.window_for_position(p, 2000, 1022)
        assert e - s == 1022 and 1 <= loc <= 1022 and s <= p - 1 < e       # p inside window


def test_combine_variant_scores_zscore_average():
    # identical models preserve order; shared-key intersection; constant model contributes 0
    c = NB.combine_variant_scores([{"A": 1, "B": 2, "C": 3}, {"A": 1, "B": 2, "C": 3}])
    assert c["A"] < c["B"] < c["C"]
    inter = NB.combine_variant_scores([{"A": 1, "B": 2, "C": 3}, {"B": 9, "C": 8, "D": 7}])
    assert set(inter) == {"B", "C"}
    mixed = NB.combine_variant_scores([{"A": 5, "B": 5}, {"A": 1, "B": 2}])  # first is constant -> z=0
    assert mixed["A"] < mixed["B"]
    assert NB.combine_variant_scores([]) == {} and NB.combine_variant_scores([{"A": 1}, {"B": 1}]) == {}


def test_ensemble_merge_recomputes_rho(tmp_path):
    # two model runs (--keep-scores) over one shared assay -> ensemble combine -> finite rho, 1 assay
    var0 = [[f"A{i}G", float(i), float(i)] for i in range(1, 26)]        # model0 perfectly ranks y
    var1 = [[f"A{i}G", float(i) * 0.9, float(i)] for i in range(1, 26)]  # model1 nearly so
    import json as _j
    f0 = tmp_path / "m0.json"; f1 = tmp_path / "m1.json"
    _j.dump({"results": [{"dms_id": "X", "n": 25, "rho": 1.0, "rho_shuf": 0.0, "var": var0}]}, f0.open("w"))
    _j.dump({"results": [{"dms_id": "X", "n": 25, "rho": 1.0, "rho_shuf": 0.0, "var": var1}]}, f1.open("w"))
    m = NB.ensemble_merge([str(f0), str(f1)])
    assert m["n_assays"] == 1 and m["median_abs_rho"] > 0.9


def test_phase2_helpers_match_canonical():
    # drift guard extension: the new shared logic must be byte-equivalent across notebook + canonical
    for args in [(1, 2000, 1022), (2000, 2000, 1022), (1000, 2000, 1022), (7, 500, 1022)]:
        assert NB.window_for_position(*args) == CANON.window_for_position(*args)
    a = NB.combine_variant_scores([{"A": 1, "B": 3, "C": 2}, {"A": 2, "B": 1, "C": 5}])
    b = CANON.combine_variant_scores([{"A": 1, "B": 3, "C": 2}, {"A": 2, "B": 1, "C": 5}])
    assert a == b


def test_self_test_report_counts():
    refs = {
        "a": ("/x/a.csv", "M" * 100),    # short, present -> scorable in skip-mode
        "b": ("/x/b.csv", "M" * 2000),   # long, present  -> only scorable with windowing
        "c": ("/x/c.csv", "M" * 50),     # short, absent  -> not scorable
    }
    present = {"/x/a.csv", "/x/b.csv"}
    rep = NB.self_test_report(refs, 1022, lambda p: p in present)
    assert rep == {"total": 3, "have_csv": 2, "long_gt_maxlen": 1,
                   "scorable_skip_mode": 1, "scorable_window_mode": 2}


def test_self_test_report_matches_canonical():
    refs = {"a": ("/x/a.csv", "M" * 100), "b": ("/x/b.csv", "M" * 3000)}
    ef = lambda p: True
    assert NB.self_test_report(refs, 1022, ef) == CANON.self_test_report(refs, 1022, ef)


def test_self_test_on_real_reference_when_online(tmp_path):
    """R3 real-surface: fetch the REAL ProteinGym reference + verify load_reference + self_test_report.
    Skips cleanly when offline so the suite stays green without network."""
    import urllib.request
    try:
        data = urllib.request.urlopen(NB._REF_URL, timeout=20).read()
    except Exception:
        pytest.skip("offline / ProteinGym reference unreachable")
    (tmp_path / "DMS_substitutions.csv").write_bytes(data)
    refs = NB.load_reference(str(tmp_path))
    assert len(refs) > 200                         # real benchmark is 217 assays
    rep = NB.self_test_report(refs, 1022, lambda p: False)   # no per-assay CSVs attached here
    assert rep["total"] == len(refs) and rep["have_csv"] == 0 and rep["long_gt_maxlen"] >= 10
