"""Tests for the forward blosum62 |Spearman| ProteinGym sweep (bounds the shipped forward default)."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts.forward_blosum_proteingym_sweep import REF, _spearman, score_one

ART = Path(__file__).resolve().parent.parent / "wiki" / "forward_blosum_proteingym_2026-07-17.json"


def test_spearman_uses_mid_ranks_so_a_constant_is_zero_not_one():
    assert _spearman([1, 2, 3], [5, 5, 5]) == 0.0            # the documented tie-order trap
    assert _spearman([1, 2, 3], [10, 20, 30]) == pytest.approx(1.0)


def test_spearman_is_scale_invariant():
    assert _spearman([1, 2, 3], [1, 100, 10000]) == pytest.approx(1.0)


@pytest.mark.skipif(not ART.exists(), reason="forward-blosum artifact not present")
def test_the_shipped_forward_default_is_median_modest_not_the_headline_035():
    """The cell LEADS with TEM-1 0.35; at N>=200 that is top-decile, not representative. Pin the honest
    scale number so a reader cannot take 0.35 as typical."""
    d = json.loads(ART.read_text(encoding="utf-8"))
    assert d["n_scored"] >= 200
    assert d["abs_spearman"]["median"] < 0.25               # 'modest' median, well below the 0.35 headline
    assert d["n_above_0.3"] < d["n_scored"] * 0.2           # 0.30+ is a minority (~top-13%)


@pytest.mark.skipif(not REF.exists(), reason="ProteinGym cache (D:) not mounted")
def test_a_real_assay_yields_a_bounded_abs_spearman():
    import csv
    from scripts.forward_blosum_proteingym_sweep import DMS_DIR
    with REF.open(encoding="utf-8") as fh:
        row = next(r for r in csv.DictReader(fh) if (DMS_DIR / f"{r['DMS_id']}.csv").exists())
    res = score_one(row)
    if res["status"] == "SCORED":
        assert 0.0 <= res["abs_spearman"] <= 1.0
