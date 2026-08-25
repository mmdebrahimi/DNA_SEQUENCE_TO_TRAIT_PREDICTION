"""The mutation-order pooling correction (scripts/epistasis_pooling_check.py).

The 2026-07-27 epistasis sweep published per-protein joint-vs-additive Deltas POOLED across mutation
orders, and headlined ParD at -0.283 as "joint can be MUCH worse". Pooling across a variable BOTH scores
track inflates the metric -- the same shape as clonality inflation on the AMR side, with mutation order in
place of the clone. Within order the ParD Delta is -0.053, five times smaller.

These tests pin the arithmetic (pure, no assay files needed) and the corrected published numbers.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from scripts.epistasis_pooling_check import pooling_report, within_order  # noqa: E402


def test_within_order_is_an_n_weighted_mean_not_a_plain_mean():
    """Orders carry different n; a plain mean would silently over-weight a small order."""
    po = {"2": {"additive": 0.0, "joint": 0.0, "n": 900}, "3": {"additive": 1.0, "joint": 1.0, "n": 100}}
    assert within_order(po, "additive") == pytest.approx(0.1)      # not 0.5


def test_pooling_report_recovers_the_inflation_and_skips_single_order_proteins():
    results = [
        # a protein whose pooled rho EXCEEDS its rho at every order -> the confounder signature
        {"dms": "ParD_like", "seqlen": 93, "additive": 0.543, "joint": 0.260,
         "per_order": {"2": {"additive": 0.301, "joint": 0.355, "n": 300},
                       "3": {"additive": 0.356, "joint": 0.270, "n": 300},
                       "4": {"additive": 0.120, "joint": -0.007, "n": 300}}},
        # single-order protein: a pooling effect is undefined, so it must be dropped, not reported as 0
        {"dms": "OneOrder", "seqlen": 200, "additive": 0.5, "joint": 0.5,
         "per_order": {"2": {"additive": 0.5, "joint": 0.5, "n": 300}}},
    ]
    rows = pooling_report(results)
    assert [r["protein"] for r in rows] == ["ParD_like"]
    r = rows[0]
    assert r["within_delta"] == pytest.approx(-0.0527, abs=1e-3)
    assert r["pooled_delta"] == pytest.approx(-0.283, abs=1e-3)
    # the anomaly is ADDITIVE's outsized pooling bonus, not joint collapsing
    assert r["additive_pooling_gain"] > 4 * r["joint_pooling_gain"]
    assert abs(r["pooled_delta"]) > 4 * abs(r["within_delta"])


def test_the_committed_sweep_still_shows_the_artifact_it_was_corrected_for():
    """Guards the correction against a silent re-run that changes the underlying numbers."""
    p = ROOT / "wiki" / "forward_epistasis_sweep_2026-07-27.json"
    if not p.exists():
        pytest.skip("sweep artifact absent")
    rows = {r["protein"]: r for r in pooling_report(json.loads(p.read_text(encoding="utf-8"))["results"])}
    pard = rows["F7YBW8_MESOW_Aakre_2015"]
    assert pard["pooled_delta"] == pytest.approx(-0.283, abs=0.002)
    assert pard["within_delta"] == pytest.approx(-0.053, abs=0.002)
    # joint's pooling gain is near-constant across proteins; additive's is not
    joint_gains = [r["joint_pooling_gain"] for r in rows.values()]
    assert max(joint_gains) - min(joint_gains) < 0.02
    assert max(r["additive_pooling_gain"] for r in rows.values()) > 0.25


def test_the_published_artifacts_no_longer_headline_the_pooled_number():
    """The report card is AUTO-GENERATED -- fixing only the .md would be wiped on the next rebuild, so
    the generator is what carries the correction."""
    card = ROOT / "wiki" / "forward_validation_report_card.md"
    gen = ROOT / "scripts" / "build_forward_report_card.py"
    if not (card.exists() and gen.exists()):
        pytest.skip("forward report card or its generator absent")
    for f in (card, gen):
        text = f.read_text(encoding="utf-8", errors="replace")
        assert "delta -0.283" not in text, f"{f.name} still headlines the pooled figure"
    assert "WITHIN-ORDER delta -0.053" in gen.read_text(encoding="utf-8", errors="replace")


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
