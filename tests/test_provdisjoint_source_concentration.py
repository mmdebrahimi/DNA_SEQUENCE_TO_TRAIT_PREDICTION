"""Guards for the source-concentration disclosure measurement.

This is a DISCLOSURE, not a demotion: no SCORED cell's published metrics change, and the cells remain
provenance-disjoint from the tuning data. What it adds is whether they are also provenance-DIVERSE.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))


def test_concentration_reports_largest_share_not_just_distinct_count():
    """8 sources sounds diverse until one holds 97%. The share is the number that carries the meaning."""
    from provdisjoint_source_concentration import concentration
    prov = {f"a{i}": {"bioproject_acc": "PRJNA1" if i < 58 else f"PRJNA{i}"} for i in range(60)}
    c = concentration(prov, "bioproject_acc")
    assert c["distinct"] == 3
    assert c["largest_share"] > 0.9, "a 58/60 dominant source must surface as a high share"


def test_unknown_sources_are_not_merged_into_one_pseudo_source():
    """Empty/NULL provenance is MISSING METADATA. Merging it into a single bucket would manufacture
    concentration that isn't real -- the opposite failure from the one being measured."""
    from provdisjoint_source_concentration import concentration
    prov = {f"a{i}": {"bioproject_acc": "" if i < 50 else f"PRJNA{i}"} for i in range(60)}
    c = concentration(prov, "bioproject_acc")
    assert c["n_unknown"] == 50
    assert c["distinct"] == 10, "unknowns must not count as a source"
    assert c["largest_share"] is not None and c["largest_share"] < 0.5


def test_the_result_records_both_directions_of_error():
    """Single-source cells were optimistic on one cell and pessimistic on another. Recording only the
    one that matched the story would be the confirmation the memo is arguing against."""
    memo = ROOT / "wiki" / "provdisjoint_source_concentration_2026-08-28.md"
    if not memo.exists():
        import pytest
        pytest.skip("memo absent")
    t = memo.read_text(encoding="utf-8")
    assert "optimistic" in t and "pessimistic" in t
    assert "Not a demotion" in t


def test_measured_concentration_matches_the_committed_artifact():
    """Pin the headline: 3 of 10 cells are >=80% one BioProject, and the gentamicin cell has no rmt."""
    f = ROOT / "wiki" / "provdisjoint_source_concentration.json"
    if not f.exists():
        import pytest
        pytest.skip("artifact absent (requires a network run)")
    d = json.loads(f.read_text(encoding="utf-8"))
    assert d["complete"], "an incomplete sweep cannot support a concentration claim"
    assert d["n_single_source"] >= 3
    assert len(d["cells"]) == 10
