"""Pins for the target-site denominator probe: the shape conditions must be MEASURED, not asserted."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
ART = ROOT / "wiki" / "doubt_target_site_denominator_probe.json"


@pytest.mark.skipif(not ART.is_file(), reason="probe artifact not present")
def test_all_three_shape_conditions_were_measured():
    d = json.loads(ART.read_text(encoding="utf-8"))
    for k in ("negative_class_present", "enough_units", "purity_discriminating"):
        assert isinstance(d["shape"][k], bool), f"{k} not measured"


@pytest.mark.skipif(not ART.is_file(), reason="probe artifact not present")
def test_purity_must_discriminate_for_the_verdict_to_mean_anything():
    """If every candidate unit were pure, 'pure' would separate nothing and a hit would be vacuous."""
    d = json.loads(ART.read_text(encoding="utf-8"))
    if d["verdict"] != "TWO_VOCABULARIES":
        assert d["n_pure_units"] < d["n_units_with_min_carriers"]
        assert d["shape"]["purity_discriminating"] is True


@pytest.mark.skipif(not ART.is_file(), reason="probe artifact not present")
def test_a_hit_survives_correction_over_the_units_actually_tested():
    d = json.loads(ART.read_text(encoding="utf-8"))
    alpha = 0.05 / max(d["n_units_with_min_carriers"], 1)
    for s in d["surviving"]:
        assert s["p"] <= alpha, f"{s['substitution']} does not survive family-wise correction"


@pytest.mark.skipif(not ART.is_file(), reason="probe artifact not present")
def test_hits_are_outside_the_deployed_catalog_positions():
    """A unit the rule can ALREADY represent is not a completeness gap."""
    d = json.loads(ART.read_text(encoding="utf-8"))
    cat = set(d["catalogued_positions"])
    for s in d["surviving"]:
        pos = int("".join(c for c in s["substitution"] if c.isdigit()))
        assert pos not in cat, f"{s['substitution']} is already representable"


@pytest.mark.skipif(not ART.is_file(), reason="probe artifact not present")
def test_the_artifact_disclaims_being_a_curation_recommendation():
    """L2 qualifies a call; it never changes one. Data-derived NNRTI curation was measured and declined."""
    d = json.loads(ART.read_text(encoding="utf-8"))
    assert any("not itself a curation recommendation" in lim for lim in d["honest_limits"])
