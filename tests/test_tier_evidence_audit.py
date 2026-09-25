"""Tests for the tier-vs-evidence audit.

The load-bearing tests are CALIBRATION guards. This detector was miscalibrated in both directions during
its own build: first it flagged 46 of 115 cells as over-claiming (a metric regex that missed this corpus's
bare-count idiom), then 39, then 29 -- each time because the reader, not the registry, was wrong. A hit rate
near 40% on a two-sided detector is a statement about the detector. So the rate is pinned from both sides.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts.tier_evidence_audit import (
    ADJUDICATED_FALSE_POSITIVE,
    MEASURED_TIERS,
    NEVER_MEASURED_TIERS,
    _measured_figures,
    audit,
)

REPO = Path(__file__).resolve().parent.parent
ARTIFACT = REPO / "wiki" / "tier_evidence_audit_2026-09-24.json"


# --- figure detection ----------------------------------------------------------------------------

def test_metric_named_figures_are_detected():
    for blob in ("accuracy 0.8125", "sens 1.00 / spec 0.939", "Spearman 0.7611", "AUROC 0.695"):
        assert _measured_figures(blob), blob


def test_bare_counts_are_detected():
    """This corpus's dominant idiom. Missing it flagged 46 of 115 cells as over-claiming."""
    for blob in ("85/85 core-comparable", "agreement 232/264", "64/64 exact", "160/161"):
        assert _measured_figures(blob), blob


def test_bare_correlations_are_detected():
    assert _measured_figures("tracks owner-reported height at r=+0.619")


def test_non_performance_numbers_are_not_figures():
    """Cohort sizes, coordinates, versions and allele frequencies are not performance."""
    for blob in ("3277 dogs x 29M SNVs", "chr2:233759924", "canFam3.1 liftover",
                 "OMIA 001199-9796"):
        assert not _measured_figures(blob), blob


# --- tier partition ------------------------------------------------------------------------------

def test_tier_sets_are_disjoint_and_cover_the_claim_directions():
    assert not (NEVER_MEASURED_TIERS & MEASURED_TIERS)
    assert "NOT_CENSUSED" in NEVER_MEASURED_TIERS          # the forward-cell under-claim shape
    assert "INDEPENDENT_MEASURED" in MEASURED_TIERS


# --- CALIBRATION: the rate is pinned from both sides ---------------------------------------------

def test_hit_rate_is_plausible_not_a_detector_failure():
    """A two-sided detector flagging a large fraction of the corpus is miscalibrated, not insightful.
    Pinned at <=20% because this one really did read 40% at first."""
    rep = audit()
    flagged = rep["n_under_claim_suspect"] + rep["n_over_claim_suspect"]
    assert flagged / rep["n_cells"] <= 0.20, f"{flagged}/{rep['n_cells']} is a detector failure"


def test_detector_is_not_vacuous_it_still_flags_something():
    """The opposite failure: widening the vocabulary until nothing flags. A clean zero here would mean
    the known under-claim candidates stopped being visible."""
    rep = audit()
    assert rep["n_under_claim_suspect"] >= 1
    assert rep["verdict"] == "ADJUDICATION_REQUIRED"


def test_a_planted_under_claim_is_caught():
    """Non-vacuity, planted: a never-measured tier next to a real figure and no denial must flag."""
    from scripts.tier_evidence_audit import DENIES_MEASUREMENT
    blob = "curated catalog; agreement 232/264 on the cohort"
    assert _measured_figures(blob)
    assert not any(p in blob.lower() for p in DENIES_MEASUREMENT)


def test_a_denial_phrase_makes_a_never_measured_tier_coherent():
    """A cell that says it has no free label must NOT be flagged for quoting a null baseline."""
    from scripts.tier_evidence_audit import DENIES_MEASUREMENT
    blob = "VALIDATION_DATA_WALL -- no free INDEPENDENT-colour cohort; literature contingency 0.95"
    assert any(p in blob.lower() for p in DENIES_MEASUREMENT)


def test_tool_comparator_makes_faithful_to_tool_coherent():
    """ktype is correctly FAITHFUL_TO_TOOL *because* Kaptive is a tool, even though it is measured
    (agreement 232/264). Flagging that would push a cell toward an over-claim."""
    rep = audit()
    ktype = [r for r in rep["cells"] if r["cell_id"] == "typing:Klebsiella:ktype"]
    assert ktype and ktype[0]["status"] == "CONSISTENT"


def test_adjudicated_false_positives_are_named_with_reasons():
    """Exemptions are named one at a time -- the alternative is widening until nothing flags."""
    assert len(ADJUDICATED_FALSE_POSITIVE) <= 5
    for cell_id, why in ADJUDICATED_FALSE_POSITIVE.items():
        assert cell_id and len(why) > 60, f"{cell_id} needs a real reason"


# --- the substantive finding ---------------------------------------------------------------------

def test_projected_cells_resolve_via_the_standing_card():
    """AMR/viral cells are projected and carry no prose figures; their evidence lives in a standing card.
    If this regresses, ~35 cells read as over-claiming and the audit is useless."""
    rep = audit()
    cipro = [r for r in rep["cells"]
             if r["cell_id"] == "amr:Escherichia_coli_Shigella:ciprofloxacin"]
    assert cipro and cipro[0]["status"] == "CONSISTENT"


def test_the_three_declined_cells_are_flagged():
    """The substantive finding: a tier asserting a measurement where the standing card explicitly
    declines to report one (UNDERPOWERED / ABSTAINS_BY_DESIGN)."""
    rep = audit()
    over = {r["cell_id"] for r in rep["suspects"] if r["status"] == "OVER_CLAIM_SUSPECT"}
    assert "amr:Salmonella:ciprofloxacin" in over
    assert "amr:Acinetobacter:meropenem" in over


@pytest.mark.skipif(not ARTIFACT.exists(), reason="artifact not built")
def test_artifact_records_its_posture_and_scope_limit():
    d = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    assert "TRIAGE" in d["posture"].upper()
    # it must say what it deliberately does NOT adjudicate
    assert "not_flagged_on_purpose" in d
