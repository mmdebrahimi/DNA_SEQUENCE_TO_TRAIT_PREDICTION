"""The serovar caller scored against the pinned REFERENCE TOOL, not an undocumented production field.

The 2026-09-04 run measured a -0.142 delta against NCBI `computed_types`, whose tool, version and
configuration are undocumented, and recorded that such a delta cannot separate "worse than the reference
METHOD" from "diverges from one NCBI pipeline". These tests pin the run that resolved it against
SeqSero2 1.3.2 -- and the accounting that keeps the resolution honest, since the answer came out
AGAINST our caller.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from salmserovar_seqsero2_validate import NON_SPECIFIC, score  # noqa: E402

ARTIFACT = ROOT / "wiki" / "salmserovar_seqsero2_2026-09-08.json"


@pytest.fixture(scope="module")
def art() -> dict:
    if not ARTIFACT.exists():
        pytest.skip("seqsero2 comparison artifact absent")
    return json.loads(ARTIFACT.read_text(encoding="utf-8"))


# --- scoring is the SAME judgement for every caller -----------------------------------------------

def test_abstention_is_not_a_miss():
    """Pooling the two would flatter whichever caller abstains most -- here, ours."""
    idx = {}
    assert score(None, "Enteritidis", idx) == "no_call"
    assert score("", "Enteritidis", idx) == "no_call"
    for token in ("-", "NA", "pending", "unknown", "untypeable"):
        assert score(token, "Enteritidis", idx) == "no_call", token


def test_non_specific_vocabulary_is_matched_case_insensitively():
    assert "pending" in NON_SPECIFIC and "untypeable" in NON_SPECIFIC
    assert score("PENDING", "Enteritidis", {}) == "no_call"


def test_a_real_call_is_scored_not_abstained():
    assert score("Enteritidis", "enteritidis", {}) == "hit"
    assert score("Typhimurium", "enteritidis", {}) == "miss"


# --- the headline resolution ----------------------------------------------------------------------

def test_the_delta_is_against_a_pinned_versioned_tool(art):
    """The whole point: a named version, not a production field of unknown provenance."""
    rt = art["reference_tool"]
    assert rt["name"] == "SeqSero2" and rt["version"] == "1.3.2"
    assert "seqsero2:1.3.2" in rt["image"]


def test_delta_is_computed_from_the_two_reported_accuracies(art):
    assert art["delta_vs_reference_tool"] == pytest.approx(
        art["ours"]["accuracy"] - art["seqsero2"]["accuracy"], abs=1e-9)


def test_the_verdict_follows_the_measured_delta(art):
    d = art["delta_vs_reference_tool"]
    expected = ("OURS_TRAILS_THE_REFERENCE_TOOL" if d < -0.02 else
                "OURS_MATCHES_THE_REFERENCE_TOOL" if abs(d) <= 0.02 else
                "OURS_LEADS_THE_REFERENCE_TOOL")
    assert art["verdict"] == expected


def test_all_three_callers_share_one_denominator(art):
    """Excluding a genome for one caller and not the others would compare different populations."""
    n = art["n_comparable"]
    for key in ("ours", "seqsero2", "ncbi_computed_types"):
        t = art[key]
        assert t["n"] == n
        assert t["hit"] + t["miss"] + t["no_call"] == n


def test_the_comparison_is_non_vacuous(art):
    assert art["n_comparable"] >= 100
    assert art["ours"]["hit"] > 0 and art["seqsero2"]["hit"] > 0


# --- honesty rails --------------------------------------------------------------------------------

def test_the_reference_tool_is_not_claimed_to_be_correct(art):
    joined = " ".join(art["honest_limits"]).lower()
    assert "tool, not the wet-lab assay" in joined
    assert "not a claim that seqsero2 is correct" in joined


def test_the_cohort_is_declared_prevalence_flattened(art):
    """SeqSero2 scores below its own published accuracy here; the cap-12 design is why, and it
    depresses every caller. Without this the LEVELS would read as population rates."""
    joined = " ".join(art["honest_limits"]).lower()
    assert "cap 12" in joined
    assert "not" in joined and "population-weighted" in joined


def test_residual_circularity_is_bounded_not_denied(art):
    joined = " ".join(art["honest_limits"]).lower()
    assert "bounded, not eliminated" in joined


def test_fairness_is_recorded_as_structural(art):
    assert "same" in art["fairness"].lower()
    assert "abstentions counted separately" in art["fairness"].lower()
