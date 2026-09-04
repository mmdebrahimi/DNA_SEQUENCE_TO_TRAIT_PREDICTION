"""The abstention diagnosis must partition by the FIRST failing axis, not by a trailing dash.

The original diagnosis counted formulas ending in '-' and concluded phase-2 flagellin was the dominant
defect. `4:H?:-` ends in '-' but failed on H1, so that count answers a different question than the one
asked. These pin the distinction and the zero headroom that killed the tempting fix.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from salmserovar_nocall_anatomy import classify  # noqa: E402

ART = ROOT / "wiki" / "salmserovar_nocall_anatomy_2026-09-04.json"
PAIRS = {("4", "i"), ("9", "g,m")}


def test_an_empty_H2_is_not_automatically_an_H2_failure():
    """THE BUG: both formulas end in '-', but only one is actually blocked by phase 2."""
    assert classify("4:H?:-", PAIRS) == "H1_phase1_flagellin_unresolved"
    assert classify("4:i:-", PAIRS) == "O_H1_valid_only_H2_blocks_it"


def test_each_axis_is_attributed_to_the_first_failure():
    assert classify("O?:i:1,2", PAIRS) == "O_antigen_unresolved"
    assert classify("O?:H?:-", PAIRS) == "both_O_and_H1_unresolved"
    assert classify("77:z:-", PAIRS) == "O_H1_called_but_pair_absent_from_table"
    assert classify(None, PAIRS) == "no_formula_at_all"
    assert classify("garbage", PAIRS) == "no_formula_at_all"


@pytest.mark.skipif(not ART.exists(), reason="anatomy artifact absent")
def test_the_committed_anatomy_contradicts_the_original_diagnosis():
    d = json.loads(ART.read_text(encoding="utf-8"))
    assert d["verdict"] == "H2_IS_NOT_THE_DOMINANT_CAUSE"
    causes = d["causes"]
    biggest = max(causes.items(), key=lambda kv: kv[1])[0]
    assert biggest == "O_antigen_unresolved", "the largest cause is the O axis, not H2"
    assert d["phase2_reachable_bucket"] < causes["O_antigen_unresolved"]
    assert "supersedes" in d, "the artifact must name the claim it corrects"


@pytest.mark.skipif(not ART.exists(), reason="anatomy artifact absent")
def test_the_zero_headroom_of_the_tempting_fix_is_recorded():
    """A fix with zero measured headroom must be documented as such so nobody builds it."""
    d = json.loads(ART.read_text(encoding="utf-8"))
    assert d["oh1_unique_fallback_headroom"]["recoverable"] == 0
    assert "zero" in d["why"]


@pytest.mark.skipif(not ART.exists(), reason="anatomy artifact absent")
def test_the_real_priority_and_its_cost_are_stated():
    d = json.loads(ART.read_text(encoding="utf-8"))
    assert d["fix_priority_by_measured_size"][0] == "O_antigen_unresolved"
    assert any("data engineering" in s for s in d["honest_limits"]), (
        "the O-antigen fix is DB coverage, not a code change -- that cost must be stated")


# --- O-antigen threshold-vs-coverage probe ------------------------------------------------------
# The second asserted cause ("DB coverage / data engineering") was also wrong. These pin the measured
# replacement AND the refusal to assert a third cause without measuring it.

O_PROBE = ROOT / "wiki" / "salmserovar_o_antigen_probe_2026-09-04.json"


@pytest.mark.skipif(not O_PROBE.exists(), reason="O-probe artifact absent")
def test_the_db_coverage_claim_was_measured_and_refuted():
    d = json.loads(O_PROBE.read_text(encoding="utf-8"))
    assert d["verdict"] == "THRESHOLD_IS_THE_BINDING_CONSTRAINT"
    v = d["verdicts"]
    assert v["sub_threshold_hit_exists"] > v.get("no_O_hit_at_any_threshold", 0)
    assert "tests_the_claim" in d, "the artifact must name the claim it tests"


@pytest.mark.skipif(not O_PROBE.exists(), reason="O-probe artifact absent")
def test_sub_threshold_hits_are_CORRECT_not_merely_present():
    """A rejected hit only matters if it names the right O group; presence alone proves nothing."""
    c = json.loads(O_PROBE.read_text(encoding="utf-8"))["sub_threshold_correctness"]
    assert c["names_wrong_O_group"] == 0 and c["names_correct_O_group"] == c["n"]


@pytest.mark.skipif(not O_PROBE.exists(), reason="O-probe artifact absent")
def test_the_finding_is_scoped_to_its_actual_shape_not_generalised():
    """11 of 14 are one O group -- calling this a uniform threshold problem would over-generalise."""
    s = json.loads(O_PROBE.read_text(encoding="utf-8"))["sub_threshold_shape"]
    assert max(s["O_group_counts"].values()) > sum(s["O_group_counts"].values()) / 2
    assert s["coverage_max"] < 80.0, "all sub-threshold hits must be below the deployed coverage cut"
    assert s["identity_median"] > 95.0, "the point is near-perfect identity at partial coverage"


@pytest.mark.skipif(not O_PROBE.exists(), reason="O-probe artifact absent")
def test_it_REFUSES_to_assert_a_third_cause():
    """Two asserted causes were measured wrong; the third must stay explicitly open."""
    d = json.loads(O_PROBE.read_text(encoding="utf-8"))
    open_qs = d["what_is_still_NOT_established"]
    assert any("NOT measured" in q for q in open_qs)
    assert any("wrong-call" in q.lower() or "wrong call" in q.lower() for q in open_qs)
    assert any("Infantis" in q for q in open_qs), "the unexplained within-serovar split must be kept"
