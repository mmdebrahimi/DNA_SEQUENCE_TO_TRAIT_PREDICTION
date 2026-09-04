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
