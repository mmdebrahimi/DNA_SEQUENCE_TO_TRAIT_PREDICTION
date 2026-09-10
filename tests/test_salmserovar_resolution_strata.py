"""Pins for the resolution-route stratification of the serovar cell's published accuracy.

The measurement answers "of the 168 hits, how many were uniquely resolved?" -- a question the headline
could not answer, which is the denominator error this repo has hit before in other costumes. These
tests pin the classifier's logic and the refusal that makes the method sound; the numbers themselves
are pinned against the committed artifact so a silent recompute cannot drift them.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.salmserovar_resolution_strata import (  # noqa: E402
    ARBITRARY, EXACT, FALLBACK, UNRESOLVED, classify_route, main, split_formula,
)

ARTIFACT = ROOT / "wiki" / "salmserovar_resolution_strata_2026-09-09.json"
TABLE = {("4", "i", "1,2"): "Typhimurium", ("9", "g,m", "-"): "Enteritidis"}


# --- the classifier ---------------------------------------------------------------------------

def test_an_abstention_is_unresolved_not_a_route():
    assert classify_route("O?:r:1,5", None, False, TABLE) == UNRESOLVED
    assert classify_route("4:i:1,2", "-", False, TABLE) == UNRESOLVED


def test_an_exact_key_hit_is_the_exact_route():
    assert classify_route("4:i:1,2", "Typhimurium", False, TABLE) == EXACT


def test_a_call_that_depends_on_the_winner_is_arbitrary_even_when_the_key_is_present():
    """The measured fact outranks the table, which records only winners and cannot say on its own
    which keys were contested -- reading it as authoritative would classify every arbitrary winner
    as a clean unique resolution."""
    assert classify_route("4:i:1,2", "Typhimurium", True, TABLE) == ARBITRARY


def test_a_missing_exact_key_falls_back():
    assert classify_route("8:r:l,w", "Goldcoast", False, TABLE) == FALLBACK


def test_an_unresolved_axis_is_never_used_as_a_lookup_key():
    """`O?`/`H?` mean the axis did not resolve. Treating them as literal keys would silently look up
    a formula the caller never actually determined."""
    assert split_formula("O?:r:1,5") is None
    assert split_formula("4:H?:1,2") is None
    assert split_formula("4:i:1,2") == ("4", "i", "1,2")
    assert split_formula(None) is None
    assert split_formula("4:i") is None


# --- the refusal that makes the method sound --------------------------------------------------

def test_it_refuses_when_the_two_policy_runs_disagree_on_any_formula(tmp_path: Path):
    """The whole method rests on the two runs differing ONLY in the name lookup. If an antigenic
    formula moved, a serovar diff is not attributable to the ambiguity policy and the strata are
    fiction -- so this must refuse rather than report, and the refusal must be non-vacuous."""
    a = tmp_path / "a.jsonl"
    b = tmp_path / "b.jsonl"
    a.write_text(json.dumps({"asm_acc": "X", "serovar": "Typhimurium", "formula": "4:i:1,2"}) + "\n",
                 encoding="utf-8")
    b.write_text(json.dumps({"asm_acc": "X", "serovar": "Typhimurium", "formula": "9:g,m:-"}) + "\n",
                 encoding="utf-8")
    cohort = tmp_path / "c.json"
    cohort.write_text(json.dumps({"isolates": [{"asm_acc": "X", "serovar_raw": "Typhimurium"}]}),
                      encoding="utf-8")
    rc = main(["--cohort", str(cohort), "--first", str(b), "--omit", str(a),
               "--out", str(tmp_path / "o.json")])
    assert rc == 2, "a formula mismatch must refuse, not silently report unattributable strata"
    assert not (tmp_path / "o.json").exists(), "a refused run must not leave an artifact behind"


# --- the committed result ----------------------------------------------------------------------

@pytest.mark.skipif(not ARTIFACT.exists(), reason="stratification artifact not present")
def test_the_strata_reconcile_with_the_published_headline():
    """The parts must sum to the whole. An arithmetic check catches a classifier that drops or
    double-counts an isolate, which no amount of reading the code would surface."""
    d = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    assert d["n"] == 200
    assert sum(sum(v.values()) for v in d["strata"].values()) == d["n"]
    assert sum(d["hits_by_route"].values()) == d["headline_reproduced"]["hit"] == 168
    assert d["headline_reproduced"]["accuracy"] == 0.84, "must reproduce the published number exactly"


@pytest.mark.skipif(not ARTIFACT.exists(), reason="stratification artifact not present")
def test_the_overwhelming_majority_of_hits_are_uniquely_resolved():
    """The headline survives the stratification: 165 of 168 hits do not depend on a tie-break."""
    d = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    assert d["hits_by_route"][EXACT] == 165
    assert d["hits_by_route"][ARBITRARY] == 3
    assert d["unique_resolution_accuracy"] == 0.825


@pytest.mark.skipif(not ARTIFACT.exists(), reason="stratification artifact not present")
def test_the_arbitrary_winner_trade_records_both_directions():
    """It gains 3 hits AND manufactures 5 confidently-wrong calls. Recording only the first is the
    framing the name-table memo used; this cell's own threshold work committed to the opposite
    asymmetry, so both numbers ship and the choice between them stays an acceptance-bar question."""
    t = json.loads(ARTIFACT.read_text(encoding="utf-8"))["arbitrary_winner_trade"]
    assert t["hits_gained_vs_omit"] == 3
    assert t["abstentions_converted_to_wrong_calls"] == 5
    assert t["isolates_whose_call_depends_on_the_winner"] == 8


@pytest.mark.skipif(not ARTIFACT.exists(), reason="stratification artifact not present")
def test_the_zero_fallback_stratum_is_a_measurement_not_an_unexercised_path():
    """A stratum of exactly zero is the shape a broken classifier produces, so the artifact must
    also show the fallback was REACHED and failed -- otherwise 'contributes nothing' would be
    indistinguishable from 'never ran'."""
    f = json.loads(ARTIFACT.read_text(encoding="utf-8"))["fallback_finding"]
    assert f["hits"] == 0
    assert f["reached_but_unresolved"] >= 1, (
        "if nothing ever reached the fallback, the zero says nothing about the fallback")
