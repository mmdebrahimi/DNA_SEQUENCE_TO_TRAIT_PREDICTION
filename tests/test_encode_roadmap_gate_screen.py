"""Tests for the ENCODE regulatory-activity gate screen.

The load-bearing one is `test_every_insufficient_data_gate_names_WHY`. An `insufficient_data` with
no reason reads as "go measure this", which is the honest answer for exactly one of the six gates
here and the WRONG answer for the other five, where the quantity is UNDEFINED for this candidate
rather than unmeasured. A screen whose gaps are unexplained is a shrug wearing a verdict.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "wiki" / "encode_roadmap_gate_screen_2026-11-24.json"

pytestmark = pytest.mark.skipif(not ARTIFACT.exists(), reason="gate-screen artifact absent")

ALL_GATES = ("G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8", "G9", "G10")


@pytest.fixture(scope="module")
def rec() -> dict:
    return json.loads(ARTIFACT.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def gates(rec) -> dict:
    return rec["per_gate_detail"]["gates"]


# ----------------------------------------------------------------- the refusal is the result
def test_the_whole_candidate_screen_REFUSED(rec):
    """Neither declared layer fits, and the honest record of that is a refusal."""
    assert rec["whole_candidate_result"]["verdict"] == "REFUSED"


def test_the_refusal_is_not_recorded_as_a_rejection(rec):
    """A refused screen neither cleared nor rejected the label question."""
    note = rec["whole_candidate_result"]["what_this_does_NOT_mean"]
    assert "NOT a rejection" in note


def test_no_false_layer_was_declared(rec):
    """Declaring L4 would flip five gates that dispatch on `variation_is_constructed`."""
    assert rec["whole_candidate_result"]["intended_layer"] not in (
        "L1_AMR_RS", "L4_forward_continuous")


def test_the_structural_declarations_are_honest(rec):
    """Natural variation must be declared natural, or the wrong gate branches run."""
    g = rec["per_gate_detail"]["gates"]
    # G2/G7 reaching a measured branch at all proves variation_is_constructed was False:
    # every one of those gates short-circuits to not_applicable when it is True.
    assert g["G2"]["verdict"] != "not_applicable"
    assert g["G7"]["verdict"] != "not_applicable"


# ----------------------------------------------------------------- gate coverage
def test_all_ten_gates_are_reported(gates):
    assert set(gates) == set(ALL_GATES), "a gate missing from the record is a silent gap"


def test_no_gate_tripped(rec):
    assert rec["per_gate_detail"]["summary"]["trip"] == []


def test_the_four_answerable_gates_passed(rec):
    """G1 and G3 by assertion, G2 and G7 from measurements taken this session."""
    assert set(rec["per_gate_detail"]["summary"]["pass"]) == {"G1", "G2", "G3", "G7"}


def test_every_insufficient_data_gate_names_WHY(rec, gates):
    """The anti-shrug test: an unexplained gap is indistinguishable from laziness."""
    for g in rec["per_gate_detail"]["summary"]["insufficient_data"]:
        d = gates[g]
        assert d.get("withheld_input"), f"{g}: no withheld input named"
        why = d.get("why_not_supplied") or ""
        assert len(why) > 80, f"{g}: reason too thin to audit"


def test_undefined_is_distinguished_from_unmeasured(gates):
    """Four gates are undefined for this candidate; one is merely unmeasured. Not the same."""
    # Case-insensitive on purpose: the assertion is about the gate declaring its quantity
    # non-existent, not about how that word is capitalised.
    undefined = [g for g in ("G4", "G9", "G10")
                 if "undefined" in gates[g]["why_not_supplied"].lower()]
    assert len(undefined) == 3, "a gate whose quantity does not exist must say so"
    assert "unmeasured" in gates["G6"]["why_not_supplied"], (
        "G6 is genuinely measurable and must not be filed with the undefined ones")


def test_the_two_tempting_false_substitutions_are_refused_by_name(gates):
    """60,055 bigWigs are not assemblies; biosample terms are not lineages."""
    assert "60,055" in gates["G5"]["why_not_supplied"]
    assert "not an assembly count" in gates["G5"]["why_not_supplied"]
    assert "not an effective" in gates["G8"]["why_not_supplied"]


# ----------------------------------------------------------------- judgments + measurements
def test_both_judgment_gates_were_asserted_with_evidence(rec):
    j = rec["judgment_gates_supplied"]
    for g in ("G1", "G3"):
        assert j[g]["assertion"] and len(j[g]["basis"]) > 60
    assert "theatre" in j["note"], "the reason defaults are refused must be recorded"


def test_each_measured_input_records_how_it_was_measured(rec):
    for name, m in rec["measured_this_session"].items():
        assert m["how"], f"{name}: no method recorded"
        assert isinstance(m["value"], (int, float))


def test_the_g2_pass_carries_its_ill_posedness_caveat(rec):
    """G2 asks whether SOURCE predicts CLASS; there is no class here."""
    note = rec["measured_this_session"]["largest_source_share"]["scope_note"]
    assert "ILL-POSED" in note
    assert "NOT as" in note, "the PASS must not read as 'the confound was tested and absent'"


# ----------------------------------------------------------------- the regime cell
def test_the_regime_cell_is_recorded_with_the_MOLECULAR_endpoint(rec):
    q = rec["regime_map"]["correct_cell"]["queried"]
    assert "endpoint=molecular" in q


def test_the_works_verdict_carries_an_explicit_do_not_over_read(rec):
    """Its evidence is HIV variant effect in one virus, not chromatin signal."""
    cell = rec["regime_map"]["correct_cell"]
    assert cell["result"]["verdict"] == "WORKS"
    warn = cell["DO_NOT_OVER_READ"]
    assert "over-read" in warn and "HIV" in warn


def test_the_plans_wrong_cell_is_recorded_as_a_correction(rec):
    c = rec["regime_map"]["cell_the_plan_claimed"]
    assert c["verdict"] == "REQUIRES_DECONFOUNDING"
    assert "WRONG one" in c["correction"]
    assert "CAUTIOUS" in c["correction"], (
        "the direction of an error is part of the correction -- this one under-claimed")


def test_regime_py_was_not_edited(rec):
    assert rec["regime_map"]["regime_py_NOT_edited"] is True


# ----------------------------------------------------------------- scope
def test_the_verdict_scope_bounds_the_label_question_only(rec):
    s = rec["verdict_scope"]
    assert "LABEL question ONLY" in s
    assert "not a build recommendation" in s


def test_honest_limits_record_that_the_screen_did_not_run(rec):
    assert any("did not run" in lim for lim in rec["honest_limits"])
