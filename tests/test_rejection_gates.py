"""The rejection-gate screen: it must REFUSE rather than guess, and must reproduce the hand verdicts."""
from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

from dna_decode.eval.rejection_gates import (
    DECODER_GATES, INSUFFICIENT_DATA, L1_AMR_RS, L4_FORWARD_CONTINUOUS, LABEL_GATES,
    MAX_MODE_SHARE, MIN_DISTINCT_VALUES, NEEDS_HUMAN_EVIDENCE, NOT_APPLICABLE, PASS, TRIP,
    g1_circular_label, g2_study_equals_class, g3_sampling_defined_label, g6_phenotype_censoring,
    screen_candidate,
)

_SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "screen_candidate_gates.py"


def _runner():
    spec = importlib.util.spec_from_file_location("screen_candidate_gates", _SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


CLEAN_L1 = {
    "label_provenance_evidence": "broth microdilution MIC read in the lab",
    "label_is_measured": True,
    "label_semantics_evidence": "an MIC reading",
    "label_is_assay_reading": True,
    "largest_source_share": 0.31,
    "non_ecosystem_min_class_n": 60,
    "n_fetchable_assemblies": 120,
    "censored_fraction": 0.1,
    "provenance_field_populated_fraction": 0.98,
    "min_effective_lineages": 8,
    "loci_without_recorded_variant_fraction": 0.0,
    "off_panel_variant_fraction": 0.05,
}


# --- the anchor: the screen must re-derive what was decided by hand -------------------------------

def test_it_reproduces_both_committed_hand_verdicts():
    """If this fails, either the screen or a committed memo is wrong. Both are load-bearing."""
    mod = _runner()
    for key, spec in mod.CANDIDATES.items():
        ok, bad = mod.check_one(spec)
        assert ok, f"{key} disagrees with {spec['memo']}: {bad}"


def test_the_reproduction_check_is_not_vacuous():
    """A check that passes no matter what proves nothing. Corrupt an expectation -> it must fail."""
    mod = _runner()
    spec = dict(mod.HBV)
    spec["expected"] = {"G1": PASS}          # HBV's label is tool-derived; G1 genuinely TRIPS
    ok, bad = mod.check_one(spec)
    assert not ok and any("G1" in b for b in bad)


# --- the two judgement gates REFUSE; they never default to pass -----------------------------------

@pytest.mark.parametrize("gate", [g1_circular_label, g3_sampling_defined_label])
def test_a_judgement_gate_refuses_when_no_human_reading_is_supplied(gate):
    assert gate({}).verdict == NEEDS_HUMAN_EVIDENCE


@pytest.mark.parametrize("gate,text_key", [(g1_circular_label, "label_provenance_evidence"),
                                           (g3_sampling_defined_label, "label_semantics_evidence")])
def test_prose_alone_does_not_satisfy_a_judgement_gate(gate, text_key):
    """Supplying narrative without the explicit assertion must not read as a pass."""
    assert gate({text_key: "some words about the dataset"}).verdict == NEEDS_HUMAN_EVIDENCE


def test_a_screen_missing_a_judgement_reading_REFUSES_and_does_not_clear():
    ev = {k: v for k, v in CLEAN_L1.items() if k != "label_provenance_evidence"}
    res = screen_candidate("x", L1_AMR_RS, ev)
    assert res.verdict == "REFUSED" and "G1" in res.reason


def test_everything_measured_and_clean_clears():
    res = screen_candidate("x", L1_AMR_RS, CLEAN_L1)
    assert res.verdict == "CLEARS"


def test_a_clearing_screen_still_says_it_is_not_a_build_recommendation():
    d = screen_candidate("x", L1_AMR_RS, CLEAN_L1).as_dict()
    assert "NOT a build recommendation" in d["contract"]


# --- a trip is decisive, even with everything else unmeasured -------------------------------------

def test_a_tripped_gate_rejects_before_missing_measurements_are_considered():
    """HBV's shape: G1 alone settles it. The screen must not demand the rest first."""
    res = screen_candidate("x", L1_AMR_RS, {"label_provenance_evidence": "a rules engine's output",
                                            "label_is_measured": False})
    assert res.verdict == "REJECTED" and "G1" in res.reason


# --- the layer decides which gates apply ----------------------------------------------------------

def test_intended_layer_is_required_and_never_inferred():
    assert screen_candidate("x", "guess", CLEAN_L1).verdict == "REFUSED"


def test_g6_dispatches_on_layer_and_the_two_forms_genuinely_differ():
    """Same evidence, different layer -> different criterion. This is the PEAR lesson; if both layers
    returned the same thing the dispatch would be decorative."""
    ev = {"censored_fraction": 0.9, "mode_share": 0.10, "n_distinct_values": 40}
    assert g6_phenotype_censoring(ev, L1_AMR_RS).verdict == TRIP           # breakpoint-censored
    assert g6_phenotype_censoring(ev, L4_FORWARD_CONTINUOUS).verdict == PASS  # well-spread continuous


def test_g6_on_a_continuous_assay_is_OPEN_not_passed_when_degeneracy_is_unscreened():
    """The trap this exists for: an unscreened censored assay posts the BEST number (CcdB, 79.3% tied)."""
    assert g6_phenotype_censoring({}, L4_FORWARD_CONTINUOUS).verdict == INSUFFICIENT_DATA


def test_g6_trips_on_a_degenerate_continuous_assay():
    r = g6_phenotype_censoring({"mode_share": 0.793, "n_distinct_values": 8}, L4_FORWARD_CONTINUOUS)
    assert r.verdict == TRIP


def test_the_l4_degeneracy_bars_match_the_shipped_gate_they_were_taken_from():
    """Drift guard: these are copies of forward_inverse_roundtrip's constants, not new thresholds."""
    src = (Path(__file__).resolve().parents[1] / "scripts" / "forward_inverse_roundtrip.py").read_text(
        encoding="utf-8")
    assert f"MAX_MODE_SHARE = {MAX_MODE_SHARE}" in src
    assert f"MIN_DISTINCT_VALUES = {MIN_DISTINCT_VALUES}" in src


# --- applicability before measurement (the bug --verify actually caught) --------------------------

def test_an_inapplicable_gate_says_so_rather_than_asking_for_a_measurement():
    """`insufficient_data` reads as 'go measure this'. For a gate that CANNOT apply that is misleading,
    and it is the defect the hand-verdict reproduction check found in G2."""
    r = g2_study_equals_class({"variation_is_constructed": True})
    assert r.verdict == NOT_APPLICABLE


def test_every_constructed_variation_gate_is_inapplicable_without_any_measurement():
    ev = {"label_provenance_evidence": "wet-lab growth", "label_is_measured": True,
          "label_semantics_evidence": "assay reading", "label_is_assay_reading": True,
          "variation_is_constructed": True, "genotype_defined_by_construction": True,
          "loci_without_recorded_variant_fraction": 0.0, "off_panel_variant_fraction": 0.0,
          "mode_share": 0.05, "n_distinct_values": 500}
    res = screen_candidate("x", L4_FORWARD_CONTINUOUS, ev)
    got = {g.gate: g.verdict for g in res.gates}
    for gate in ("G2", "G4", "G5", "G7", "G8"):
        assert got[gate] == NOT_APPLICABLE, f"{gate} -> {got[gate]}"
    assert res.verdict == "CLEARS"


def test_all_ten_gates_run_every_time():
    res = screen_candidate("x", L1_AMR_RS, CLEAN_L1)
    assert [g.gate for g in res.gates] == list(LABEL_GATES) + list(DECODER_GATES)
