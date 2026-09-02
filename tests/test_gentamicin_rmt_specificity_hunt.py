"""Pins for the `rmt` specificity hunt: the field parser, the imported pattern, and the artifact's logic.

What is pinned is the machinery that decides whether 60 counter-examples retract a deployed rule. The
correlation-free parts each had a way to be silently wrong.
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
HUNT = ROOT / "wiki" / "gentamicin_rmt_specificity_hunt.json"
CONTROL = ROOT / "wiki" / "gentamicin_rmt_project_control.json"


def _mod(name: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / f"{name}.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


# --- the AMR_genotypes parser: two real traps in the live field ------------------------------------

def test_amr_genotype_parser_strips_the_wrapping_quotes():
    """The cell is quote-wrapped, so the quote rides on the FIRST and LAST tokens -- the same trap that
    made `AST_phenotypes` unmatchable at either end before it was fixed."""
    p = _mod("gentamicin_rmt_specificity_hunt").parse_amr_genotypes
    syms = p('"aac(3)-IId,rmtB1,npmA"')
    assert "aac(3)-IId" in syms, "leading quote not stripped"
    assert "npmA" in syms, "trailing quote not stripped"


def test_amr_genotype_parser_strips_the_equals_annotation():
    """Symbols carry =PARTIAL / =POINT / =MISTRANSLATION / =PARTIAL_END_OF_CONTIG suffixes."""
    p = _mod("gentamicin_rmt_specificity_hunt").parse_amr_genotypes
    assert p("rmtB=PARTIAL,armA=MISTRANSLATION") == {"rmtB", "armA"}


def test_amr_genotype_parser_treats_null_as_empty():
    p = _mod("gentamicin_rmt_specificity_hunt").parse_amr_genotypes
    assert p("NULL") == set() and p(None) == set() and p("") == set()


# --- the pattern must be the DEPLOYED one, not a retyped copy --------------------------------------

def test_the_hunt_imports_the_deployed_rescue_pattern():
    from dna_decode.eval.amr_rules import DRUG_RULE
    assert _mod("gentamicin_rmt_specificity_hunt").RESCUE == DRUG_RULE["gentamicin"]["symbol_rescue"]


def test_arma_is_not_rescued_because_the_frozen_rule_already_counts_it():
    """armA is filed by AMRFinder under Subclass GENTAMICIN, so it was never the gap; rescuing it would
    double-count and would misattribute the v2 sensitivity gain."""
    m = _mod("gentamicin_rmt_specificity_hunt")
    assert not m.RESCUE_RE.match("armA")
    assert m.RESCUE_RE.match("rmtB1") and m.RESCUE_RE.match("npmA")


# --- the committed artifacts stay self-consistent --------------------------------------------------

@pytest.mark.skipif(not HUNT.is_file(), reason="hunt artifact not present")
def test_the_hunt_run_was_complete_or_says_it_was_not():
    """A capped or error-bearing run cannot be quoted as a specificity bound."""
    d = json.loads(HUNT.read_text(encoding="utf-8"))
    assert d["complete"] is True
    assert not d["errors"]
    assert d["n_labelled_isolates"] > 10000


@pytest.mark.skipif(not HUNT.is_file(), reason="hunt artifact not present")
def test_the_hunt_searched_a_population_not_our_cache():
    """The whole point of the inversion: the previous hunt was cache-bounded and could not reach these."""
    d = json.loads(HUNT.read_text(encoding="utf-8"))
    assert "NOT restricted to the local" in d["population"]
    assert len(d["rmt_S_accessions"]) == d["rmt_carrier_counts"]["S"]


@pytest.mark.skipif(not CONTROL.is_file(), reason="control artifact not present")
def test_the_control_uses_an_undisputed_determinant_as_its_yardstick():
    """aac(3) is what the FROZEN rule already counted, so it tests the project's labels without
    reference to the rescue under scrutiny."""
    d = json.loads(CONTROL.read_text(encoding="utf-8"))
    assert d["aac3_R_rate_inside"] is not None and d["aac3_R_rate_outside"] is not None
    assert d["verdict"] in ("LABEL_ARTIFACT", "SPECIFIC_TO_RMT", "INCONCLUSIVE")


@pytest.mark.skipif(not CONTROL.is_file(), reason="control artifact not present")
def test_the_label_artifact_verdict_rests_on_a_real_divergence():
    """Non-vacuity: the verdict must be backed by the numbers, not just asserted."""
    d = json.loads(CONTROL.read_text(encoding="utf-8"))
    if d["verdict"] == "LABEL_ARTIFACT":
        assert d["aac3_R_rate_inside"] < 0.5 < 0.8 < d["aac3_R_rate_outside"]


@pytest.mark.skipif(not (HUNT.is_file() and CONTROL.is_file()), reason="artifacts not present")
def test_the_counter_examples_do_not_retract_the_rule_only_because_they_are_one_submission():
    """The load-bearing arithmetic: every S carrier is inside the artifact project, and every carrier
    outside it is R. If that ever stops holding, the conclusion must be revisited."""
    h = json.loads(HUNT.read_text(encoding="utf-8"))
    c = json.loads(CONTROL.read_text(encoding="utf-8"))
    assert c["inside_totals"]["rmt"].get("R", 0) == 0
    assert c["outside_totals"]["rmt"].get("S", 0) == 0
    assert h["rmt_carrier_counts"]["R"] > 100


@pytest.mark.skipif(not HUNT.is_file(), reason="hunt artifact not present")
def test_the_artifact_records_that_specificity_stays_untested():
    """A label artifact means the counter-examples cannot TEST the rule -- never that it is vindicated."""
    d = json.loads(HUNT.read_text(encoding="utf-8"))
    assert any("tool-derived" in lim for lim in d["honest_limits"])
