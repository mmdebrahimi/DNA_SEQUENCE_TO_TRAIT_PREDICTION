"""Tests for the CYP4F2 concordance harness + the panel-derived AF audit.

Two load-bearing pins:
  * the MECHANISTIC IDENTITY (star-name agreement == complement of "carries *2"), which is what makes the
    naming gap a complete explanation rather than a partial one;
  * the AF audit's REFUSAL behaviour -- a claim with no local panel must never count as a pass, and an
    audit that measured nothing must not print a clean report.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts.cyp4f2_getrm_concordance import (
    ASSERTED_AF,
    INVISIBLE_TO_CALLER,
    allele_frequencies,
    carries_invisible_allele,
    star3_count,
)
from scripts.pgx_af_panel_audit import CLAIMS, audit

REPO = Path(__file__).resolve().parent.parent
ARTIFACT = REPO / "wiki" / "cyp4f2_getrm_concordance_2026-09-24.json"
AF_ARTIFACT = REPO / "wiki" / "pgx_af_panel_audit_2026-09-24.json"
BAR = REPO / "wiki" / "cyp4f2_star_acceptance_bar.json"


# --- truth parsing -------------------------------------------------------------------------------

@pytest.mark.parametrize("dip,expected", [
    ("*1/*1", 0),
    ("*1/*3", 1),
    ("*3/*3", 2),
    ("(*2)/*3", 1),          # *2 haplotype carries no *3
    ("*1/(*2)", 0),
    ("*1/*2", 0),
])
def test_star3_count(dip, expected):
    assert star3_count(dip) == expected


def test_unparseable_truth_returns_none_not_zero():
    """A parse failure defaulted to 0 would score an unparsed row as a *3 non-carrier, inflating
    specificity. It must return None so the caller can exclude and COUNT it."""
    assert star3_count("") is None
    assert star3_count("*1/*2/*3") is None
    assert star3_count("nothing") is None


def test_star2_is_the_only_invisible_allele():
    """The caller reads exactly one SNP (rs2108622 -> *3). *1 and *3 are nameable; *2 is not."""
    assert set(INVISIBLE_TO_CALLER) == {"*2"}
    assert carries_invisible_allele("*1/*3") is False
    assert carries_invisible_allele("*3/*3") is False
    assert carries_invisible_allele("(*2)/*3") is True
    assert carries_invisible_allele("*1/*2") is True


# --- allele-frequency math -----------------------------------------------------------------------

def test_allele_frequencies_basic():
    dos = {"a": 0, "b": 1, "c": 2, "d": 2}
    pop = {"a": "EUR", "b": "EUR", "c": "EAS", "d": "EAS"}
    af = allele_frequencies(dos, pop)
    assert af["EUR"]["alt_freq"] == pytest.approx(0.25)      # 1 of 4 alleles
    assert af["EAS"]["alt_freq"] == pytest.approx(1.0)       # 4 of 4
    assert af["ALL"]["alt_freq"] == pytest.approx(0.625)     # 5 of 8


def test_allele_frequencies_ignores_samples_without_a_population():
    """A sample missing from the ped must not silently land in some bucket."""
    af = allele_frequencies({"a": 2, "unknown": 0}, {"a": "EUR"})
    assert af["EUR"]["n_samples"] == 1
    assert af["ALL"]["n_samples"] == 1


# --- artifact pins -------------------------------------------------------------------------------

@pytest.mark.skipif(not ARTIFACT.exists(), reason="artifact not built")
def test_artifact_headline():
    d = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    assert d["n_scored"] == 64
    assert d["n_parse_failures"] == 0
    assert d["star_verdict"] == "STAR_PROXY_PARTIAL"
    assert d["af_verdict"] == "AF_CLAIM_WRONG"


@pytest.mark.skipif(not ARTIFACT.exists(), reason="artifact not built")
def test_the_mechanistic_identity_is_exact():
    """star_name_agrees == NOT carries *2, for EVERY sample -- zero exceptions. Same identity that held
    on SLCO1B1; if it breaks, the naming gap has a second cause and the memo is incomplete."""
    d = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    exceptions = [r for r in d["rows"]
                  if r["star_name_agrees"] == r["carries_allele_invisible_to_caller"]]
    assert exceptions == []


@pytest.mark.skipif(not ARTIFACT.exists(), reason="artifact not built")
def test_the_registered_prediction_held():
    """The bar predicted ~0.81 with exactly 12 disagreements -- a band deliberately different from
    SLCO1B1's 0.398, so it could have been falsified in either direction."""
    d = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    assert d["registered_prediction_held"] is True
    assert d["n_scored"] - d["star_name_axis"]["exact_name_agreement"] == 12


@pytest.mark.skipif(not ARTIFACT.exists(), reason="artifact not built")
def test_function_axis_is_definitional_and_labelled_as_such():
    d = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    c = d["function_axis_control"]
    assert c["zygosity_exact"] == c["n"] == 64
    assert c["IS_A_CONTROL_NOT_VALIDATION"] is True


@pytest.mark.skipif(not ARTIFACT.exists(), reason="artifact not built")
def test_no_safety_miss_and_scope_is_stated():
    d = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    assert d["safety_axis"]["n_called_normal_function_anyway"] == 0
    # the *2-functional-status scope limit must travel WITH the safety number
    assert "*2" in d["safety_axis"]["scope"]


@pytest.mark.skipif(not ARTIFACT.exists(), reason="artifact not built")
def test_the_run_is_not_vacuous():
    d = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    assert d["safety_axis"]["n_truth_star3_carriers"] >= 10
    assert d["star_name_axis"]["n_carrying_allele_invisible_to_caller"] >= 10
    assert 0 < d["star_name_axis"]["exact_name_agreement"] < d["n_scored"]


@pytest.mark.skipif(not ARTIFACT.exists(), reason="artifact not built")
def test_the_eas_af_claim_is_recorded_as_wrong_with_the_complement_note():
    """The measured EAS *3 frequency (~0.213) and the asserted 0.79 differ by ~0.58; 1-0.213 = 0.787
    matches the assertion to ~0.003. Pin the numbers so the correction cannot silently drift back."""
    d = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    eas = d["af_corroboration_axis"]["checks"]["EAS"]
    assert eas["within_tolerance"] is False
    assert eas["measured"] < 0.30
    assert d["af_corroboration_axis"]["checks"]["EUR"]["within_tolerance"] is True


@pytest.mark.skipif(not BAR.exists(), reason="bar not present")
def test_bar_frozen_with_a_falsifiable_band():
    d = json.loads(BAR.read_text(encoding="utf-8"))
    assert d["frozen_before_reading_any_result"] is True
    # all three star verdict bands must exist, else the prediction could not have been wrong
    for k in ("star_proxy_not_reference_agreeing", "star_proxy_partial", "star_proxy_agrees"):
        assert k in d["verdicts"]
    assert d["registered_predictions"]["predicted_disagreement_count"] == 12


# --- the AF panel audit: refusal behaviour is the point ------------------------------------------

def test_af_audit_reports_unmeasurable_rather_than_passing():
    """A claim whose variant has no local panel must be UNMEASURABLE, never OK. Silently passing it is
    how an unchecked wrong number survives."""
    rep = audit()
    unmeasurable = [c for c in rep["claims"] if c["status"] == "UNMEASURABLE_NO_LOCAL_PANEL"]
    assert unmeasurable, "expected ABCG2/NUDT15 to have no local panel"
    for c in unmeasurable:
        assert c["checks"] == {}          # no verdict may be manufactured
        assert c["reason"]


def test_af_audit_catches_the_cyp4f2_eas_claim():
    rep = audit()
    c4 = [c for c in rep["claims"] if c["gene"] == "CYP4F2"][0]
    assert c4["status"] == "CLAIM_OUT_OF_BAND"
    assert c4["checks"]["EAS"]["within_tolerance"] is False
    assert c4["checks"]["EUR"]["within_tolerance"] is True
    assert rep["verdict"] == "CLAIMS_OUT_OF_BAND"


def test_af_audit_flags_the_complement_signature():
    """1-ALT matching the assertion better than ALT is the allele-polarity signature; it must be
    surfaced as a field, not left for a reader to notice."""
    rep = audit()
    c4 = [c for c in rep["claims"] if c["gene"] == "CYP4F2"][0]
    assert c4["checks"]["EAS"]["matches_complement_better"] is True
    assert c4["checks"]["EUR"]["matches_complement_better"] is False


def test_af_audit_is_non_vacuous_the_other_claims_are_accurate():
    """3 of 4 measurable claims land within a hair. This is the CONTROL: if the AF computation were
    broken, all four would be off, and 'CYP4F2 is wrong' would be a pipeline artifact instead."""
    rep = audit()
    ok = [c for c in rep["claims"] if c["status"] == "OK"]
    assert len(ok) >= 3
    for c in ok:
        for p, chk in c["checks"].items():
            assert abs(chk["delta"]) < 0.02, f"{c['gene']} {p} delta {chk['delta']}"


def test_af_audit_claim_table_points_at_where_each_claim_lives():
    """A failing claim is only actionable if the audit says which file asserts it."""
    for c in CLAIMS:
        assert c["asserted_in"] and ("/" in c["asserted_in"] or ".py" in c["asserted_in"])
        assert c["asserted"], "a claim with no asserted frequency cannot be audited"
