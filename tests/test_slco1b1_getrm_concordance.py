"""Tests for the SLCO1B1 star-name vs function concordance harness.

The load-bearing test is the MECHANISTIC IDENTITY pin: star-name agreement is exactly the complement of
"carries an allele the caller cannot see". If that ever stops holding, the naming gap has a second cause
and the memo's explanation is incomplete.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts.slco1b1_getrm_concordance import (
    INVISIBLE_TO_CALLER,
    REDUCED_FUNCTION,
    carries_invisible_allele,
    haplotypes,
    norm_star,
    reduced_count,
)

REPO = Path(__file__).resolve().parent.parent
ARTIFACT = REPO / "wiki" / "slco1b1_getrm_concordance_2026-09-22.json"
BAR = REPO / "wiki" / "slco1b1_star_acceptance_bar.json"


@pytest.mark.parametrize("dip,expected", [
    ("*1/*1", 0),
    ("*1/*15", 1),          # *15 = 388 + 521C -> reduced
    ("*15/*15", 2),
    ("*5/*15", 2),
    ("*1/*17", 1),          # *17 = *15 + promoter -> reduced
    ("*1B/*1B", 0),         # 388A>G ONLY -> NOT reduced
    ("*1A/*1B", 0),
    ("*1/*14", 0),
    ("*1/*21", 0),
])
def test_reduced_function_count(dip, expected):
    assert reduced_count(dip) == expected


def test_only_521c_bearing_alleles_count_as_reduced():
    """*1B/*14/*21 carry 388A>G or other variants but NOT 521T>C. Counting them reduced would invert the
    biology and manufacture a safety signal that is not there."""
    assert set(REDUCED_FUNCTION) == {"*5", "*15", "*17"}
    for a in ("*1", "*1A", "*1B", "*14", "*21"):
        assert reduced_count(f"*1/{a}") == 0


def test_unparseable_truth_returns_none_not_zero():
    assert reduced_count("") is None
    assert reduced_count("*1/*5/*15") is None
    assert reduced_count("nothing") is None
    assert haplotypes("*1/") is None


def test_norm_star_is_order_insensitive():
    """*1/*5 and *5/*1 are the same diplotype; a name comparison must not fail on ordering alone."""
    assert norm_star("*1/*5") == norm_star("*5/*1")
    assert norm_star("*1/*15") != norm_star("*1/*5")


def test_carries_invisible_allele():
    assert carries_invisible_allele("*1/*1") is False
    assert carries_invisible_allele("*1/*5") is False      # *5 IS 521C -- visible to the caller
    for a in INVISIBLE_TO_CALLER:
        assert carries_invisible_allele(f"*1/{a}") is True


def test_star5_is_not_invisible_but_star15_is():
    """The distinction the whole memo rests on: *5 is defined by 521T>C (which we genotype), *15 is
    *5 PLUS 388A>G (which we do not)."""
    assert "*5" not in INVISIBLE_TO_CALLER
    assert "*15" in INVISIBLE_TO_CALLER


# --- artifact pins -------------------------------------------------------------------------------

@pytest.mark.skipif(not ARTIFACT.exists(), reason="artifact not built")
def test_artifact_headline():
    d = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    assert d["n_scored"] == 88
    assert d["n_parse_failures"] == 0
    assert d["verdict"] == "STAR_PROXY_NOT_REFERENCE_AGREEING"


@pytest.mark.skipif(not ARTIFACT.exists(), reason="artifact not built")
def test_the_mechanistic_identity_is_exact():
    """star_name_agrees == NOT carries_allele_invisible_to_caller, for EVERY sample. This is what makes
    the 388A>G blindness a COMPLETE explanation of the naming gap rather than a partial one."""
    d = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    exceptions = [r for r in d["rows"]
                  if r["star_name_agrees"] == r["carries_allele_invisible_to_caller"]]
    assert exceptions == []


@pytest.mark.skipif(not ARTIFACT.exists(), reason="artifact not built")
def test_function_axis_is_exactly_definitional():
    """Perfect 1:1 dosage agreement. Pinned NOT as a quality claim but because any deviation would mean
    the axis is not the self-comparison the memo says it is."""
    d = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    c = d["function_axis_control"]
    assert c["zygosity_exact"] == c["n"] == 88
    assert c["IS_A_CONTROL_NOT_VALIDATION"] is True


@pytest.mark.skipif(not ARTIFACT.exists(), reason="artifact not built")
def test_no_safety_miss():
    """A reduced-function carrier called Normal Function is the dangerous direction for statin myopathy."""
    d = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    assert d["safety_axis"]["n_called_normal_function_anyway"] == 0
    assert d["safety_axis"]["n_truth_reduced_function_carriers"] >= 10


@pytest.mark.skipif(not ARTIFACT.exists(), reason="artifact not built")
def test_the_run_is_not_vacuous():
    """Both axes need real positives: reduced-function carriers for safety, invisible-allele carriers for
    the naming axis. A run with neither would report flattering numbers about nothing."""
    d = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    assert d["star_name_axis"]["n_carrying_allele_invisible_to_caller"] >= 10
    # and the naming gap must not be trivially 0 or 1 -- it has to discriminate
    assert 0 < d["star_name_axis"]["exact_name_agreement"] < d["n_scored"]


@pytest.mark.skipif(not BAR.exists(), reason="bar not present")
def test_bar_frozen_with_predictions_and_both_directions():
    d = json.loads(BAR.read_text(encoding="utf-8"))
    assert d["frozen_before_reading_any_result"] is True
    assert "star_proxy_agrees" in d["verdicts"]                     # a bar that can only fail is not a bar
    assert "star_proxy_not_reference_agreeing" in d["verdicts"]
    assert "registered_predictions" in d
