"""Tests for the UGT1A1 tag-vs-GeT-RM concordance harness.

Pure-logic tests (parsing, confusion math, repeat-allele mapping) plus a NON-VACUITY pin that the
committed artifact's headline numbers still reproduce -- so the measurement cannot silently drift.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts.ugt1a1_getrm_concordance import (
    REPEAT_ALLELE_TO_STAR,
    _dose,
    confusion,
    truth_reduced_count,
    truth_star_count,
)

REPO = Path(__file__).resolve().parent.parent
ARTIFACT = REPO / "wiki" / "ugt1a1_getrm_concordance_2026-09-22.json"
BAR = REPO / "wiki" / "ugt1a1_tag_acceptance_bar.json"


# --- truth parsing -------------------------------------------------------------------------------

@pytest.mark.parametrize("diplotype,expected", [
    ("*1/*1", 0),
    ("*1/*28", 1),
    ("*28/*28", 2),
    ("*60/*28", 1),
    ("(*28 + *60)/(*28 + *60)", 2),          # compound haplotypes, both carry *28
    ("(*28 + *60)/(*37 + *60)", 1),          # one *28 haplotype, one *37
    ("*28/(*60)", 1),                        # parenthesised single allele
    ("(*28)/(*28 + *60)", 2),                # the NA20509 truth string
    ("*1/*36", 0),
])
def test_truth_star28_count(diplotype, expected):
    assert truth_star_count(diplotype, "*28") == expected


def test_unparseable_truth_returns_none_not_zero():
    """A parse failure must NOT be defaulted to 0 -- that would score an unparsed row as a true
    negative, silently inflating specificity. This is the accounting trap the harness guards."""
    assert truth_star_count("", "*28") is None
    assert truth_star_count("no-alleles-here", "*28") is None
    assert truth_star_count("*1/*28/*36", "*28") is None      # not two haplotypes
    assert truth_star_count("*1/", "*28") is None             # empty haplotype


def test_star_match_is_exact_not_substring():
    """*28 must not match *2 or *280-style tokens, and *3 must not match *37."""
    assert truth_star_count("*2/*2", "*28") == 0
    assert truth_star_count("*1/*37", "*3") == 0


# --- the post-hoc reduced-function framing -------------------------------------------------------

@pytest.mark.parametrize("diplotype,expected", [
    ("*1/*1", 0),
    ("*1/*28", 1),                           # TA7, reduced
    ("*1/*37", 1),                           # TA8, also reduced -- the whole point of the arm
    ("(*37)/*60", 1),
    ("(*28 + *60)/(*37 + *60)", 2),          # both haplotypes reduced-function
    ("*1/*36", 0),                           # TA5 is INCREASED function, must NOT count
    ("*36/*36", 0),
])
def test_truth_reduced_function_count(diplotype, expected):
    assert truth_reduced_count(diplotype) == expected


def test_reduced_function_excludes_star36():
    """*36 (TA5) is a SHORTER repeat = increased expression. Counting it as reduced would invert the
    biology and manufacture agreement."""
    assert truth_reduced_count("*36/*36") == 0
    assert truth_reduced_count("*28/*36") == 1


# --- repeat-allele mapping -----------------------------------------------------------------------

def test_repeat_indel_to_star_mapping():
    """Reference is A(TA)6TAA = *1. +1 TA -> *28 (TA7); +2 TA -> *37 (TA8); -1 TA -> *36 (TA5)."""
    assert REPEAT_ALLELE_TO_STAR[("C", "CAT")] == "*28"
    assert REPEAT_ALLELE_TO_STAR[("C", "CATAT")] == "*37"
    assert REPEAT_ALLELE_TO_STAR[("CAT", "C")] == "*36"


def test_repeat_insert_lengths_are_consistent_with_ta_dinucleotide_steps():
    """Each mapped indel must change length by a whole number of TA dinucleotides (2 bp)."""
    for (ref, alt), star in REPEAT_ALLELE_TO_STAR.items():
        assert abs(len(alt) - len(ref)) % 2 == 0, f"{star} is not a whole-TA step"


# --- genotype dosage -----------------------------------------------------------------------------

@pytest.mark.parametrize("gt,expected", [
    ("0|0", 0), ("0|1", 1), ("1|0", 1), ("1|1", 2),
    ("0/1", 1), ("1/1", 2),
    ("0|1:99", 1),          # trailing FORMAT fields ignored
    (".|.", None), ("./.", None), (".", None),
])
def test_dose(gt, expected):
    assert _dose(gt) == expected


# --- confusion math ------------------------------------------------------------------------------

def test_confusion_counts_and_rates():
    pairs = [(1, 1), (2, 2), (0, 0), (1, 0), (0, 1)]
    c = confusion(pairs)
    assert (c["tp"], c["fn"], c["fp"], c["tn"]) == (2, 1, 1, 1)
    assert c["sensitivity"] == pytest.approx(2 / 3, abs=1e-4)
    assert c["ppv"] == pytest.approx(2 / 3, abs=1e-4)
    assert c["zygosity_exact"] == 3          # (1,1) (2,2) (0,0)


def test_confusion_zygosity_is_stricter_than_carriership():
    """A caller that gets carriership right but zygosity wrong must NOT score as exact -- the CPIC
    activity score turns on copy number, so 1-vs-2 copies is a real Poor/Intermediate mis-assignment."""
    c = confusion([(2, 1)])
    assert c["tp"] == 1 and c["sensitivity"] == 1.0    # carriership: correct
    assert c["zygosity_exact"] == 0                    # zygosity: wrong


def test_confusion_rates_are_none_when_undefined():
    c = confusion([(0, 0)])
    assert c["sensitivity"] is None and c["ppv"] is None


# --- artifact non-vacuity pins -------------------------------------------------------------------

@pytest.mark.skipif(not ARTIFACT.exists(), reason="concordance artifact not built")
def test_committed_artifact_headline_numbers():
    d = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    assert d["n_scored"] == 65
    assert d["n_parse_failures"] == 0
    # both arms must be scored on the FULL set -- an arm scored on an easier subset is not comparable
    assert d["tag_arm"]["n"] == 65
    assert d["direct_repeat_arm"]["n"] == 65
    assert d["post_hoc_reduced_function_arm"]["n"] == 65


@pytest.mark.skipif(not ARTIFACT.exists(), reason="concordance artifact not built")
def test_the_run_is_not_vacuous_it_contains_real_positives():
    """A tag cannot be tested against absent positives. Pin that the scored set really does carry
    *28 carriers AND non-carriers, so sensitivity and specificity are both informative."""
    d = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    assert d["n_truth_star28_carriers"] >= 20
    assert d["tag_arm"]["tn"] >= 20


@pytest.mark.skipif(not ARTIFACT.exists(), reason="concordance artifact not built")
def test_star6_control_is_perfect_so_the_plumbing_is_validated():
    """*6 is a DIRECT SNP, not a tag. If this control were not exact, the tag/repeat numbers would be
    uninterpretable (a parsing artifact rather than a biological result)."""
    d = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    c = d["star6_control"]
    assert c["n"] == 65 and c["zygosity_exact"] == 65


@pytest.mark.skipif(not ARTIFACT.exists(), reason="concordance artifact not built")
def test_post_hoc_arm_is_labelled_post_hoc():
    """The reduced-function framing was added AFTER the bar was scored. It must never be presented as
    the frozen primary."""
    d = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    rf = d["post_hoc_reduced_function_arm"]
    assert rf["post_hoc"] is True and rf["not_in_frozen_bar"] is True


@pytest.mark.skipif(not BAR.exists(), reason="bar not present")
def test_bar_was_frozen_before_results_and_names_both_verdict_directions():
    d = json.loads(BAR.read_text(encoding="utf-8"))
    assert d["frozen_before_reading_any_result"] is True
    # a bar that can only confirm is not a bar
    assert "tag_demote_and_disclose" in d["verdicts"]
    assert "tag_confirmed" in d["verdicts"]
