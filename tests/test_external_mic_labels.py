"""Unit tests for external-cohort MIC -> tiered R/S labeling."""
from __future__ import annotations

import json

import pytest

from dna_decode.data import external_mic_labels as eml


# --------------------------------------------------------------------------- #
# canonical_drug
# --------------------------------------------------------------------------- #
def test_canonical_drug_aliases():
    assert eml.canonical_drug("CIP") == "ciprofloxacin"
    assert eml.canonical_drug("cipro") == "ciprofloxacin"
    assert eml.canonical_drug("CRO") == "ceftriaxone"
    assert eml.canonical_drug("CN") == "gentamicin"
    assert eml.canonical_drug(" Gentamicin ") == "gentamicin"


def test_canonical_drug_reject_unknown():
    assert eml.canonical_drug("meropenem") is None
    assert eml.canonical_drug("amikacin") is None
    assert eml.canonical_drug("") is None


# --------------------------------------------------------------------------- #
# parse_mic_token
# --------------------------------------------------------------------------- #
def test_parse_mic_token_plain():
    assert eml.parse_mic_token("8") == 8.0
    assert eml.parse_mic_token(0.5) == 0.5


def test_parse_mic_token_censored():
    assert eml.parse_mic_token(">8") == 8.0
    assert eml.parse_mic_token("<=0.25") == 0.25
    assert eml.parse_mic_token(">=16") == 16.0


def test_parse_mic_token_units():
    assert eml.parse_mic_token("4 mg/L") == 4.0
    assert eml.parse_mic_token(">2 ug/mL") == 2.0


def test_parse_mic_token_na():
    assert eml.parse_mic_token("NA") is None
    assert eml.parse_mic_token("") is None
    assert eml.parse_mic_token(None) is None
    assert eml.parse_mic_token("garbage") is None


# --------------------------------------------------------------------------- #
# tier_for_isolate (real classify_tier; cipro clsi_r=2 clsi_s=0.5)
# --------------------------------------------------------------------------- #
def test_tier_high_r_plain():
    assert eml.tier_for_isolate(["8"], set(), "ciprofloxacin") == "HIGH_R"   # plain 8 >= 4*2


# --- operator-aware censoring (Step 2) ---
def test_tier_censored_lower_bound_high_r():
    # ">8" is a lower bound whose bound itself lands HIGH_R -> CENSORED_HIGH_R
    assert eml.tier_for_isolate([">8"], set(), "ciprofloxacin") == "CENSORED_HIGH_R"


def test_tier_censored_upper_bound_high_s():
    # "<=0.06" upper bound, bound <= clsi_s/4 -> CENSORED_HIGH_S
    assert eml.tier_for_isolate(["<=0.06"], set(), "ciprofloxacin") == "CENSORED_HIGH_S"


def test_tier_censored_lower_bound_midrange_excluded():
    # ">2" lower bound, bound is BORDERLINE not HIGH_R -> excluded (conservative)
    assert eml.tier_for_isolate([">2"], set(), "ciprofloxacin") == "CENSORED_EXCLUDED"


def test_tier_upper_bound_never_calls_r():
    # "<=8" upper bound at a high value must NOT be called R (the strip-to-bound bug)
    assert eml.tier_for_isolate(["<=8"], set(), "ciprofloxacin") == "CENSORED_EXCLUDED"


def test_tier_lower_bound_never_calls_s():
    # ">=0.06" lower bound at a low value must NOT be called S
    assert eml.tier_for_isolate([">=0.06"], set(), "ciprofloxacin") == "CENSORED_EXCLUDED"


def test_parse_mic_value_preserves_operator():
    assert eml.parse_mic_value(">8") == eml.MicValue(8.0, ">", ">8")
    assert eml.parse_mic_value("4") == eml.MicValue(4.0, "=", "4")
    assert eml.parse_mic_value("<=0.5 mg/L").operator == "<="
    assert eml.parse_mic_value("NA") is None


def test_censored_high_r_is_strict_label():
    res = eml.build_drug_labels({"GCA_x": [">16"]}, "ciprofloxacin")
    assert res["strict"]["GCA_x"] == "R"
    assert res["buckets"].get("CENSORED_HIGH_R") == 1


def test_tier_high_s():
    assert eml.tier_for_isolate(["0.06"], set(), "ciprofloxacin") == "HIGH_S"  # <= 0.5/4


def test_tier_borderline_excluded():
    # MIC 0.75: CLSI + EUCAST both intermediate (no AMBIGUOUS), inside [clsi_s/2, 2*clsi_r].
    assert eml.tier_for_isolate(["0.75"], set(), "ciprofloxacin") == "BORDERLINE"


def test_tier_ambiguous_excluded():
    # MIC 1.0: CLSI intermediate vs EUCAST R -> AMBIGUOUS (also excluded).
    assert eml.tier_for_isolate(["1"], set(), "ciprofloxacin") == "AMBIGUOUS"


def test_tier_conflict_from_calls():
    # contradictory categorical calls -> CONFLICT (excluded)
    assert eml.tier_for_isolate(["8"], {"R", "S"}, "ciprofloxacin") == "CONFLICT"


# --------------------------------------------------------------------------- #
# build_drug_labels
# --------------------------------------------------------------------------- #
def test_build_drug_labels_strict_relaxed_buckets():
    isolate_to_mics = {
        "GCA_hr": [">16"],    # HIGH_R   -> strict R + relaxed R
        "GCA_hs": ["0.03"],   # HIGH_S   -> strict S + relaxed S
        "GCA_dr": ["5"],      # DECISIVE_R -> relaxed R only
        "GCA_bl": ["0.75"],   # BORDERLINE -> excluded
    }
    res = eml.build_drug_labels(isolate_to_mics, "cipro")
    assert res["drug"] == "ciprofloxacin"
    assert res["strict"] == {"GCA_hr": "R", "GCA_hs": "S"}
    assert res["relaxed"] == {"GCA_hr": "R", "GCA_hs": "S", "GCA_dr": "R"}
    assert res["buckets"]["BORDERLINE"] == 1
    assert res["n_total"] == 4
    assert res["n_strict"] == 2
    assert res["n_relaxed"] == 3


def test_build_drug_labels_rejects_non_pilot_drug():
    with pytest.raises(ValueError):
        eml.build_drug_labels({"GCA_1": ["8"]}, "meropenem")


def test_bucket_counts_sum_to_total():
    isolate_to_mics = {f"GCA_{i}": ["8"] for i in range(5)}
    res = eml.build_drug_labels(isolate_to_mics, "ciprofloxacin")
    assert sum(res["buckets"].values()) == res["n_total"] == 5


# --------------------------------------------------------------------------- #
# write_labels round-trip
# --------------------------------------------------------------------------- #
def test_write_labels(tmp_path):
    isolate_to_mics = {"GCA_hr": [">16"], "GCA_hs": ["0.03"], "GCA_dr": ["5"]}
    res = eml.build_drug_labels(isolate_to_mics, "ciprofloxacin")
    paths = eml.write_labels(tmp_path, res)
    strict = (tmp_path / "selected_strict.tsv").read_text().splitlines()
    relaxed = (tmp_path / "selected_relaxed.tsv").read_text().splitlines()
    assert "GCA_hr\tR" in strict and "GCA_hs\tS" in strict
    assert len(strict) == 2 and len(relaxed) == 3
    bj = json.loads((tmp_path / "buckets_ciprofloxacin.json").read_text())
    assert bj["drug"] == "ciprofloxacin"
    assert bj["breakpoint_version"] == eml.BREAKPOINT_VERSION
    assert paths["strict"].endswith("selected_strict.tsv")
