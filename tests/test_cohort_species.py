"""Tests for the pure species-composition logic (Step 2). Offline, synthetic, no pandas, no Docker."""
from __future__ import annotations

import pytest

from dna_decode.eval.cohort_species import (
    COMPOSITION_CLEAN,
    COMPOSITION_MIXED,
    INSUFFICIENT_RESOLUTION,
    MATCHES,
    OTHER_GENUS,
    SAME_GENUS_OTHER_SPECIES,
    UNRESOLVED,
    CrosstabMismatch,
    classify_species,
    compose,
    composition_verdict,
    vocabulary_status,
)

KP = "Klebsiella pneumoniae BIDMC 53"
KA = "Klebsiella aerogenes MGH 77"
EC = "Escherichia coli K-12"


# --------------------------------------------------------------------------- classification

def test_a_pure_cohort_has_no_off_species():
    c = classify_species({"a": KP, "b": "Klebsiella pneumoniae"}, "Klebsiella_pneumoniae")
    assert c.count(MATCHES) == 2
    assert c.count(SAME_GENUS_OTHER_SPECIES) == 0
    assert c.count(OTHER_GENUS) == 0


def test_aerogenes_is_same_genus_other_species_NOT_a_match():
    """THE discriminating case. `Klebsiella aerogenes` shares the genus with the scored organism, which
    is exactly why a genus-named cohort hides it -- and it must not be counted as expected."""
    c = classify_species({"a": KP, "b": KA}, "Klebsiella_pneumoniae")
    assert c.by_bucket[MATCHES] == ["a"]
    assert c.by_bucket[SAME_GENUS_OTHER_SPECIES] == ["b"]
    assert c.count(OTHER_GENUS) == 0


def test_a_different_genus_is_its_own_bucket():
    c = classify_species({"a": KP, "b": EC}, "Klebsiella_pneumoniae")
    assert c.by_bucket[OTHER_GENUS] == ["b"]
    assert c.count(SAME_GENUS_OTHER_SPECIES) == 0


def test_a_GENUS_level_amrfinder_flag_compares_at_genus_only():
    """`-O Escherichia` is genus-level, so it cannot pin a species. Comparing against an invented
    species would manufacture mismatches for every E. coli strain in the cohort."""
    c = classify_species({"a": EC, "b": "Escherichia fergusonii", "c": KP}, "Escherichia")
    assert c.expected_species is None
    assert sorted(c.by_bucket[MATCHES]) == ["a", "b"]   # both Escherichia
    assert c.by_bucket[OTHER_GENUS] == ["c"]


def test_the_four_buckets_PARTITION_the_input():
    """A denominator that does not sum is the failure mode this module exists to avoid."""
    resolved = {"a": KP, "b": KA, "c": EC, "d": None, "e": "", "f": "Klebsiella"}
    c = classify_species(resolved, "Klebsiella_pneumoniae")
    total = sum(c.count(b) for b in (MATCHES, SAME_GENUS_OTHER_SPECIES, OTHER_GENUS, UNRESOLVED))
    assert total == len(resolved) == c.n_total


def test_a_genus_only_name_is_not_credited_as_the_expected_species():
    """`Klebsiella` alone does not establish `Klebsiella pneumoniae`."""
    c = classify_species({"a": "Klebsiella"}, "Klebsiella_pneumoniae")
    assert c.count(MATCHES) == 0
    assert c.by_bucket[SAME_GENUS_OTHER_SPECIES] == ["a"]


# --------------------------------------------------------------------------- the refusal

def test_too_many_unresolved_REFUSES_and_carries_no_numbers():
    """Asserted by KEY ABSENCE, so a future change cannot leak a composition figure past the refusal."""
    resolved = {f"a{i}": None for i in range(5)} | {"b": KP}
    v = composition_verdict(classify_species(resolved, "Klebsiella_pneumoniae"))
    assert v["verdict"] == INSUFFICIENT_RESOLUTION
    for leaked in ("n_matches_expected", "n_same_genus_other_species", "n_other_genus",
                   "species_counts", "off_species_accessions"):
        assert leaked not in v, f"refused verdict leaked {leaked}"
    assert "max_unresolved_fraction" in v and "reason" in v


def test_an_empty_cohort_refuses_rather_than_reporting_clean():
    v = composition_verdict(classify_species({}, "Klebsiella_pneumoniae"))
    assert v["verdict"] == INSUFFICIENT_RESOLUTION
    assert "n_matches_expected" not in v


def test_a_little_unresolved_still_reports():
    resolved = {f"k{i}": KP for i in range(19)} | {"bad": None}   # 1/20 = 5% <= 10% bar
    v = composition_verdict(classify_species(resolved, "Klebsiella_pneumoniae"))
    assert v["verdict"] == COMPOSITION_CLEAN
    assert v["n_unresolved"] == 1


def test_mixed_vs_clean_verdicts():
    clean = composition_verdict(classify_species({"a": KP}, "Klebsiella_pneumoniae"))
    mixed = composition_verdict(classify_species({"a": KP, "b": KA}, "Klebsiella_pneumoniae"))
    assert clean["verdict"] == COMPOSITION_CLEAN
    assert mixed["verdict"] == COMPOSITION_MIXED
    assert mixed["off_species_accessions"] == ["b"]


# --------------------------------------------------------------------------- the cross-tab

def test_crosstab_splits_outcomes_by_species():
    resolved = {"r1": KP, "r2": KA, "s1": KP}
    preds = {"r1": "R", "r2": "S", "s1": "S"}
    labels = {"r1": 1, "r2": 1, "s1": 0}
    got = compose(resolved, preds, labels)["by_species"]
    assert got["Klebsiella pneumoniae"]["tp"] == 1
    assert got["Klebsiella pneumoniae"]["tn"] == 1
    assert got["Klebsiella aerogenes"]["fn"] == 1


def test_crosstab_RECONCILES_against_the_supplied_confusion():
    resolved = {"a": KP, "b": KA}
    preds, labels = {"a": "R", "b": "S"}, {"a": 1, "b": 1}
    ok = compose(resolved, preds, labels, confusion={"tp": 1, "fp": 0, "tn": 0, "fn": 1})
    assert ok["totals"]["tp"] == 1 and ok["totals"]["fn"] == 1


def test_a_crosstab_that_does_not_reconcile_RAISES():
    """A decomposition that disagrees with the matrix it decomposes is not evidence."""
    resolved = {"a": KP}
    with pytest.raises(CrosstabMismatch, match="refusing"):
        compose(resolved, {"a": "R"}, {"a": 1}, confusion={"tp": 99, "fp": 0, "tn": 0, "fn": 0})


def test_abstain_is_its_own_cell_not_folded_into_S():
    got = compose({"a": KP}, {"a": "INDETERMINATE"}, {"a": 1})
    assert got["by_species"]["Klebsiella pneumoniae"]["abstain"] == 1
    assert got["by_species"]["Klebsiella pneumoniae"]["fn"] == 0


def test_an_unresolved_species_gets_its_own_row_not_dropped():
    got = compose({"a": None}, {"a": "S"}, {"a": 1})["by_species"]
    assert got["(unresolved)"]["fn"] == 1


def test_a_prediction_without_a_label_is_skipped_not_guessed():
    got = compose({"a": KP, "b": KP}, {"a": "R", "b": "R"}, {"a": 1})
    assert got["totals"]["tp"] == 1


# --------------------------------------------------------------------------- vocabulary

def test_aerogenes_is_reported_as_outside_the_vocabulary_rather_than_crashing():
    """K. aerogenes is not one of the 14 supported organisms, and it is the finding -- it must survive
    to the artifact instead of raising UnknownOrganism mid-audit."""
    assert vocabulary_status("Klebsiella aerogenes") == "not_in_organism_vocab"
    assert vocabulary_status("Klebsiella pneumoniae") == "in_organism_vocab"
    assert vocabulary_status(None) == "unresolved"
