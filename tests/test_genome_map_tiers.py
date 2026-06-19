"""Step 3 tests — tier classifier (table-driven over real + synthetic wording)."""
from __future__ import annotations

import pytest

from dna_decode.genome_map import (
    TIER_CURATED_FUNCTION,
    TIER_HOMOLOGY_HYPOTHESIS,
    TIER_UNKNOWN,
)
from dna_decode.genome_map.tiers import classify_feature_tier
from dna_decode.genome_map.tier_vocab import matches_any, LOW_CONFIDENCE_PATTERNS


# (product, gene_symbol, expected_tier)
_CASES = [
    # named gene symbol -> curated regardless of product
    ("DNA gyrase subunit A", "gyrA", TIER_CURATED_FUNCTION),
    ("", "parC", TIER_CURATED_FUNCTION),
    # specific product, no gene symbol -> curated
    ("DNA gyrase subunit A", "", TIER_CURATED_FUNCTION),
    ("beta-lactamase TEM", "", TIER_CURATED_FUNCTION),
    # low-confidence wording, no gene symbol -> homology
    ("putative oxidoreductase", "", TIER_HOMOLOGY_HYPOTHESIS),
    ("ABC transporter, by similarity", "", TIER_HOMOLOGY_HYPOTHESIS),
    ("DUF1234 domain-containing protein", "", TIER_HOMOLOGY_HYPOTHESIS),
    ("uncharacterized membrane protein", "", TIER_HOMOLOGY_HYPOTHESIS),  # "uncharacterized" = low-conf
    ("uncharacterized protein", "", TIER_UNKNOWN),  # full "uncharacterized protein" = unknown
    # unknown wording / empty -> unknown
    ("hypothetical protein", "", TIER_UNKNOWN),
    ("", "", TIER_UNKNOWN),
    ("protein of unknown function", "", TIER_UNKNOWN),
]


@pytest.mark.parametrize("product,gene,expected", _CASES)
def test_classify_feature_tier_table(product, gene, expected):
    tier, reason = classify_feature_tier(product, gene)
    assert tier == expected, f"{product!r}/{gene!r} -> {tier} (reason: {reason})"
    assert reason  # reason always populated


def test_reason_names_matched_pattern():
    _, reason = classify_feature_tier("putative kinase", "")
    assert "putative" in reason


def test_reason_names_gene_symbol():
    _, reason = classify_feature_tier("anything", "gyrA")
    assert "gyrA" in reason


def test_low_confidence_wording_demotes_even_specific_sounding():
    # G1 type (a): "putative DNA gyrase" with NO firm gene call is demoted.
    tier, reason = classify_feature_tier("putative DNA gyrase subunit", "")
    assert tier == TIER_HOMOLOGY_HYPOTHESIS


def test_named_gene_outranks_weak_product_wording():
    # A firm gene assignment is curated; the weak wording is the caller's secondary.
    tier, _ = classify_feature_tier("putative DNA gyrase subunit", "gyrA")
    assert tier == TIER_CURATED_FUNCTION


def test_matches_any_helper():
    assert matches_any("a putative thing", LOW_CONFIDENCE_PATTERNS) == "putative"
    assert matches_any("DNA gyrase", LOW_CONFIDENCE_PATTERNS) is None
