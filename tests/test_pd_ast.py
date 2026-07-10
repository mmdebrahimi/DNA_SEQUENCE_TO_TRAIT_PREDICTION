"""Tests for the canonical NCBI-PD `AST_phenotypes` parser (shared by the census + the accrual fetch)."""
from __future__ import annotations

import pytest

from dna_decode.data.pd_ast import ast_label_for, parse_ast_phenotypes

# Verbatim shape of the live PD field: comma-separated, wrapped in literal double quotes.
REAL = '"ampicillin=ND,cefazolin=ND,ceftriaxone=R,ciprofloxacin=S,gentamicin=ND"'


def test_parses_the_real_quoted_comma_separated_field():
    assert parse_ast_phenotypes(REAL, {"ceftriaxone", "ciprofloxacin", "gentamicin"}) == {
        "ceftriaxone": "R", "ciprofloxacin": "S"}          # ND dropped


def test_matches_a_drug_at_the_first_position():
    """Regression: the opening quote rides on the first token."""
    assert parse_ast_phenotypes('"ciprofloxacin=R,gentamicin=S"', {"ciprofloxacin"}) == {"ciprofloxacin": "R"}


def test_matches_a_drug_at_the_last_position():
    """Regression: the closing quote rides on the last token."""
    assert parse_ast_phenotypes('"gentamicin=S,ciprofloxacin=R"', {"ciprofloxacin"}) == {"ciprofloxacin": "R"}


def test_single_drug_field_is_fully_quoted():
    assert parse_ast_phenotypes('"ciprofloxacin=R"', {"ciprofloxacin"}) == {"ciprofloxacin": "R"}


def test_the_naive_membership_test_would_miss_the_end_positions():
    """Pins WHY this module exists — the old bug was silent, never an exception."""
    field = '"ciprofloxacin=R,gentamicin=S"'
    assert "ciprofloxacin=R" not in field.split(",")        # the old census test — misses it
    assert parse_ast_phenotypes(field, {"ciprofloxacin"})   # the correct parser finds it


def test_drops_non_binary_calls():
    got = parse_ast_phenotypes('"a=ND,b=I,c=NS,d=R"', {"a", "b", "c", "d"})
    assert got == {"d": "R"}


def test_drops_null_and_empty():
    assert parse_ast_phenotypes("NULL", {"a"}) == {}
    assert parse_ast_phenotypes("", {"a"}) == {}
    assert parse_ast_phenotypes(None, {"a"}) == {}


def test_ignores_offpanel_drugs():
    assert parse_ast_phenotypes('"amoxicillin=R"', {"ciprofloxacin"}) == {}


def test_tolerates_the_semicolon_form_of_the_derived_candidates_tsv():
    assert parse_ast_phenotypes("ceftriaxone=R;ciprofloxacin=R", {"ceftriaxone", "ciprofloxacin"}) == {
        "ceftriaxone": "R", "ciprofloxacin": "R"}


def test_is_case_insensitive_on_drug_and_label():
    assert parse_ast_phenotypes('"Ciprofloxacin=r"', {"ciprofloxacin"}) == {"ciprofloxacin": "R"}


def test_malformed_tokens_are_skipped():
    assert parse_ast_phenotypes('"garbage,=,cipro,a=R"', {"a"}) == {"a": "R"}


@pytest.mark.parametrize("drug,expected", [("ceftriaxone", "R"), ("ciprofloxacin", "S"), ("gentamicin", None)])
def test_ast_label_for(drug, expected):
    assert ast_label_for(REAL, drug) == expected
