"""The equivalence rule decides who wins, so it is pinned in both directions.

A rule too LENIENT flatters both callers and hides real misses; a rule too STRICT punishes notation and
understates them. Both failure modes are pinned, along with the parsing trap that manufactured fake
disagreements before it was fixed.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from dna_decode.salmserovar.equivalence import (  # noqa: E402
    canonical, equivalent, is_formula, load_formula_index,
)

TABLE = ROOT / "data" / "salmserovar_db" / "serovar_table.tsv"
RESULT = ROOT / "wiki" / "salmserovar_validation_2026-09-04.json"


@pytest.fixture(scope="module")
def idx():
    if not TABLE.exists():
        pytest.skip("White-Kauffmann table not built on this host")
    return load_formula_index(TABLE)


# --- leniency must not leak -------------------------------------------------------------------
def test_distinct_serovars_are_never_equivalent(idx):
    """The failure mode that would silently inflate every number in the artifact."""
    for a, b in [("Johannesburg", "Cubana"), ("Newport", "Newbrunswick"),
                 ("Typhimurium", "Typhi"), ("Enteritidis", "Dublin"),
                 ("Infantis", "I -:r:1,5")]:
        ok, _ = equivalent(a, b, idx)
        assert ok is False, f"{a!r} and {b!r} must NOT be equivalent"


def test_no_fuzzy_matching_on_near_identical_names(idx):
    ok, _ = equivalent("Newport", "Newpor", idx)
    assert ok is False, "a one-character difference must not be forgiven"


# --- notation must be forgiven ------------------------------------------------------------------
def test_notation_variants_are_equivalent(idx):
    for a, b in [("Typhimurium", "Typhimurium var. 5-"),
                 ("typhimurium - monophasic", "I 4,[5],12:i:-"),
                 ("4,[5],12:i:-", "i 4,5,12:i:-"),
                 ("Salmonella enterica subsp. enterica serovar Newport", "Newport")]:
        ok, why = equivalent(a, b, idx)
        assert ok is True, f"{a!r} vs {b!r} should be equivalent ({why})"


def test_formula_resolves_to_name_through_the_committed_table(idx):
    ok, why = equivalent("Enteritidis", "9:g,m:-", idx)
    assert ok is True and "white-kauffmann" in why


def test_without_the_table_a_formula_does_not_silently_match_a_name():
    """Degrading to no table must LOSE matches, never invent them."""
    ok, _ = equivalent("Enteritidis", "9:g,m:-", None)
    assert ok is False


# --- the parsing trap ---------------------------------------------------------------------------
def test_computed_types_split_respects_quotes_because_formulas_contain_commas():
    """The bug that manufactured fake disagreements: naive comma-split shreds `I 4,[5],12:i:-`."""
    from build_salmserovar_cohort import parse_computed_serotype, parse_computed_types
    raw = '"serotype=I 4,[5],12:i:-","antigen_formula=4:i:-"'
    assert parse_computed_serotype(raw) == "I 4,[5],12:i:-"
    assert parse_computed_types('"serotype=Heidelberg","antigen_formula=4:r:1,2"')[
        "antigen_formula"] == "4:r:1,2"


def test_placeholder_labels_are_rejected_as_serovars():
    """`pending` reached the first cohort and would have scored as a miss."""
    from build_salmserovar_cohort import NON_SPECIFIC, norm_serovar
    for junk in ("pending", "Pending", "unknown", "not determined", "NULL"):
        assert norm_serovar(junk) in NON_SPECIFIC


def test_is_formula_distinguishes_formula_from_name():
    assert is_formula("4,5,12:i:1,2") and is_formula("9:g,m:-")
    assert not is_formula("Typhimurium") and not is_formula("Newport")


def test_canonical_never_empties_a_real_serovar():
    for s in ("Typhimurium", "I 4,[5],12:i:-", "Salmonella enterica serovar Newport"):
        assert canonical(s), f"{s!r} canonicalised to empty"


# --- the committed result ------------------------------------------------------------------------
@pytest.mark.skipif(not RESULT.exists(), reason="validation artifact absent")
def test_the_committed_result_reports_the_incumbent_beside_us():
    """A wrapper number without its incumbent is uninterpretable -- both must ship."""
    d = json.loads(RESULT.read_text(encoding="utf-8"))
    assert d["ours_accuracy"] is not None and d["tool_accuracy"] is not None
    assert d["delta_ours_minus_incumbent"] < 0, "the measured result is an underperformance"
    # Abstention must stay visible and separate from error.
    assert d["ours_vs_wetlab_label"]["no_call"] > 0
    assert "no_call is reported separately from miss" in " ".join(d["honest_limits"])
    # The cohort must have cleared the project's own diversity bar.
    assert d["cohort"]["passes_source_diversity_bar"] is True
