"""PointFinder vs an independent curated caller on catalogued QRDR positions.

Two things these tests protect, both of which were live errors during the build:

1. ABSTENTION IS NOT AGREEMENT. On a fragmented draft assembly no reference gene aligns, the caller
   returns status 'ok' with zero mutations, and a naive set comparison scores that against an equally
   empty comparator set as a PERFECT MATCH. That counts a genome the caller did no work on as a
   success.
2. A FAIRNESS RESTRICTION THAT REMOVES NOTHING IS NOT A CONTROL. The sibling ResFinder run shipped
   claiming a POINT-row exclusion that turned out to exclude nothing. Here the position restriction is
   required to demonstrate that it actually removes comparator calls.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from pointfinder_amrfinder_concordance import (  # noqa: E402
    amrfinder_point_calls, codon_of, split_symbol,
)

ARTIFACT = ROOT / "wiki" / "pointfinder_amrfinder_concordance_2026-09-05.json"


# --- symbol parsing -------------------------------------------------------------------------------

def test_symbol_splits_on_the_LAST_underscore():
    """Gene names contain underscores; splitting on the FIRST one turns '16S_rrsC_A226G' into '16S'."""
    assert split_symbol("gyrA_S83L") == ("gyrA", "S83L")
    assert split_symbol("16S_rrsC_A226G") == ("16S_rrsC", "A226G")
    assert split_symbol("noUnderscore") is None


def test_codon_parses_single_and_multi_residue_forms():
    assert codon_of("S83L") == 83
    assert codon_of("DVADD871EPEEE") == 871      # real AMRFinder symbol shape
    assert codon_of("not-a-mutation") is None


# --- the WILDTYPE filter --------------------------------------------------------------------------

def test_wildtype_rows_are_not_counted_as_calls(tmp_path):
    """mutations.tsv is --mutation_all output: it lists every SCREENED position, including clean ones.
    Keeping them would turn 'screened and found nothing' into a called mutation."""
    p = tmp_path / "mutations.tsv"
    p.write_text(
        "Element symbol\tElement name\n"
        "gyrA_S83L\tgyrA quinolone resistance\n"
        "gyrA_D87D\tEscherichia gyrA [WILDTYPE]\n",
        encoding="utf-8")
    calls = amrfinder_point_calls(p, {"gyrA"})
    assert calls == [("gyrA", "S83L", 83)]


def test_genes_outside_the_reference_set_are_ignored(tmp_path):
    p = tmp_path / "mutations.tsv"
    p.write_text("Element symbol\tElement name\n"
                 "gyrA_S83L\tx\n23S_G2802A\ty\n", encoding="utf-8")
    assert [c[0] for c in amrfinder_point_calls(p, {"gyrA"})] == ["gyrA"]


def test_a_missing_mutations_file_yields_no_calls_not_an_error(tmp_path):
    assert amrfinder_point_calls(tmp_path / "absent.tsv", {"gyrA"}) == []


# --- the committed artifact -----------------------------------------------------------------------

@pytest.fixture(scope="module")
def art() -> dict:
    if not ARTIFACT.exists():
        pytest.skip("pointfinder concordance artifact absent")
    return json.loads(ARTIFACT.read_text(encoding="utf-8"))


def test_abstentions_are_excluded_from_agreement_not_counted_as_matches(art):
    """The load-bearing accounting check: scored + abstained must partition the genomes that ran."""
    ab = art["abstention"]
    assert ab["n_scored"] + ab["n_abstained"] == ab["n_genomes_ran"]
    assert art["n_genomes"] == ab["n_scored"]
    assert art["n_genomes_exact_set_match"] <= ab["n_scored"]


def test_the_position_restriction_actually_removes_calls(art):
    """A restriction that removes nothing is not a fairness control -- the sibling ResFinder run
    shipped exactly that over-claim and had to be corrected."""
    pr = art["position_restriction"]
    assert pr["restriction_is_live"] is True
    assert pr["n_amrfinder_calls_removed"] > 0
    assert art["n_amrfinder_calls_catalogued"] < art["n_amrfinder_calls_all"]


def test_the_comparison_is_non_vacuous(art):
    assert art["n_pointfinder_calls"] > 0
    assert art["n_amrfinder_calls_catalogued"] > 0


def test_abstention_cost_is_reported_not_hidden(art):
    """An abstention is only harmless when the comparator also found nothing; the run must report how
    many real comparator calls fell inside abstained genomes."""
    ab = art["abstention"]
    assert "n_calls_missed_via_abstention" in ab
    assert ab["n_calls_missed_via_abstention"] >= 0
    assert ab["n_abstained_where_amrfinder_found_catalogued_mutations"] <= ab["n_abstained"]


def test_the_tier_does_not_move_and_the_limits_say_why(art):
    joined = " ".join(art["honest_limits"]).lower()
    assert "not a wet-lab label" in joined or "not correctness" in joined
    assert "faithful_to_tool" in joined
    assert "catalogue difference" in joined      # a discordance need not be a caller defect


def test_per_mutation_breakdown_is_present_rather_than_a_bare_rate(art):
    """Averaging discordance away is how the earlier defects stayed hidden; the split must survive."""
    assert art["per_mutation"]
    for _k, v in art["per_mutation"].items():
        assert set(v) == {"both", "pf_only", "amr_only"}
