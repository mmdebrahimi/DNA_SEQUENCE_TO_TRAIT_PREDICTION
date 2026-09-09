"""The ktype caller measured against full-locus Kaptive — and the accounting that makes it honest.

`ktype` was the last cell on the shared "never measured" default, believed cohort-blocked. It was not:
the caller types the capsule from ONE gene (`wzi`) while Kaptive types the FULL K locus, so Kaptive is
an independent comparator and needs no wet-lab cohort.

The bar is not 100%. wzi -> K-type is ~94% predictive and NOT one-to-one (Brisse 2013 JCM), so these
tests pin the things that keep the measurement from over- or under-reading that: an uncertain reference
cannot adjudicate, an abstention is not an error, and an ambiguous multi-KL call is a miss in the strict
figure but is reported separately rather than silently punished.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from ktype_kaptive_concordance import TYPEABLE, klebsiella_genomes  # noqa: E402

ARTIFACT = ROOT / "wiki" / "ktype_kaptive_concordance_2026-09-08.json"


@pytest.fixture(scope="module")
def art() -> dict:
    if not ARTIFACT.exists():
        pytest.skip("ktype concordance artifact absent")
    return json.loads(ARTIFACT.read_text(encoding="utf-8"))


# --- cohort discovery -----------------------------------------------------------------------------

def test_cohort_discovery_reads_headerless_selected_tsv():
    """selected.tsv is bare `<accession>\\t<label>` with NO header. A DictReader silently eats the
    first genome -- which returned a vacuous 0 accessions on the first attempt at this run."""
    g = klebsiella_genomes()
    if not g:
        pytest.skip("Klebsiella cohorts/assemblies not present (external drive offline)")
    assert len(g) > 100
    assert all(a.startswith(("GCA_", "GCF_")) for a in g)


# --- the accounting that makes the number honest --------------------------------------------------

def test_an_uncertain_reference_is_excluded_not_counted(art):
    """Kaptive rows it does not call typeable cannot adjudicate; they leave the denominator."""
    assert art["n_confident_reference"] + art["n_kaptive_untypeable_excluded"] \
        + art["n_kaptive_error"] == art["n_genomes"]
    assert TYPEABLE == "typeable"


def test_abstentions_are_separated_from_disagreements(art):
    """An abstention is not a wrong answer; pooling the two would understate the caller."""
    assert art["n_comparable"] + art["n_our_abstentions"] == art["n_confident_reference"]
    assert art["n_agree"] + art["n_disagree"] == art["n_comparable"]
    assert art["n_our_abstentions"] > 0          # non-vacuous: the caller does abstain here


def test_the_strict_figure_is_the_headline_and_lenient_is_reported_beside_it(art):
    """An ambiguous multi-KL string is not a resolved call. Crediting it silently would inflate."""
    assert art["agreement"] <= art["agreement_lenient_crediting_ambiguous"]
    assert art["n_ambiguous_containing_kaptive_answer"] <= art["n_ambiguous_calls_scored_as_miss"]
    assert "strict" in art["which_figure_is_the_headline"]


def test_agreement_is_computed_on_the_comparable_denominator(art):
    assert art["agreement"] == pytest.approx(art["n_agree"] / art["n_comparable"], abs=1e-9)


# --- the verdict ----------------------------------------------------------------------------------

def test_verdict_matches_the_measured_agreement(art):
    """The verdict must follow the number, not be asserted alongside it."""
    a, ceil = art["agreement"], art["published_ceiling"]
    expected = ("REACHES_THE_WZI_METHOD_CEILING" if a >= ceil else
                "NEAR_THE_WZI_METHOD_CEILING" if a >= ceil - 0.10 else
                "FALLS_SHORT_OF_THE_WZI_METHOD_CEILING")
    assert art["verdict"] == expected


def test_the_ceiling_is_sourced_not_asserted(art):
    """~94% is a published property with a citation, not a number chosen to make this pass."""
    assert art["published_ceiling"] == 0.94
    assert "Brisse" in art["ceiling_source"]
    assert "not one-to-one" in art["ceiling_source"].lower()


def test_the_comparison_is_non_vacuous(art):
    assert art["n_comparable"] > 100
    assert art["n_agree"] > 0 and art["n_disagree"] > 0


# --- tier discipline ------------------------------------------------------------------------------

def test_the_tier_does_not_move_and_the_limits_say_why(art):
    joined = " ".join(art["honest_limits"]).lower()
    assert "tool, not a wet-lab label" in joined
    assert "faithful_to_tool" in joined
    assert "enriched for resistance" in joined


def test_the_ceiling_is_not_treated_as_an_exact_bar(art):
    """It is a published property of OTHER cohorts; a single hard threshold would over-read it."""
    joined = " ".join(art["honest_limits"]).lower()
    assert "not a constant" in joined


# --- the KL15 lead, resolved ----------------------------------------------------------------------

def test_the_kl15_lead_was_resolved_as_a_method_limit(art):
    """It was recorded as "lead, not diagnosis". One allele, wzi_50, maps to three K loci -- correct on
    10, KL52 on 6, KL51 on 4 -- which is the published non-one-to-one behaviour, not a defect."""
    r = art.get("kl15_lead_resolved")
    assert r, "the KL15 lead must be resolved or still marked open"
    assert r["verdict"] == "CONFIRMED_METHOD_LIMIT_NOT_AN_IMPLEMENTATION_DEFECT"
    c = r["counts"]
    assert c["wzi_50_correct_KL15"] == 10
    assert c["wzi_50_to_KL52"] + c["wzi_50_to_KL51"] == 10   # right half the time, wrong half
    assert r["action"].startswith("none")


def test_the_resolution_names_why_no_code_change_helps(art):
    why = art["kl15_lead_resolved"]["why"].lower()
    assert "single-gene" in why and "cannot separate" in why
