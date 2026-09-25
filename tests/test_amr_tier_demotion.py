"""Pins the AMR tier demotion: 'a label source exists' must never be projected as 'a measurement exists'.

`_AMR_STATUS_MAP` maps phenotype_source_status 'ncbi_pd' -> NEAR_INDEPENDENT + SCORED. That status means a
free NCBI-PD label source EXISTS for the cell, not that a measurement was obtained. For 10 cells those
coincide; for 3 they did not, and those inherited a measured tier while the standing report card reported no
metric for them at all (found 2026-09-24 by scripts/tier_evidence_audit.py).

THIS FILE IS THE TRIPWIRE for a deliberately FAIL-SOFT mechanism. `_amr_unscored()` returns no demotions if
its evidence sources cannot be read, which preserves old behaviour rather than mass-demoting on a missing
file -- but that direction favours the OVER-claim. So these tests assert the demotions actually FIRE: if a
source is moved, renamed or emptied, the suite breaks loudly instead of quietly restoring the over-claim.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from dna_decode.data.cell_registry import AbstentionVocab, EvidenceTier, cells
from dna_decode.data.cell_registry import _amr_unscored

REPO = Path(__file__).resolve().parent.parent

DEMOTED = {
    "amr:Salmonella:ciprofloxacin": "UNDERPOWERED",
    "amr:Acinetobacter:meropenem": "ABSTAINS_BY_DESIGN",
    "amr:Pseudomonas_aeruginosa:meropenem": "ABSTAINS_BY_DESIGN",
}
# the 10 cells the card reports as SCORED must be untouched by the demotion
STILL_SCORED = (
    "amr:Escherichia_coli_Shigella:ciprofloxacin",
    "amr:Escherichia_coli_Shigella:gentamicin",
    "amr:Klebsiella:meropenem",
    "amr:Campylobacter:ciprofloxacin",
)


def _by_id():
    return {c.cell_id: c for c in cells()}


# --- the tripwire --------------------------------------------------------------------------------

def test_the_demotion_sources_resolve_and_the_demotion_fires():
    """THE tripwire. _amr_unscored() is fail-soft; if its sources vanish it returns {} and the over-claim
    silently returns. This asserts it is non-empty and names the expected keys."""
    unscored = _amr_unscored()
    assert unscored, ("no demotions resolved -- calibrated_amr_rules.json or "
                      "wiki/provdisjoint_census_results.json is missing/moved/empty")
    keys = set(unscored)
    assert ("Acinetobacter", "meropenem") in keys
    assert ("Pseudomonas_aeruginosa", "meropenem") in keys
    assert ("Salmonella", "ciprofloxacin") in keys


@pytest.mark.parametrize("cell_id,native", sorted(DEMOTED.items()))
def test_declined_cells_do_not_hold_a_measured_tier(cell_id, native):
    c = _by_id()[cell_id]
    assert c.evidence_tier is EvidenceTier.KNOWLEDGE_BASELINE
    assert c.native_abstention == native
    assert c.abstention_vocab in (AbstentionVocab.UNDERPOWERED, AbstentionVocab.ABSTAIN_BY_DESIGN)


@pytest.mark.parametrize("cell_id,_native", sorted(DEMOTED.items()))
def test_a_demoted_cell_says_why_and_is_not_mislabelled_as_sourceless(cell_id, _native):
    """A demoted cell is NOT a 'no free source' cell -- the NCBI-PD source exists and must stay named.
    Routing it through that branch would assert something false in the opposite direction."""
    c = _by_id()[cell_id]
    assert c.validation_slice.startswith("NOT SCORED --")
    assert "NCBI Pathogen Detection" in c.label_provenance
    assert "no free isolate-level phenotype source" not in c.validation_slice


@pytest.mark.parametrize("cell_id", STILL_SCORED)
def test_genuinely_scored_cells_are_untouched(cell_id):
    """The demotion must be surgical. If it caught a SCORED cell it would be an under-claim of exactly the
    kind this repo treats as equally serious."""
    c = _by_id()[cell_id]
    assert c.evidence_tier is EvidenceTier.NEAR_INDEPENDENT
    assert c.abstention_vocab is AbstentionVocab.SCORED
    assert c.validation_slice.startswith("NCBI-PD provenance-disjoint")


def test_exactly_three_ncbi_pd_cells_are_demoted():
    """Pinned as an exact count: a 4th demotion means a SCORED cell got caught, and a 2nd means one
    stopped being caught. Either is a regression worth failing on."""
    demoted = [c for c in cells()
               if c.track == "amr" and c.claim_status == "ncbi_pd"
               and c.evidence_tier is not EvidenceTier.NEAR_INDEPENDENT]
    assert {c.cell_id for c in demoted} == set(DEMOTED)


def test_the_registry_agrees_with_the_standing_card():
    """The point of reading the card's own sources: the two surfaces must not disagree about which cells
    were scored. Derived from the card, not hardcoded."""
    card = json.loads((REPO / "wiki" / "decoder_validation_report_card.json")
                      .read_text(encoding="utf-8"))
    rows = card.get("cells", card)
    by_key = {(str(c.organism).lower(), str(c.target).lower()): c for c in cells() if c.track == "amr"}
    for r in rows:
        cell = by_key.get((str(r.get("organism", "")).lower(), str(r.get("drug", "")).lower()))
        if cell is None:
            continue
        measured = cell.evidence_tier in (EvidenceTier.NEAR_INDEPENDENT,
                                          EvidenceTier.INDEPENDENT_MEASURED)
        card_scored = r.get("state") == "SCORED"
        assert measured == card_scored, (
            f"{cell.cell_id}: registry measured={measured} but card state={r.get('state')}")


def test_frozen_surface_is_only_read_never_rewritten():
    """The demotion changes what the REGISTRY derives, not the frozen surface it derives from -- which is
    one of the five files the prospective lock pins."""
    from dna_decode.data.shipped_decoder_surface import shipped_decoder_rows
    statuses = {r["phenotype_source_status"] for r in shipped_decoder_rows()}
    assert "ncbi_pd" in statuses          # the surface still says a label source exists
    rows = [r for r in shipped_decoder_rows()
            if (r["organism"], r["drug"]) == ("Salmonella", "ciprofloxacin")]
    assert rows and rows[0]["phenotype_source_status"] == "ncbi_pd"
