"""Tests for the HCMV report card builder + the trust-surface HCMV tier (UNKNOWN -> IN_DISTRIBUTION)."""
from __future__ import annotations

from dna_decode.data.trust_surface import lookup_trust
from scripts.build_hcmv_report_card import build


def test_card_is_in_distribution_not_independent():
    rc = build()
    assert rc["tier"] == "IN_DISTRIBUTION"                 # honest: catalog curated from the same fold-change
    assert "closed" in rc["independence"].lower()          # independence is a CLOSED negative for free data
    assert rc["n_cells"] == 5                              # GCV/valGCV/CDV/FOS/letermovir


def test_card_catalog_counts_per_gene():
    genes = build()["genes"]
    assert set(genes) == {"UL97", "UL54", "UL56"}
    for g in genes.values():
        assert g["n_resistance"] >= 5 and g["n_benign"] >= 2


def test_trust_tier_hcmv_is_in_distribution():
    # the CLI trust surface must render a REAL tier for HCMV drugs, not UNKNOWN
    for drug in ("ganciclovir", "valganciclovir", "cidofovir", "foscarnet", "letermovir"):
        r = lookup_trust(drug, "HCMV")
        assert r["tier"] == "IN_DISTRIBUTION", f"{drug} tier {r['tier']} != IN_DISTRIBUTION"


def test_trust_tier_guarded_by_genus():
    # a wrong genus for an HCMV drug must NOT silently borrow the HCMV tier (namespace guard)
    r = lookup_trust("ganciclovir", "Escherichia")
    assert r["tier"] != "IN_DISTRIBUTION"
