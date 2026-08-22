"""Tests for the demand-completion probe (`scripts/fba_demand_completion_probe.py`).

The probe's answer depends entirely on which metabolites it is willing to consider, so the currency
exclusion list is the load-bearing piece: an over-broad entry silently HIDES a real missing demand, and a
missing entry floods every gene with trivially-shared cofactors. The last test checks that against the
real artifact rather than in the abstract.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.fba_demand_completion_probe import (  # noqa: E402
    CURRENCY_BASES,
    base_name,
    is_currency,
)


def test_base_name_strips_only_the_compartment_tag():
    assert base_name("adphep_LD_c") == "adphep_LD"
    assert base_name("gmhep7p_c") == "gmhep7p"
    assert base_name("feenter_p") == "feenter"
    assert base_name("atp_c") == "atp"


def test_base_name_leaves_an_untagged_id_alone():
    assert base_name("weird") == "weird"


def test_currency_is_detected_in_every_compartment():
    assert is_currency("atp_c")
    assert is_currency("h2o_p")
    assert is_currency("h_e")


def test_real_pathway_intermediates_are_not_currency():
    """These are the actual answers the probe returned; if any were classified as currency the result
    would silently disappear."""
    for m in ("adphep_LD_c", "gmhep7p_c", "gmhep1p_c", "hlipa_c", "hhlipa_c", "phhlipa_c",
              "gicolipa_c", "feenter_c", "fe3dhbzs3_c", "dtdp4addg_c", "hemeO_c", "udpg_c",
              "camp_c", "gdptp_c"):
        assert not is_currency(m), f"{m} must not be treated as a currency metabolite"


def test_currency_list_does_not_swallow_any_reported_answer():
    """INTEGRITY CHECK against the real artifact: no metabolite the probe named as a sole-route product
    may appear in the currency list. If it did, re-running with a corrected list would surface results
    the committed artifact claims do not exist."""
    art = Path("wiki/fba_demand_completion_2026-08-22.json")
    if not art.exists():                      # artifact is regenerable; skip rather than fail CI
        return
    d = json.loads(art.read_text(encoding="utf-8"))
    named = {x["metabolite"] for g in d["genes"] for x in g["sole_route_metabolites"]}
    assert named, "artifact reports no sole-route metabolites -- the probe found nothing to guard"
    swallowed = sorted(m for m in named if base_name(m) in CURRENCY_BASES)
    assert not swallowed, f"currency list would hide reported answers: {swallowed}"


def test_every_verified_gene_actually_flipped():
    """The artifact's headline is '21/21 flip'. Pin that the per-gene records agree with it, so a future
    run that quietly stopped flipping cannot keep the same summary line."""
    art = Path("wiki/fba_demand_completion_2026-08-22.json")
    if not art.exists():
        return
    d = json.loads(art.read_text(encoding="utf-8"))
    verified = [g for g in d["genes"] if g.get("verification")]
    assert verified, "no verification records present"
    assert all(g["verification"]["flips_to_essential"] for g in verified)
    assert len(verified) == d["n_verified_flip_to_essential_under_demand"]
