"""Unit tests for the HCMV antiviral-resistance cell (dna_decode/data/hcmv_amr.py).

Pins: catalog self-integrity (WT consistency), drug->gene routing, mutant-level R/S semantics, per-mutation
CROSS-RESISTANCE drug specificity, and the phenotyped-BENIGN specificity anchor (a benign polymorphism must
return S, never R — the HCMV analogue of SARS-CoV-2 Mpro P132H). Pure; no network / no BLAST.
"""
from __future__ import annotations

from dna_decode.data.hcmv_amr import (
    BENIGN_BY_GENE,
    RESISTANCE_BY_GENE,
    all_supported_hcmv_drugs,
    call_hcmv_observed,
    genes_for_hcmv_drug,
    is_hcmv_resistance_mutation,
    wt_consistency_ok,
)


def test_supported_drugs():
    assert set(all_supported_hcmv_drugs()) == {
        "ganciclovir", "valganciclovir", "cidofovir", "foscarnet", "letermovir"}


def test_drug_gene_routing():
    assert genes_for_hcmv_drug("ganciclovir") == ("UL97", "UL54")
    assert genes_for_hcmv_drug("valganciclovir") == ("UL97", "UL54")
    assert genes_for_hcmv_drug("cidofovir") == ("UL54",)
    assert genes_for_hcmv_drug("foscarnet") == ("UL54",)
    assert genes_for_hcmv_drug("letermovir") == ("UL56",)
    assert genes_for_hcmv_drug("aspirin") == ()


def test_catalog_wt_consistency_integrity_gate():
    # every catalogued POINT mutation at the same (gene, pos) must assert the same WT residue
    assert wt_consistency_ok() is True


def test_catalog_nonempty_per_gene():
    for gene in ("UL97", "UL54", "UL56"):
        assert len(RESISTANCE_BY_GENE[gene]) >= 5
        assert len(BENIGN_BY_GENE[gene]) >= 2


def test_r_on_canonical_gcv_ul97():
    call = call_hcmv_observed("ganciclovir", {"UL97": {"M460V"}})
    assert call.prediction == "R" and call.determinants == ["UL97:M460V"]


def test_s_on_benign_ul54_polymorphism():
    # C304S is a recombinant-phenotyped no-effect polymorphism -> must be S (specificity anchor)
    call = call_hcmv_observed("ganciclovir", {"UL54": {"C304S"}})
    assert call.prediction == "S"
    assert call.undetectable_mechanisms  # S surfaces the honest blind-spot list


def test_every_benign_returns_s_across_its_drugs():
    drug_for_gene = {"UL97": "ganciclovir", "UL54": "cidofovir", "UL56": "letermovir"}
    for gene, benign in BENIGN_BY_GENE.items():
        drug = drug_for_gene[gene]
        for b in benign:
            call = call_hcmv_observed(drug, {gene: {b}})
            assert call.prediction == "S", f"benign {gene}:{b} wrongly called {call.prediction} for {drug}"


def test_letermovir_high_grade_ul56():
    assert call_hcmv_observed("letermovir", {"UL56": {"C325W"}}).prediction == "R"
    assert call_hcmv_observed("letermovir", {"UL56": {"R369S"}}).prediction == "R"


def test_cross_resistance_is_per_mutation_drug_specific():
    # Q578H is catalogued GCV+CDV+FOS -> R for all three
    for drug in ("ganciclovir", "cidofovir", "foscarnet"):
        assert call_hcmv_observed(drug, {"UL54": {"Q578H"}}).prediction == "R"
    # S585A is FOS-only -> R for foscarnet, S for GCV/CDV
    assert call_hcmv_observed("foscarnet", {"UL54": {"S585A"}}).prediction == "R"
    assert call_hcmv_observed("ganciclovir", {"UL54": {"S585A"}}).prediction == "S"
    assert call_hcmv_observed("cidofovir", {"UL54": {"S585A"}}).prediction == "S"


def test_drug_only_scored_against_its_target_genes():
    # M460V is a UL97 GCV marker; cidofovir routes to UL54 only -> S (not scored against UL97)
    assert call_hcmv_observed("cidofovir", {"UL97": {"M460V"}}).prediction == "S"
    # a GCV UL54 marker in a UL97-only observed set is still found for GCV (GCV covers both genes)
    assert call_hcmv_observed("ganciclovir", {"UL54": {"F412L"}}).prediction == "R"


def test_multi_gene_observed_gcv():
    call = call_hcmv_observed("ganciclovir", {"UL97": {"M460V"}, "UL54": {"F412L"}})
    assert call.prediction == "R"
    assert set(call.determinants) == {"UL97:M460V", "UL54:F412L"}


def test_indeterminate_on_unknown_drug():
    assert call_hcmv_observed("rifampicin", {"UL97": {"M460V"}}).prediction == "INDETERMINATE"


def test_is_hcmv_resistance_mutation_drug_filter():
    assert is_hcmv_resistance_mutation("UL54", "S585A", "foscarnet") is True
    assert is_hcmv_resistance_mutation("UL54", "S585A", "ganciclovir") is False
    assert is_hcmv_resistance_mutation("UL54", "C304S") is False   # benign, not in resistance catalog


def test_valganciclovir_equals_ganciclovir_mechanism():
    assert call_hcmv_observed("valganciclovir", {"UL97": {"C603W"}}).prediction == "R"
