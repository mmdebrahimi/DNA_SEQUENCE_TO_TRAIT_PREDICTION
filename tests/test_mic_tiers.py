"""Tests for `dna_decode/data/mic_tiers.py` — the shared MIC-tier classifier.

Parametrized over the 4 currently-supported drugs (cipro, cef, tet, gent).
Pins:
  - Breakpoint catalog values (per CLSI 2024 + EUCAST 14.0 E. coli)
  - HIGH_R / HIGH_S thresholds (4x safety margin from breakpoints)
  - DECISIVE_R / DECISIVE_S split (between 2x breakpoint and 4x breakpoint)
  - BORDERLINE gray zone (between CLSI_S/2 and 2*CLSI_R)
  - AMBIGUOUS skip when EUCAST=None (tet)
  - CONFLICT detection from distinct R + S labels
  - NO_MIC when no numeric values
  - AMRFinder Class filter per drug (cipro keeps MULTIDRUG; cef keeps β-lactam family)
"""
from __future__ import annotations

import pytest

from dna_decode.data.mic_tiers import (
    DRUG_BREAKPOINTS,
    DRUG_AMRFINDER_CLASSES,
    RELAXED_MIC_TIERS,
    STRICT_MIC_TIERS,
    UnknownDrugError,
    amrfinder_classes_for,
    breakpoints_for,
    classify_tier,
    is_relaxed_r,
    is_relaxed_s,
    is_strict_r,
    is_strict_s,
    supported_drugs,
)


# ---- breakpoint catalog --------------------------------------------------


def test_supported_drugs_returns_6_drugs():
    # meropenem added 2026-06-07 (carbapenem); oxacillin added 2026-06-07 (1st Gram-positive, S. aureus mecA).
    assert set(supported_drugs()) == {"ciprofloxacin", "ceftriaxone", "tetracycline",
                                      "gentamicin", "meropenem", "oxacillin"}


def test_breakpoints_for_unknown_drug_raises():
    with pytest.raises(UnknownDrugError):
        breakpoints_for("unobtainium")


def test_cipro_breakpoints_match_clsi_2024_eucast_14():
    bp = breakpoints_for("ciprofloxacin")
    assert bp["clsi_r"] == 2.0
    assert bp["clsi_s"] == 0.5
    assert bp["eucast_r"] == 1.0
    assert bp["eucast_s"] == 0.25


def test_cef_breakpoints_match_clsi_2024_eucast_14():
    bp = breakpoints_for("ceftriaxone")
    assert bp["clsi_r"] == 4.0
    assert bp["clsi_s"] == 1.0
    assert bp["eucast_r"] == 2.0
    assert bp["eucast_s"] == 1.0


def test_tet_breakpoints_have_no_eucast():
    bp = breakpoints_for("tetracycline")
    assert bp["clsi_r"] == 16.0
    assert bp["clsi_s"] == 4.0
    assert bp["eucast_r"] is None
    assert bp["eucast_s"] is None


def test_gent_breakpoints_match_clsi_2024_eucast_14():
    bp = breakpoints_for("gentamicin")
    assert bp["clsi_r"] == 16.0
    assert bp["clsi_s"] == 4.0
    assert bp["eucast_r"] == 4.0
    assert bp["eucast_s"] == 2.0


# ---- classify_tier ------------------------------------------------------


@pytest.fixture(params=["ciprofloxacin", "ceftriaxone", "tetracycline", "gentamicin"])
def drug(request):
    return request.param


def test_conflict_when_both_r_and_s_labels(drug):
    bp = breakpoints_for(drug)
    assert classify_tier([1.0], {"R", "S"}, bp) == "CONFLICT"


def test_no_mic_when_empty(drug):
    bp = breakpoints_for(drug)
    assert classify_tier([], set(), bp) == "NO_MIC"


def test_no_mic_when_all_nan(drug):
    from math import nan
    bp = breakpoints_for(drug)
    assert classify_tier([nan, nan], set(), bp) == "NO_MIC"


def test_high_r_at_4x_clsi_r(drug):
    bp = breakpoints_for(drug)
    mic_high = bp["clsi_r"] * 4
    assert classify_tier([mic_high], {"R"}, bp) == "HIGH_R"


def test_high_s_at_clsi_s_div_4(drug):
    bp = breakpoints_for(drug)
    mic_low = bp["clsi_s"] / 4
    assert classify_tier([mic_low], {"S"}, bp) == "HIGH_S"


def test_borderline_at_clsi_s_div_2(drug):
    bp = breakpoints_for(drug)
    mic = bp["clsi_s"] / 2
    assert classify_tier([mic], set(), bp) == "BORDERLINE"


def test_borderline_at_2x_clsi_r(drug):
    bp = breakpoints_for(drug)
    mic = bp["clsi_r"] * 2
    assert classify_tier([mic], set(), bp) == "BORDERLINE"


def test_decisive_r_between_2x_and_4x_clsi_r(drug):
    bp = breakpoints_for(drug)
    # Pick MIC strictly between 2*CLSI_R and 4*CLSI_R (e.g., 3*CLSI_R)
    mic = bp["clsi_r"] * 3
    assert classify_tier([mic], set(), bp) == "DECISIVE_R"


def test_decisive_s_between_clsi_s_div_4_and_clsi_s_div_2(drug):
    bp = breakpoints_for(drug)
    # Strictly between CLSI_S/4 and CLSI_S/2: use 3*CLSI_S/8
    mic = bp["clsi_s"] * 3 / 8
    assert classify_tier([mic], set(), bp) == "DECISIVE_S"


def test_tet_skips_ambiguous_check_when_eucast_none():
    """Tet has no EUCAST breakpoints; the AMBIGUOUS branch must short-circuit."""
    bp = breakpoints_for("tetracycline")
    # MIC = 100 is far above CLSI_R (16); without EUCAST check this is HIGH_R, not AMBIGUOUS
    assert classify_tier([100.0], {"R"}, bp) == "HIGH_R"


def test_cipro_ambiguous_when_clsi_eucast_disagree():
    """Cipro CLSI_S=0.5, EUCAST_S=0.25. MIC=0.4 is S by CLSI, I by EUCAST → AMBIGUOUS."""
    bp = breakpoints_for("ciprofloxacin")
    # MIC = 0.4: CLSI call = S (<=0.5), EUCAST call = I (between 0.25 and 1.0) → disagree
    assert classify_tier([0.4], set(), bp) == "AMBIGUOUS"


def test_median_used_not_mean(drug):
    bp = breakpoints_for(drug)
    # 4 values: 1, 1, 100, 100 — median = 50.5 (between 1 and 100). For cipro that's HIGH_R.
    # For tet, median between 1 and 100 is 50.5 which is > 4*16=64? No, 50.5 < 64; would be DECISIVE_R.
    # Use a cleaner test: 3 values where median != mean.
    bp_low = bp["clsi_s"] / 4
    bp_high = bp["clsi_r"] * 4
    # Median of [bp_low, bp_high, bp_high] = bp_high → HIGH_R
    assert classify_tier([bp_low, bp_high, bp_high], {"R"}, bp) == "HIGH_R"


# ---- AMRFinder Class filter ---------------------------------------------


def test_cipro_amrfinder_classes_include_multidrug():
    cls = amrfinder_classes_for("ciprofloxacin")
    assert "QUINOLONE" in cls
    assert "FLUOROQUINOLONE" in cls
    assert "MULTIDRUG" in cls


def test_cef_amrfinder_classes_include_beta_lactam_family():
    cls = amrfinder_classes_for("ceftriaxone")
    assert "BETA-LACTAM" in cls
    assert "CARBAPENEM" in cls
    assert "CEPHALOSPORIN" in cls
    assert "MULTIDRUG" in cls


def test_tet_amrfinder_classes():
    cls = amrfinder_classes_for("tetracycline")
    assert "TETRACYCLINE" in cls
    assert "MULTIDRUG" in cls


def test_gent_amrfinder_classes():
    cls = amrfinder_classes_for("gentamicin")
    assert "AMINOGLYCOSIDE" in cls
    assert "MULTIDRUG" in cls


def test_amrfinder_classes_unknown_drug_raises():
    with pytest.raises(UnknownDrugError):
        amrfinder_classes_for("unobtainium")


# ---- Strict vs relaxed pass helpers --------------------------------------


def test_strict_constants_only_high():
    assert STRICT_MIC_TIERS == frozenset({"HIGH_R", "HIGH_S"})


def test_relaxed_constants_include_decisive():
    assert RELAXED_MIC_TIERS == frozenset({"HIGH_R", "HIGH_S", "DECISIVE_R", "DECISIVE_S"})


def test_is_strict_r_only_high_r():
    assert is_strict_r("HIGH_R")
    assert not is_strict_r("DECISIVE_R")
    assert not is_strict_r("HIGH_S")


def test_is_relaxed_r_accepts_high_and_decisive():
    assert is_relaxed_r("HIGH_R")
    assert is_relaxed_r("DECISIVE_R")
    assert not is_relaxed_r("HIGH_S")
    assert not is_relaxed_r("BORDERLINE")


def test_is_relaxed_s_accepts_high_and_decisive():
    assert is_relaxed_s("HIGH_S")
    assert is_relaxed_s("DECISIVE_S")
    assert not is_relaxed_s("HIGH_R")
    assert not is_relaxed_s("BORDERLINE")


# ---- Per-drug mechanism catalog -----------------------------------------


def test_cipro_loci_catalog_matches_existing_audit_script():
    """Catalog must match the canonical Phase 1 cipro catalog in
    scripts/cipro_mechanism_audit.py — same loci, same mechanism names."""
    from dna_decode.data.mic_tiers import loci_by_mechanism_for
    catalog = loci_by_mechanism_for("ciprofloxacin")
    # QRDR has the 4 canonical target genes
    assert catalog["QRDR_target_alteration"] == {"gyrA", "gyrB", "parC", "parE"}
    # Plasmid-protect includes qnr family + aac(6')-Ib-cr variants
    assert "qnrB" in catalog["plasmid_protect_modify"]
    assert "qnrS" in catalog["plasmid_protect_modify"]
    assert "aac(6')-Ib-cr" in catalog["plasmid_protect_modify"]
    # Co-resistance modifiers
    assert "acrB" in catalog["efflux"]
    assert "tolC" in catalog["efflux"]
    assert "ompC" in catalog["porin_loss"]
    assert "marR" in catalog["regulatory"]
    assert "acrR" in catalog["regulatory"]


def test_cef_catalog_includes_textbook_beta_lactamases():
    from dna_decode.data.mic_tiers import loci_by_mechanism_for
    catalog = loci_by_mechanism_for("ceftriaxone")
    # ESBL families that hydrolyze 3rd-gen cephalosporins
    assert "blaCTX-M" in catalog["acquired_beta_lactamase"]
    assert "blaCMY" in catalog["acquired_beta_lactamase"]
    assert "blaTEM" in catalog["acquired_beta_lactamase"]
    # Carbapenemases (also affect cef)
    assert "blaKPC" in catalog["acquired_beta_lactamase"]
    assert "blaNDM" in catalog["acquired_beta_lactamase"]
    # Chromosomal ampC hyperproduction
    assert "ampC" in catalog["ampC_hyperproduction"]


def test_tet_catalog_includes_efflux_and_rpps():
    from dna_decode.data.mic_tiers import loci_by_mechanism_for
    catalog = loci_by_mechanism_for("tetracycline")
    # Tet efflux family — the canonical distributed-mechanism dataset
    assert "tetA" in catalog["tet_efflux"]
    assert "tetB" in catalog["tet_efflux"]
    # Ribosomal protection proteins
    assert "tetM" in catalog["tet_ribosomal_protection"]
    assert "tetO" in catalog["tet_ribosomal_protection"]
    # Enzymatic inactivation
    assert "tetX" in catalog["tet_enzymatic"]


def test_gent_catalog_includes_modifying_enzymes_and_rmt():
    from dna_decode.data.mic_tiers import loci_by_mechanism_for
    catalog = loci_by_mechanism_for("gentamicin")
    # Acetyl + phospho + nucleotidyl transferases
    assert "aac(3)-IIa" in catalog["aminoglycoside_modifying_enzymes"]
    assert "aph(2'')-Ia" in catalog["aminoglycoside_modifying_enzymes"]
    assert "ant(2'')-Ia" in catalog["aminoglycoside_modifying_enzymes"]
    # 16S methyltransferases — high-level aminoglycoside resistance
    assert "armA" in catalog["16S_rRNA_methyltransferase"]
    assert "rmtB" in catalog["16S_rRNA_methyltransferase"]


def test_loci_unknown_drug_raises():
    from dna_decode.data.mic_tiers import loci_by_mechanism_for
    with pytest.raises(UnknownDrugError):
        loci_by_mechanism_for("unobtainium")


# ---- Per-drug primary mechanisms ----------------------------------------


def test_cipro_primary_mechanisms_match_phenotype_merge_script():
    """Primary mechanisms must match scripts/cipro_mechanism_phenotype_merge.py."""
    from dna_decode.data.mic_tiers import primary_mechanisms_for
    prim = primary_mechanisms_for("ciprofloxacin")
    assert prim == frozenset({"QRDR_target_alteration", "plasmid_protect_modify"})


def test_cef_primary_mechanisms():
    from dna_decode.data.mic_tiers import primary_mechanisms_for
    prim = primary_mechanisms_for("ceftriaxone")
    assert "acquired_beta_lactamase" in prim
    assert "ampC_hyperproduction" in prim


def test_tet_primary_mechanisms_include_efflux_and_rpps():
    from dna_decode.data.mic_tiers import primary_mechanisms_for
    prim = primary_mechanisms_for("tetracycline")
    assert "tet_efflux" in prim
    assert "tet_ribosomal_protection" in prim


def test_gent_primary_mechanisms():
    from dna_decode.data.mic_tiers import primary_mechanisms_for
    prim = primary_mechanisms_for("gentamicin")
    assert "aminoglycoside_modifying_enzymes" in prim
    assert "16S_rRNA_methyltransferase" in prim


def test_co_resistance_mechanisms_shared_across_drugs():
    """Efflux + regulatory + porin_loss are co-resistance modifiers, not primary."""
    from dna_decode.data.mic_tiers import (
        CO_RESISTANCE_MECHANISMS,
        DRUG_PRIMARY_MECHANISMS,
    )
    assert CO_RESISTANCE_MECHANISMS == frozenset({"efflux", "regulatory", "porin_loss"})
    # Primary + co-resistance must be disjoint for every drug
    for drug, primary in DRUG_PRIMARY_MECHANISMS.items():
        assert primary.isdisjoint(CO_RESISTANCE_MECHANISMS), \
            f"{drug}: primary mechanisms overlap co-resistance"


# ---- classify_gene_symbol -----------------------------------------------


def test_classify_gene_symbol_cipro_qrdr():
    from dna_decode.data.mic_tiers import classify_gene_symbol
    assert classify_gene_symbol("ciprofloxacin", "gyrA") == "QRDR_target_alteration"
    assert classify_gene_symbol("ciprofloxacin", "parC") == "QRDR_target_alteration"


def test_classify_gene_symbol_cipro_strips_mutation_suffix():
    """gyrA_S83L should classify as QRDR (prefix split on underscore)."""
    from dna_decode.data.mic_tiers import classify_gene_symbol
    assert classify_gene_symbol("ciprofloxacin", "gyrA_S83L") == "QRDR_target_alteration"


def test_classify_gene_symbol_cipro_prefix_match():
    """qnrB19 should match qnrB as a prefix."""
    from dna_decode.data.mic_tiers import classify_gene_symbol
    assert classify_gene_symbol("ciprofloxacin", "qnrB19") == "plasmid_protect_modify"


def test_classify_gene_symbol_cef_beta_lactamase_prefix():
    """blaCTX-M-15 should match blaCTX-M as a prefix."""
    from dna_decode.data.mic_tiers import classify_gene_symbol
    assert classify_gene_symbol("ceftriaxone", "blaCTX-M-15") == "acquired_beta_lactamase"


def test_classify_gene_symbol_tet_efflux():
    from dna_decode.data.mic_tiers import classify_gene_symbol
    assert classify_gene_symbol("tetracycline", "tetA") == "tet_efflux"
    assert classify_gene_symbol("tetracycline", "tetM") == "tet_ribosomal_protection"
    assert classify_gene_symbol("tetracycline", "tetX") == "tet_enzymatic"


def test_classify_gene_symbol_gent_amg_enzymes():
    from dna_decode.data.mic_tiers import classify_gene_symbol
    assert classify_gene_symbol("gentamicin", "aac(3)-IIa") == "aminoglycoside_modifying_enzymes"
    assert classify_gene_symbol("gentamicin", "armA") == "16S_rRNA_methyltransferase"


def test_classify_gene_symbol_unknown_returns_empty():
    from dna_decode.data.mic_tiers import classify_gene_symbol
    assert classify_gene_symbol("ciprofloxacin", "unknownGene") == ""
    assert classify_gene_symbol("ciprofloxacin", "") == ""


def test_classify_gene_symbol_unknown_drug_raises():
    from dna_decode.data.mic_tiers import classify_gene_symbol
    with pytest.raises(UnknownDrugError):
        classify_gene_symbol("unobtainium", "gyrA")
