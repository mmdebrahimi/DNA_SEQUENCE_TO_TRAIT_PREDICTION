"""Tests for the HIV-1 NNRTI determinant catalog (dna_decode/data/hiv_amr) — Wave B v0.

Pure-logic, no BLAST / no network / no dataset. Pins the catalog (sourced verbatim from the Stanford
HIVDB genotype-phenotype dataset page's NNRTI Major DR Positions) + the deterministic call shape. The
dataset-based VALIDATION (PhenoSense fold-change, within-subtype) is a separate step that needs the
downloaded file.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from dna_decode.data.hiv_amr import (  # noqa: E402
    HIV_UNDETECTABLE_MECHANISMS, NNRTI_RT_MAJOR_DRMS, NRTI_MAJOR_POSITIONS, NRTI_UNDETECTABLE_MECHANISMS,
    _RT_WT, call_from_observed_substitutions, call_nrti_from_observed, gene_for_drug,
    is_nrti_major_position_mutation, is_resistance_mutation, resistance_mutations_for,
    supported_hiv_drugs, supported_nrti_drugs,
)


def test_supported_drugs():
    drugs = set(supported_hiv_drugs())
    # the major NNRTIs the dataset carries fold-change for
    assert {"efavirenz", "nevirapine", "etravirine", "rilpivirine"} <= drugs


def test_catalog_matches_stanford_nnrti_major_positions():
    # VERBATIM from the Stanford HQ-filtered dataset page (2026-04-22): 100I, 101P, 103N, 106A/M,
    # 181C/I/A, 188C/L/H, 190A/E/S/Q, 230L (consensus-B wild-types).
    expected = {
        "L100I", "K101P", "K103N", "V106A", "V106M",
        "Y181C", "Y181I", "Y181A", "Y188C", "Y188L", "Y188H",
        "G190A", "G190E", "G190S", "G190Q", "M230L",
    }
    assert NNRTI_RT_MAJOR_DRMS == expected


def test_position_188_is_tyrosine_not_glycine():
    # The provenance catch: consensus-B wild-type at RT-188 is Y (Tyrosine), so the DRMs are Y188C/L/H,
    # NOT the "G188" a web summarizer garbled. Pin it so a future edit can't silently regress.
    assert _RT_WT[188] == "Y"
    for mut in ("Y188C", "Y188L", "Y188H"):
        assert mut in NNRTI_RT_MAJOR_DRMS
    assert not any(m.startswith("G188") for m in NNRTI_RT_MAJOR_DRMS)


def test_is_resistance_mutation():
    assert is_resistance_mutation("efavirenz", "RT", "K103N") is True
    assert is_resistance_mutation("efavirenz", "RT", "Y188L") is True
    assert is_resistance_mutation("efavirenz", "RT", "Q999Z") is False
    assert is_resistance_mutation("efavirenz", "PR", "K103N") is False  # wrong gene


def test_call_resistant_on_major_drm():
    c = call_from_observed_substitutions("efavirenz", {"RT": {"K103N"}})
    assert c.prediction == "R" and c.determinants == ["RT:K103N"]
    assert c.undetectable_mechanisms == []
    assert "PhenoSense" in c.caveat and "Sierra" in c.caveat  # the circularity-safe instruction


def test_call_resistant_multiple_determinants_sorted():
    c = call_from_observed_substitutions("nevirapine", {"RT": {"Y181C", "G190A"}})
    assert c.prediction == "R"
    assert c.determinants == ["RT:G190A", "RT:Y181C"]  # sorted


def test_call_susceptible_surfaces_blind_spots():
    # no major DRM -> S, but the blind spots surface (S != definitely NNRTI-susceptible)
    c = call_from_observed_substitutions("efavirenz", {"RT": {"A98G"}})
    assert c.prediction == "S"
    assert c.undetectable_mechanisms == HIV_UNDETECTABLE_MECHANISMS
    assert "minor" in c.caveat.lower()


def test_class_level_v0_scope_documented_in_blind_spots():
    # v0 is class-level; per-drug differential resistance is an explicit named blind spot.
    assert "per_drug_differential_resistance" in HIV_UNDETECTABLE_MECHANISMS


def test_unknown_drug_indeterminate():
    c = call_from_observed_substitutions("notadrug", {"RT": {"K103N"}})
    assert c.prediction == "INDETERMINATE"


def test_gene_for_drug():
    assert gene_for_drug("efavirenz") == "RT"
    assert gene_for_drug("notadrug") is None


def test_resistance_mutations_for_unknown_raises():
    try:
        resistance_mutations_for("notadrug")
    except KeyError:
        return
    raise AssertionError("expected KeyError for unknown drug")


# ---------- NRTI (v0 position-based) ----------

def test_nrti_supported_drugs_and_positions():
    assert {"lamivudine", "abacavir", "zidovudine", "tenofovir"} <= set(supported_nrti_drugs())
    # Stanford NRTI major positions (2026-04-22)
    assert set(NRTI_MAJOR_POSITIONS) == {41, 65, 70, 74, 75, 151, 184, 210, 215}


def test_is_nrti_major_position_mutation():
    assert is_nrti_major_position_mutation("M184V") is True
    assert is_nrti_major_position_mutation("K65R") is True
    assert is_nrti_major_position_mutation("T215Y") is True
    assert is_nrti_major_position_mutation("A98G") is False   # 98 is not an NRTI major position


def test_call_nrti_resistant_on_major_position():
    c = call_nrti_from_observed("lamivudine", {"RT": {"M184V"}})
    assert c.prediction == "R" and c.determinants == ["RT:M184V"]
    c2 = call_nrti_from_observed("zidovudine", {"RT": {"T215Y", "M41L"}})
    assert c2.prediction == "R" and c2.determinants == ["RT:M41L", "RT:T215Y"]


def test_call_nrti_overcalls_t215_revertant_by_design():
    # POSITION-BASED v0 DELIBERATELY calls a T215 revertant R (it's non-consensus at a major position).
    # This documents the known over-call the validation quantifies — NOT a bug.
    c = call_nrti_from_observed("zidovudine", {"RT": {"T215S"}})
    assert c.prediction == "R"
    assert "over-calls" in c.caveat


def test_call_nrti_susceptible_and_unknown():
    c = call_nrti_from_observed("lamivudine", {"RT": {"A98G"}})  # not at a major position
    assert c.prediction == "S" and c.undetectable_mechanisms == NRTI_UNDETECTABLE_MECHANISMS
    assert call_nrti_from_observed("notadrug", {"RT": {"M184V"}}).prediction == "INDETERMINATE"


# ---------- PI / INSTI / CAI (protease / integrase / capsid target-site classes) ----------

from dna_decode.data.hiv_amr import (  # noqa: E402
    CAI_CLASS, INSTI_CLASS, PI_CLASS, all_supported_hiv_drugs, call_hiv_observed,
    call_target_class, gene_for_hiv_drug, hiv_target_class_for, supported_cai_drugs,
    supported_insti_drugs, supported_pi_drugs,
)


def test_new_class_supported_drugs():
    assert "darunavir" in supported_pi_drugs() and "lopinavir" in supported_pi_drugs()
    assert "dolutegravir" in supported_insti_drugs() and "raltegravir" in supported_insti_drugs()
    assert supported_cai_drugs() == ["lenacapavir"]


def test_catalog_positions_sourced():
    # Stanford HIVDB major positions (PI/INSTI) + CAPELLA capsid emergent set (CAI). Pin so an edit can't
    # silently drift the sourced catalog.
    assert set(PI_CLASS.positions) == {30, 32, 33, 46, 47, 48, 50, 54, 76, 82, 84, 88, 90}
    assert set(INSTI_CLASS.positions) == {66, 92, 118, 138, 140, 143, 147, 148, 155, 263}
    assert set(CAI_CLASS.positions) == {56, 66, 67, 70, 74, 105, 107}


def test_gene_for_hiv_drug_routes_all_five_classes():
    assert gene_for_hiv_drug("efavirenz") == "RT"     # NNRTI
    assert gene_for_hiv_drug("lamivudine") == "RT"    # NRTI
    assert gene_for_hiv_drug("darunavir") == "PR"     # PI
    assert gene_for_hiv_drug("dolutegravir") == "IN"  # INSTI
    assert gene_for_hiv_drug("lenacapavir") == "CA"   # CAI
    assert gene_for_hiv_drug("notadrug") is None


def test_pi_position_based_call():
    c = call_hiv_observed("lopinavir", {"PR": {"V82A"}})
    assert c.prediction == "R" and c.determinants == ["PR:V82A"]
    # POSITION-BASED: any non-consensus residue at a major position counts (deliberate over-call)
    assert call_hiv_observed("lopinavir", {"PR": {"V82Z"}}).prediction == "R"
    s = call_hiv_observed("lopinavir", {"PR": {"L10I"}})  # 10 is not a major PI position
    assert s.prediction == "S" and s.undetectable_mechanisms


def test_insti_position_based_call():
    c = call_hiv_observed("dolutegravir", {"IN": {"R263K"}})
    assert c.prediction == "R" and c.determinants == ["IN:R263K"]
    assert call_hiv_observed("raltegravir", {"IN": {"N155H", "Q148H"}}).prediction == "R"


def test_cai_is_mutant_level_not_position_based():
    # CAI is MUTANT-LEVEL: a catalogued substitution -> R ...
    assert call_hiv_observed("lenacapavir", {"CA": {"M66I"}}).prediction == "R"
    assert call_hiv_observed("lenacapavir", {"CA": {"Q67H"}}).prediction == "R"
    # ... but a NON-catalogued residue at a major capsid position (capsid polymorphism) -> S
    # (this is the whole reason CAI is mutant-level: position-based over-called the resistance-enriched set).
    assert CAI_CLASS.major_drms is not None
    assert "A105V" not in CAI_CLASS.major_drms        # A105T/S are catalogued, A105V is not
    assert call_hiv_observed("lenacapavir", {"CA": {"A105V"}}).prediction == "S"
    assert {"M66I", "Q67H", "K70R", "N74D", "A105T", "T107N"} <= CAI_CLASS.major_drms


def test_call_target_class_unknown_drug_indeterminate():
    assert call_target_class("notadrug", {"PR": {"V82A"}}, PI_CLASS).prediction == "INDETERMINATE"


def test_hiv_target_class_for_only_new_classes():
    assert hiv_target_class_for("darunavir") is PI_CLASS
    assert hiv_target_class_for("efavirenz") is None     # NNRTI not in the PI/INSTI/CAI registry
    assert hiv_target_class_for("lenacapavir") is CAI_CLASS


def test_all_supported_hiv_drugs_covers_five_classes():
    drugs = set(all_supported_hiv_drugs())
    assert {"efavirenz", "lamivudine", "darunavir", "dolutegravir", "lenacapavir"} <= drugs


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for fn in fns:
        fn(); print(f"PASS {fn.__name__}")
    print(f"\n{len(fns)} passed")
