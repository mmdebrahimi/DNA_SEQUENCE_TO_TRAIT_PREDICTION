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
    HIV_UNDETECTABLE_MECHANISMS, NNRTI_RT_MAJOR_DRMS, _RT_WT, call_from_observed_substitutions,
    gene_for_drug, is_resistance_mutation, resistance_mutations_for, supported_hiv_drugs,
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


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for fn in fns:
        fn(); print(f"PASS {fn.__name__}")
    print(f"\n{len(fns)} passed")
