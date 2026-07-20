"""Tests for the fungal AMR determinant catalog (dna_decode/data/fungal_amr) — EP-7 step 1.

Pure-logic, no BLAST / no network. Pins the catalog + the deterministic call shape.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from dna_decode.data.fungal_amr import (  # noqa: E402
    FUNGAL_UNDETECTABLE_MECHANISMS, call_from_observed_substitutions, is_resistance_mutation,
    resistance_mutations_for, supported_fungal_drugs,
)


def test_supported_drugs():
    drugs = set(supported_fungal_drugs())
    assert {"fluconazole", "voriconazole", "caspofungin", "micafungin"} <= drugs


def test_erg11_azole_hotspots_present():
    erg11 = resistance_mutations_for("fluconazole")["ERG11"]
    for s in ("Y132F", "K143R", "F126T"):   # C. auris clade-specific Lockhart substitutions
        assert s in erg11


def test_is_resistance_mutation():
    assert is_resistance_mutation("fluconazole", "ERG11", "Y132F") is True
    assert is_resistance_mutation("fluconazole", "ERG11", "Q999Z") is False
    assert is_resistance_mutation("fluconazole", "FKS1", "S639F") is False   # FKS1 is echinocandin, not azole


def test_call_resistant_on_erg11_hotspot():
    c = call_from_observed_substitutions("fluconazole", {"ERG11": {"Y132F"}})
    assert c.prediction == "R" and c.determinants == ["ERG11:Y132F"]
    assert c.undetectable_mechanisms == []


def test_call_susceptible_surfaces_blind_spots():
    # no catalogued mutation -> S, but efflux/aneuploidy blind spots surfaced (S != definitely susceptible)
    c = call_from_observed_substitutions("fluconazole", {"ERG11": {"A1B"}})
    assert c.prediction == "S"
    assert c.undetectable_mechanisms == FUNGAL_UNDETECTABLE_MECHANISMS
    assert "efflux" in c.caveat.lower()


def test_call_echinocandin_fks1():
    c = call_from_observed_substitutions("caspofungin", {"FKS1": {"S639F"}})
    assert c.prediction == "R" and c.determinants == ["FKS1:S639F"]


def test_unknown_drug_indeterminate():
    c = call_from_observed_substitutions("notadrug", {"ERG11": {"Y132F"}})
    assert c.prediction == "INDETERMINATE"


def test_causal_marker_is_high_confidence():
    c = call_from_observed_substitutions("fluconazole", {"ERG11": {"Y132F"}})
    assert c.prediction == "R" and c.confidence == "HIGH"
    assert c.lineage_only_determinants == ()


def test_clade_iv_haplotype_only_is_low_confidence():
    # K177R/N335S/E343D is the C. auris clade IV ERG11 background haplotype: non-discriminative
    # (identical genotype in an R + an S isolate on the AR Bank). R call preserved (sensitivity) but LOW.
    c = call_from_observed_substitutions("fluconazole", {"ERG11": {"K177R", "N335S", "E343D"}})
    assert c.prediction == "R"                       # sensitivity preserved -- no missed resistance
    assert c.confidence == "LOW_LINEAGE_ONLY"
    assert set(c.lineage_only_determinants) == {"ERG11:K177R", "ERG11:N335S", "ERG11:E343D"}
    assert "lineage-associated" in c.caveat.lower()


def test_causal_plus_lineage_stays_high():
    # a genuine causal marker present -> HIGH even if clade-background markers also present
    c = call_from_observed_substitutions("fluconazole", {"ERG11": {"Y132F", "E343D"}})
    assert c.prediction == "R" and c.confidence == "HIGH"


def test_susceptible_confidence_is_na():
    c = call_from_observed_substitutions("fluconazole", {"ERG11": set()})
    assert c.prediction == "S" and c.confidence == "NA"


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for fn in fns:
        fn(); print(f"PASS {fn.__name__}")
    print(f"\n{len(fns)} passed")
