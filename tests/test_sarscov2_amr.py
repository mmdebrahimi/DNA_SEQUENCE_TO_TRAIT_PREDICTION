"""Offline tests for the SARS-CoV-2 Mpro resistance catalog (dna_decode/data/sarscov2_amr.py)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.data.sarscov2_amr import (  # noqa: E402
    MPRO_MAJOR_DRMS, MPRO_WT, SARSCOV2_MPRO_DRUGS, all_supported_sarscov2_drugs,
    call_sarscov2_observed, gene_for_sarscov2_drug, is_mpro_major_drm, supported_sarscov2_drugs,
)


def test_catalog_has_canonical_nirmatrelvir_set():
    for s in ("E166V", "E166A", "L50F", "S144A", "A173V", "H172Y", "L167F", "T21I"):
        assert is_mpro_major_drm(s), f"{s} should be a catalogued Mpro major DRM"
    # WT/pos/mut self-consistency: every catalogued mutation's WT matches MPRO_WT at that position
    for s in MPRO_MAJOR_DRMS:
        wt, mut = s[0], s[-1]
        pos = int("".join(c for c in s[1:-1] if c.isdigit()))
        assert MPRO_WT[pos] == wt, f"{s}: WT {wt} != MPRO_WT[{pos}]={MPRO_WT[pos]}"
        assert mut != wt, f"{s}: mutant equals WT"


def test_call_R_on_major_substitution():
    c = call_sarscov2_observed("nirmatrelvir", {"Mpro": {"E166V"}})
    assert c.prediction == "R"
    assert c.determinants == ["Mpro:E166V"]
    assert c.undetectable_mechanisms == []        # R call surfaces no blind spots


def test_call_S_on_benign_omicron_polymorphism():
    """P132H is fixed in Omicron and is NOT nirmatrelvir-resistant -> mutant-level catalog must NOT call it R
    (the reason Mpro is mutant-level, not position-based)."""
    c = call_sarscov2_observed("nirmatrelvir", {"Mpro": {"P132H"}})
    assert c.prediction == "S"
    assert len(c.undetectable_mechanisms) == 6    # S call surfaces the blind spots


def test_call_S_on_empty_and_indeterminate_on_unknown_drug():
    assert call_sarscov2_observed("ensitrelvir", {"Mpro": set()}).prediction == "S"
    assert call_sarscov2_observed("aspirin", {"Mpro": {"E166V"}}).prediction == "INDETERMINATE"


def test_class_level_all_mpro_drugs_share_the_set():
    for drug in SARSCOV2_MPRO_DRUGS:
        assert call_sarscov2_observed(drug, {"Mpro": {"E166V"}}).prediction == "R"
        assert gene_for_sarscov2_drug(drug) == "Mpro"
    assert gene_for_sarscov2_drug("ciprofloxacin") is None
    assert set(supported_sarscov2_drugs()) == set(all_supported_sarscov2_drugs()) == set(SARSCOV2_MPRO_DRUGS)


def test_caveat_names_independent_validation_and_circularity():
    c = call_sarscov2_observed("nirmatrelvir", {"Mpro": {"E166V"}})
    assert "CoV-RDB measured fold-change" in c.caveat
    assert "circular" in c.caveat.lower()
    assert c.rule == "sarscov2_mpro_major_drm_v0"


if __name__ == "__main__":
    import pytest
    sys.exit(pytest.main([__file__, "-q"]))
