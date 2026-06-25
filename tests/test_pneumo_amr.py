"""Pins the S. pneumoniae gene-presence AMR rule (non-frozen organism_rules cell). Pure-logic, always runs."""
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
from dna_decode.organism_rules.pneumo_amr import (  # noqa: E402
    GENE_PRESENCE_RULES, RULE_STATUS, call_drug, resolve_drug,
)


def test_branding_non_frozen():
    assert RULE_STATUS == "KNOWLEDGE_BASELINE"   # never confused with the frozen deployed surface


def test_macrolide_ermB_calls_R():
    c = call_drug("erythromycin", ["ermBL_HG799494; ermB_16_X82819"])
    assert c.prediction == "R" and "ermb" in c.matched_tokens


def test_macrolide_mefA_calls_R():
    c = call_drug("erythromycin", ["mefA_10_AF376746; msrD_2_AF274302"])
    assert c.prediction == "R" and ("mefa" in c.matched_tokens or "msrd" in c.matched_tokens)


def test_macrolide_no_determinant_calls_S():
    c = call_drug("erythromycin", ["_"])
    assert c.prediction == "S" and c.matched_tokens == ()


def test_tetracycline_tetM_calls_R():
    c = call_drug("tetracycline", ["tetM_12_FR671418"])
    assert c.prediction == "R" and "tetm" in c.matched_tokens


def test_tetracycline_no_determinant_calls_S():
    c = call_drug("tetracycline", ["_"])
    assert c.prediction == "S"


def test_macrolide_alias_routes():
    assert resolve_drug("azithromycin") == "erythromycin"
    assert call_drug("macrolide", ["ermB_16_X82819"]).prediction == "R"


def test_betalactam_out_of_scope():
    # β-lactams are deferred (PBP engine) — the gene-presence rule returns None, never a bogus call.
    assert call_drug("penicillin", ["pbp2x_..."]) is None
    assert resolve_drug("ceftriaxone") is None


def test_rules_cover_expected_drugs():
    assert set(GENE_PRESENCE_RULES) == {"erythromycin", "tetracycline"}


if __name__ == "__main__":
    import pytest
    raise SystemExit(pytest.main([__file__, "-q"]))
