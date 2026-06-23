"""Offline tests for the AMR Portal feasibility pure logic (SIR binarize + alias-aware leakage + summarize)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.amr_portal_feasibility import binarize_sir, iso_leaked, summarize  # noqa: E402


def test_binarize_sir():
    assert binarize_sir("resistant") == "R"
    assert binarize_sir("non-susceptible") == "R"
    assert binarize_sir("susceptible") == "S"
    assert binarize_sir("susceptible-dose dependent") == "S"
    for v in ("intermediate", "", None, "not tested"):
        assert binarize_sir(v) is None


def test_iso_leaked_alias_aware_case_insensitive():
    leak = {"GCA_001.1", "SAMEA999", "ERS5"}
    assert iso_leaked({"samea999", "ERSX"}, leak) is True       # biosample matches (case-insensitive)
    assert iso_leaked({"gca_001.1"}, leak) is True              # assembly matches
    assert iso_leaked({"SAMN_NEW", "ERS_NEW", "GCA_NEW"}, leak) is False
    assert iso_leaked({"", None}, leak) is False                # empty aliases never leak


def test_summarize_disjoint_powering_and_dedupe():
    leak = {"SAMEA_LEAK"}
    rows = [
        # E. coli cipro: one disjoint R, one disjoint S, one leaked (excluded from disjoint)
        ("Escherichia coli", "ciprofloxacin", "SAMEA_A", "ERS_A", "GCA_A", "resistant"),
        ("Escherichia coli", "ciprofloxacin", "SAMEA_B", "ERS_B", "GCA_B", "susceptible"),
        ("Escherichia coli", "ciprofloxacin", "SAMEA_LEAK", "ERS_L", "GCA_L", "resistant"),
        # duplicate of A (same isolate+drug+sir) -> deduped, not double-counted
        ("Escherichia coli", "ciprofloxacin", "SAMEA_A", "ERS_A", "GCA_A", "resistant"),
        # an intermediate -> counted as intermediate_or_blank, not in total
        ("Escherichia coli", "ciprofloxacin", "SAMEA_C", "ERS_C", "GCA_C", "intermediate"),
    ]
    agg = summarize(rows, leak)
    c = agg[("Escherichia coli", "ciprofloxacin")]
    assert c["total"] == 3                       # A, B, LEAK (intermediate excluded)
    assert c["disjoint"] == 2 and c["disjoint_R"] == 1 and c["disjoint_S"] == 1
    assert c["leaked"] == 1
    assert c["intermediate_or_blank"] == 1


if __name__ == "__main__":
    import pytest
    sys.exit(pytest.main([__file__, "-q"]))
