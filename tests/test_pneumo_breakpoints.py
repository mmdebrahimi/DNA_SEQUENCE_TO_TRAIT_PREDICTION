"""Pins the non-frozen pneumo β-lactam breakpoint set + the breakpoint-context-aware classify."""
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
from dna_decode.data.pneumo_breakpoints import classify, breakpoints_for  # noqa: E402


def test_penicillin_breakpoint_ambiguity():
    # The SAME MIC (1.0) is R under meningitis (S<=0.06) but S under non-meningitis (S<=2) -- the
    # ambiguity that forces explicit breakpoint context.
    assert classify("penicillin", "meningitis", 1.0) == "R"
    assert classify("penicillin", "non_meningitis", 1.0) == "S"


def test_intermediate_band():
    assert classify("penicillin", "non_meningitis", 4.0) == "I"   # between S<=2 and R>=8


def test_unknown_context_or_drug_is_none():
    assert classify("penicillin", "no_such_context", 1.0) is None
    assert classify("not_a_drug", "meningitis", 1.0) is None
    assert classify("penicillin", "meningitis", None) is None


def test_breakpoints_present():
    assert breakpoints_for("ceftriaxone", "meningitis") == {"S": 0.5, "R": 2.0}
    assert breakpoints_for("meropenem", "non_meningitis") == {"S": 0.25, "R": 1.0}


if __name__ == "__main__":
    import pytest
    raise SystemExit(pytest.main([__file__, "-q"]))
