"""The guard: a doubt block can NEVER contain a resistance call.

This is the constraint that keeps L2 out of the learned-predictor regime (0-for-5, de-confounded).
It is enforced inside `DoubtBlock.as_dict` rather than only asserted here, so a bug raises at emit
time instead of shipping a product-surface falsehood.

A guard that cannot fail proves nothing, so every test below either constructs a violation and
requires the raise, or pins that a legitimate block survives it.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from dna_decode.eval.doubt import (DoubtBlock, DoubtSignal, assert_no_call,  # noqa: E402
                                   completeness_signal, position_novelty_signal)


# --- the guard actually catches things (non-vacuity) ---

def test_a_call_shaped_key_is_refused():
    bad = DoubtBlock([DoubtSignal("k", "weak", "a reason", {"prediction": "resistant"})])
    with pytest.raises(ValueError, match="call-shaped key"):
        bad.as_dict()


def test_a_call_shaped_value_is_refused_even_under_an_innocent_key():
    """The dangerous case: no suspicious key name, just an R sitting in the evidence."""
    bad = DoubtBlock([DoubtSignal("k", "weak", "a reason", {"observed_outcome": "R"})])
    with pytest.raises(ValueError, match="call-shaped value"):
        bad.as_dict()


def test_the_guard_reaches_inside_lists():
    bad = DoubtBlock([DoubtSignal("k", "weak", "a reason", {"per_isolate": [{"phenotype": "S"}]})])
    with pytest.raises(ValueError, match="call-shaped key"):
        bad.as_dict()


def test_the_guard_reaches_arbitrarily_deep():
    with pytest.raises(ValueError, match="call-shaped value"):
        assert_no_call({"a": {"b": [{"c": ["SUSCEPTIBLE"]}]}})


def test_key_match_is_exact_not_substring():
    """`carriers_labelled_r` is a cohort statistic, not a call about this sample -- it must pass,
    or the guard would forbid the evidence the doubt layer exists to carry."""
    assert_no_call({"carriers_labelled_r": 36, "carriers_labelled_s": 0})


def test_a_substitution_containing_R_is_not_mistaken_for_a_call():
    """`K103R` contains an R; only an exact match is a call."""
    assert_no_call({"novel_substitutions": ["K103R", "V179D"]})


def test_a_bare_amino_acid_letter_would_be_refused_and_that_is_deliberate():
    """Known, accepted tradeoff: a lone `S`/`R`/`I` is indistinguishable from a call, so the guard
    refuses it. No shipped evidence field carries a bare residue letter; erring toward refusal is
    the safe direction, and pinning it here keeps the tradeoff explicit rather than a surprise."""
    with pytest.raises(ValueError, match="call-shaped value"):
        assert_no_call({"wildtype_residue": "S"})


# --- legitimate blocks survive ---

def test_a_real_completeness_block_emits_cleanly_and_carries_no_prediction():
    block = DoubtBlock([completeness_signal("rmtE1", "AMINOGLYCOSIDE", 36, 0, 0.517, 131)])
    d = block.as_dict()
    assert d["schema"] == "decoder-doubt-v1"
    assert d["any_doubt"] is True and d["max_tier"] == "strong"
    assert "prediction" not in d and "call" not in d
    assert "NEVER" in d["contract"]


def test_a_real_position_novelty_block_emits_cleanly():
    DoubtBlock([position_novelty_signal(["K103W"], "hiv-nnrti-rt")]).as_dict()


def test_an_empty_block_is_honest_about_having_no_doubt():
    d = DoubtBlock([]).as_dict()
    assert d["any_doubt"] is False and d["max_tier"] == "none" and d["signals"] == []


def test_max_tier_reports_the_strongest_signal_present():
    block = DoubtBlock([
        completeness_signal("weak-one", "C", 4, 0, 0.583, 125),
        completeness_signal("rmtE1", "AMINOGLYCOSIDE", 36, 0, 0.517, 131),
    ])
    assert block.max_tier == "strong"
