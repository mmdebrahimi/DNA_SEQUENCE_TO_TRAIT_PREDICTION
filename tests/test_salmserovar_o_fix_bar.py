"""The O-fix acceptance bar, registered BEFORE the fix exists — and the reading that grounds it.

Two artifacts are pinned here. (1) `wiki/salmserovar_o_fix_acceptance_bar.json`: frozen 2026-09-09 with
row-level baseline counts, so a differential-O-typing change cannot be scored against a bar tuned to its
own result. (2) `wiki/seqsero2_o_algorithm_reading_2026-09-09.md`: the source reading that replaced
inferring rules from header names, which had already produced three wrong causal claims in this repo.

The specific trap the bar exists to close: an earlier draft guarded the H axes "where both resolve", a
denominator that SHRINKS when abstentions rise. A fix that converted wrong O calls into no-calls would
have scored as an improvement.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
BAR = ROOT / "wiki" / "salmserovar_o_fix_acceptance_bar.json"
READING = ROOT / "wiki" / "seqsero2_o_algorithm_reading_2026-09-09.md"
SS2 = ROOT / "wiki" / "salmserovar_seqsero2_2026-09-08.json"


@pytest.fixture(scope="module")
def bar() -> dict:
    if not BAR.exists():
        pytest.skip("acceptance bar absent")
    return json.loads(BAR.read_text(encoding="utf-8"))


def test_the_bar_is_registered_before_any_implementation(bar):
    assert bar["status"] == "REGISTERED_BEFORE_IMPLEMENTATION"
    assert bar["frozen"] == "2026-09-09"


def test_the_baseline_matches_the_measured_artifact(bar):
    """A bar whose baseline drifts from the artifact it was derived from is not a baseline."""
    if not SS2.exists():
        pytest.skip("seqsero2 artifact absent")
    art = json.loads(SS2.read_text(encoding="utf-8"))
    b = bar["baseline_2026-09-08"]
    assert b["n_cohort"] == art["n_comparable"] == 200
    for k in ("hit", "miss", "no_call", "accuracy"):
        assert b["ours_overall"][k] == art["ours"][k]


def test_the_accuracy_denominator_is_wet_lab_not_tool_agreement(bar):
    """Porting SeqSero2's logic makes SeqSero2-agreement circular as a performance metric."""
    d = bar["accuracy_denominator"].lower()
    assert "wet-lab" in d and "never" in d
    assert "circular" in d


def test_the_abstention_conversion_trap_is_explicitly_disqualified(bar):
    joined = " ".join(bar["explicitly_not_a_success_signal"]).lower()
    assert "converting wrong calls into abstentions" in joined
    assert "restricted to post-change resolved rows" in joined


def test_the_bar_guards_both_H_axes_on_resolution_not_just_agreement(bar):
    """H2 has 54 unresolved rows at baseline; guarding agreement alone would ignore that."""
    conds = " ".join(bar["adopt_if_ALL_hold"]).lower()
    for axis in ("h1", "h2"):
        assert f"{axis} ours_unresolved does not rise" in conds
    assert bar["baseline_2026-09-08"]["per_axis_vs_seqsero2"]["H2"]["ours_unresolved"] == 54


def test_the_o_target_is_absolute_isolates_on_a_fixed_denominator(bar):
    conds = " ".join(bar["adopt_if_ALL_hold"])
    assert "+15" in conds and "159" in conds
    assert "not on rows where both resolve" in conds


def test_any_failure_rejects(bar):
    assert bar["reject_if_ANY_fails"] is True


# --- the reading that grounds it ------------------------------------------------------------------

def test_the_reading_records_both_corrections():
    if not READING.exists():
        pytest.skip("reading artifact absent")
    t = " ".join(READING.read_text(encoding="utf-8").lower().split())
    assert "correction 1" in t and "correction 2" in t
    assert "overstated" in t              # my own claim, corrected
    assert "not unimplementable" in t     # the review's C1, corrected


def test_the_reading_states_the_score_is_coverage_of_reference():
    """The translatability answer: the k-mer score is a coverage-like quantity, so the branch
    thresholds map onto BLAST coverage rather than being untranslatable."""
    t = " ".join(READING.read_text(encoding="utf-8").lower().split())
    assert "coverage-of-reference" in t
    assert "not a measured one" in t      # the mapping is argued, not validated


def test_the_reading_keeps_the_two_not_in_markers_asymmetric():
    t = " ".join(READING.read_text(encoding="utf-8").lower().split())
    assert "opposite" in t
    assert "never" in t and "positive call" in t
