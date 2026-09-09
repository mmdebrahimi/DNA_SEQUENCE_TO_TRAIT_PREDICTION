"""The name-table rebuild, and the SECOND mis-specified bar it exposed.

Two independent design choices were conflated in one change and are now measured apart: canonical
NAMING (remove_list / rename_dict / bracket expansion / subspecies tie-break) versus the AMBIGUITY
POLICY for a key that is still contested afterwards. The naming half pays; the policy half costs.

These tests pin the shipped policy, the measured numbers, and the fact that a second frozen bar
returned REJECT on a strictly-dominant change -- because two overrides in one session is a signal about
the bar-writing step, not a licence to ignore bars.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.build_salmserovar_db import (  # noqa: E402
    PREFERRED_SUBSPECIES,
    build_table,
    expand_optional_factors,
)

ART = ROOT / "wiki" / "salmserovar_name_table_result_2026-09-09.json"
BAR = ROOT / "wiki" / "salmserovar_name_table_acceptance_bar.json"
MEMO = ROOT / "wiki" / "salmserovar_name_table_2026-09-09.md"


@pytest.fixture(scope="module")
def art() -> dict:
    if not ART.exists():
        pytest.skip("name-table result artifact absent")
    return json.loads(ART.read_text(encoding="utf-8"))


# --- the bracket expansion, which is why a bare `r` can reach Bovismorbificans -----------------------

def test_optional_factors_expand_to_every_reachable_spelling():
    """`r,[i]` brackets an OPTIONAL factor, so a genome presenting only `r` IS that serovar. Without
    expansion `('8','r','1,5')` and `('8','r,[i]','1,5')` are unrelated dict keys."""
    got = expand_optional_factors("r,[i]")
    assert "r" in got and "i,r" in got and "r,[i]" in got


def test_expansion_sorts_factors_because_the_scheme_is_not_alphabetical():
    """SeqSero2 sorts before comparing precisely because entries like `r,[i]` are out of order."""
    assert "i,r" in expand_optional_factors("r,[i]")
    assert "r,i" not in expand_optional_factors("r,[i]")


def test_expansion_is_identity_for_an_unbracketed_antigen():
    assert expand_optional_factors("1,5") == {"1,5"}
    assert expand_optional_factors("-") == {"-"}


def test_a_multi_optional_antigen_expands_to_every_subset():
    got = expand_optional_factors("g,[m],[s]")
    assert {"g", "g,m", "g,s", "g,m,s"} <= got


# --- the shipped policy -----------------------------------------------------------------------------

def test_the_default_ambiguity_policy_is_the_measured_dominant_one():
    """`omit` is the PRINCIPLED option and it measured strictly worse (165 hits, 3 regressions) than
    `first` (168 hits, 0 regressions). The default must be the one the evidence supports, and the
    docstring must carry both numbers so the trade is not re-litigated from memory."""
    import inspect
    assert inspect.signature(build_table).parameters["ambiguous_policy"].default == "first"
    doc = " ".join((build_table.__doc__ or "").split())
    assert "165 hits (+0), 3 hit->no_call REGRESSIONS" in doc
    assert "168 hits (+3), ZERO regressions" in doc


def test_the_subspecies_tiebreak_is_declared_as_ours_not_seqsero2s():
    """SeqSero2 filters by a subspecies it determines independently; we have no such signal. Presenting
    our weak I-preference as SeqSero2's algorithm would be an over-claim."""
    assert PREFERRED_SUBSPECIES == "I"
    src = (ROOT / "scripts" / "build_salmserovar_db.py").read_text(encoding="utf-8")
    assert "This is OUR heuristic, NOT SeqSero2's algorithm" in src


# --- the measured result ----------------------------------------------------------------------------

def test_the_change_is_strictly_dominant(art):
    """Zero regressions in EITHER direction, fewer confident errors, higher accuracy. This is what
    licenses shipping something the bar rejected -- without it there is no argument."""
    assert art["n_regressions"] == 0
    assert art["transitions"].get("hit->miss", 0) == 0
    assert art["transitions"].get("hit->no_call", 0) == 0
    assert art["ours_after"]["hit"] == 168
    assert art["ours_after"]["miss"] < 14
    assert art["ours_after"]["accuracy"] > 0.825


def test_the_verdict_still_says_reject(art):
    """The bar was not edited to fit the result and the verdict was not rewritten."""
    assert art["verdict"] == "REJECT"
    failed = [c for c in art["bar_checks"] if not c["pass"]]
    assert len(failed) == 1 and "hits rise" in failed[0]["check"]


def test_the_index_matched_control_separates_scoring_from_calling(art):
    """This change rebuilds the table that `load_formula_index` reads, so it moves the SCORING
    function too. +1 hit is equivalence, not calling, and must not be credited to the fix."""
    ctl = art["index_matched_control"]
    assert ctl["index_only_delta"] == 1
    assert ctl["previous_run_names_rescored_under_current_index"]["hit"] == 166
    caller_attributable = art["ours_after"]["hit"] - ctl["previous_run_names_rescored_under_current_index"]["hit"]
    assert caller_attributable == 2


def test_the_antigen_axes_did_not_move(art):
    """A name-table change must not touch the antigen calling at all."""
    for axis in ("H1", "H2"):
        assert art["per_axis_vs_seqsero2"][axis]["before"] == art["per_axis_vs_seqsero2"][axis]["after"]


def test_the_cell_still_trails_the_reference_tool(art):
    assert art["delta_vs_reference_tool_after"] < 0


# --- the meta-finding, which is the part worth keeping ----------------------------------------------

def test_the_memo_records_that_this_is_the_SECOND_misspecified_bar():
    """Two overrides in one session is a signal about the bar-writing step. If that observation is
    lost, the next session learns 'override when inconvenient' instead of 'analyse before you freeze'."""
    if not MEMO.exists():
        pytest.skip("memo absent")
    t = " ".join(MEMO.read_text(encoding="utf-8").lower().split())
    assert "second mis-specified bar" in t
    assert "unachievable by construction" in t
    assert "analysis has to come first" in t
    assert "not a habit to normalise" in t


def test_the_memo_quotes_the_caller_attributable_gain_not_just_the_raw_one():
    if not MEMO.exists():
        pytest.skip("memo absent")
    t = " ".join(MEMO.read_text(encoding="utf-8").lower().split())
    assert "quote +2 when attributing to the fix" in t


def test_the_bar_file_was_not_edited_to_fit_the_result():
    if not BAR.exists():
        pytest.skip("bar absent")
    bar = json.loads(BAR.read_text(encoding="utf-8"))
    assert bar["status"] == "REGISTERED_BEFORE_IMPLEMENTATION"
    assert bar["reject_if_ANY_fails"] is True
    assert any("+4" in c for c in bar["adopt_if_ALL_hold"])
