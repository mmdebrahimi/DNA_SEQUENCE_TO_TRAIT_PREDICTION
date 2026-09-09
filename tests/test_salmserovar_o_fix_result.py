"""The O-port result, and the OVERRIDE it was shipped under.

The frozen bar returned REJECT on one clause of seven, by one isolate, while the change was strictly
dominant on the measured evidence. The code shipped anyway. These tests exist so that override can
never be quietly forgotten or reread as a pass: the artifact must keep saying REJECT, the memo must
keep calling the adoption an override, and the evidence that made the override defensible (zero
regressions) must keep being true of the committed numbers.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

ART = ROOT / "wiki" / "salmserovar_o_fix_result_2026-09-09.json"
MEMO = ROOT / "wiki" / "salmserovar_o_fix_result_2026-09-09.md"
BAR = ROOT / "wiki" / "salmserovar_o_fix_acceptance_bar.json"
READING = ROOT / "wiki" / "seqsero2_o_algorithm_reading_2026-09-09.md"


@pytest.fixture(scope="module")
def art() -> dict:
    if not ART.exists():
        pytest.skip("result artifact absent")
    return json.loads(ART.read_text(encoding="utf-8"))


def test_the_verdict_field_still_says_reject(art):
    """The bar was NOT edited to fit the result and the verdict was NOT rewritten. If either had
    happened the pre-registration would have been theatre, so this is the load-bearing assertion."""
    assert art["verdict"] == "REJECT"
    assert any(not c["pass"] for c in art["bar_checks"])
    assert "override" in art["verdict_is_mechanical"].lower()


def test_exactly_one_clause_failed_and_it_is_the_no_call_ceiling(art):
    failed = [c for c in art["bar_checks"] if not c["pass"]]
    assert len(failed) == 1, failed
    assert "no_call" in failed[0]["check"]


def test_the_change_is_strictly_dominant_which_is_what_licensed_the_override(art):
    """Zero `hit -> anything else`. Without this the override has no argument at all, so it is pinned
    against the committed numbers rather than left in prose."""
    assert art["n_regressions"] == 0
    for k in art["transitions"]:
        assert not (k.startswith("hit->") and k != "hit->hit"), k


def test_the_no_call_rise_is_net_and_its_parts_are_recorded(art):
    """The clause is a NET ceiling policing a DIRECTIONAL failure mode -- which is the defect. The
    parts must stay visible: wrong calls becoming abstentions, AND abstentions becoming hits."""
    t = art["transitions"]
    assert t.get("miss->no_call", 0) > 0
    assert t.get("no_call->hit", 0) > 0
    delta = art["ours_after"]["no_call"] - art["ours_before"]["no_call"]
    assert delta == t.get("miss->no_call", 0) - t.get("no_call->hit", 0)


def test_the_gain_is_correct_calls_not_abstentions(art):
    """The bar's own `explicitly_not_a_success_signal` forbids buying agreement with abstentions.
    The rescue must dominate the abstention movement by a wide margin, or the override is unsound."""
    t = art["transitions"]
    assert t.get("miss->hit", 0) >= 4 * t.get("miss->no_call", 0)
    assert art["ours_after"]["accuracy"] > art["ours_before"]["accuracy"]
    assert art["ours_after"]["miss"] < art["ours_before"]["miss"]


def test_the_H_axes_did_not_move_at_all(art):
    """H1/H2 headers were left byte-identical, so the change CANNOT reach them. A move here would mean
    the DB rebuild did something it was not supposed to."""
    for axis in ("H1", "H2"):
        assert art["per_axis_vs_seqsero2"][axis]["before"] == art["per_axis_vs_seqsero2"][axis]["after"]


def test_the_o_axis_improved_far_past_the_bar(art):
    o = art["per_axis_vs_seqsero2"]["O"]
    assert o["after"]["agree"] - o["before"]["agree"] >= 15
    assert o["after"]["differ"] < o["before"]["differ"]


def test_the_cell_is_not_claimed_to_have_caught_the_reference_tool(art):
    """Better, not fixed. An artifact that lost this would be over-claiming."""
    assert art["delta_vs_reference_tool_after"] < 0
    assert art["delta_vs_reference_tool_after"] > art["delta_vs_reference_tool_before"]


def test_the_memo_calls_the_adoption_an_override_not_a_pass():
    if not MEMO.exists():
        pytest.skip("memo absent")
    t = " ".join(MEMO.read_text(encoding="utf-8").lower().split())
    assert "override" in t
    assert "the bar is **not edited**" in t or "the bar is not edited" in t
    assert "still reads `reject`" in t or "still reads reject" in t


def test_the_memo_keeps_the_named_cost_of_refusing_the_tyr_refinement():
    """A true O-2 genome will be called O-9. That is a real clinical limitation (Paratyphi A) and must
    stay in the memo, not only in a code comment."""
    if not MEMO.exists():
        pytest.skip("memo absent")
    t = " ".join(MEMO.read_text(encoding="utf-8").lower().split())
    assert "paratyphi a" in t
    assert "will be called o-9" in t
    assert "no o-2 serovar" in t          # why the positive direction was untestable here


def test_the_memo_records_that_the_port_failed_first_time():
    """Run 1 scored 0.7300 and was rejected on two clauses. Presenting only the second run would make
    a two-attempt result look like a first-attempt one."""
    if not MEMO.exists():
        pytest.skip("memo absent")
    t = " ".join(MEMO.read_text(encoding="utf-8").lower().split())
    assert "did **not** work first time" in t or "did not work first time" in t
    assert "0.7300" in t


def test_the_reading_memo_is_marked_partly_superseded():
    """Its threshold-translation argument was measured and failed for the paralog case. The memo is
    kept verbatim as the pre-implementation record, so the warning has to be a header on it."""
    if not READING.exists():
        pytest.skip("reading memo absent")
    t = " ".join(READING.read_text(encoding="utf-8").lower().split())
    assert "partly superseded" in t
    assert "salmserovar_o_fix_result_2026-09-09.md" in t
    # the original text must SURVIVE -- the bar was frozen against it
    assert "not a measured one" in t
    assert "correction 1" in t and "correction 2" in t


def test_the_frozen_bar_file_was_not_edited_to_fit_the_result():
    if not BAR.exists():
        pytest.skip("bar absent")
    bar = json.loads(BAR.read_text(encoding="utf-8"))
    assert bar["status"] == "REGISTERED_BEFORE_IMPLEMENTATION"
    assert bar["reject_if_ANY_fails"] is True
    assert bar["baseline_2026-09-08"]["ours_overall"]["no_call"] == 20
