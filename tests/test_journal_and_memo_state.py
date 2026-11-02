"""No published surface may claim the ceiling criterion is `met` while an arm is untried.

The positional arm's memo opened with *"The ceiling criterion now reads `met`, and whether to accept
that is yours."* That was true when it was rendered and became FALSE the moment three successor arms
were declared on 2026-11-02. A memo is generated FROM its artifacts plus the live verdict, so the fix
is a re-render and never a hand edit -- and these tests are what make the staleness detectable rather
than something a reader has to notice.

The journal half pins the two structural properties that a hand-appended row gets wrong: the row must
sit INSIDE the Action-Log anchor (a row appended after it would be picked up as the numbering anchor
by the next writer) and its class must come from the ledger's own 7-value enum.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
MEMO = REPO / "wiki/positional_arm_result_2026-11-01.md"
RESULT = REPO / "wiki/positional_arm_result_2026-11-01.json"
LEDGER = REPO / "soraya_runs/arms-ledger.json"
JOURNAL = REPO / "project_state/glm-learned-representation-2026-10-07.md"

ANCHOR = "<!-- project-state:end:action-log -->"
#: The enum the ledger's own "Allowed Action Classes" section defines.
ACTION_CLASSES = {"propose", "research", "write-plan", "edit-local-code", "run-tests",
                  "ask-user", "stop"}


def _planned() -> set[str]:
    d = json.loads(LEDGER.read_text(encoding="utf-8"))
    crit = next(k for k in d["criteria"] if k.startswith("ceiling:"))
    return {n for n, a in d["criteria"][crit].items() if a.get("status") == "planned"}


# --------------------------------------------------------------------------------------------------
# the memo
# --------------------------------------------------------------------------------------------------
@pytest.mark.skipif(not MEMO.is_file(), reason="memo not rendered")
def test_the_memo_does_not_claim_met_while_an_arm_is_untried():
    md = MEMO.read_text(encoding="utf-8")
    if not _planned():
        pytest.skip("every arm is measured; the criterion may legitimately read met")
    assert "now reads `met`" not in md, (
        "the memo still carries the stale `met` claim; re-render it rather than editing it")
    assert "blocked:untried-arm" in md


@pytest.mark.skipif(not MEMO.is_file(), reason="memo not rendered")
def test_the_memo_names_the_unrun_arms_rather_than_leaving_them_to_be_inferred():
    planned = _planned()
    if not planned:
        pytest.skip("nothing is planned")
    decisions = MEMO.read_text(encoding="utf-8").split("## CORRECTIONS", 1)[0]
    for name in sorted(planned):
        assert name in decisions, (
            "%s is declared-but-unrun and the DECISIONS section does not name it" % name)


@pytest.mark.skipif(not (MEMO.is_file() and RESULT.is_file()), reason="memo or result absent")
def test_the_memo_on_disk_matches_a_fresh_render():
    """No figure may be hand-edited in after generation."""
    from scripts.freeze_positional_arm_bar import OUT as BAR_OUT
    from scripts.positional_arm_memo import render

    res = json.loads(RESULT.read_text(encoding="utf-8"))
    bar = json.loads(BAR_OUT.read_text(encoding="utf-8"))
    on_disk = MEMO.read_text(encoding="utf-8")
    # the DECISIONS section depends on the LIVE verdict, which `render(..., None)` cannot know, so
    # compare from `## The result` onward -- the half that is a pure function of the artifacts.
    fresh_tail = render(res, bar, None).split("## The result", 1)[1]
    assert on_disk.split("## The result", 1)[1].strip() == fresh_tail.strip()


@pytest.mark.skipif(not MEMO.is_file(), reason="memo not rendered")
def test_the_memo_still_carries_the_arms_own_verdict_and_corrections():
    """The re-render must not have cost the content the memo exists for."""
    md = MEMO.read_text(encoding="utf-8")
    res = json.loads(RESULT.read_text(encoding="utf-8"))
    assert res["verdict"] in md
    assert "0.5086" in md and "0.4893" in md, "the orientation finding must survive a re-render"
    i_dec, i_cor, i_res = (md.index("## DECISIONS NEEDED"), md.index("## CORRECTIONS"),
                           md.index("## The result"))
    assert i_dec < i_cor < i_res, "ordered by required reader action"


# --------------------------------------------------------------------------------------------------
# the journal
# --------------------------------------------------------------------------------------------------
def _action_rows(text: str) -> list[tuple[int, str]]:
    rows = []
    for m in re.finditer(r"^\| (\d+) \| \d{4}-\d{2}-\d{2} \| ([a-z-]+) \|", text, re.M):
        rows.append((int(m.group(1)), m.group(2)))
    return rows


def test_every_action_log_row_sits_INSIDE_the_anchor():
    """A row appended AFTER the anchor would be read as the numbering anchor by the next writer."""
    text = JOURNAL.read_text(encoding="utf-8")
    assert text.count(ANCHOR) == 1
    before, after = text.split(ANCHOR, 1)
    assert _action_rows(before), "no Action-Log rows found before the anchor"
    assert not _action_rows(after), (
        "Action-Log row(s) found AFTER the anchor: %r" % _action_rows(after))


def test_every_action_class_is_from_the_seven_value_enum():
    rows = _action_rows(JOURNAL.read_text(encoding="utf-8").split(ANCHOR, 1)[0])
    bad = [(n, c) for n, c in rows if c not in ACTION_CLASSES]
    assert not bad, "row(s) with an off-enum action class: %r" % bad


def test_the_action_log_numbering_is_contiguous_from_one():
    nums = [n for n, _ in _action_rows(JOURNAL.read_text(encoding="utf-8").split(ANCHOR, 1)[0])]
    assert nums == sorted(nums), "rows are out of order"
    assert nums == list(range(1, len(nums) + 1)), (
        "numbering has a gap or a duplicate: %r" % nums)


def test_the_answered_authority_call_is_AMENDED_rather_than_deleted():
    """The answer plus its date IS the record; deleting the question would erase that a decision
    was made."""
    text = JOURNAL.read_text(encoding="utf-8")
    assert "ACCEPT `met` as the ceiling criterion's terminal" in text, (
        "the authority call was deleted rather than amended")
    assert "[ANSWERED 2026-11-02]" in text
    # and the genuinely-new second decision is flagged rather than folded into the answered one
    assert "SECOND, GENUINELY NEW decision" in text


def test_the_journal_records_the_decision_that_restored_the_verdict():
    text = JOURNAL.read_text(encoding="utf-8")
    row21 = [b for b in text.split("\n") if b.startswith("| 21 | ")]
    assert len(row21) == 1
    body = row21[0]
    assert "blocked:untried-arm" in body and "met" in body
    assert "13 -> 16" in body, "the arm-count change must be recorded"
    assert "FOUR of the eight" in body, "the measured guard outcome must be recorded"
    assert "ZERO" in body, "the byte-identical re-derivation result must be recorded"
    assert "11008" in body, "the measured heads geometry lock must be recorded"
