"""CLAUDE.md's report-card figures must match the artifact they describe.

WHY THIS EXISTS
`CLAUDE.md` is loaded into every session, so a stale number there is repeated to every future reader.
On 2026-08-23 its headline trust claim read *"25 cells (6 SCORED / 4 NOT_CENSUSED / 1 UNDERPOWERED …)"*
while the artifact said **27 cells / 10 SCORED / 3 UNDERPOWERED / 0 NOT_CENSUSED** — and the SAME file
said "10 provenance-disjoint SCORED cells" two bullets away. It contradicted itself on its own headline,
and it UNDER-stated the validated surface by 4 SCORED cells.

Under-claiming is as much a trust-surface falsehood as over-claiming, and it was caught only because a
rebuild happened to run. This converts that luck into a test.

The repo already guards documentation drift this way (`test_every_wiki_artifact_a_contract_cites_actually_exists`,
`test_no_fba_memo_cites_a_json_artifact_that_does_not_exist`); this is the same pattern applied to the
number that matters most.

SCOPE: this asserts CLAUDE.md agrees with the COMMITTED artifact. It does not re-run the builder — that
would make a unit test depend on the whole scoring surface. If the artifact itself is stale, rebuild with
`uv run python scripts/build_validation_report_card.py` (read-only, exit 0 always).
"""
from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

import pytest

ARTIFACT = Path("wiki/decoder_validation_report_card.json")
CLAUDE_MD = Path("CLAUDE.md")


def _states() -> Counter:
    d = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    rows = d.get("cells") or d.get("rows") or []
    if not rows:
        pytest.skip("report card artifact carries no cell rows")
    return Counter(r.get("state") or r.get("cell_state") for r in rows)


@pytest.mark.skipif(not ARTIFACT.exists() or not CLAUDE_MD.exists(),
                    reason="report card artifact or CLAUDE.md absent")
def test_claude_md_states_the_correct_total_cell_count():
    total = sum(_states().values())
    text = CLAUDE_MD.read_text(encoding="utf-8")
    assert f"{total} cells" in text, (
        f"CLAUDE.md does not state the artifact's total of {total} cells. "
        f"Re-derive with `uv run python scripts/build_validation_report_card.py` and update the "
        f"report-card bullet in CLAUDE.md."
    )


@pytest.mark.skipif(not ARTIFACT.exists() or not CLAUDE_MD.exists(),
                    reason="report card artifact or CLAUDE.md absent")
def test_claude_md_states_the_correct_count_for_every_populated_state():
    """Every state the artifact actually has must appear in CLAUDE.md with the right number.

    Only populated states are checked: an empty bucket need not be enumerated (though the current text
    does mention `0 NOT_CENSUSED`, because an EMPTY NOT_CENSUSED bucket is itself the meaningful claim —
    no shipped decoder is rendering invisibly).
    """
    text = CLAUDE_MD.read_text(encoding="utf-8")
    wrong = []
    for state, n in sorted(_states().items()):
        if n and f"{n} {state}" not in text:
            wrong.append(f"{n} {state}")
    assert not wrong, (
        "CLAUDE.md's report-card figures disagree with wiki/decoder_validation_report_card.json. "
        f"Missing or wrong: {wrong}. This is the exact drift that made CLAUDE.md claim 6 SCORED when "
        "the artifact said 10. Re-derive and update the report-card bullet."
    )


@pytest.mark.skipif(not CLAUDE_MD.exists(), reason="CLAUDE.md absent")
def test_claude_md_does_not_contradict_itself_on_the_scored_count():
    """The 2026-08-23 defect was not only staleness — the file asserted BOTH 6 and 10 SCORED.

    Any two different `<n> SCORED` claims in one file means at least one is wrong, whichever the
    artifact happens to agree with.
    """
    import re
    text = CLAUDE_MD.read_text(encoding="utf-8")
    # Ignore the explanatory correction note, which deliberately quotes the old wrong figure.
    text = re.sub(r"\*\*re-derived from the artifact.*?\*\*", "", text, flags=re.S)
    counts = set(re.findall(r"(\d+) SCORED\b", text))
    assert len(counts) <= 1, (
        f"CLAUDE.md states conflicting SCORED-cell counts: {sorted(counts)}. "
        "Exactly one can be right; re-derive from wiki/decoder_validation_report_card.json."
    )
