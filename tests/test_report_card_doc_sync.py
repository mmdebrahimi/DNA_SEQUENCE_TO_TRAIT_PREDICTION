"""Every LIVE doc's report-card figures must match the artifact they describe.

WHY THIS EXISTS
`CLAUDE.md` is loaded into every session, so a stale number there is repeated to every future reader.
On 2026-08-23 its headline trust claim read *"25 cells (6 SCORED / 4 NOT_CENSUSED / 1 UNDERPOWERED …)"*
while the artifact said **27 cells / 10 SCORED / 3 UNDERPOWERED / 0 NOT_CENSUSED** — and the SAME file
said "10 provenance-disjoint SCORED cells" two bullets away. It contradicted itself on its own headline,
and it UNDER-stated the validated surface by 4 SCORED cells.

Under-claiming is as much a trust-surface falsehood as over-claiming, and it was caught only because a
rebuild happened to run. This converts that luck into a test.

WHY IT NOW COVERS MORE THAN CLAUDE.md (2026-08-26)
The guard hand-listed `CLAUDE.md` as a single constant -- so `docs/ARCHITECTURE.md` went on carrying the
SAME stale "25 cells (6 SCORED / 4 NOT_CENSUSED)" figure for three more days, invisible, and was caught
only by a human reading it. A hand-listed target set beside the data that defines it always drifts; this
is the fourth instance of that bug class in this repo. The doc set is now DERIVED from the tree.

WHY DATED wiki/ ARTIFACTS ARE EXEMPT (load-bearing, not laziness)
A dated memo -- `wiki/reproducibility_freeze_2026-06-13.md`, every `*_2026-*.md` result packet -- RECORDS
WHAT WAS TRUE THEN. Forcing it to today's numbers would make it lie about its own moment and destroy the
audit trail. The exemption is scoped by SHAPE (a dated artifact), never by an allowlist of files someone
has to remember to extend, and `test_the_exemption_does_not_swallow_the_live_docs` pins that the live set
stays non-empty so the exemption can never quietly become "check nothing".

SCOPE: this asserts the docs agree with the COMMITTED artifact. It does not re-run the builder -- that
would make a unit test depend on the whole scoring surface. If the artifact itself is stale, rebuild with
`uv run python scripts/build_validation_report_card.py` (read-only, exit 0 always).
"""
from __future__ import annotations

import json
import re
from collections import Counter
from pathlib import Path

import pytest

ARTIFACT = Path("wiki/decoder_validation_report_card.json")
CLAUDE_MD = Path("CLAUDE.md")

# A doc is LIVE if it describes the project as it stands now. Derived from the tree rather than listed,
# so a new top-level or docs/ markdown file is covered the day it lands.
_LIVE_DOC_DIRS = (Path("."), Path("docs"))
# A dated artifact records a moment; see the module docstring. Matched by SHAPE, not by name.
_DATED = re.compile(r"_?\d{4}-\d{2}-\d{2}")

# APPEND-ONLY DATED-ENTRY LOGS are exempt for the SAME reason dated artifacts are: each entry is scoped
# to the moment it records. CHANGELOG.md line ~533 records "Current card: 25 cells (6 SCORED /
# 4 NOT_CENSUSED ...)" INSIDE a past release entry, where it was TRUE; LESSONS_LEARNED.md quotes the same
# superseded figure inside the dated lesson ABOUT that correction. Forcing either to today's numbers
# would falsify the record -- the guard making a file lie so it could go green.
#
# BOTH were found by this guard rather than anticipated, and BOTH were verified by reading the file
# before exempting. LESSONS_LEARNED.md in particular caught me writing the lesson that describes this
# very guard: the blanket rule demanded a HISTORY log enumerate every CURRENT count, which is the
# "blanket guard forces a false claim on the exempt" failure -- scope to the records the rule applies to.
_HISTORY_FILES = {"CHANGELOG.md", "LESSONS_LEARNED.md"}


def live_docs() -> list[Path]:
    """Every live project doc, derived. Dated artifacts and release-history files are excluded."""
    out = []
    for d in _LIVE_DOC_DIRS:
        for p in sorted(d.glob("*.md")):
            if _DATED.search(p.name) or p.name in _HISTORY_FILES:
                continue
            out.append(p)
    return out


_STATE_COUNT = re.compile(r"\b\d+\s+(SCORED|UNDERPOWERED|NOT_CENSUSED|LABEL_CONFOUNDED|"
                          r"ABSTAINS_BY_DESIGN|NO_FREE_PHENOTYPE_SOURCE)\b")


def cited_states(text: str) -> set[str]:
    """The cell-states this text actually states a count for. Pure."""
    return {m.group(1) for m in _STATE_COUNT.finditer(text)}


def docs_citing_counts() -> list[Path]:
    """The live docs that state at least one report-card count -- where a stale number can hide."""
    return [p for p in live_docs() if cited_states(p.read_text(encoding="utf-8", errors="replace"))]


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
    text = CLAUDE_MD.read_text(encoding="utf-8")
    # Ignore the explanatory correction note, which deliberately quotes the old wrong figure.
    text = re.sub(r"\*\*re-derived from the artifact.*?\*\*", "", text, flags=re.S)
    counts = set(re.findall(r"(\d+) SCORED\b", text))
    assert len(counts) <= 1, (
        f"CLAUDE.md states conflicting SCORED-cell counts: {sorted(counts)}. "
        "Exactly one can be right; re-derive from wiki/decoder_validation_report_card.json."
    )


# ------------------------------------------------------------------ every LIVE doc, not just CLAUDE.md

@pytest.mark.skipif(not ARTIFACT.exists(), reason="report card artifact absent")
def test_every_live_doc_that_cites_a_count_agrees_with_the_artifact():
    """THE GAP THIS CLOSES. `docs/ARCHITECTURE.md` carried the same stale "25 cells / 6 SCORED /
    4 NOT_CENSUSED" figure for three days after CLAUDE.md was corrected, because the guard hand-listed
    CLAUDE.md. Any live doc stating a count is now checked against the artifact.
    """
    states = _states()
    wrong = {}
    for p in docs_citing_counts():
        text = p.read_text(encoding="utf-8", errors="replace")
        # a doc may quote a superseded figure while EXPLAINING the correction; that is honest, so a
        # stale number is only a failure when the CORRECT one is absent from the same file.
        # SCOPED to the states the doc actually cites (narrowed 2026-08-31). The blanket rule -- every
        # doc mentioning any count must enumerate ALL of them -- failed NEXT.md's "3 of 10 SCORED AMR
        # cells rest on one BioProject" for not also stating `2 ABSTAINS_BY_DESIGN` and `11
        # NO_FREE_PHENOTYPE_SOURCE`, forcing a scoped source-concentration sentence to assert a census
        # it never made. That is this guard's OWN documented failure mode (see the LESSONS_LEARNED note
        # above). Every count a doc DOES state is still checked, so a stale figure anywhere still fails;
        # CLAUDE.md's completeness is separately enforced by
        # `test_claude_md_states_the_correct_count_for_every_populated_state`.
        cited = cited_states(text)
        missing = [f"{n} {s}" for s, n in sorted(states.items())
                   if n and s in cited and f"{n} {s}" not in text]
        if missing:
            wrong[str(p)] = missing
    assert not wrong, (
        f"live docs disagree with wiki/decoder_validation_report_card.json: {wrong}. "
        "Re-derive with `uv run python scripts/build_validation_report_card.py` and update them."
    )


@pytest.mark.skipif(not ARTIFACT.exists(), reason="report card artifact absent")
def test_no_live_doc_contradicts_itself_on_the_scored_count():
    """Generalises the self-contradiction check past CLAUDE.md. Two different `<n> SCORED` claims in one
    live file means at least one is wrong, whichever the artifact happens to agree with."""
    bad = {}
    for p in docs_citing_counts():
        text = re.sub(r"\*\*re-derived from the artifact.*?\*\*", "",
                      p.read_text(encoding="utf-8", errors="replace"), flags=re.S)
        text = re.sub(r"supersed\w+ a stale.*?\*\*", "", text, flags=re.S)
        counts = set(re.findall(r"(\d+) SCORED\b", text))
        if len(counts) > 1:
            bad[str(p)] = sorted(counts)
    assert not bad, f"live docs state conflicting SCORED-cell counts: {bad}"


def test_the_exemption_does_not_swallow_the_live_docs():
    """NON-VACUITY. The dated-artifact exemption is matched by SHAPE, so a naming change could in
    principle exempt everything and leave the guard checking nothing. Pin that it does not."""
    live = live_docs()
    assert CLAUDE_MD in live, "CLAUDE.md fell out of the live set — the exemption is too broad"
    assert Path("docs/ARCHITECTURE.md") in live, "ARCHITECTURE.md is the file this guard was widened for"
    assert len(live) >= 8, f"only {len(live)} live docs discovered: {[str(p) for p in live]}"


def test_dated_artifacts_are_exempt_because_they_record_a_moment():
    """A dated memo must NOT be forced to today's numbers — that would make it lie about its own moment
    and destroy the audit trail. Pins the exemption's REASON, not just its effect."""
    assert _DATED.search("reproducibility_freeze_2026-06-13.md")
    assert _DATED.search("decoder_report_card_2026-01-02.md")
    assert not _DATED.search("ARCHITECTURE.md")
    assert not _DATED.search("CLAUDE.md")
    # and no dated file may sneak into the live set
    assert not [p for p in live_docs() if _DATED.search(p.name)]


@pytest.mark.parametrize("name", sorted(_HISTORY_FILES))
def test_each_history_log_is_exempt_for_a_verified_reason_not_a_convenience(name):
    """Both exemptions were found BY this guard, not anticipated, and both were verified by reading the
    file. Each is an append-only dated-entry log: a past entry recording the then-current card is TRUE in
    place, and rewriting it to today's numbers would falsify the record.

    This pins the FACT that justifies each exemption. If a log ever stops carrying superseded figures
    inside past entries, the carve-out should be revisited rather than assumed to still apply.
    """
    p = Path(name)
    if not p.exists():
        pytest.skip(f"{name} absent")
    assert p not in live_docs()
    text = p.read_text(encoding="utf-8", errors="replace")
    assert "6 SCORED" in text, (
        f"{name} no longer records a superseded card figure — re-check whether it still needs the "
        f"history exemption, rather than leaving a stale carve-out in place")


def test_the_history_exemption_is_narrow_enough_to_still_catch_the_real_defect():
    """NON-VACUITY for the exemption itself. Two files are now exempt; the guard is worthless if the
    headline-carrying docs drifted into that set."""
    live = {p.name for p in live_docs()}
    assert "CLAUDE.md" in live and "ARCHITECTURE.md" in live
    assert not (_HISTORY_FILES & live), "a history file leaked into the live set"
    assert len(_HISTORY_FILES) <= 3, (
        "the history exemption is growing — each addition must be a verified append-only dated-entry "
        "log, not a file that was easier to exempt than to fix")


# --- the count check was SCOPED to cited states (2026-08-31); pin it still catches the defect ---

def test_a_stale_count_still_fails_even_when_only_one_state_is_cited():
    """The narrowing must not create a hiding place. A doc citing ONE state is still checked on it."""
    states = _states()
    n_scored = states.get("SCORED", 0)
    assert n_scored, "precondition: the artifact has SCORED cells"
    stale = f"the card holds {n_scored + 1} SCORED cells"
    cited = cited_states(stale)
    assert cited == {"SCORED"}
    assert f"{n_scored} SCORED" not in stale        # exactly what the assertion checks -> would fail


def test_a_doc_is_not_forced_to_state_counts_it_never_claimed():
    """The false positive that motivated the narrowing: NEXT.md's source-concentration sentence.

    Demanding it also state `2 ABSTAINS_BY_DESIGN` and `11 NO_FREE_PHENOTYPE_SOURCE` would force the
    file to assert a census it never made -- the blanket-guard-forces-a-false-claim failure.
    """
    cited = cited_states("3 of 10 SCORED AMR cells rest on one BioProject; 110 cells are registered")
    assert cited == {"SCORED"}, "a passing SCORED mention must not drag in the other states"
    assert cited_states("no counts here at all") == set()


def test_the_defect_this_guard_was_built_for_is_still_in_scope():
    """`docs/ARCHITECTURE.md` carried "25 cells / 6 SCORED / 4 NOT_CENSUSED" stale for three days."""
    assert cited_states("25 cells -- 6 SCORED / 4 NOT_CENSUSED") == {"SCORED", "NOT_CENSUSED"}


def test_the_scoped_check_is_not_vacuous():
    """A predicate that selected nothing would make every check above pass while guarding nothing."""
    scoped = docs_citing_counts()
    assert scoped, "no live doc is checked -- the guard is vacuous"
    assert CLAUDE_MD in scoped, "CLAUDE.md carries the full census and must stay in scope"
