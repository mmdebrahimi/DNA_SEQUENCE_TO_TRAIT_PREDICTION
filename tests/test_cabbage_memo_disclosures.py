"""Every limit a CABBAGE cell's ARTIFACT declares must also appear in the memo a human reads.

WHY THIS EXISTS
This repo's own recorded rule is *"a block carried only in JSON is not a disclosure"* — learned when two
of five `trust_surface` disclosure layers reached the record and never a renderer, so the gentamicin
organism-scope warning and the whole prospective tier were invisible to anyone running the CLI. The same
split exists one level out: each cell's `honest_limits` is machine-readable, and the `wiki/*.md` memo is
what anyone actually reads. Nothing checked that the second carried the first.

Four cells were scored on 2026-10-10/11 and each memo carries its caveats BY HAND. A fifth cell written
in a hurry — or a memo edited later — can silently drop "not clonality-corrected" or quote a
single-country cohort as source-diverse, and no existing guard would notice: `test_report_card_doc_sync`
checks that live docs match artifact NUMBERS and deliberately EXEMPTS dated `wiki/*.md`.

THAT EXEMPTION MUST NOT EXTEND HERE, and the distinction is the point. A dated memo is exempt from
number-sync because it records what was true then. A DISCLOSURE is not a number that drifts — "these
counts are isolate-level, not clonality-corrected" describes how the cell was measured and is as true
today as the day it was written. So dated memos are exempt from one check and bound by this one.

EXEMPTIONS OVER wiki/*.md, ENUMERATED — because an exemption phrased as a property of the FILE is
inherited silently by every later obligation, and two obligations can relate OPPOSITELY to the same
property. Add a third check over these files and state explicitly which row it joins.

  | check                                   | dated memo exempt? | why |
  |-----------------------------------------|--------------------|-----|
  | test_report_card_doc_sync (numbers)     | YES                | it records what was true THEN; forcing today's figures would make it lie about its own moment |
  | this file (disclosures)                 | **NO**             | a disclosure describes how the cell was MEASURED and does not drift |

WHAT IS DERIVED AND WHAT IS ASSERTED (stated so a reader can judge it)
DERIVED: the artifact→memo mapping (each memo cites its own `cabbage_scored_cell_*.json`), and WHICH
topics are required for a given cell (read off that cell's own artifact fields). A hand-listed mapping
beside the data that defines it has drifted five times in this repo.
ASSERTED: the keyword alternatives that count as disclosing a topic. Wording varies between memos, so
exact-sentence matching would fail; these are deliberately generous, which makes this a FLOOR — it
catches an ABSENT disclosure, not a weakly-worded one.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

WIKI = Path("wiki")
CELL_GLOB = "cabbage_scored_cell_*.json"

# ASSERTED: generous keyword alternatives per topic. A floor, not a quality bar.
TOPIC_PATTERNS = {
    "clonality": r"clonalit|clonally|clonality-corrected|lineage",
    "source_concentrated": r"concentrat",
    "source_diverse": r"source[- ]diverse|source_diverse",
    "fetch_failures": r"fetch fail|fetch_fail|fetch-fail",
    "not_prospective": r"prospective",
}


def _cells() -> list[Path]:
    return sorted(WIKI.glob(CELL_GLOB))


def _memos() -> list[Path]:
    return sorted(p for p in WIKI.glob("cabbage_*.md"))


def _memo_for(artifact: Path) -> list[Path]:
    """DERIVED: the memo(s) citing this artifact by filename."""
    return [m for m in _memos() if artifact.name in m.read_text(encoding="utf-8")]


def _required_topics(d: dict) -> set[str]:
    """DERIVED from the artifact's own fields: which disclosures this cell owes a reader."""
    need = set()
    limits = " ".join(d.get("honest_limits") or [])
    if re.search(r"clonalit", limits, re.I):
        need.add("clonality")
    sc = d.get("source_concentration_of_scored_set") or {}
    if sc.get("verdict") == "source_concentrated":
        need.add("source_concentrated")
    elif sc.get("verdict") == "source_diverse":
        need.add("source_diverse")
    if (d.get("failures") or {}).get("n_fetch_failed", 0) > 0:
        need.add("fetch_failures")
    if re.search(r"prospective", limits, re.I):
        need.add("not_prospective")
    return need


def test_there_are_cells_and_memos_to_check():
    """Guards the guard: an empty glob would make every test below pass vacuously."""
    assert _cells(), "no cabbage_scored_cell_*.json artifacts found"
    assert _memos(), "no cabbage_*.md memos found"


@pytest.mark.parametrize("artifact", _cells(), ids=lambda p: p.stem)
def test_every_scored_cell_artifact_is_described_by_a_memo(artifact):
    """A cell scored with no human-readable memo discloses nothing to anyone not reading JSON."""
    assert _memo_for(artifact), (
        f"{artifact.name} is cited by no wiki/cabbage_*.md memo. A cell whose only record is JSON has "
        "no disclosure surface for a human reader.")


@pytest.mark.parametrize("artifact", _cells(), ids=lambda p: p.stem)
def test_the_memo_discloses_every_limit_its_artifact_declares(artifact):
    d = json.loads(artifact.read_text(encoding="utf-8"))
    need = _required_topics(d)
    assert need, f"{artifact.name} declares no limits at all, which is itself suspect"
    memos = _memo_for(artifact)
    text = " ".join(m.read_text(encoding="utf-8") for m in memos).lower()
    missing = [t for t in sorted(need)
               if not re.search(TOPIC_PATTERNS[t], text, re.I)]
    assert not missing, (
        f"{artifact.name}: its own artifact declares these and the memo(s) "
        f"{[m.name for m in memos]} do not disclose them: {missing}")


@pytest.mark.parametrize("artifact", _cells(), ids=lambda p: p.stem)
def test_a_concentrated_cohort_is_never_described_as_source_diverse(artifact):
    """The dangerous direction. Only ONE of the four cells clears this repo's diversity bar, and a
    memo calling a 96%-one-country cohort source-diverse would be a trust-surface falsehood."""
    d = json.loads(artifact.read_text(encoding="utf-8"))
    sc = d.get("source_concentration_of_scored_set") or {}
    if sc.get("verdict") != "source_concentrated":
        pytest.skip("cell is not source_concentrated")
    for m in _memo_for(artifact):
        body = m.read_text(encoding="utf-8")
        for hit in re.finditer(TOPIC_PATTERNS["source_diverse"], body, re.I):
            window = body[max(0, hit.start() - 220): hit.end() + 220].lower()
            # A memo may legitimately mention the term while NEGATING it, or while contrasting with
            # the one cell that does clear the bar.
            negated = re.search(r"not |fails|never|do not|does not|only .*e\.? ?coli|unlike", window)
            assert negated, (
                f"{m.name} uses a source-diverse phrase near offset {hit.start()} with no negation "
                f"or contrast in context, but {artifact.name} is source_concentrated "
                f"(largest country share {sc.get('largest_country_share')})")


def test_the_required_topic_set_is_derived_not_hardcoded():
    """If `_required_topics` ignored the artifact it would demand the same topics of every cell, and a
    cell with a genuinely different limit set would pass by accident."""
    diverse = {"source_concentration_of_scored_set": {"verdict": "source_diverse"},
               "honest_limits": ["clonality note"], "failures": {"n_fetch_failed": 0}}
    conc = {"source_concentration_of_scored_set": {"verdict": "source_concentrated"},
            "honest_limits": ["clonality note", "not prospective"],
            "failures": {"n_fetch_failed": 5}}
    a, b = _required_topics(diverse), _required_topics(conc)
    assert a != b, "the topic set does not vary with the artifact"
    assert "source_diverse" in a and "source_concentrated" in b
    assert "fetch_failures" in b and "fetch_failures" not in a
    assert _required_topics({}) == set(), "an empty artifact should demand nothing"


def test_the_check_FIRES_on_a_memo_that_drops_a_disclosure(tmp_path):
    """Non-vacuity, proven rather than promised: a guard never shown to fail is not known to work.
    Rebuilds the real comparison against a memo with the clonality disclosure removed."""
    art = {"source_concentration_of_scored_set": {"verdict": "source_concentrated",
                                                  "largest_country_share": 0.96},
           "honest_limits": ["NOT CLONALITY-CORRECTED: isolate-level counts."],
           "failures": {"n_fetch_failed": 12}}
    need = _required_topics(art)
    assert "clonality" in need and "fetch_failures" in need

    good = "This cell is source_concentrated. Not clonality-corrected. 12 fetch failures excluded."
    bad = "This cell is source_concentrated. 12 fetch failures excluded."

    def missing(text):
        return [t for t in sorted(need) if not re.search(TOPIC_PATTERNS[t], text, re.I)]

    assert missing(good) == [], f"the complete memo should pass, got {missing(good)}"
    assert missing(bad) == ["clonality"], (
        "dropping the clonality sentence must be detected; otherwise this guard is decoration")
