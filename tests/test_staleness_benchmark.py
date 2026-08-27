"""The staleness benchmark + auditor: the properties that make its verdict trustworthy.

WHY THIS EXISTS. ~91% of the substantive errors in the 2026-08-25/27 sessions were stale/unverified
semantic claims, and `test_claude_md_citations.py` records in its own docstring that it CANNOT catch that
class. A mechanical proximity screen was tried, scored 0/5 TP and 5/5 FP, and was not shipped. Before
spending GPU on a semantic auditor, the benchmark that would judge it has to be sound -- otherwise a PASS
means nothing.

These tests pin the four things that could make the verdict a lie: a gameable score, a parser that
fabricates hits, a baseline that is asserted rather than computed, and a negative set built from memory.
"""
from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from kaggle_staleness_auditor import (  # noqa: E402
    ARTIFACT_HEAD_CHARS, build_prompt, load_items_with_artifacts, parse_verdict, render_kernel,
    score_results,
)
from staleness_benchmark import (  # noqa: E402
    PASS_MAX_FP, PASS_MIN_TP, POSITIVES, baseline_verdicts, build, mechanical_screen,
    regenerate_negatives, score,
)


# ------------------------------------------------------------------ the benchmark's own soundness

def test_the_benchmark_has_both_classes():
    """A one-class benchmark cannot distinguish a good auditor from one that flags everything."""
    items = build()
    pos = [i for i in items if i.label == "stale"]
    neg = [i for i in items if i.label == "supported"]
    assert len(pos) == 5, f"expected 5 curated positives, got {len(pos)}"
    assert len(neg) >= 3, (
        f"only {len(neg)} negatives regenerated — the mechanical screen found too few hits on the live "
        f"CLAUDE.md, so the false-positive half of the benchmark is too weak to bind")


def test_the_negatives_are_REGENERATED_not_hand_typed():
    """The 5 false positives were described in prose but never persisted. Hand-copying them from a comment
    would make the negative set my own recollection — the exact thing this benchmark exists to distrust."""
    neg = regenerate_negatives()
    assert neg, "no negatives regenerated"
    assert all(n.source == "regenerated" for n in neg)
    # and they must come from the LIVE file, not a fixture
    assert all((ROOT / n.artifact).exists() for n in neg)


def test_every_positive_names_the_artifact_that_refutes_it():
    """A positive without a resolvable artifact cannot be judged by any auditor — it would be unfair."""
    for p in POSITIVES:
        assert (ROOT / p.artifact).exists(), f"{p.item_id} cites a missing artifact: {p.artifact}"
        assert p.why, f"{p.item_id} records no ground-truth reason"


def test_every_item_ships_a_nonempty_artifact_excerpt():
    """The model is asked to judge a claim against an excerpt. An empty excerpt makes the item unanswerable
    and would silently depress the score."""
    empty = [i["item_id"] for i in load_items_with_artifacts() if not i["artifact_text"].strip()]
    assert not empty, f"items with no artifact text: {empty}"


# ------------------------------------------------------------------ the score cannot be gamed

def test_flagging_everything_FAILS():
    """THE LOAD-BEARING PROPERTY. An auditor that shouts 'stale' at everything gets 5/5 true positives.
    If the score rewarded that, the benchmark would certify a useless model."""
    items = build()
    s = score(items, {i.item_id: "stale" for i in items})
    assert s.tp == 5
    assert not s.passed, "an all-stale auditor must FAIL — it has 5 false positives"


def test_flagging_nothing_FAILS():
    """The mirror: an auditor that never flags has zero false positives but catches nothing."""
    items = build()
    s = score(items, {i.item_id: "supported" for i in items})
    assert s.fp == 0
    assert not s.passed, "a never-flag auditor must FAIL — it catches no real staleness"


def test_hedging_everything_FAILS():
    """`unclear` counts as NOT flagged, so a model cannot pass by refusing to commit."""
    items = build()
    s = score(items, {i.item_id: "unclear" for i in items})
    assert not s.passed


def test_a_perfect_run_passes_and_beats_the_baseline():
    items = build()
    s = score(items, {i.item_id: i.label for i in items})
    assert s.tp == 5 and s.fp == 0 and s.passed


def test_an_unanswered_item_is_recorded_not_silently_dropped():
    """A kernel that crashes halfway must not look like a cautious model."""
    items = build()
    s = score(items, {})
    assert len(s.missing) == len(items)
    assert not s.passed


# ------------------------------------------------------------------ the baseline is computed, not asserted

def test_the_mechanical_baseline_reproduces_its_documented_failure():
    """The bar to beat is 0/5 TP and 5/5 FP. That figure is quoted in several places from prose; this
    RE-DERIVES it from the live screen, so the comparison rests on a computation."""
    items = build()
    base = score(items, baseline_verdicts(items))
    assert base.tp == 0, f"the mechanical screen should catch NO real staleness, got tp={base.tp}"
    assert base.fp == len([i for i in items if i.label == "supported"])
    assert not base.passed


def test_the_mechanical_screen_still_finds_its_false_positives():
    """If CLAUDE.md is ever rewritten such that the screen finds nothing, the negative set silently
    empties and the benchmark stops binding. Fail loudly instead."""
    text = (ROOT / "CLAUDE.md").read_text(encoding="utf-8", errors="replace")
    hits = mechanical_screen(text)
    assert len(hits) >= 3, f"only {len(hits)} mechanical hits — the negative half is eroding"


# ------------------------------------------------------------------ the parser cannot fabricate a hit

@pytest.mark.parametrize("raw,expect,ok", [
    ('{"verdict": "stale", "evidence": "ran on CPU"}', "stale", True),
    ('<think>reasoning</think>{"verdict":"supported","evidence":"consistent"}', "supported", True),
    ("I think it is stale because the file shipped", "unclear", False),   # prose, no JSON
    ('{"verdict": "stal', "unclear", False),                              # truncated
    ('{"verdict": "maybe", "evidence": "x"}', "unclear", False),          # out-of-vocabulary
    ("", "unclear", False),
])
def test_parse_verdict_fails_closed(raw, expect, ok):
    """FAIL-CLOSED is the point. A parser that guessed `stale` on malformed output would manufacture true
    positives, and because `unclear` scores as not-flagged, the bug would look like a cautious model."""
    v = parse_verdict(raw)
    assert v["verdict"] == expect
    assert v["parse_ok"] is ok


def test_the_word_stale_in_prose_does_not_become_a_flag():
    """A reasoning model narrating 'this is not stale' must not be read as flagging stale."""
    assert parse_verdict("The claim is not stale at all.")["verdict"] == "unclear"


# ------------------------------------------------------------------ prompt + kernel shape

def test_the_prompt_carries_the_claim_and_the_artifact():
    p = build_prompt("some claim", "wiki/x.md", "artifact body here")
    assert "some claim" in p and "wiki/x.md" in p and "artifact body here" in p


def test_the_prompt_truncates_a_huge_artifact():
    p = build_prompt("c", "wiki/x.md", "A" * (ARTIFACT_HEAD_CHARS * 3))
    assert len(p) < ARTIFACT_HEAD_CHARS * 2, "excerpt not truncated — would blow the context window"


def test_the_kernel_pins_the_kaggle_gotchas_this_repo_already_paid_for():
    """PYTHONUTF8 (or logs return 0-byte) and a T4-appropriate dtype. Recorded gotchas, applied."""
    k = render_kernel()
    assert "PYTHONUTF8" in k
    assert "bfloat16" in k or "float16" in k
    assert "/kaggle/working/results.json" in k


def test_scoring_a_finished_run_reports_both_the_model_and_the_baseline():
    """A model number without the baseline beside it is not a verdict — it is a number."""
    items = build()
    fake = [{"item_id": i.item_id, "raw": json.dumps({"verdict": i.label, "evidence": "f"})}
            for i in items]
    p = Path(tempfile.mkdtemp()) / "r.json"
    p.write_text(json.dumps(fake), encoding="utf-8")
    r = score_results(p)
    assert r["model_score"]["passed"] is True
    assert r["mechanical_baseline"]["passed"] is False
    assert r["beat_baseline"] is True
    assert r["n_unparseable"] == 0


def test_the_pass_condition_is_pinned_where_a_reader_can_see_it():
    """Pre-registration only counts if the threshold cannot drift after results land."""
    assert (PASS_MIN_TP, PASS_MAX_FP) == (3, 1)


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
