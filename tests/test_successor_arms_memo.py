"""Tests for the three-arm memo.

EVERY RAIL GETS A FIRING TEST. The memo's three guards each re-read the rendered text, so each one
could silently be a no-op; a test that only renders the real artifacts would report all three as
exercised while never having run any. So for each rail there is a hand-built text that violates it
and must raise, beside the real render that must not.

THE ORDERING IS ASSERTED, NOT EYEBALLED. Decisions-before-corrections-before-findings is the whole
design -- a reader who stops after one screen must have seen the decisions -- so the section offsets
are compared rather than the headings merely being present.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

import scripts.successor_arms_memo as M
from scripts.successor_arms_memo import (
    MemoRefused,
    assert_gains_carry_context,
    assert_limits_present,
    assert_no_ineligible_headline,
    build,
    headings,
    honest_limits_of,
    ineligible_numbers,
    load_arms,
    render_memo,
)

REPO = Path(__file__).resolve().parents[1]
ASGIVEN = "asgiven-position-resolved-promoteratlas-150bp"
HEADS = "shipped-expression-heads-zeroshot-150bp"
FINETUNE = "finetuned-promoteratlas-encoder-150bp"

HEADS_ART = REPO / "wiki/expression_heads_result_2026-11-02.json"
ASGIVEN_ART = REPO / "wiki/asgiven_arm_result_2026-11-02.json"

pytestmark = pytest.mark.skipif(
    not (HEADS_ART.is_file() and ASGIVEN_ART.is_file()),
    reason="the committed arm artifacts are absent here")


def _report(collected=(ASGIVEN, HEADS), refused=(FINETUNE,)) -> dict:
    arms = {}
    for n in collected:
        arms[n] = {"arm": n, "collected": True, "metric": 0.5, "verdict": "X",
                   "artifact": "wiki/%s.json" % n, "interval": None}
    for n in refused:
        arms[n] = {"arm": n, "collected": False, "artifact": "wiki/%s.json" % n,
                   "refusal": "job is NOT_STARTED, not COMPLETE"}
    return {"arms": arms, "n_collected": len(collected), "n_arms": len(arms)}


VERDICT = {"verdict": "blocked:untried-arm", "detail": "1 arm(s) declared but never run",
           "measured": [], "planned": [FINETUNE], "killed": [], "best": ASGIVEN}


def _real_arts() -> dict:
    return load_arms()


def _real_text() -> str:
    return build(_real_arts(), _report(), VERDICT)


# ----------------------------------------------------------------- limits completeness


def test_every_declared_limit_appears_verbatim_in_the_real_memo():
    arts = _real_arts()
    text = _real_text()
    total = 0
    for art in arts.values():
        for lim in honest_limits_of(art):
            total += 1
            assert lim in text
    assert total >= 12, "both committed arms declare six limits each; found %d" % total


def test_the_limits_rail_FIRES_when_one_is_dropped():
    """Non-vacuity proof: without this the completeness check could be a no-op."""
    arts = _real_arts()
    text = _real_text()
    tampered = {ASGIVEN: {"honest_limits": ["A LIMIT NOBODY RENDERED."]}}
    with pytest.raises(MemoRefused, match="absent from the memo"):
        assert_limits_present(text, tampered)


def test_honest_limits_handles_BOTH_artifact_shapes():
    """The as-given and heads artifacts carry a list; the fine-tune scorer emits a dict."""
    assert honest_limits_of({"honest_limits": ["a", "b"]}) == ["a", "b"]
    out = honest_limits_of({"honest_limits": {"what_a_null_bounds": "only shuffling"}})
    assert len(out) == 1 and "only shuffling" in out[0]
    assert honest_limits_of({}) == []
    assert honest_limits_of(None) == []


def test_a_dict_shaped_limit_with_a_nested_value_is_still_rendered():
    out = honest_limits_of({"honest_limits": {"multiplicity": {"correction": "none"}}})
    assert len(out) == 1 and "correction" in out[0]


# ----------------------------------------------------------------- the ineligible rail


def test_the_heads_ineligible_secondary_is_detected():
    art = json.loads(HEADS_ART.read_text(encoding="utf-8"))
    vals = ineligible_numbers(art)
    assert "-0.2795" in vals, "the ineligible secondary head was not detected: %r" % vals


def test_an_ELIGIBLE_block_is_not_treated_as_ineligible():
    art = {"secondaries": {"h": {"eligible_as_the_arms_result": True,
                                 "median_spearman": 0.9}}}
    assert ineligible_numbers(art) == []


def test_no_ineligible_value_reaches_a_heading_in_the_real_memo():
    arts = _real_arts()
    text = _real_text()
    assert_no_ineligible_headline(text, arts)
    heads_line = " ".join(headings(text))
    assert "-0.2795" not in heads_line


def test_the_ineligible_rail_FIRES_on_a_headlined_secondary():
    """Non-vacuity proof, driven with the real ineligible value."""
    arts = {HEADS: json.loads(HEADS_ART.read_text(encoding="utf-8"))}
    bad = "# Memo\n\n## The arm scored -0.2795\n\nbody\n"
    with pytest.raises(MemoRefused, match="INELIGIBLE"):
        assert_no_ineligible_headline(bad, arts)


def test_an_ineligible_value_in_the_BODY_is_allowed():
    """Disclosing that ineligible secondaries exist is required; headlining them is not."""
    arts = {HEADS: json.loads(HEADS_ART.read_text(encoding="utf-8"))}
    assert_no_ineligible_headline("# Memo\n\nbody mentions -0.2795 as ineligible\n", arts)


def test_the_real_memo_DISCLOSES_that_ineligible_secondaries_exist():
    assert "INELIGIBLE" in _real_text()


# ----------------------------------------------------------------- gains carry context


def test_the_gain_context_rail_FIRES_without_a_margin():
    with pytest.raises(MemoRefused, match="margin"):
        assert_gains_carry_context("the arm shows a gain of +0.0010 and nothing else")


def test_the_gain_context_rail_FIRES_without_the_ceiling():
    with pytest.raises(MemoRefused, match="0.941"):
        assert_gains_carry_context("a gain of +0.0010 against a margin of 0.02, incumbent 0.5076")


def test_the_gain_context_rail_is_silent_when_no_gain_is_stated():
    assert_gains_carry_context("this memo makes no such claim at all")


def test_the_gain_rail_does_not_trip_on_the_word_AGAINST():
    """Regression: the first version substring-matched "gain", whose letters sit inside "against",
    so a memo with no gain claim tripped the rail. Whole-word matching is the fix."""
    assert_gains_carry_context("the arm was scored against the incumbent and nothing more")


def test_the_real_memo_states_margin_incumbent_and_ceiling():
    text = _real_text()
    assert "margin" in text and "incumbent" in text and "0.941" in text


def test_no_metric_is_rendered_as_a_BARE_number():
    """Every metric line must carry the incumbent and the ceiling fraction beside it."""
    for line in _real_text().splitlines():
        if line.strip().startswith("- metric:"):
            assert "incumbent" in line and "replicate ceiling" in line, line


# ----------------------------------------------------------------- ordering by reader action


def test_decisions_come_before_corrections_which_come_before_findings():
    text = _real_text()
    i_dec = text.index("## Decisions needed")
    i_cor = text.index("## Corrections")
    i_fin = text.index("## What each arm measured")
    assert i_dec < i_cor < i_fin


def test_the_decisions_section_is_within_the_first_screen():
    """Truncation, not attention: content below a cut is absent, not merely less-read."""
    text = _real_text()
    assert text.index("## Decisions needed") < len("\n".join(text.splitlines()[:8])) + 200


# ----------------------------------------------------------------- the verdict is reported only


def test_the_memo_states_the_verdict_and_that_accepting_it_is_the_users_call():
    text = _real_text()
    assert "blocked:untried-arm" in text
    assert "Nothing here acts on the verdict" in text


def test_the_memo_distinguishes_arm_set_exhausted_from_ceiling_reached():
    """The central over-claim this arm structure exists to prevent."""
    text = _real_text()
    assert "ARM SET is exhausted" in text
    assert "would **not** say" in text or "not** say the" in text


def test_an_outstanding_arm_is_named_as_a_decision_not_buried():
    text = _real_text()
    dec = text[text.index("## Decisions needed"):text.index("## Corrections")]
    assert FINETUNE in dec


def test_an_unrun_arm_may_not_be_killed_for_redundancy_and_the_memo_says_so():
    assert "dies only by its own executed evidence" in _real_text()


# ----------------------------------------------------------------- unreadable arms


def test_an_unreadable_arm_is_REPORTED_not_omitted():
    """Silently dropping an arm would make the memo look complete when it is not."""
    arts = dict(_real_arts())
    arts[FINETUNE] = None
    text = build(arts, _report(), VERDICT)
    assert FINETUNE in text
    assert "Not readable" in text


def test_load_arms_returns_None_for_an_absent_artifact(tmp_path):
    arts = load_arms(tmp_path)
    assert set(arts) >= {ASGIVEN, HEADS, FINETUNE}
    assert all(v is None for v in arts.values())


def test_a_memo_over_no_readable_arms_still_renders_the_decisions():
    """The worst case must still tell the reader what is outstanding."""
    arts = {ASGIVEN: None, HEADS: None, FINETUNE: None}
    text = build(arts, _report(collected=(), refused=(ASGIVEN, HEADS, FINETUNE)), VERDICT)
    assert "## Decisions needed" in text
    assert text.count("Not readable") == 3


# ----------------------------------------------------------------- end to end


def test_the_real_memo_renders_end_to_end_and_main_writes_it(tmp_path):
    out = tmp_path / "memo.md"
    assert M.main(["--out", str(out)]) == 0
    text = out.read_text(encoding="utf-8")
    assert text.startswith("# The three successor arms")
    assert "ASGIVEN_TIES" in text and "HEADS_NO_GAIN" in text


def test_render_is_pure_and_does_not_touch_disk(tmp_path, monkeypatch):
    """`render_memo` must be callable on synthetic data, or the rails cannot be tested."""
    monkeypatch.chdir(tmp_path)
    text = render_memo({ASGIVEN: None}, _report(collected=(), refused=(ASGIVEN,)), VERDICT)
    assert "## Decisions needed" in text
