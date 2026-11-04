"""Tests for the K562 eukaryotic-arm memo renderer.

THE CENTRAL RAIL IS VACUOUS ON THE REAL ARTIFACT, AND THAT IS THE REASON THESE TESTS EXIST IN THIS
SHAPE. The committed artifact's `provisional_verdict_if_complete` is `None`, so "the provisional
verdict never reaches a heading" is trivially true of it -- nothing could leak because there is
nothing to leak. A suite that only rendered the real artifact would report a guard as exercised
while never having run it. So the rail is driven three ways: against a SYNTHETIC artifact that
carries a provisional value, against a hand-built text that puts that value in a heading (proving
the guard FIRES), and against the real artifact (proving it does not false-positive on `None`).
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

import scripts.eukaryotic_arm_memo as M
from scripts.eukaryotic_arm_memo import (
    MemoRefused,
    assert_limits_complete,
    assert_provisional_not_headlined,
    build,
    headings,
    render_memo,
)

REPO = Path(__file__).resolve().parents[1]
REAL = REPO / "wiki/eukaryotic_arm_result_K562_2026-11-01.json"


def _real() -> dict:
    return json.loads(REAL.read_text(encoding="utf-8"))


def synthetic(**over) -> dict:
    """A minimal valid artifact. Defaults carry a PROVISIONAL verdict, which the real one lacks."""
    art = {
        "schema": "eukaryotic-arm-result-v1",
        "cell": "SYNTH",
        "date": "2026-01-01",
        "bar": "wiki/some_bar.json",
        "criterion": "ceiling: a criterion",
        "verdict": "SUBSTRATE_LIMITED",
        "verdict_detail": "the CI lies above the anchor",
        "provisional_verdict_if_complete": "EUKARYOTIC_ADDS",
        "honest_limits": ["LIMIT ONE is load bearing.", "LIMIT TWO is also load bearing."],
        "refused": {"a_comparison": {"why": "the metric is unstated",
                                     "what_is_done_instead": "recorded verbatim"}},
        "gain_bootstrap": {"ci_lo": 0.1, "ci_hi": 0.2, "mean": 0.15, "n_folds": 10},
        "medians": {"conv": 0.5, "gc": 0.2},
        "bacterial_anchor": {"conv_minus_gc": 0.0122, "why_this_is_the_anchor": "within-substrate",
                             "best_arm": 0.5076, "conv": 0.3234, "gc": 0.3112,
                             "replicate_ceiling": 0.941,
                             "fraction_of_ceiling_captured_by_best": 0.5394,
                             "headroom_best_minus_gc": 0.1964,
                             "corroborating_artifact": "wiki/x.md"},
        "substrate": {"n_measured": 10, "n_designed": 20, "join_rate": 0.5,
                      "sequence_length": 230, "sequence_column": "seq", "n_folds_published": 10},
        "target_reparameterization": {"applied": "y = x + 1", "computed_over": "the whole cell",
                                      "why": "log1p is undefined below -1",
                                      "evaluation_unaffected": "Spearman is invariant"},
        "per_fold": [{"fold": 1, "conv_spearman": 0.5, "gc_spearman": 0.2, "gain": 0.3,
                      "n_train": 9, "n_test": 1}],
        "wall_seconds": 1.0,
    }
    art.update(over)
    return art


# ----------------------------------------------------------------- the enforced provisional rail


def test_the_guard_FIRES_when_the_provisional_verdict_reaches_a_heading():
    """Non-vacuity proof. Without this, the rail could be a no-op and every other test would pass."""
    art = synthetic()
    text = "# A title\n\n## Verdict: EUKARYOTIC_ADDS\n\nbody\n"
    with pytest.raises(MemoRefused, match="appears in the heading"):
        assert_provisional_not_headlined(text, art)


def test_the_renderer_does_not_headline_a_provisional_verdict_that_EXISTS():
    """Driven on a synthetic artifact precisely because the real one's field is None."""
    art = synthetic()
    text = build(art)
    prov = art["provisional_verdict_if_complete"]
    assert prov not in " ".join(headings(text))
    assert any("SUBSTRATE_LIMITED" in h for h in headings(text))


def test_a_provisional_verdict_is_still_DISCLOSED_in_the_body_not_hidden():
    """Not headlining it is not the same as concealing it; the reader must know it exists."""
    text = build(synthetic())
    assert "provisional_verdict_if_complete" in text
    assert "not** this arm's result" in text or "not this arm's result" in text


def test_the_guard_does_not_false_positive_on_a_None_provisional():
    assert_provisional_not_headlined("## Verdict: SUBSTRATE_LIMITED\n",
                                     synthetic(provisional_verdict_if_complete=None))


def test_an_empty_string_provisional_is_treated_as_absent():
    assert_provisional_not_headlined("## Verdict: SUBSTRATE_LIMITED\n",
                                     synthetic(provisional_verdict_if_complete="   "))


def test_the_guard_checks_EVERY_heading_level():
    art = synthetic()
    for hashes in ("#", "##", "###", "####", "#####", "######"):
        with pytest.raises(MemoRefused):
            assert_provisional_not_headlined("%s Verdict EUKARYOTIC_ADDS\n" % hashes, art)


def test_a_provisional_value_in_the_BODY_is_not_a_violation():
    """The rail is about headings. Disclosure in prose is required, not forbidden."""
    assert_provisional_not_headlined("# Title\n\nbody mentions EUKARYOTIC_ADDS here\n",
                                     synthetic())


# ----------------------------------------------------------------- the limits-completeness rail


def test_every_honest_limit_appears_verbatim():
    art = synthetic()
    text = build(art)
    for lim in art["honest_limits"]:
        assert lim in text


def test_the_limits_guard_FIRES_when_one_is_dropped():
    """Non-vacuity proof for the completeness check."""
    art = synthetic()
    text = build(art)
    art_with_extra = synthetic(
        honest_limits=art["honest_limits"] + ["A THIRD LIMIT nobody rendered."])
    with pytest.raises(MemoRefused, match="honest limits are absent"):
        assert_limits_complete(text, art_with_extra)


def test_the_limits_guard_names_how_many_are_missing():
    art = synthetic()
    with pytest.raises(MemoRefused, match=r"1 of 3"):
        assert_limits_complete(build(art),
                               synthetic(honest_limits=art["honest_limits"] + ["MISSING ONE."]))


# ----------------------------------------------------------------- refusals to render


def test_an_unrecognised_schema_is_REFUSED():
    with pytest.raises(MemoRefused, match="schema"):
        render_memo(synthetic(schema="something-else-v9"))


@pytest.mark.parametrize("field", ["verdict", "verdict_detail", "honest_limits"])
def test_a_missing_required_field_is_REFUSED(field):
    """Rendering a memo with no verdict would state a result the artifact does not have."""
    with pytest.raises(MemoRefused, match=field):
        render_memo(synthetic(**{field: None}))


def test_main_REFUSES_an_absent_artifact(tmp_path):
    assert M.main(["--artifact", str(tmp_path / "nope.json"),
                   "--out", str(tmp_path / "o.md")]) == 2


def test_main_REFUSES_rather_than_writing_a_bad_memo(tmp_path):
    bad = tmp_path / "bad.json"
    bad.write_text(json.dumps(synthetic(verdict=None)), encoding="utf-8")
    out = tmp_path / "o.md"
    assert M.main(["--artifact", str(bad), "--out", str(out)]) == 2
    assert not out.exists(), "a refused render must leave no partial memo behind"


# ----------------------------------------------------------------- the refused block, verbatim


def test_the_refused_block_is_carried_VERBATIM_not_summarised():
    art = synthetic()
    text = build(art)
    blk = art["refused"]["a_comparison"]
    assert blk["why"] in text
    assert blk["what_is_done_instead"] in text


# ----------------------------------------------------------------- the real committed artifact


@pytest.mark.skipif(not REAL.is_file(), reason="the committed K562 artifact is absent here")
def test_the_real_artifact_renders_end_to_end():
    art = _real()
    text = build(art)
    assert text.splitlines()[0].startswith("# ")
    assert "SUBSTRATE_LIMITED" in text.splitlines()[0]


@pytest.mark.skipif(not REAL.is_file(), reason="the committed K562 artifact is absent here")
def test_the_real_headline_carries_the_verdict_and_the_cell():
    text = build(_real())
    head = text.splitlines()[0]
    assert "SUBSTRATE_LIMITED" in head and "K562" in head


@pytest.mark.skipif(not REAL.is_file(), reason="the committed K562 artifact is absent here")
def test_all_five_real_honest_limits_are_present():
    art = _real()
    text = build(art)
    assert len(art["honest_limits"]) == 5
    for lim in art["honest_limits"]:
        assert lim in text


@pytest.mark.skipif(not REAL.is_file(), reason="the committed K562 artifact is absent here")
def test_the_real_refusal_about_the_unstated_metric_is_carried_verbatim():
    art = _real()
    text = build(art)
    blk = art["refused"]["published_incumbents_head_to_head"]
    assert blk["why"] in text


@pytest.mark.skipif(not REAL.is_file(), reason="the committed K562 artifact is absent here")
def test_the_real_memo_states_that_no_provisional_verdict_exists():
    """The real artifact's field is None; the memo must say so rather than stay silent, or a
    reader cannot tell 'none recorded' from 'not rendered'."""
    art = _real()
    assert art["provisional_verdict_if_complete"] is None
    assert "no** provisional verdict" in build(art)


@pytest.mark.skipif(not REAL.is_file(), reason="the committed K562 artifact is absent here")
def test_the_real_memo_reports_the_within_substrate_anchor_not_raw_spearman():
    """The comparison that crosses a substrate boundary is the gain; saying so is load-bearing."""
    text = build(_real())
    assert "within-substrate" in text
    assert "0.0122" in text


@pytest.mark.skipif(not REAL.is_file(), reason="the committed K562 artifact is absent here")
def test_main_writes_the_real_memo(tmp_path):
    out = tmp_path / "memo.md"
    assert M.main(["--artifact", str(REAL), "--out", str(out)]) == 0
    assert "SUBSTRATE_LIMITED" in out.read_text(encoding="utf-8")


@pytest.mark.skipif(not REAL.is_file(), reason="the committed K562 artifact is absent here")
def test_every_per_fold_row_reaches_the_table():
    art = _real()
    text = build(art)
    assert len(art["per_fold"]) == 10
    for f in art["per_fold"]:
        assert "| %s |" % f["fold"] in text
