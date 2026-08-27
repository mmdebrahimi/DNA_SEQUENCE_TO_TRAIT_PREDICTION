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


# ------------------------------------------------------- the two v2 fixes (2026-08-27), pinned

def test_structural_facts_state_existence_and_shape():
    """THE P3 FIX. The genome-map miss was structural, not textual: browser.py calls ITSELF 'the deferred
    v1 graphical browser' (describing what it implements) and contains ZERO occurrences of 'SHIPPED', so a
    bigger excerpt would only have supplied more of the same misleading prose. What refutes the claim is
    that the module EXISTS and defines functions."""
    from kaggle_staleness_auditor import structural_facts
    body = "def a():\n    pass\ndef b():\n    pass\nclass C:\n    pass\n"
    f = structural_facts("dna_decode/x.py", body)
    assert "exists: yes" in f
    assert "2 function(s)" in f and "1 class(es)" in f


def test_structural_facts_do_not_leak_a_verdict():
    """Honest context, not a thumb on the scale: the facts must never name a verdict, or the benchmark
    would be measuring the prompt rather than the model."""
    from kaggle_staleness_auditor import structural_facts
    f = structural_facts("dna_decode/genome_map/browser.py",
                         (ROOT / "dna_decode/genome_map/browser.py").read_text(encoding="utf-8"))
    low = f.lower()
    assert "stale" not in low and "supported" not in low and "deferred" not in low


def test_the_p3_item_now_carries_its_structural_refutation():
    """The item the model got wrong must now ship the fact that decides it."""
    from kaggle_staleness_auditor import load_items_with_artifacts
    p3 = {i["item_id"]: i for i in load_items_with_artifacts()}["P3_genome_map_browser"]
    assert "exists: yes" in p3["facts"]
    assert "function(s)" in p3["facts"], "a .py artifact must report its implemented-code count"


def test_every_item_carries_nonempty_facts():
    """A blank facts line would silently un-do the fix for that item."""
    from kaggle_staleness_auditor import load_items_with_artifacts
    blank = [i["item_id"] for i in load_items_with_artifacts() if not i["facts"].strip()]
    assert not blank, f"items shipping empty ARTIFACT FACTS: {blank}"


def test_the_generation_budget_was_raised_off_the_value_that_truncated():
    """THE P4 FIX, and it is evidence-driven: at 1200 tokens exactly ONE item blew the cap (5207 chars,
    `</think>` never closed) while the other nine closed at ~1600-3000. The cap bound on one hard item,
    not on the task."""
    from kaggle_staleness_auditor import MAX_NEW_TOKENS, render_kernel
    assert MAX_NEW_TOKENS >= 2000, "P4 truncated at 1200; the budget must clear it with margin"
    k = render_kernel()
    assert f"max_new_tokens={MAX_NEW_TOKENS}" in k
    assert "max_new_tokens=1200" not in k, "the kernel still carries the cap that truncated P4"


def test_the_kernel_prompt_matches_the_local_prompt_shape():
    """The kernel builds its prompt inline; if it drifts from build_prompt, the locally-tested prompt is
    not the one that runs on the GPU — and the tests would be measuring a different system."""
    from kaggle_staleness_auditor import build_prompt, render_kernel
    k = render_kernel()
    local = build_prompt("c", "x.py", "body")
    for field in ("CLAIM (from project documentation)", "ARTIFACT the claim is about",
                  "ARTIFACT FACTS", "ARTIFACT EXCERPT"):
        assert field in local, f"local prompt lost {field!r}"
        assert field in k, f"kernel prompt lost {field!r}"


# ------------------------------------------------------- the corpus extractor (full-run input)

def test_every_extracted_pair_mentions_its_own_artifact():
    """THE PAIRING BUG THIS CAUGHT. The first extractor paired a bullet with EVERY path it cited, so a
    claim about `viz/browser.py` got paired with `config/datasources.yaml` -- sending the model to judge
    a claim against the wrong file. The local-window rule fixed it; this pins that it stays fixed."""
    from staleness_corpus import extract_pairs
    from pathlib import Path as _P
    text = (ROOT / "CLAUDE.md").read_text(encoding="utf-8", errors="replace")
    bad = [p.pair_id for p in extract_pairs(text) if p.artifact not in p.claim]
    assert not bad, f"pairs whose claim does not mention their artifact: {bad[:5]}"


def test_a_pointer_without_a_status_word_is_not_extracted():
    """"see X for the schema" asserts nothing about X, so it has nothing to be stale about. Auditing it
    would spend GPU to learn nothing."""
    from staleness_corpus import extract_pairs
    pointer = "- The format is described in `docs/schema.md` for reference.\n"
    assert extract_pairs(pointer) == []


def test_a_status_claim_near_a_path_IS_extracted():
    from staleness_corpus import extract_pairs
    claim = "- The browser is deferred; see `dna_decode/genome_map/browser.py`.\n"
    pairs = extract_pairs(claim)
    assert len(pairs) == 1
    assert pairs[0].artifact == "dna_decode/genome_map/browser.py"


def test_one_bullet_citing_two_files_yields_two_separately_judged_pairs():
    """A bullet can be accurate about one artifact and stale about another -- the ProSST case was exactly
    that shape, so they must be judged separately rather than collapsed."""
    from staleness_corpus import extract_pairs
    b = ("- Work shipped in `dna_decode/a.py` but the run is still deferred per "
         "`wiki/b.md` pending review.\n")
    pairs = extract_pairs(b)
    assert {p.artifact for p in pairs} == {"dna_decode/a.py", "wiki/b.md"}


def test_the_corpus_kernel_uses_the_BENCHMARKED_prompt_verbatim():
    """The full run inherits the benchmark's measured performance ONLY if it runs the same prompt. If
    these drift, the 5/5-TP result does not transfer and the run means nothing."""
    from kaggle_staleness_auditor import MAX_NEW_TOKENS, SYSTEM_PROMPT
    k = (ROOT / "scripts/kaggle/staleness_corpus_kernel.py").read_text(encoding="utf-8")
    assert SYSTEM_PROMPT in k, "corpus kernel's system prompt drifted from the benchmarked one"
    assert f"MAX_NEW_TOKENS = {MAX_NEW_TOKENS}" in k
    for field in ("CLAIM (from project documentation)", "ARTIFACT FACTS", "ARTIFACT EXCERPT"):
        assert field in k


# --- the adjudication record + the base-rate arithmetic that explains the benchmark-vs-field gap ---

def test_adjudication_covers_every_flag_the_pass_raised():
    from scripts.staleness_adjudication import ADJUDICATIONS, tally
    assert tally()["flags"] == len(ADJUDICATIONS) == 9
    assert {a.verdict for a in ADJUDICATIONS} <= {"true_positive", "false_positive"}
    # every adjudication must carry a REASON -- a bare verdict is not an adjudication
    for a in ADJUDICATIONS:
        assert len(a.why) > 60, f"{a.artifact} has no substantive reasoning"


def test_base_rate_explains_the_precision_gap_not_a_regression():
    """The load-bearing claim: the SAME detector posts 0.8 on a 50/50 benchmark and 0.22 in the field.

    If this ever fails, the "it's the base rate" explanation is wrong and the detector really did regress.
    """
    from scripts.staleness_adjudication import precision_from_base_rate, specificity
    spec = specificity()
    at_benchmark = precision_from_base_rate(0.50, sens=1.0, spec=spec)
    at_field = precision_from_base_rate(0.025, sens=1.0, spec=spec)
    assert at_benchmark > 0.85, "a 50/50 set should flatter precision"
    assert at_field < 0.35, "a 2.5%-positive corpus should crush precision at the same specificity"
    # and the field specificity is NOT worse than the benchmark's 4/5 -- the detector did not degrade
    assert spec >= 0.80


def test_field_specificity_beats_the_benchmark_so_precision_drop_is_not_degradation():
    from scripts.staleness_adjudication import specificity, precision
    assert specificity() > 0.80 > precision(), (
        "the whole point: specificity went UP while precision went DOWN")
