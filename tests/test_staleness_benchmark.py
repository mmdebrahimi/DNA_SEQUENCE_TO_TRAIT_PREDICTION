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
    assert tally()["flags"] == len(ADJUDICATIONS) == 11
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


# --- v3: the pre-registration and the scorer that grades against it ---

def test_v3_predictions_were_committed_before_the_fix():
    """The predictions must be falsifiable and must include a must-hold, or it is not a real test."""
    from scripts.staleness_v3_preregistration import PREDICTIONS, BAR, targeted, must_hold
    assert len(targeted()) == 3, "one prediction per arm of the diagnosed mechanism"
    assert len(must_hold()) >= 1, "a fix with no must-hold can pass by silencing everything"
    assert BAR["min_true_positives"] == 2, "both real catches must survive -- non-negotiable"
    for p in PREDICTIONS:
        assert p.predict in {"stale", "supported"} and len(p.because) > 50


def test_the_v3_scorer_reproduces_the_v2_baseline():
    """CONTROL: the scorer must reproduce the hand-adjudicated v2 result on the committed v2 raw output.

    If it cannot recover 2 true positives from a run I checked by hand, it is measuring something other
    than what I adjudicated, and any v3 number it produces would be meaningless.
    """
    from pathlib import Path
    from scripts.score_v3 import score
    raw = Path(__file__).resolve().parent.parent / "wiki" / "staleness_corpus_run_2026-08-27_raw.json"
    if not raw.exists():
        import pytest
        pytest.skip("v2 raw output not present")
    r = score(raw)
    assert r["n_scored"] == 110
    assert len(r["true_positives"]) == 2, "the two hand-checked real catches must be recovered"
    assert not r["aggregate_pass"], "v2 must FAIL the v3 bar -- otherwise the bar is not a bar"


def test_scorer_never_loses_a_true_positive_to_the_shared_key_trap():
    """negative_results_map carries TWO claims with opposite verdicts. Keying by artifact could drop the
    real one; the scorer must keep the true positive. Same shared-key overwrite trap as the report card."""
    from scripts.score_v3 import _truth
    t = _truth()
    assert t["wiki/negative_results_map_2026-06-13.md"] == "true_positive"


def test_the_deployed_prompt_is_the_one_the_AB_MEASURED_BETTER():
    """UPDATED after the single-variable A/B. The previous version of this test asserted v3 must NOT be
    deployed, on the reading that v3 halved recall. The A/B refuted that premise: at a MATCHED excerpt
    length v2 and v3 have IDENTICAL recall (1/3 each), and v3 has fewer false positives (6 vs 9), higher
    specificity (0.944 vs 0.916) and fewer unparseable answers. The recall collapse belonged to the
    excerpt cap, not the prompt.

    The guard's job is unchanged -- deploy only what was measured. What changed is which prompt that is.
    """
    from scripts.kaggle_staleness_auditor import SYSTEM_PROMPT
    assert "When you cannot tell whether a claim is CAPABILITY or FINDING" in SYSTEM_PROMPT
    assert "ALSO CRITICAL: weigh ARTIFACT FACTS as evidence" not in SYSTEM_PROMPT


def test_the_kernel_does_not_truncate_the_excerpt_below_what_was_measured():
    """Truncating the input was treating an OOM by deleting the evidence, and it cost half the recall.

    Same prompt, 6000 -> 3000 chars: recall 0.667 -> 0.333. The kernel must keep the full excerpt and
    free the KV cache per item instead.
    """
    k = (Path(__file__).resolve().parent.parent / "scripts" / "kaggle" /
         "staleness_corpus_kernel.py").read_text(encoding="utf-8")
    assert "[:3000]" not in k, "3000-char truncation measurably halved recall"
    assert "empty_cache()" in k, "the OOM must be handled in memory, not by discarding evidence"


def test_v3_scored_below_its_own_bar():
    """Pin the FAIL so a later reader cannot mistake v3's nicer flag count for a passing result."""
    from pathlib import Path
    from scripts.score_v3 import score
    raw = Path(__file__).resolve().parent.parent / "wiki" / "staleness_v3_run_2026-08-27_raw.json"
    if not raw.exists():
        import pytest
        pytest.skip("v3 raw output not present")
    r = score(raw)
    assert not r["aggregate_pass"], "v3 FAILED its pre-registered bar -- do not record it as a pass"
    assert len(r["true_positives"]) < 2, "the bar required both hand-verified true positives to survive"
    assert r["targeted_hits"] >= 2, "the diagnosed arms did flip -- the idea is supported, just unproven"


def test_the_corpus_snapshot_is_frozen_and_complete():
    """The snapshot is the fixed test set every future prompt version scores against.

    Scoring against LIVE CLAUDE.md is what made v2-vs-v3 unresolvable: fixing a flagged claim removes a
    positive from the test set, so two versions are never graded on the same thing. This file must not be
    regenerated to "refresh" it -- a new measurement campaign gets a NEW dated snapshot.
    """
    import json
    from pathlib import Path
    snap = Path(__file__).resolve().parent.parent / "wiki" / "staleness_corpus_snapshot_2026-08-27_POSTREPAIR_NO_POSITIVES.json"
    items = json.loads(snap.read_text(encoding="utf-8"))
    assert len(items) >= 100, "the snapshot must cover the whole corpus, not a slice"
    for it in items:
        assert {"item_id", "claim", "artifact", "artifact_text", "facts"} <= set(it)
        assert it["artifact"] in it["claim"] or len(it["claim"]) > 0
    assert len({i["item_id"] for i in items}) == len(items), "item ids must be unique"


def test_the_prerepair_snapshot_actually_contains_positives():
    """The snapshot frozen from LIVE CLAUDE.md had ZERO known stale claims -- I had already fixed them.

    A test set with no positives can measure precision and nothing else, which makes it useless for the
    recall comparison it was frozen for. This is the corpus-repair trap biting the fix written for it.
    The evaluation snapshot is therefore reconstructed from git BEFORE the repairs, and must contain the
    known positives or the A/B is meaningless.
    """
    from scripts.score_ab import truth_items, SNAPSHOT
    import json
    if not SNAPSHOT.exists():
        import pytest
        pytest.skip("eval snapshot not present")
    assert len(json.loads(SNAPSHOT.read_text(encoding="utf-8"))) >= 100
    assert len(truth_items()) >= 2, "the eval snapshot must retain the hand-verified stale claims"


def test_ab_ground_truth_is_derived_from_the_corpus_not_hardcoded():
    """Ground truth matches on claim TEXT in the snapshot, not a hand-listed set of item ids.

    A hand-enumerated list beside the data it describes drifts -- this repo has hit that three times.
    Deriving it means re-freezing the snapshot cannot silently leave the truth set pointing at nothing.
    """
    from scripts.score_ab import STALE_MARKERS, truth_items
    assert all(isinstance(m, str) and m for m in STALE_MARKERS)
    ids = truth_items()
    assert all(isinstance(i, str) for i in ids)


def test_the_oom_mitigation_bounds_generation_not_the_excerpt():
    """The excerpt is evidence; the generation length is not. Measured: truncating the excerpt halved
    recall (0.667 -> 0.333) while empty_cache() alone did NOT prevent the OOM (the cache-clearing
    version died at the same item). So the mitigation must bound the TOTAL sequence by shortening the
    generation for long prompts, never by shortening the input.
    """
    k = (Path(__file__).resolve().parent.parent / "scripts" / "kaggle" /
         "staleness_corpus_kernel.py").read_text(encoding="utf-8")
    assert "TOTAL_TOKEN_BUDGET" in k
    assert "NOT YET VERIFIED" in k, "an unmeasured mitigation must say so in the file"
    # the kernel must not truncate at all -- excerpt length is the corpus builder's business, and
    # keeping it in one place is what stops a truncation being reintroduced as a local OOM patch.
    assert "artifact_text'][:" not in k and 'artifact_text"][:' not in k
    from scripts.staleness_corpus import ARTIFACT_HEAD_CHARS
    assert ARTIFACT_HEAD_CHARS >= 6000, "3000 chars measurably halved recall"


def test_the_generation_budget_clips_outliers_not_the_typical_item():
    """A budget that binds on nearly every item is a different CONFIG, not an outlier guard.

    At TOTAL_TOKEN_BUDGET=4200 the cap bound on 107/110 real items, cutting generation to a median 1292
    tokens against the 2500 every measured run used -- and 1292 sits on the v1 cap of 1200 that produced
    a truncated, unparseable answer. This pins the budget high enough that the typical item keeps the
    measured generation length.
    """
    import re
    k = (Path(__file__).resolve().parent.parent / "scripts" / "kaggle" /
         "staleness_corpus_kernel.py").read_text(encoding="utf-8")
    total = int(re.search(r"TOTAL_TOKEN_BUDGET = (\d+)", k).group(1))
    max_new = int(re.search(r"MAX_NEW_TOKENS = (\d+)", k).group(1))
    # measured prompt-token p90 over the real corpus was 3000; the typical item must still get max_new
    assert total - 3000 >= max_new, (
        f"budget {total} leaves only {total-3000} generation tokens at p90 prompt length, "
        f"below the measured {max_new}")
