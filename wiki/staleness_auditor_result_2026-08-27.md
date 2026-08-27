# F2 PASSES: a semantic auditor beats the mechanical screen that failed (2026-08-27)

**Verdict: PASS on the pre-registered bar, and it beat the baseline.** Qwen3-8B on a free Kaggle T4 scored
**3/5 true positives, 0/5 false positives** against a pre-registered `>=3 TP AND <=1 FP`. The mechanical
proximity screen this replaces scored **0/5 TP and 5/5 FP**.

The threshold and the benchmark were committed **before** the run (`aeb1114`), so this is a falsification
that could have failed and didn't.

## The numbers

| auditor | TP | FN | FP | TN | passes bar? |
|---|---:|---:|---:|---:|---|
| mechanical proximity screen (the incumbent) | 0 | 5 | **5** | 0 | no |
| **Qwen3-8B, free Kaggle T4** | **3** | 2 | **0** | **5** | **yes** |

Artifact: `wiki/staleness_auditor_result_2026-08-27.json`. Kernel:
`emanueleebrahimi/dna-staleness-auditor` (private, free tier, deletable).

## The result that matters most is the zero

**0 false positives — on the exact 5 regions the mechanical screen flagged.** That is the whole reason the
mechanical version was never shipped: a ~100% FP rate gets a guard disabled. Every one of those regions is
prose where a deferral marker sits near a file path but refers to something else, and in three cases the
marker lives inside *correction* text explaining a past fix. The model declined all five.

So the failure mode that killed the mechanical approach is the one the semantic approach handles best.

## Both misses are diagnosed, and neither is arbitrary

**P3 (genome-map browser) — a genuine trap, and a fair one.** The claim was "a visual browser deferred".
The artifact is `dna_decode/genome_map/browser.py`, whose own module docstring reads *"Standalone graphical
browser for a genome-map JSON — the deferred v1 'graphical browser'"*. The model quoted exactly that and
concluded the claim was supported. It is describing what the module *implements*; the file's existence is
the refutation. A careful human reading only the excerpt could make the same mistake — which says the
excerpt, not the model, is the weak link.

**P4 (BV-BRC census) — a budget artifact, not a judgment failure.** The output was truncated mid-`<think>`
at the 1200-token generation cap, so the parser fell through to `unclear` — the fail-closed path working as
designed. It is scored as a miss (correctly: hedging must not earn credit), but it is not evidence the
model cannot judge the item. It never finished reading.

**The honest read: 3/5 understates the ceiling.** One miss is a fixable excerpt/prompt problem, one is a
fixable token budget. Neither is "the model cannot do this". But the number recorded is the number measured,
and the fixes must be tested rather than assumed.

## What this licenses, and what it does not

**Licenses:** running the full corpus (542 memos / ~525k tokens + CLAUDE.md) as a *flagging* pass, with
every flag adjudicated against its artifact by hand before any doc changes.

**Does NOT license:** editing any documentation on the model's say-so. A 3/5-recall flagger is a search
aid, not an authority. The adjudication step is the product; the model only decides where to look.

Also note what the benchmark *is*: 10 items, one model, one prompt. It shows the semantic approach clears a
bar the mechanical one could not. It does not measure recall on the real corpus, where claims are messier
and the artifact excerpt may not contain the deciding text (exactly the P3 failure at scale).

## Method notes worth keeping

- **The negatives were regenerated, not transcribed.** The mechanical screen's 5 false positives were
  described in a code comment but never persisted, so `scripts/staleness_benchmark.py` re-runs that screen
  to rebuild them. A hand-typed negative set would have been my recollection of what the screen found.
- **The baseline is computed, not quoted.** `baseline_verdicts()` re-derives 0/5 and 5/5 from the live
  screen — independent confirmation that the prose figure was accurate.
- **The score is non-gameable in both directions,** pinned by tests: flag-everything gets 5/5 TP and still
  fails on 5 FP; flag-nothing and hedge-everything both fail. `unclear` counts as not-flagged, which is
  what makes the hedge strategy lose.
- **The parser fails closed.** A parser that guessed `stale` on malformed output would manufacture true
  positives — and because hedging scores as not-flagged, that bug would have looked like a cautious model
  rather than a defect. P4 is the live proof it works.

## Next, in order

1. **Raise `max_new_tokens`** (1200 → ~2500) and re-run the benchmark. P4 alone may move 3/5 → 4/5, and it
   is one number in one file.
2. **Widen the artifact excerpt** or select the relevant section rather than the head — P3's deciding
   evidence was the file's *existence and purpose*, not its opening lines.
3. **Only then** run the full corpus, adjudicating every flag by hand.

Re-run: `uv run python scripts/kaggle_staleness_auditor.py --emit-kernel` then push via
`scripts/kaggle_push_poll.py`. Score with `--score <results.json>`.
