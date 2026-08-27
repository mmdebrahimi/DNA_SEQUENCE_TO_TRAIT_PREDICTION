# Adjudicating all 9 flags: 2 real, 7 false — and my benchmark could never have predicted that

Follow-up to `wiki/staleness_corpus_run_2026-08-27.md`. Every flag from the 80-item pass has now been
checked by hand against its artifact. The result is worth more than the two fixes it produced, because it
says something about the *benchmark* rather than the model.

## The tally

| | |
|---|---:|
| flags raised | 9 |
| **true positives** | **2** |
| false positives | 7 |
| **precision in the field** | **0.22** |
| specificity | **0.910** |
| base rate of stale claims | 0.025 (2 in 80) |

Re-derive with `uv run python scripts/staleness_adjudication.py`.

## The headline is not "it got worse" — specificity went UP

The benchmark measured **4/5 precision (spec 0.80)** on a set I built 50% positive. The field run posts
**0.22 precision at spec 0.910**. The detector did not degrade — it is *better* at not flagging clean items
than the benchmark suggested. Precision collapsed for one reason:

```
base rate 0.500 -> expected precision 0.92    [the benchmark]
base rate 0.025 -> expected precision 0.22    [the real corpus]
```

Same specificity, same sensitivity, both numbers computed by `precision_from_base_rate`. **A screen with
excellent specificity still returns mostly false positives when the thing it screens for is rare** — and a
documentation corpus is exactly that regime: 2 stale claims in 80.

**The design error is mine, and it is the reusable lesson: a 50/50 benchmark is structurally incapable of
predicting field precision.** I built the benchmark balanced because that is how one measures a
classifier — but the number a *user* experiences is precision, and precision is a function of the base
rate the benchmark deliberately destroyed. The benchmark validated the right quantity (specificity) and I
read the wrong one off it.

Three tests now pin this so the explanation can't quietly become an excuse — including one that fails if
field specificity ever drops below the benchmark's, which is what a genuine regression would look like.

## The two real catches

1. **`wiki/negative_results_map_2026-06-13.md`** — CLAUDE.md said *"8 reusable rejection GATES … G1–G8"*;
   I added G9/G10 the day before and never updated the bullet. Fixed.
2. **`tests/test_models_cache.py`** — claim said *"8 unit tests"*; the file has **42**. Fixed by removing
   the hardcoded count (the same drifting-count class already fixed in README).

**Honesty about #2: right file, wrong reason.** The model's stated evidence — *"the tests do not reference
`verify_complete`"* — is simply false; it appears 19 times. The claim is genuinely stale on the count, so
it scores as a true positive, but the *evidence* would not have survived adjudication on its own. That is
an argument for adjudication, not against the flag.

## The seven false positives, and what they teach

One is the **capability-vs-finding** category error, now confirmed on an unseen item (`bvbrc_strict_mic_4drug_census.py`)
after producing the benchmark's only FP. Code existing refutes a *capability* claim ("X is not built"); it
**supports** a *finding* claim, because the script is the instrument that produced the finding.

The other six are ordinary misreadings, and three are notable because the model **asserted the opposite of
the text in front of it**:

- `prosst_scorer.py` — claimed the artifact says it ran on Kaggle. The wiki has a section headed
  *"Local, not Kaggle"* and the line *"Ran entirely LOCALLY on CPU — no Kaggle"*.
- `decoder_v0_ux_and_success_criterion.md` — read the heading *"Explicit non-criteria for v0"* correctly
  and then drew the inverted conclusion.
- `provdisjoint_census_results.json` — conflated two different fields (`powered: false` → UNDERPOWERED, 3
  cells, which the claim itself reports) with the NOT_CENSUSED bucket (verified 0 of 27).

These are not subtle-judgment failures that a better prompt fixes cheaply. They are the reason the
adjudication step is the product.

## What this changes

**The framing.** Not "1-in-5 flags is wrong" (the benchmark's number) but **"roughly 1-in-4 flags is
right"**. That is still worth running: 2 real stale claims surfaced from 80 bullets nobody re-reads, at
the cost of checking 9 candidates — about ten minutes of adjudication. The tool earns its keep as a
**triage funnel**, not an oracle.

**It strengthens, not weakens, the standing rule:** editing a doc on the model's say-so was never
licensed, and at 0.22 precision it is further from licensed than the benchmark implied.

## Scope

- 80 of 110 pairs. The remaining 30 are running with a capped excerpt (the OOM fix) and **inlined** rather
  than dataset-attached, which removes the failure mode that killed the first corpus attempt.
- Base rate 0.025 is measured on this corpus after two of its stale claims were already fixed by earlier
  sweeps — a fresher corpus would likely run slightly higher.
- Sensitivity is **not** measured in the field: I can count the flags that were wrong, not the stale claims
  that were never flagged. The benchmark's 5/5 recall is the only recall number, and it comes from a set
  where the positives were curated to be findable.
