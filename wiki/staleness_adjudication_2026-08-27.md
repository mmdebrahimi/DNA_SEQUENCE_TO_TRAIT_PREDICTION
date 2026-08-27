# Adjudicating all 11 flags: 2 real, 9 false — and my benchmark could never have predicted that

Follow-up to `wiki/staleness_corpus_run_2026-08-27.md`. The run is COMPLETE at 110/110 and every flag
has now been checked by hand against its artifact. The result is worth more than the two fixes it produced, because it
says something about the *benchmark* rather than the model.

## The tally

**UPDATED — the run is now COMPLETE at 110/110.** The remaining 30 finished with the capped excerpt
(no OOM, 0 unparseable) and raised 2 more flags, both false positives.

| | 80-item | **full 110** |
|---|---:|---:|
| flags raised | 9 | **11** |
| **true positives** | 2 | **2** |
| false positives | 7 | **9** |
| **precision in the field** | 0.22 | **0.18** |
| specificity | 0.910 | **0.917** |
| base rate of stale claims | 0.025 | **0.018** |

Re-derive with `uv run python scripts/staleness_adjudication.py`.

## The headline is not "it got worse" — specificity went UP

The benchmark measured **4/5 precision (spec 0.80)** on a set I built 50% positive. The field run posts
**0.18 precision at spec 0.917**. The detector did not degrade — it is *better* at not flagging clean items
than the benchmark suggested. Precision collapsed for one reason:

```
base rate 0.500 -> expected precision 0.92    [the benchmark]
base rate 0.018 -> expected precision 0.18    [the real corpus]
```

Same specificity, same sensitivity, both numbers computed by `precision_from_base_rate`. **A screen with
excellent specificity still returns mostly false positives when the thing it screens for is rare** — and a
documentation corpus is exactly that regime: 2 stale claims in 110.

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

## The false positives, and what they teach

One is the **capability-vs-finding** category error, now confirmed on an unseen item (`bvbrc_strict_mic_4drug_census.py`)
after producing the benchmark's only FP. Code existing refutes a *capability* claim ("X is not built"); it
**supports** a *finding* claim, because the script is the instrument that produced the finding.

The others are ordinary misreadings, and three are notable because the model **asserted the opposite of
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
right"** — and on the full run, closer to 1-in-6. That is still worth running: 2 real stale claims
surfaced from 110 bullets nobody re-reads, at the cost of checking 11 candidates — about ten minutes. The tool earns its keep as a
**triage funnel**, not an oracle.

**It strengthens, not weakens, the standing rule:** editing a doc on the model's say-so was never
licensed, and at 0.18 precision it is further from licensed than the benchmark implied.

## Scope

- Base rate 0.018 is measured on this corpus after two of its stale claims were already fixed by earlier
  sweeps — a fresher corpus would likely run slightly higher.
- Sensitivity is **not** measured in the field: I can count the flags that were wrong, not the stale claims
  that were never flagged. The benchmark's 5/5 recall is the only recall number, and it comes from a set
  where the positives were curated to be findable.


## The last 30 items, and the single most diagnostic flag of the run

The rerun (capped excerpt, items inlined rather than dataset-attached) completed cleanly: **26 supported,
2 flagged, 2 unclear, 0 unparseable.** Both flags are false positives, and one of them is the sharpest
result the whole exercise produced.

**`dna_decode/genome_map/browser.py`** — flagged *stale* because "FACTS show browser.py exists with 9
implemented functions, contradicting the claim it was 'deferred'". But the only "deferred"+browser text in
CLAUDE.md **is the correction recording that the browser shipped on 2026-07-11**. My system prompt names
this case explicitly and says it is `supported`:

> *"A sentence that says 'X was deferred, and here is the correction recording that it shipped' is
> `supported` — it is accurate text about a past error, not a stale claim."*

**The v2 `ARTIFACT FACTS` rule overrode the prompt's own stated exception.** That is the third instance of
the same unscoped instruction (after the benchmark's N5 and the field census flag), and it is now
characterised precisely: *"implemented code exists → stale"* silently outranks every exception written
above it, including the ones written to catch exactly this.

The second flag re-flagged `negative_results_map` on the claim that the colour family *"needed and did not
have"* G9/G10 — a **past-tense** statement about what was missing, not refuted by the gates now existing.
Notably it is the same file as the run's one clean true positive: the pass flagged both the bullet that
was genuinely stale and an accurate historical one about the same artifact.

## Final scope

- **110/110 pairs scored.** Nothing is unrun.
- The **capability-vs-finding** fix is now confirmed necessary on three independent instances, and its
  scope is sharper than when it was recorded: it must not merely be *scoped* to capability claims, it must
  be subordinated to the correction-text exception it currently overrides.
- Sensitivity remains unmeasured in the field — flags that were wrong are countable; stale claims never
  flagged are not.
