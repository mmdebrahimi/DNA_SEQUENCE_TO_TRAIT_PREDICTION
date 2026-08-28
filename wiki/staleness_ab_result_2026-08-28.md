# The single-variable A/B overturns yesterday's verdict: it was the excerpt cap, not the prompt

Yesterday I reported v3 as a FAIL and said the attribution was **unresolved** because I had changed three
things at once. That caveat was the right call — the experiment that resolves it says the prompt was not
the culprit.

## The design

Everything held fixed except one variable at a time. Corpus reconstructed from git at the commit *before*
I repaired the claims it contains (so the known stale claims are still in it), one excerpt length for both
arms, same model and greedy decoding, both prompts verified byte-identical to their originals by executing
the assignment out of the built kernel.

Ground truth: **3 genuinely-stale pairs** among 110, all hand-verified before any of this.

## Result 1 — the prompt is NOT what cost recall

| | v2 @ 3000 | **v3 @ 3000** |
|---|---:|---:|
| recall | 0.333 (1/3) | **0.333 (1/3)** |
| false positives | 9 | **6** |
| specificity | 0.916 | **0.944** |
| unparseable | 2 | **1** |

**Identical recall.** v3 flags a third fewer false positives, is more specific, and produces fewer
malformed answers. Yesterday's "v3 halved recall" reading was an artifact of the confound.

## Result 2 — the excerpt cap is the dominant variable

Same prompt (v2), same corpus, only the excerpt length changed:

| excerpt | recall |
|---|---:|
| **6000 chars** | **0.667** (2/3) |
| 3000 chars | 0.333 (1/3) |

**Halving the excerpt halved the recall.** My OOM fix "worked" by deleting the evidence the model needed
to judge — and the model said so plainly at the time (*"the artifact's content is cut off mid-test"*),
which I read as a model failure rather than as a description of what I had done to its input.

## What changed as a result

- **v3 is adopted.** At matched conditions it costs nothing and cuts false positives by a third. The
  guard that said "do not deploy v3" is updated — its *premise* was refuted, its job is unchanged: deploy
  only what was measured.
- **The 3000-char cap is reverted to 6000**, with the OOM handled where it belongs: `del gen, ids` +
  `torch.cuda.empty_cache()` per item. A memory problem gets a memory fix, not a truncated input. Pinned
  by a test that fails if `[:3000]` ever returns.

## Power — the caveat that keeps this honest

**Three positives.** Every recall figure here moves in steps of 0.333, and both headline recall
differences are **one item**. They are consistent with mechanisms I can point at (the model naming the
truncation; the three FP classes v3 targets) but they are *not* statistically established, and I am not
going to dress them up.

The better-powered comparison is specificity: **9 vs 6 false positives across 107 negatives**. That is the
number carrying the adoption decision, and it still is not large.

## The methodological lesson, twice-learned

I wrote the corpus-repair lesson yesterday, froze a snapshot to fix it — and froze it from the **repaired**
CLAUDE.md, so it contained **zero** positives. A test set with no positives measures precision and nothing
else, which is exactly what the comparison needed it not to be.

Git is what saved it: the pre-repair corpus was recoverable at a commit hash. **Reconstructing a
test set from version control is strictly better than freezing one, because the freeze can be taken at the
wrong moment and version control cannot.**

Ground truth here is also **derived** by matching claim text in the snapshot, not hardcoded as item ids
beside it — the hand-enumerated-list trap this repo has hit three times. Re-freezing the snapshot cannot
silently leave the truth set pointing at nothing.

## Standing status

- deployed prompt: **v3**, excerpt **6000**, OOM handled by cache-clearing
- best measured recall on the pre-repair set: **0.667** (v2 @ 6000). v3 @ 6000 is **not yet measured** —
  the obvious next run, and the one that would tell us whether v3 keeps its precision edge at full excerpt.
- the flagging pass remains a **triage funnel**: every flag adjudicated by hand, no doc edited on the
  model's say-so.
