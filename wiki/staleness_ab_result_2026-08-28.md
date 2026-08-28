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

## Result 3 — the deployed config measured, and it is the best cell

I adopted v3 @ 6000 without having measured that exact combination, which is the gap my own guards exist
to prevent, so I ran it. It OOM'd at item ~80 (below), but all three ground-truth items completed, so the
80 that finished are a valid like-for-like set:

| config (same 80 items) | recall | false positives |
|---|---:|---:|
| v2 @ 6000 | 2/3 | 7 |
| v2 @ 3000 | 1/3 | 7 |
| v3 @ 3000 | 1/3 | 5 |
| **v3 @ 6000** | **2/3** | **3** |

**Stated at the strength the evidence supports:** on this completed 80-item prefix, the longer excerpt
recovered one additional positive in both prompt arms, and v3 reduced false positives in both excerpt
arms. That is a descriptive interaction check, **not** an independence result — with 3 positives, recall
moves in single-item jumps. The deployed cell more than halves false positives at equal recall.

The missing 30 are a **contiguous document-order tail** (completed items are exactly positions 0-79), not
a length-selected subset — measured, median prompt 6522 chars completed vs 6794 missing. So the cut is not
OOM-length bias; it is a recency gradient (CLAUDE.md is roughly chronological, so the tail is the newest
and least-rotted content). All 3 known positives fall in the prefix, and the tail is unadjudicated.

## CORRECTION — the memory fix I claimed does not work

I wrote, in the previous commit, that the OOM was handled "where it belongs" with a per-item
`empty_cache()`. **That is wrong and the run disproved it**: the cache-clearing version died at the same
item. The OOM is a single-item PEAK (one 3.94 GiB allocation on a 14.56 GiB T4), not fragmentation
accumulating between items, so freeing between items cannot help.

The mitigation now bounds the TOTAL sequence — capping the *generation* for long prompts rather than the
input — because the excerpt is evidence and the reasoning length is not (the P4 fix showed ~1600-3000
tokens is typical, so the cap rarely binds). **It is NOT YET VERIFIED at full-corpus scale**, is labelled
so in the file, and a test requires that label to stay until a clean 110/110 run lands.

## What changed as a result

- **v3 is adopted.** At matched conditions it costs nothing and cuts false positives by a third. The
  guard that said "do not deploy v3" is updated — its *premise* was refuted, its job is unchanged: deploy
  only what was measured.
- **The 3000-char cap is reverted to 6000.** The excerpt is evidence; the OOM is bounded by capping the
  GENERATION for long prompts instead (see the correction above — the `empty_cache()` idea did not work).
  Pinned by a test asserting the kernel truncates nothing and that the corpus builder keeps >= 6000.

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

- deployed prompt: **v3**, excerpt **6000**, OOM bounded by an adaptive generation budget (unverified)
- **v3 @ 6000 is measured for PROMPT and EXCERPT only** — recall 2/3 with 3 false positives against v2 @
  6000's 2/3 with 7 — and it was measured under a FIXED 2500-token generation. The deployed kernel adds an
  adaptive generation budget that no measured run used, so "deployed = measured" holds for the prompt and
  excerpt, not for the runtime policy.
- the OOM mitigation is a **hypothesis**, not a result: a clean 110/110 run has not happened.
- the flagging pass remains a **triage funnel**: every flag adjudicated by hand, no doc edited on the
  model's say-so.
