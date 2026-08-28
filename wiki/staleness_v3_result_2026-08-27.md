# v3 FAILS its pre-registered bar — and I destroyed my own baseline getting there

The pre-registration (committed at `868b401`, before the fix was written) did its job twice: it caught an
over-correction the aggregate score would have hidden, and then it exposed that my attribution was
confounded by variables I changed myself.

## The verdict

**FAIL.** The bar required both hand-verified true positives to survive. One did.

| | v2 | **v3** |
|---|---:|---:|
| flags | 10 | **7** |
| **true positives** (bar: ≥2) | 2 | **1** |
| false positives | 8 | 6 |
| unparseable | 4 | **1** |
| **aggregate bar** | FAIL | **FAIL** |
| targeted flips (the diagnosis) | 1/4 | **2/4** |

Precision improved and recall got worse. **That is the trade the must-hold prediction existed to make
visible** — without it, "flags 10 → 7, all three predicted flips landed" reads like a clean success while
the tool quietly stops finding things.

## What is solid

**The three targeted arms flipped exactly as predicted**, in the full 110-item run:

- **finding** — `bvbrc_strict_mic_4drug_census.py` → `supported`. Holds at both excerpt lengths.
- **correction** — `browser.py` → `supported` at 3000 chars.
- **historical** — `negative_results_map` → `supported` on the historical claim, *while still flagging the
  genuinely stale gate-count claim on the same file*. That discrimination between two claims sharing one
  artifact is the single best thing v3 did.

## What is NOT solid, and why that is my fault

**I changed three variables between v2 and v3**: the prompt, the excerpt cap (6000 → 3000, my OOM fix),
and — decisively — **the corpus text itself, by fixing the two stale claims the v2 run had found.**

A single-variable isolation run (same v3 prompt, 6000-char excerpts) shows how much that matters:

- `test_models_cache.py` → **`stale` again**. So the must-hold MISS was driven by the excerpt cap, not the
  prompt: at 3000 chars the model said *"the artifact's content is cut off mid-test"* and answered
  `unclear`, which the fail-closed parser correctly refuses to treat as a flag.
- `browser.py` → **`stale` again**. So the correction arm did **not** hold at 6000. Its apparent success at
  3000 was at least partly the excerpt length, not the subordination fix.

But that isolation run rebuilt the corpus from *live* CLAUDE.md — which I had already repaired — so its
claim text differs from what v3 scored. **The attribution is therefore unresolved, and I am not going to
assert one.**

## The real lesson: fixing what you measure destroys the baseline

This is the methodological finding, and it generalizes past this tool:

> **A detector evaluated on a corpus it is also used to repair cannot be re-measured against that corpus.**
> Every true positive I fixed removed a positive from the test set. The base rate — already 0.018 — moves
> under the instrument, so v2 and v3 were never scored on the same thing.

The fix is cheap and I did not do it: **snapshot the corpus** at first measurement and score every later
version against the frozen snapshot, repairing only a copy. The v3 dataset happened to be a pre-fix
upload, which is the only reason any comparison was possible at all — luck, not design.

## A bug in my own scorer, found by disbelieving the result

The first v3 score read **0 true positives**, which would have meant the fix destroyed the tool. It was
wrong. `score()` built `{artifact: verdict}` as a dict, and `negative_results_map` carries two claims — so
the second silently **overwrote** the first and deleted a correctly-caught true positive.

This is the shared-key overwrite trap. I had written a test for it in `_truth()` in the same file, and left
the identical bug ten lines below. It is now keyed by item, with the control test (reproduce v2's 2 TP)
that catches exactly this class.

## Status

- **v3 is NOT adopted.** It fails its pre-registered bar, and its wins are confounded.
- **v2 remains the measured configuration** (11 flags / 2 TP / 9 FP by hand adjudication).
- The subordination idea is *supported but unproven*: the finding arm is robust across excerpt lengths; the
  correction arm is not; the historical arm produced a genuine two-claims-one-artifact discrimination.

**Next, in order:** freeze a corpus snapshot; re-run v2 and v3 against it at a single fixed excerpt length;
only then attribute. Roughly one kernel-hour, and it converts three confounded variables into one.
