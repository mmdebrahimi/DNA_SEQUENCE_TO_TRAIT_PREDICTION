# Pre-registration — the 270-seed null resolution for the fragment-pyramid arm

<!-- doc-status: LIVE as of 2026-10-30. This is a PRE-REGISTRATION: it is written BEFORE any seed
     in the range 50..269 exists, and its job is to fix the reading of the outcome in advance.
     If you are reading it after the batches have landed, its value is as a FROZEN record of what
     was predicted; do not edit the predictions. -->

**DECISION THIS DISCHARGES:** the user's call of 2026-10-30 on §2.1 of
`wiki/session_handoff_2026-10-30.md` — resolve the arm's null-median validity gate by **spending
seeds**, chosen over re-specifying the bar and over accepting no verdict.

**WHAT MUST NOT HAPPEN:** the result being read after the fact in whichever direction suits. The
most likely single outcome, forecast below, **confirms a small real positive bias rather than
clearing the gate**, so the reading is fixed here first.

---

## 1. The state being resolved

| quantity | value | source |
|---|---|---|
| bar | `\|null median\| < 0.02` | `scripts/glm_fragment_pyramid_arm.py` `validity_controls.null_median_bar` |
| seeds measured to date | 50 (0..9 twice-published + 10..49 fresh) | `wiki/fragment_null_ext_summary_2026-10-29.json` |
| pooled null median | **+0.0242** | pooled over those 50 |
| pooled sd | 0.0654 | ibid. |
| measured bootstrap CI95 | **`[-0.0126, +0.0466]` — CONTAINS ZERO** | ibid. |
| consequence | `verdict: INVALID_NULL_SCORED`; nothing may be graded | arm's verdict ladder |

Underneath the veto, the pre-registered comparison on seeds 0..9 is already measured:
`paired.median −0.0127`, `wins 2/10`, `ci95 [-0.0205, -0.0040]` → the verdict that would land is
**`NO_GAIN`**. The gate is suppressing a *negative*, not a positive.

## 2. Target and sizing

n = 270, where `SE(median) ≈ 1.2533·sd/√n` reaches **0.0050**, a quarter of the bar.
Four batches of 55 seeds over 50..269 (contiguous, verified by assertion in the kernel), at the
measured 4.945 min/seed ≈ **4.5 h actual per batch, ~18.1 h total T4**.

## 3. THE SIZING MODEL IS NOT THE DECISION STATISTIC — and they disagree

Normal theory at n=50 gives `[+0.0015, +0.0469]`, which **excludes zero**. The measured bootstrap at
the same n gives `[-0.0126, +0.0466]`, which **contains zero**. They disagree on the one question
that matters, so the median's sampling distribution here is not well described by `1.2533·sd/√n` —
unsurprising with sd inflated by a heavy-tailed per-seed spread (individual seeds reach ±0.13).

**Therefore:** the normal-theory numbers in §2 and §4 are for **SIZING ONLY**. The decision
statistic is the **bootstrap CI recomputed on the pooled per-seed values**, by the same method that
produced `[-0.0126, +0.0466]`. Anyone who reports a normal-theory interval as the verdict has used
the wrong instrument.

## 4. FORECAST — recorded before the data

Assuming the true median equals today's point estimate (0.0242):

| pooled n | SE (normal theory) | CI straddles the 0.02 bar | P(point estimate < bar) |
|---|---|---|---|
| 50 (now) | 0.0116 | yes | 0.359 |
| 105 (after b1) | 0.0080 | yes | 0.300 |
| 160 (after b2) | 0.0065 | yes | 0.258 |
| 215 (after b3) | 0.0056 | yes | 0.226 |
| **270 (after b4)** | **0.0050** | **yes** | **0.200** |

**Two predictions, both falsifiable:**

- **P1 — the gate most likely still fails.** ~80% chance the n=270 point estimate is ≥ 0.02.
- **P2 — n=270 does NOT resolve the bar with confidence.** The CI at 270 is forecast to *straddle*
  0.02 (`[+0.0144, +0.0340]`). Excluding the bar from above would need **~1,463 seeds ≈ 121 h of
  T4** — four weeks of the 30 h/week Kaggle quota. **So 270 seeds buys a point estimate, not a
  confident resolution.** This is a property of the chosen option, stated before it ran, not a
  reason to abandon it: the point estimate is what the arm's rule actually consumes.

P1 and P2 are independent: P1 is about the arm's implemented point test, P2 about confidence.

## 5. THE READING, fixed in advance

Let `M` be the pooled median over all 270 seeds and `CI` its recomputed bootstrap interval.

- **A — CLEAN NULL.** `|M| < 0.02`. → The validity gate PASSES. The arm's pre-registered verdict on
  seeds 0..9 becomes computable and is **`NO_GAIN`**; the arm closes as a measured negative. Report
  `CI` beside `M`; if `CI` still straddles 0.02, say so — a passing point test on a straddling
  interval is a *weak* pass and must be labelled one.
- **B — REAL SMALL POSITIVE BIAS.** `M ≥ 0.02` and `CI` excludes zero. → The bar is failing
  CORRECTLY: the shuffled-train-label null genuinely scores. `INVALID_NULL_SCORED` stands, and the
  frontier becomes a **MECHANISM question — what makes a shuffled-label null positive?** That is a
  validity defect in the pipeline and a real finding, not a shrug. Candidate mechanisms to test, in
  the order they are cheap: the train/eval block split leaking composition; the early-stopping rule
  selecting on eval; an asymmetry in how the null's ridge is fit. **Do not** re-specify the bar at
  this point — a bar that is correctly detecting a defect is doing its job.
- **C — STILL UNRESOLVED.** `M ≥ 0.02` but `CI` contains zero. → The seed lever is **retired by
  measurement**: 270 seeds was the prescribed n and it did not settle it. The decision reverts to the
  two options the user declined on 2026-10-30, now with the third closed. Per §4-P2 this is a
  *likely* outcome and is NOT a failed run.

**In all three cases the pyramid-vs-kmer4 science is unchanged** and does not depend on this gate:
paired CI excludes zero on the negative side at both seed ranges (−0.0127 `[-0.0205,-0.0040]` at
0..9; −0.0061 `[-0.0133,-0.0037]` at 10..49). The pyramid loses to a 4-mer baseline either way. This
gate decides what the LEDGER may record, not what is true.

## 6. Validity conditions for the run itself

- **`PINNED_SHA` must stay `74b7a28df74d46974f36c5f55ed601b9cca3b156`** — it matches the `head` of
  the 10..49 artifact, so all 270 seeds come from byte-identical code. A bumped pin makes the pooled
  median a mixture of two estimators and voids the pooling.
- **Each batch pushes to its OWN slug.** Kaggle keeps one output per kernel; re-pushing
  `glm-fragment-null-ext` would destroy the retrieved 10..49 artifacts.
- **No batch may emit a pre-registered verdict.** Each reports a shifted seed range and so returns
  `PARTIAL_SEED_RANGE` / `verdict: null` by design. A batch that emitted a verdict would be
  re-adjudicating an answered comparison on a range chosen after the fact.
- **The `--seed-start` shift is confirmed, not trusted.** The kernel refuses (exit 4) if the arm
  reports a `seed_range.start` other than the batch's own — a silently ignored flag would redraw
  already-measured seeds and resolve nothing while looking healthy.
