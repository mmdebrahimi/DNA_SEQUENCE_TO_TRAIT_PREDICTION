# The fragment arm's null gate PASSES at 160 seeds — interim, branch A

<!-- doc-status: LIVE as of 2026-10-31, and INTERIM BY CONSTRUCTION: batches b3/b4 are in flight and
     will take the pool to 270. The headline below is expected to move slightly, not to flip. The
     reading of each outcome was fixed in advance in wiki/fragment_null_270seed_prereg_2026-10-30.md
     (sections 3a and 5); read those before quoting anything here. -->

## The result

| pool | n | null median | bar `\|median\| < 0.02` |
|---|---|---|---|
| seeds 0..49 (all that existed on 2026-10-30) | 50 | **+0.0242** | OUTSIDE — every verdict void |
| **seeds 0..159 (with b1 + b2)** | **160** | **+0.0112** | **WITHIN — branch A, CLEAN NULL** |

Exact order-statistic interval for the median: **`[-0.0126, +0.0354]`**, achieved coverage 0.9522.

**The gate passes. This is pre-registered branch A, and it is a WEAK PASS — which the
pre-registration required be labelled as such**, because the interval still contains the 0.02 bar.
The point estimate is what the arm's rule consumes (`abs(null_median) >= bar -> INVALID`), and the
point estimate is now comfortably inside.

## Why it moved, and why the decision to spend seeds was right

The two never-drawn batches are centred essentially on zero:

| sample | n | null median | sd |
|---|---|---|---|
| seeds 0..9 (twice-published) | 10 | +0.0517 | 0.0661 |
| seeds 10..49 | 40 | +0.0191 | 0.0649 |
| **seeds 50..104 (b1, new)** | **55** | **+0.0021** | 0.0723 |
| **seeds 105..159 (b2, new)** | **55** | **−0.0047** | 0.0656 |

**The pre-registered forecast was WRONG, and in the user's favour.** It gave ~80% odds the gate would
still fail at 270, reasoning that the spread is initialisation-driven and the point estimate sat
above the bar, so more seeds would most likely *confirm a small real positive bias*. That reasoning
assumed the true median equalled the then-current point estimate of 0.0242. It did not: the first
50 seeds were the high samples, and 110 fresh ones land at ~0. **The call to spend seeds, taken over
re-specifying the bar and over accepting no verdict, is vindicated by measurement.**

It also vindicates the narrower claim the forecast was attached to: at n=50 the interval already
*contained zero*, so the honest statement was always "too blurry to tell", not "the control is
dirty". The extra seeds resolved the blur in the direction that statement allowed for.

## What passing the gate does and does not buy

**It does NOT favour the pyramid.** Underneath the voided stamp the measured result is that the
learned pyramid **loses** to the 4-mer baseline — paired median −0.0127, 2/10 seeds winning, CI
`[-0.0205, -0.0040]` on the pre-registered range. So passing validity does not rescue a result; it
**permits recording a negative one**. Branch A says the arm's verdict on seeds 0..9 becomes
computable and is **`NO_GAIN`**.

**Scope, stated plainly.** The 160-seed median is a property of the **PIPELINE**, not of any one
run. The committed 10-seed artifact still carries its own internal `null.median = 0.0517` and its
own `INVALID_NULL_SCORED` stamp; nothing here rewrites it. What the pooled estimate supports is that
**that 0.0517 was a sampling fluctuation rather than leakage** — which is exactly the question the
gate exists to answer, and exactly why pooling across batches is the right instrument for it.

**Pooling is valid because the code is identical.** All 160 seeds come from `PINNED_SHA`
`74b7a28df74d46974f36c5f55ed601b9cca3b156`, verified equal to the `head` field of the seeds-10..49
artifact before b1 was launched. A bumped pin would have made this median a mixture of two
estimators.

## Still in flight

b3 (seeds 160..214) and b4 (215..269) take the pool to **270**, plus the attested sampling sweep
(`glm-sampling-sweep-vk`). Carried unattended by `scripts/fragment_null_overnight.py`; resume at any
time with `uv run python scripts/fragment_null_batch_driver.py step`, which is idempotent.

**Expect 270 to tighten this, not to flip it — and expect it to stay a WEAK pass.** With the median
near 0.0112 the interval's upper end only clears 0.02 at roughly **n ≈ 1,190** (~98 h of T4), so the
distinction between "passes on the point estimate" and "passes with the interval clear of the bar"
will not be settled by this campaign. That was the pre-registration's second forecast (§4-P2) and it
still holds, even though its first forecast did not.

## The standing lesson

**A bar placed below its own measurement noise does not detect anything — it samples.** This gate
voided two complete runs on a reading whose interval contained zero throughout, and the reason the
control sits anywhere near ±0.13 per seed is structural and expected: scrambling the training labels
halts training almost immediately, leaving the network near its random initialisation, where it acts
as a random functional of base composition — and composition genuinely predicts this rung at
`spearman(GC, expression) = -0.1169`, the same magnitude as the whole kmer4 comparator. Real leakage
would drive the control toward that 0.11, not to 0.02.

So the bar was mis-placed from the start: set where noise lives rather than where leakage lands. The
seeds have now bought a point estimate inside it, but the placement issue is unchanged and will
recur on the next arm that inherits this bar.
