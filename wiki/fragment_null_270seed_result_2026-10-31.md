# The 270-seed null resolution — Branch A, a WEAK pass, and the arm closes NO_GAIN

<!-- doc-status: LIVE as of 2026-10-31. This is the measured resolution of the campaign
     pre-registered in wiki/fragment_null_270seed_prereg_2026-10-30.md. Read that file's
     sections 3a and 5 before quoting anything here — the reading was fixed in advance. -->

## The result

```
pooled null over 270 seeds
  median +0.0092   sd 0.0671   bar |median| < 0.02 -> WITHIN
  median CI95 [-0.0133, +0.0211]  (exact order statistics, achieved coverage 0.9556)
  contains zero: YES     contains the 0.02 bar: YES
  pre-registered branch: A (CLEAN NULL -> validity passes, arm closes NO_GAIN)
```

**Branch A fired. The validity gate PASSES**, so the fragment-pyramid arm's pre-registered verdict
on seeds 0..9 becomes computable, and it is **`NO_GAIN`**. The arm closes as a measured negative.

**It is a WEAK pass and the pre-registration requires saying so.** Branch A's own text: *"if `CI`
still straddles 0.02, say so — a passing point test on a straddling interval is a weak pass and must
be labelled one."* It straddles. The point estimate clears the bar; the interval does not exclude it.

## Pool composition

| source | seeds |
|---|---|
| `wiki/glm_fragment_pyramid_arm_kaggle_full_2026-10-28.json` | 10 |
| `wiki/fragment_null_ext_summary_2026-10-29.json` | 40 |
| `fragment_null_b1_summary.json` | 55 |
| `fragment_null_b2_summary.json` | 55 |
| `fragment_null_b3_summary.json` | 55 |
| `fragment_null_b4_summary.json` | 55 |
| **total** | **270** |

Four Kaggle T4 batches, ~4.7 h each, two concurrent — about 9 h of compute, funded by an explicit
user decision on 2026-10-30 to "spend the seeds (~270)" over three alternatives.

## The last 55 seeds moved the answer, and that is why the caveat existed

| pooled n | median | median CI95 | interval clear of the 0.02 bar? |
|---|---|---|---|
| 50 | +0.0242 | `[-0.0168, +0.0493]` | no |
| 215 (b1+b2+b4) | +0.0021 | `[-0.0149, +0.0189]` | **yes, by 0.0011** |
| **270 (all four)** | **+0.0092** | **`[-0.0133, +0.0211]`** | **NO — back across it** |

At 215 the interval cleared the bar by a hair, and the standing caveat was *"it clears by 0.0011, a
hair, and b3's final 55 seeds can still move it — do not quote 'interval clear of the bar' as settled
until 270."* b3 moved the upper bound from +0.0189 to +0.0211. **Had the 215-seed figure been banked
as a clean exclusion, it would now be a retraction.** The caveat was load-bearing, not hedging.

## Forecasts, scored

Both predictions were recorded in §4 before any batch ran.

- **P1 — FALSIFIED.** It gave ~80% that the n=270 point estimate would still be >= 0.02. It came in
  at **+0.0092**, well under. This is the *second* falsification of this same prediction: an earlier
  version was already falsified at 160 seeds, and restating it was itself an error.
- **P2 — its CLAIM holds, its PREDICTED INTERVAL does not.** The claim *"270 seeds buys a point
  estimate, not a confident resolution"* is correct: the interval straddles 0.02. But the forecast
  interval was `[+0.0144, +0.0340]`, entirely above zero, and the measured one is
  `[-0.0133, +0.0211]`, which contains zero. So the sizing model mispredicted the LOCATION while
  getting the resolution claim right.

Honest tally across this project: the two fragment-null forecasts were falsified or half-right; the
G-D fine-tune's `P(bar met) = 0.2` hit. Do not defer to my priors here.

## What this does NOT change

**The pyramid-vs-kmer4 science is untouched and does not depend on this gate.** The paired CI
excludes zero on the negative side at both seed ranges — **−0.0127 `[-0.0205, -0.0040]`** on seeds
0..9 and **−0.0061 `[-0.0133, -0.0037]`** on 10..49. The learned multi-resolution pyramid **loses to
a 4-mer-count ridge baseline** either way.

This gate decided **what the ledger may record**, not what is true. Before it, the arm's verdict was
void under `INVALID_NULL_SCORED`; now the negative is a *valid* negative rather than an unknown.

## Two traps this leaves behind

1. **The arm emits two intervals for two different estimands.** `ci95` is a normal approximation
   around the **MEAN**; `ci95_median` is the exact order-statistic interval for the **MEDIAN**. The
   headline and `clears_margin` are the MEDIAN, so reading the bare `ci95` beside them silently
   reports an interval for a different quantity. The arm is bootstrap-free by design.
2. **Do not cite `[-0.0126, +0.0466]`.** It appears in an earlier handoff labelled "bootstrap" and is
   **not reproducible by any in-repo method** (see prereg §3a). Use `median_ci95`, which reports its
   ACHIEVED coverage rather than claiming 95%.

## Reproduce

```
uv run python scripts/fragment_null_batch_driver.py status
uv run python scripts/fragment_null_batch_driver.py pool
```

Read `wiki/fragment_null_270seed_prereg_2026-10-30.md` §§3a and 5 first — the driver prints that
instruction for a reason.
