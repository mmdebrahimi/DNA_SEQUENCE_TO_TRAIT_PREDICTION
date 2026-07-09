# ESM2-650M on the full ProteinGym v1.1 substitution benchmark

**Date:** 2026-07-09 · **Kaggle GPU (Tesla T4, fp16)** · `scripts/kaggle_proteingym_sweep.py`
**Kernel:** `emanueleebrahimi/proteingym-esm2-650m-full-sweep` · **217 / 217 assays scored, 0 errors**

This replaces the n=7 humsavar anecdote with a field-comparable world-model quality number.

## Result

| | our run | ProteinGym official (per-assay) | delta |
|---|---|---|---|
| **median Spearman** | **0.490** | 0.484 | **+0.006** |
| mean Spearman | 0.429 | 0.438 | −0.009 |

**The harness reproduces the published ESM2-650M zero-shot result on the full benchmark.** Both
aggregations land within ~0.01 of ProteinGym's own per-assay numbers, with every one of the 217
substitution assays scored and no failures. Zero-shot masked-marginal, single substitutions only,
sliding window for proteins over 1022 aa.

Residual differences that explain the ~0.01 (all expected, none alarming): we skip multi-mutant rows,
we require ≥10 scorable variants per assay (no assay actually tripped this), and we use our own
sliding-window and Spearman implementations rather than ProteinGym's scorer.

## Scale regresses past 650M

Recomputed from ProteinGym's own per-assay CSV (vendored at `wiki/refs/`):

| ESM2 | per-assay median | per-assay mean | published `Average_Spearman` |
|---|---|---|---|
| 150M | 0.451 | 0.401 | 0.387 |
| **650M** | **0.484** | **0.438** | **0.414** |
| 3B | 0.467 | 0.432 | 0.406 |
| 15B | 0.438 | 0.425 | 0.400 |

**650M is the peak of the ESM2 family.** 3B and 15B are monotonically worse on all three
aggregations. This independently corroborates the humsavar finding
(`wiki/esm_am_ensemble_paired_2026-07-09.md`), where 650M → 3B regressed on 4/7 proteins with a mean
paired delta of −0.021.

The script previously cited "650M ~0.47 | 3B ~0.48 | 15B ~0.48-0.49 (curve flattens after 650M)".
That was wrong in direction — corrected in `91027f9`.

## The aggregation trap

The **published `Average_Spearman` is not a per-assay mean.** For 650M it reads 0.414, while the plain
per-assay mean of the same 217 values is 0.438. The published figure is the mean of the five
function-category averages (Activity 0.425 / Binding 0.337 / Expression 0.415 / OrganismalFitness 0.368
/ Stability 0.523), which reweights the benchmark. This script prints a plain per-assay median and
mean, so it must be compared against the median/mean columns above — never against 0.414.

Note also that ProteinGym's headline leaderboard metric is the **average**, not the median. Our script
reports the median as its primary number, so 0.484 is the right target, not 0.414.

## What this settles

The two "quality levers beyond raw model size" proposed in `2378217` are now both measured, and both
are dead ends on this substrate:

1. **Bigger model** — 3B is worse than 650M, on both humsavar (n=7, paired) and ProteinGym (n=217).
2. **Ensembling with AlphaMissense** — no paired lift at either scale (see the companion memo).

650M is the checkpoint to build on. Further gains have to come from method, not scale or naive
score-averaging.

## Provenance

- `wiki/proteingym_esm2_650m_full_2026-07-09.json` — machine-readable.
- `wiki/refs/proteingym_v1.3_DMS_level_Spearman.csv` — vendored official per-assay Spearman table.
- Future runs of the script now also emit `per_assay_spearman.csv`, enabling an assay-by-assay paired
  comparison against the official table instead of median-to-median only.
