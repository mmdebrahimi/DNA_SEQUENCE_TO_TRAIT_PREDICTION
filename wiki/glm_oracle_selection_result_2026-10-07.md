# The oracle converts its correlation into a CHOICE: 79% of the available selection value at k=10

**Date:** 2026-10-07 · **Artifact:** `wiki/glm_oracle_selection_2026-10-07.json`
**Reproduce:** `uv run python scripts/glm_oracle_selection.py` (CPU, ~67 s, no network after the first run)
**Substrate:** GEO **GSE108535** σ70 promoter MPRA, *E. coli*, 10,898 measured promoters, leave-element-out.
**Verdict:** `ORACLE_READY_AS_SCORING_HEAD` — **5/5 axes pass both pre-registered bars.**

## Why this and not another correlation

The oracle reports ρ ≈ 0.59 on held-out elements. **The generative loop does not consume a correlation — it
consumes a choice:** propose k candidates, score them, keep the best. Whether ρ 0.59 converts into a useful
choice is *not derivable from ρ*; it depends on the tail shape, on how ties fall, and on how much selection
value even exists in a pool of k. A correlation of 0.59 is compatible with capturing nearly all the
available headroom or almost none of it.

This is the move AlphaGenome's own authors make with self-distillation: the distillation metric proves
nothing, so the claim is validated **downstream** on real measured benchmarks. Here the downstream task is
selection, and the measurement is against real measured expression on promoters built from element variants
the oracle never saw.

## The result

`headroom_captured = (picked − random) / (ceiling − random)`, where all three arms score the **same** k
candidates. The ceiling is the true best of those k by measured expression — the most any scorer could
achieve — so the number is the **fraction of achievable selection value actually realised**.

| axis | ρ (one-hot) | k=2 | k=5 | **k=10** | k=20 | GC @ k=10 | null stdev @ k=10 |
|---|---|---|---|---|---|---|---|
| UP_element | +0.6843 | 0.689 | 0.791 | **0.851** | 0.877 | 0.052 | 0.053 |
| Minus35 | +0.3919 | 0.393 | 0.564 | **0.700** | 0.819 | 0.117 | 0.043 |
| Spacer | +0.6213 | 0.622 | 0.752 | **0.844** | 0.894 | 0.164 | 0.056 |
| Minus10 | +0.5963 | 0.604 | 0.691 | **0.745** | 0.795 | 0.160 | 0.053 |
| Background | +0.6632 | 0.656 | 0.783 | **0.828** | 0.877 | 0.171 | 0.058 |
| **mean** | | 0.593 | 0.716 | **0.794** | 0.853 | **0.133** | |

**The oracle's pick sits at the 78th–85th percentile of its candidate pool** (random ≈ 0.50, ceiling ≈ 0.91)
and **beats a random pick in 74–81% of trials.**

## Four readings, in order of how much they change what we do

### 1. The learned features earn their keep IN THE LOOP, not only on correlation

GC captures **0.133** of the headroom against the learned **0.794** — a 6× gap. This repo's standing lesson
is that crude baselines win more often than not (a single integer recovers 75% of a protein LM's headline; a
1992 substitution matrix ties ESM2). **It does not reproduce here, on either metric.** And on `UP_element`
GC's 0.052 does not clear its own permutation null's p95 (0.055) — on that axis GC is **indistinguishable
from a random scorer**, not weakly useful.

### 2. Propose-more-score-more genuinely pays — the same shape the inverse cell measured

Headroom rises monotonically with k on every axis (mean 0.593 → 0.716 → 0.794 → 0.853 for k = 2, 5, 10, 20).
The inverse cell's deployability rule was *propose-k / assay-k / keep-best*, with top-1 about 4× worse than
best-of-5. The loop's scoring half shows the same gradient, independently.

### 3. A weak correlation still converts — which is the thing ρ could not tell us

`Minus35` is the hardest axis on both metrics, and at **ρ = 0.392 it still captures 0.700** of the headroom.
Normal theory predicts 0.574 there. So selection value is not a simple read-off from correlation, and the
axis that looks worst by ρ is far from useless in the loop.

### 4. Correlation does NOT reliably order two scorers for loop use

`kmer6` wins mean headroom (**0.806** vs one-hot's 0.794) while **losing** Spearman on 4 of 5 axes. On
`Minus10` and `Background` it loses on ρ and wins in the loop. **Read this as "the ordering is unreliable",
NOT as "kmer6 is better"** — the 0.012 mean gap sits inside the ~0.05 null spread, so the two are
effectively tied as selectors. Same shape as the inverse cell's finding that utility does not track forward
rank, which required a per-protein check.

## What makes the number trustworthy

- **k=1 is degenerate by construction and asserted.** With one candidate the pick *is* the random draw *is*
  the best, so any lift at k=1 localises a bug to the harness. `assert_k1_is_degenerate` runs every time.
- **Permutation null per (axis, k).** Measured means at k=10: −0.014 / +0.002 / +0.004 / −0.001 / −0.002,
  stdev 0.043–0.058. The headline 0.70–0.85 sits 12–16 null-stdevs out.
- **A single headroom reading carries scorer-realisation noise that `trials` cannot reduce.** Measured on
  synthetic data: 60 independent useless scorers gave mean headroom **0.0006** with **stdev 0.096**, 31 of 60
  negative. This was found by chasing a +0.133 reading on a useless scorer that I first mistook for a
  harness bias — it is the finite-sample correlation of one scorer realisation. The null band exists because
  of it, and **bar B's margin was made noise-aware (`max(0.10, 2σ)`) at instrument-build time, before any
  real number** — a flat 0.10 is the same order as that spread and could have been cleared by noise.
- **Ties break randomly, not by index order.** `random.sample` returns its subset in random order, so a
  heavily-tied scorer (GC over 150 bp takes few distinct values) cannot be credited with whatever the index
  order gave it. Pinned by a test asserting an all-tied scorer captures ≈ 0.
- **Magnitudes are in normal-theory range:** predicted 0.574–0.879 against measured 0.700–0.851.
- **A partial run cannot emit a headline.** The first smoke printed `SELECTION_IS_THE_BOTTLENECK` on one axis
  against a bar needing three — a confident substantive verdict produced by the run's *scope*. It now returns
  `INDETERMINATE_INSUFFICIENT_AXES` when the bar is unreachable by construction.

## Honest limits

- **UPPER BOUND on the loop's realisable value, and this is the severe one.** The candidate pool here is
  **measured grid sequences**. The real loop's candidates come from a *generator* and will sit off this
  distribution, where the oracle is less reliable and where no measured label exists to check it. This
  bounds the scoring head on its home turf; **it does not simulate the loop.**
- Candidates within one held-out element set share that element, so they vary only on the other axes — a
  narrower pool than the full grid (per-axis measured spread is in the artifact: 119×–531×).
- **`fold_over_random` is unstable and should not be quoted** (2.1× to 39.6× across axes): it is a ratio of
  medians whose denominator can sit near the assay floor, so it tracks pool skew rather than scorer quality.
  Percentile and headroom are the comparable metrics.
- Both pre-registered thresholds are **asserted**, not derived; the ceiling normaliser *is* derived.
- A fully-crossed designed grid is not a natural promoter distribution.
- Ridge on fixed features is the **trivial-head half** of the design principle — a floor, not a ceiling.

## What it unblocks, and what it does not

The loop's **scoring** half is adequate: it converts its correlation into a choice worth ~79% of what a
perfect scorer could achieve, and the learned features are doing the work. The loop's **generating** half is
the measured weak point — GENERator-1.2B at default sampling is *more* distinguishable from natural
promoters than a 3-mer Markov chain. **The bottleneck is the generator, not the scorer.** Nothing here
licenses a claim about off-distribution behaviour, which is exactly the gap self-distillation was invented
for — see `wiki/alphagenome_lessons_actionable_2026-10-07.md` for why that build must wait.
