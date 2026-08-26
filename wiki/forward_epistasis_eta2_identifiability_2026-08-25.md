# The H2 memo credits the wrong η² — and the right one is not identifiable here (2026-08-25)

**Amends `wiki/forward_epistasis_h2_confirmed_2026-08-25.md`, which over-claimed.** That memo said H2 was
CONFIRMED: that η²(k) — the share of **fitness** variance between mutation orders — governs the additive
score's pooling gain. It correlated the gain against η² computed on **fitness only**.

A pooling gain needs **both** the label and the predictor to separate by group. Had the additive score been
distributed identically across orders, pooling groups that differ in fitness *mean* would add no rank
correlation at all — the gain would stay zero however large η²(fitness) grew. So η²(fitness) is NECESSARY,
not obviously SUFFICIENT, and the score-side term was never measured.

Measuring it changes the conclusion. No new compute: cached ESM2 log-prob matrices + cached assay CSVs.

## What the score side does

| protein | subsets | ρ(gain, η²_**fitness**) | ρ(gain, η²_**score**) | collinearity | partial fit \| score | partial score \| fit |
|---|---:|---:|---:|---:|---:|---:|
| GFP | 26 | +0.973 | +0.991 | +0.982 | **+0.020** | **+0.815** |
| HIS7 | 26 | +0.995 | +0.993 | +0.999 | **+0.559** | **−0.135** |
| ParD | 4 | +0.600 | +1.000 | +0.600 | n/a (degenerate) | +1.000 |

**The two eta² are near-collinear (+0.982 / +0.999), and the two well-powered proteins point in OPPOSITE
directions when partialled.** GFP attributes the gain to the score side (+0.815) with the fitness side
collapsing to nothing (+0.020); HIS7 does the reverse (fitness +0.559, score **−0.135**). That is the
signature of two variables that cannot be separated on this data — not of a winner.

So the honest claim is the **JOINT** one:

> The pooling gain is governed by **aligned between-order separation in BOTH the label and the predictor**.
> Which side carries it is **not separably identifiable** in these assays.

**Swapping "η² of the score governs it" in for "η² of fitness governs it" would just be the next overclaim** —
HIS7's −0.135 refutes the general version of exactly that statement, which an earlier reading of GFP+ParD
alone had suggested.

## What survives, and it is not nothing

**Leave-one-order-out — the honest resampling unit.** The 26 subsets are nested and share most of their
variants, so a p-value over them assumes an independence they do not have; the effective sample size is the
number of **orders** (5), not 26. Dropping each order in turn:

| protein | ρ(gain, η²_fitness) range | ρ(gain, η²_score) range | drops |
|---|---|---|---:|
| GFP | +0.955 .. +0.982 | +0.964 .. +1.000 | 5 |
| HIS7 | +0.982 .. +1.000 | +0.964 .. +1.000 | 5 |

The relationship between the gain and between-order structure is **robust to order-level resampling**. Only
the *attribution* between the two sides fails. The H2 memo's underlying finding stands; its causal wording
does not.

## The mechanism, now grounded rather than asserted

The additive score is a **sum of k per-mutation log-ratios**, so its mean scales with k *by construction*.
Measured directly:

| protein | mean additive score by order | slope |
|---|---|---|
| GFP | k2 −0.329 · k3 −0.497 · k4 −0.674 · k5 −0.854 · k6 −1.002 | ≈ −0.17 / order |
| HIS7 | k4 −6.618 · k5 −7.960 · k6 −9.288 · k7 −10.643 · k8 −11.895 | ≈ −1.32 / order |
| ParD | k2 −5.751 · k3 −8.702 · k4 −11.853 | ≈ −3.05 / order |

Linear in k on all three. **This is why the additive score separates by order at all** — and it explains the
asymmetry the correction memo flagged but could not account for: the additive score's pooling gain swings 4×
across proteins while the **joint** score's stays near-constant (+0.054 / +0.058 / +0.065). A joint score is
one forward pass on the multi-mutant; it carries no such summation, so it has no structural reason to
separate by k.

It also identifies what H1 got wrong. H1 measured the **fitness** slope vs k and was falsified (GFP has the
steeper fitness slope and the smaller gain). The quantity that mattered was the **score's** separation by k.

**LIMITATION, stated:** the joint score's own η² is **not computable here**. The 2026-07-27 sweep persisted
only aggregate per-order Spearman values, not per-variant joint scores. The mechanism above is *grounded* on
the additive side and *inferred* on the joint side.

## Corrections this forces on the H2 memo

| claim | status |
|---|---|
| "H2 CONFIRMED" | → **refined and partly superseded**; the governing quantity is the joint condition, not fitness η² |
| p = 7.4e-17 / 9.3e-26 | → **not inferential** — 26 nested subsets sharing variants; n_eff ≈ 5 orders. Descriptive only |
| "made it causal" / "entirely" / "100% between-order" | → **overstated**. Control A re-deals the same group sizes from one shuffled pool, so every group becomes an i.i.d. draw and gain → 0 *by construction*. That checks the arithmetic, not a mechanism |
| Control B "leaves it intact / untouched" | → **wrong on 2 of 3**: HIS7 +0.043→+0.052, ParD +0.148→+0.185. Direction is explainable (B zeroes within-ρ, and gain = pooled − within, so removing a positive diluting term raises it) — but it needed explaining, not asserting past |
| "η² is the cheap way to say by how much" | → **overclaims the AMR analogy**. `dna_decode/eval/clonality.py` computes a corrected estimator with uncertainty (one lineage vote, `wilson_ci`, `effective_lineage_n`, DISCORDANT clusters excluded not majority-voted). η² is a scalar exposure flag with no CI and no corrected estimate |
| deployment rule "report η² beside a pooled figure" | → **under-specified**. The check is whether **both** label and predictor separate by the grouping variable |

## Deployment rule (revised)

- Report the **within-group** metric whenever a sweep spans a grouping variable. Unchanged, and it is the
  right default for the deployed question — `predict_multi_effect` scores a user's specific mutation set,
  which is a fixed-burden question.
- State the **estimand first**. Within-group is right for *"which scorer ranks at fixed mutation burden?"*
  It is **not** right for *"can this score rank a realistic mixed-order candidate pool?"*, where cross-order
  ranking is part of the target rather than confounding.
- If a pooled figure is reported anyway, disclose **both** η²(label) and η²(predictor) — a pooled metric is
  inflated only where the two are aligned, and one alone cannot tell you.

## Reproduce

```bash
uv run python scripts/epistasis_eta2_identifiability.py   # no GPU, no network
```
Artifact: `wiki/forward_epistasis_eta2_identifiability_2026-08-25.json`. Frozen AMR surface untouched.
