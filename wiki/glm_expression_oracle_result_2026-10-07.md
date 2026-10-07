# The expression oracle works: Spearman 0.59 on held-out promoter elements

**Date:** 2026-10-07 · **Artifact:** `wiki/glm_expression_oracle_2026-10-07.json`
**Reproduce:** `uv run python scripts/glm_expression_oracle.py` (CPU, seconds, no network after first run)
**Substrate:** GEO **GSE108535** — Urtecho G, Tripp AD, Insigne KD, Kim H, Kosuri S, *Biochemistry*
2019;**58**(11):1539–1551, DOI `10.1021/acs.biochem.7b01069`. Genomically-integrated MPRA, fully-crossed
σ70 promoter grid in *E. coli*.

**This is the scoring head the generative loop was missing** — Family 3, blocked for weeks on "is there a
free measured bacterial regulatory→expression dataset?". There is, and it is 227 KB.

## The substrate, verified by download rather than from a summary

**10,898 sequence→expression pairs.** Every sequence exactly **150 bp**, ACGT-only, **10,898 distinct**
sequences, expression spanning **615×** (min 0.0569, median 0.1970, max 34.98). The table also ships the
**element identities** — `UP_element` (3 values), `Minus35` (8), `Spacer` (8), `Minus10` (8),
`Background` (8) — which is what makes the honest split below possible.

## The result

Ridge regression per feature set. **The leave-element-out column is the headline**; the random column is
reported only so the inflation is visible.

| feature | n features | random split | **leave-element-out (mean)** | inflation |
|---|---|---|---|---|
| GC content only | 1 | +0.1812 | **+0.1414** | +0.0398 |
| 4-mer | 256 | +0.7010 | **+0.4137** | +0.2873 |
| 6-mer | 4,096 | +0.6986 | **+0.5872** | +0.1114 |
| **positional one-hot** | 600 | +0.6953 | **+0.5914** | **+0.1039** |

**A working supervised sequence→expression predictor for *E. coli* promoters at ρ ≈ 0.59 on held-out
elements.** For scale, AlphaGenome reports ρ ≈ 0.5 on human expression effect sizes — so this is a credible
number for the task class, on a designed grid, from a 227 KB download and a ridge regression.

## Finding 1 — the crude baseline does NOT win here, and that is worth saying

This repo's standing lesson from the MLDE literature is that dumb baselines win more often than not: a single
integer recovers 75% of a protein LM's headline, a 1992 substitution matrix ties ESM2, plain ridge ties ESM
on TAPE fluorescence. **That pattern does not reproduce on this task.** GC content alone reaches **0.1414**
while positional one-hot reaches **0.5914** — a 4× gap. The sequence features are doing real work.

The honest reading: the lesson was never "baselines always win", it was "always *run* the baseline". Run
here, it loses, and that is informative rather than disappointing.

## Finding 2 — random splits inflate, measured a FOURTH independent time

Inflation is positive for **every** feature set: +0.04 / +0.29 / +0.11 / +0.10. This repo had already
measured ~3× inflation three separate times (GB1 1-vs-rest 0.32 → 0.92 random; one-hot 0.579 → 0.027 across
positional splits; R² 0.90 → 0.23 stratified). This is a fourth instance in a fourth setting, expressed in
Spearman rather than as a ratio.

**And 4-mers inflate worst (+0.2873).** A 4-mer histogram is the feature most able to memorise which grid
cells it has seen, so it gains most from a random split and loses most when whole elements are withheld.

## Finding 3 — per-axis difficulty, and it is not a held-out-size artifact

| held-out axis | GC | 4-mer | 6-mer | one-hot |
|---|---|---|---|---|
| UP_element (1 of 3) | 0.0387 | 0.6759 | 0.7044 | 0.6843 |
| **Minus35** (2 of 8) | 0.1514 | **0.3773** | **0.3882** | **0.3919** |
| Spacer (2 of 8) | 0.1352 | 0.5179 | 0.6095 | 0.6213 |
| Minus10 (2 of 8) | 0.1736 | 0.2978 | 0.5884 | 0.5963 |
| Background (2 of 8) | 0.2080 | 0.1994 | 0.6457 | 0.6632 |

**Generalising to unseen −35 elements is by far the hardest axis** (0.39 vs 0.62–0.70 elsewhere). This is
**not** an artifact of how much was withheld: Minus35, Spacer, Minus10 and Background all hold out **2 of 8**
values, so the comparison is matched. UP_element is the easy one partly because only 1 of 3 is withheld.

**And 4-mers collapse where position matters.** On `Background` the 4-mer falls to **0.1994** — barely above
GC — while 6-mer and one-hot hold at 0.65–0.66. Same on `Minus10` (0.2978 vs 0.59). So local composition
alone is insufficient for those axes, and **positional/longer-context features are what carry it** — the
same conclusion the generator falsifier reached independently, where adding positional features was what
broke a k-mer ceiling effect.

## Two traps encoded in code, because both are silent

1. **The join loses 65% of the data.** The expression table writes UP-element names with **hyphens**
   (`gourse-326fold-up`), the barcode map with **underscores** (`gourse_136fold_up`). A naive join on `name`
   returns **3,821** pairs and looks perfectly healthy — the model trains, the metric computes, nothing
   complains. Normalising gives all **10,898**. `load_pairs` now **enforces a ≥0.95 join-rate floor** and
   reports what the naive join would have given, so this cannot go quiet again.
   **A refinement found by its own test:** a blanket `replace("-","_")` also rewrites the mutation arrow
   `35T->A` → `35T_>A`. That still joins (it is applied to both sides) but is lossy and could in principle
   collide two distinct variants, so the arrow is now protected. Verified to still give a 100% join.
2. **The files are SPACE-delimited, not tab.** A `split('\t')` parses **zero** rows — and a zero-row parse
   that is then averaged reads as "no signal" rather than "no data".

## Honest limits

- **One organism, one assay, one designed grid.** A fully-crossed grid is not a natural promoter
  distribution, so this bounds *grid* generalisation, not genome-wide performance. The genome-wide substrate
  (GSE144621: 321,123 fragments in LB *and* M9) is downloaded-ready and untested here.
- **Ridge on fixed features is deliberately the trivial-head half** of AlphaGenome's principle — no learned
  representation is trained. **These numbers are a FLOOR, not a ceiling.**
- **No confidence intervals** on the Spearman differences.
- Spearman is rank-based and so invariant to the log transform; the *fit* is not, which is why `log_target`
  defaults on (expression spans 615×).
- Held out by **element identity**, which is stronger than random but is not a phylogenetic or
  sequence-distance split.

## What this unblocks

The generative loop had a generator and a *discriminator* but no **function predictor**. It now has one, on
real measured data, with a defensible split. That is the component AlphaGenome's transferable logic said was
required, and the one whose absence made the earlier falsifier a naturalness test rather than a fitness test.
