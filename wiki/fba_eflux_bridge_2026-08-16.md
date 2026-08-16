# The bridge experiment ran — and expression constraints changed *nothing*

**Date:** 2026-08-16 · **Verdict:** `MACHADO_PRIOR_CONFIRMED_ON_ESSENTIALITY`
**Pre-registration:** `wiki/fba_eflux_bridge_prereg_2026-08-16.md` (frozen before the first solve)
**Artifact:** `wiki/fba_eflux_bridge_2026-08-16.json` · **Runner:** `scripts/fba_eflux_bridge.py`

## Result

| arm | per-cell agreement | exact-set match |
|---|---|---|
| A — plain FBA | **0.6967** | 8/131 |
| B — E-Flux (expression-constrained) | **0.6967** | 8/131 |
| **Δ** | **+0.0000** | 0 |

Pre-registered bar was `Δ ≥ +0.02` for GO, `Δ ≤ 0` for the Machado prior. It landed on **exactly zero**.

**Zero of 11 conditions had even one confusion cell change.** Not a small effect — a null so complete
that every TP, FP, FN and TN is identical across all 1,441 scored cells.

## This is not a no-op run

E-Flux demonstrably bit:

- **2,210 reactions constrained** per condition (56 skipped as unmeasured, left unconstrained by design).
- **Wildtype growth fell substantially** in every condition — D-Fructose 0.877 → 0.643, D-Galactose
  0.868 → 0.534, N-Acetyl-D-Glucosamine 1.131 → 0.754.

So the expression layer changed the flux space materially, and the essentiality **classification** did
not move at all.

## Why — measured, and it killed my first explanation

My initial hypothesis was that E-Flux's multiplicative scaling **preserves** knockout/wildtype ratios,
so a ratio test cancels it. **That is falsified.** The per-gene ratio diagnostic (40 genes × 2
conditions, `processes=1`) shows ratios move *substantially*:

| condition | max \|Δratio\| | mean \|Δratio\| | genes with \|Δ\| > 0.01 |
|---|---|---|---|
| D-Galactose | 0.2710 | 0.054 | 10/40 |
| Potassium acetate | 0.4337 | 0.059 | 9/40 |

The actual mechanism is different and more useful: **the growth-ratio distribution is bimodal** —
genes sit at ≈1.0 (dispensable) or ≈0.0 (lethal), with almost nothing in between. Sample ratios move
`0.978 → 0.881`, `1.0 → 1.0`, `−0.0 → −0.0`. So E-Flux shifts magnitudes by up to 0.43 while
essentially never moving a gene across the **0.01** line, because that line sits in a nearly empty
region of the distribution.

**The classification is robust to large flux-magnitude changes because the readout is a threshold on a
bimodal quantity.** That is a property of the *readout*, not of E-Flux — and it predicts the null
generalizes to any constraint layer that rescales flux without flipping a knockout between
lethal and viable.

*(Recorded honestly: this replaces a hypothesis I had labelled `unfalsified` in the first version of
this memo. Measurement contradicted the expectation; the expectation was wrong.)*

## What this does and does not close

**Does:** on the 11 conditions where expression data actually exists today, expression-constrained FBA
adds nothing to conditional essentiality. It extends Machado & Herrgård 2014 (plain FBA + parsimony as
good or better, for *flux*) to *essentiality* — the exact prior the 5th advance cited for its original
NO-GO, and which the 6th advance's CONDITIONAL GO had to argue around.

**Does not:** close the bridge on the full 25-condition panel, and does not test NP881's 16-condition
panel (see below). But the measured mechanism above — a threshold on a bimodal readout — now
*predicts* the null holds on any panel, because it is a property of the essentiality readout rather
than of this substrate. That prediction is testable cheaply on the full 25 panel (the baseline
deletions already exist); it does not require NP881.

**The actionable consequence:** if expression constraints are to help at all, the readout has to
change — a continuous or graded essentiality measure rather than a 1% threshold on a bimodal ratio.
That is a redesign of the metric, not more constraint engineering.

## Correction to this morning's artifact

`wiki/fba_bridge_precise_overlap_2026-08-16.json` concluded the blocker was "fully a CODE wall (~1 day
of E-Flux work)." **That was wrong.** The 346 new NP881 profiles are published only as **raw SRA
reads** — absent from `SBRG/precise1k` (only `master`, no NP881 directory), absent from GEO (0 hits for
either BioProject), no Zenodo deposition. Using NP881 requires downloading and quantifying ~32–82
RNA-seq runs (~50–250 GB plus an alignment pipeline).

The **16/25 substrate overlap stands** — it came from SRA run metadata. What was wrong was the
availability/cost inference stacked on top of it. This run therefore used the 11 conditions covered by
PRECISE-1K, whose processed matrix *is* public (60.9 MB).

## Named caveats

1. **Strain seam.** PRECISE-1K is K-12 MG1655; the labels are Fitness Browser orgId=Keio (BW25113).
2. **Sample imbalance.** Glucose 621 samples, glycerol 111, the other nine 2–8. Nine of eleven
   conditions rest on a mean of 2–8 profiles.
3. **11 < 25**, and the 131 conditionally-essential genes here are a different subset from the
   217 on the full panel.
4. **The 0.7368 headline is not the comparator.** 0.6967 is the *baseline arm recomputed on these 11
   conditions* — comparing against the 25-condition figure would have repeated the regulatory-arm
   verdict bug. The two numbers are not commensurable and neither is a regression of the other.

## Next

The authority fork is unchanged in shape but better informed: the bridge now has a **measured null on
the available panel**, so pursuing NP881 means paying ~50–250 GB of SRA processing to re-test a
hypothesis that just returned exactly zero on 11 conditions. The cheap move is instead to finish the
ratio diagnostic and decide whether the null is mechanistic (predicting no panel will help) or
panel-specific.
