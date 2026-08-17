# E-Flux under a GRADED readout — PRE-REGISTRATION (frozen before computing)

**Date:** 2026-08-17 · **Status:** frozen before execution

## Why this test exists

Two results are now in:

1. **The binary null** (`fba_eflux_bridge_2026-08-16`): E-Flux changed 0 of 1,441 scored cells.
2. **The mechanism, at scale** (`a880001`): the growth-ratio distribution is bimodal — **0 of 16,676
   gene×condition cells** lie between 0.001 and 0.05, so the `FRAC = 0.01` threshold sits in an empty
   region and flux rescaling has nothing to flip.

Both artifacts then conclude: *"expression constraints can only help if the readout changes to something
graded."* **That conclusion has never actually been tested.** It is an inference from the shape of the
readout, not a measurement of whether E-Flux carries any essentiality signal at all.

This closes that gap. It is the decisive test of whether the null is a **binarization artifact** or a
**genuine absence of signal**.

## The question

Does E-Flux improve the **continuous** readout — ranking gene×condition cells by raw growth ratio —
even though it provably changes no binary call?

The instrument already exists and is tested: `conditional_essentiality.continuous_readout(records,
ratios, conditions)` → AUROC over the two-sided subset, plus deployed-vs-oracle confusion. It was built
for exactly this question ("is the conditional signal ABSENT, or is the binary CUTOFF discarding it?").

## Arms

Identical to the bridge run — same 11 PRECISE-1K conditions, same iML1515, same Fitness Browser labels
(orgId=Keio, `fit < -2.0`), same conditionally-essential subset. The only change is that the runner now
retains the **raw growth ratio** per cell instead of discarding it after thresholding.

| arm | readout |
|---|---|
| A — plain FBA | continuous ratio → AUROC |
| B — E-Flux | continuous ratio → AUROC |

## Frozen decision rule

- **Primary metric:** `auroc` from `continuous_readout`, computed per arm on the same cells.
- **`EFLUX_CARRIES_GRADED_SIGNAL`:** `AUROC(B) − AUROC(A) ≥ +0.02`. The binary null is then a
  **binarization artifact**, the metric-redesign recommendation has a measured payoff, and expression
  constraints are live again under a graded metric.
- **`NO_SIGNAL_EVEN_GRADED`:** `Δ ≤ 0`. Expression constraints carry **no** essentiality signal at this
  granularity, which **closes the graded-readout escape hatch** and retires the metric-redesign
  recommendation for this purpose. This would be the stronger, more useful result: it converts an open
  "maybe if we change the metric" into a closed negative.
- **Ambiguous:** `0 < Δ < +0.02`.
- **Committed:** reported in whichever direction it lands. Secondary readouts (oracle threshold, MCC,
  per-condition) are descriptive and cannot be promoted to the headline post hoc.

## Traps to avoid (recorded before running, from this repo's own history)

1. **The `conditions` argument is load-bearing.** `continuous_readout` guards with `ratios.get(c, {})`,
   so a conditions/ratios key mismatch accumulates zero cells and returns
   `{"auroc": None, "note": "degenerate: one class only"}` — which reads as a statement about the DATA,
   not about the caller. The 11 condition keys are passed explicitly and `n_cells` is asserted > 0.
2. **`oracle_*` is fitted on the evaluation set** — an upper bound, never quoted as deployable.
3. **NaN growth → ratio 0.0**, consistent with the shipped genuine-essentiality coding (an infeasible
   solve is real lethality, not a solver failure).
4. **Both arms must be scored on the same cell set**, or the AUROCs are not commensurable.
