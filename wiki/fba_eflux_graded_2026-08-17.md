> # RETRACTED 2026-08-17
>
> **The claim below is NOT supported.** The deterministic re-run (`processes=1`, the reproducible mode)
> returns delta **-0.0133** with bootstrap CI **[-0.0352, +0.0087]**, P(d<=0)=0.902 -> `NO_SIGNAL_EVEN_GRADED`.
> Three runs span -0.013 to +0.031: the run-to-run variance exceeds the effect and crosses zero, so no
> direction was ever established. Offline corroboration: outside the lethal spike the baseline is at
> CHANCE (AUROC 0.4963), and on cells that actually moved E-Flux is WORSE (0.4542 -> 0.4313).
>
> See `wiki/fba_eflux_graded_RETRACTION_2026-08-17.md`. Retained unedited below for the audit trail.

# The binary null was a BINARIZATION ARTIFACT — E-Flux does carry graded signal

**Date:** 2026-08-17 · **Verdict:** `EFLUX_CARRIES_GRADED_SIGNAL`
**Pre-registration:** `wiki/fba_eflux_continuous_prereg_2026-08-17.md` (frozen before computing)
**Artifact:** `wiki/fba_eflux_bridge_2026-08-17.json` · **Runner:** `scripts/fba_eflux_bridge.py`

## What was open

Two prior results said E-Flux does nothing, and both concluded *"expression constraints can only help if
the readout changes to something graded"*:

- the binary null (0 of 1,441 cells changed), and
- the mechanism at scale (`a880001`: 0 of 16,676 cells lie between 0.001 and 0.05, so the `FRAC = 0.01`
  threshold sits in an empty region).

**That conclusion had never been tested.** It was an inference from the readout's *shape*, not a
measurement of whether E-Flux carries any essentiality signal at all. This closes the gap.

## Result — same cells, graded instead of thresholded

| arm | binary per-cell | **continuous AUROC** |
|---|---|---|
| A — plain FBA | 0.6967 | 0.6216 |
| B — E-Flux | 0.6967 | **0.6428** |
| **Δ** | **+0.0000** | **+0.0212** |

Pre-registered bar was `Δ ≥ +0.02` = carries signal, `Δ ≤ 0` = no signal even graded. It clears.

**Gene-level paired bootstrap** (B=1000, resampling the 131 **genes** — the independent unit, since the
1,441 cells are 131 genes × 11 conditions and same-gene cells are correlated):

> **95% CI [+0.0019, +0.0397], P(Δ ≤ 0) = 0.015 → `CI_EXCLUDES_ZERO`**

So: **the same knockouts, the same conditions, the same 542 essential cells — E-Flux ranks them better,
while changing not one binary call.** The binary null was an artifact of thresholding, not an absence of
signal.

## What this overturns

The 2026-08-16 memo's "expression-constrained FBA adds nothing to conditional essentiality" is now
**scoped**: it adds nothing to the **binary call** — that part stands and is exactly reproducible — but
it does add a small, statistically-supported amount to the **graded ranking**.

The metric-redesign recommendation therefore has a **measured payoff** rather than being a hopeful
inference. The bridge is not dead; the readout was hiding it.

## Honest limits — the effect is real in DIRECTION, unstable in MAGNITUDE

1. **Run-to-run variance is comparable to the effect.** An earlier run of the identical script returned
   base 0.6193 / eflux 0.6498 → **Δ +0.0305**, versus this run's **+0.0212**. Nothing between the runs
   touched the solves. Both clear the bar and agree in sign, but **the magnitude should be read as
   "roughly +0.02 to +0.03", never as a precise figure.**
2. **Source of that variance, partially characterised.** A within-process re-solve at `processes=1` is
   bit-identical (0 of 131 cells differ), so the solver itself is deterministic; the variance comes from
   the default multiprocessing path. I could **not** test that path directly — it times out on this host
   (the same spawn storm that produced WinError 5 earlier). `processes=1` is now **pinned in the runner**
   so future runs are reproducible; this run predates the pin.
3. **The lift is modest.** AUROC 0.62 → 0.64, and the CI's lower bound is +0.0019 — barely clear of zero.
   This is a real signal, not a strong one.
4. **Inherited caveats stand:** strain seam (PRECISE-1K MG1655 vs Fitness Browser Keio/BW25113); sample
   imbalance (glucose 621, glycerol 111, the other nine 2–8); 11 of 25 conditions.
5. **`oracle_*` is fitted on the evaluation set** — an upper bound, never deployable. Noted because the
   E-Flux arm's oracle MCC (0.3252) is actually *lower* than baseline's (0.3410) even as its AUROC is
   higher: AUROC measures ranking, oracle-MCC measures one fitted operating point. They are not in
   conflict, and neither is the headline.

## Verification performed

- **My `_auroc` was wrong on first write** — inverted, and wrong on ties. Caught by unit tests before any
  result was trusted (one of my own test *expectations* was also wrong: the half-tie case is 0.875, not
  0.75). Rewritten, all 5 cases pass.
- **Cross-validated against the tested implementation:** the bootstrap's `point_delta` (0.0212) matches
  `continuous_readout`'s independently-computed delta (0.0212) exactly.
- **The silent-degeneracy trap was guarded, not assumed:** `continuous_readout` returns
  `{"auroc": None, "note": "degenerate: one class only"}` on a conditions/ratios key mismatch — a message
  that reads as a claim about the DATA rather than about the caller. The runner passes the 11 keys
  explicitly and **raises** if `n_cells` is 0.

## Next

- **Cheap:** re-run at the pinned `processes=1` to get a reproducible point estimate and settle the
  magnitude.
- **The real question this opens:** if a graded readout recovers signal, what is the *best* graded
  readout? `continuous_readout` already reports an oracle threshold; a deployable graded metric
  (calibrated, not fitted on the evaluation set) is now worth designing — that is the acceptance-metric
  decision, still the user's call, but it is no longer speculative.
