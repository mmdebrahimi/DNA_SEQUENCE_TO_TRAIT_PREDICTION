# RETRACTION — "E-Flux carries graded signal" is NOT supported

**Date:** 2026-08-17 · **Retracts:** `wiki/fba_eflux_graded_2026-08-17.md` (commit `095a986`, pushed)
**Verdict now:** `NO_SIGNAL_EVEN_GRADED` · **Artifact:** `wiki/fba_eflux_bridge_repro_2026-08-17.json`

## What I claimed, and why it was wrong

Commit `095a986` reported that the binary null was a *binarization artifact* — that under a graded
readout E-Flux improved continuous AUROC by **+0.0212**, with a gene-level bootstrap 95% CI of
**[+0.0019, +0.0397]**, P(Δ≤0) = 0.015, and concluded the metric-redesign recommendation had "a measured
payoff."

**The deterministic re-run reverses the sign.**

| run | mode | Δ AUROC | gene-level bootstrap |
|---|---|---|---|
| 1 | default multiprocessing | **+0.0305** | not run |
| 2 | default multiprocessing | **+0.0212** | CI [+0.0019, +0.0397], P(Δ≤0)=0.015 |
| **3** | **`processes=1`, deterministic** | **−0.0133** | CI [−0.0352, +0.0087], P(Δ≤0)=**0.902** |

Run 3 is the reproducible one — `processes=1` was pinned precisely because a within-process re-solve is
bit-identical while the parallel path is not. It returns `NO_SIGNAL_EVEN_GRADED`, the pre-registration's
other branch.

**The run-to-run variance is larger than the effect and spans zero.** I flagged that variance last turn
and still framed the *direction* as established. That was the error: with a spread of −0.013 to +0.031
across three runs, no direction was ever established. The bootstrap CI in run 2 was real but it measures
only gene-resampling variance — it is blind to the between-run variance that actually dominates here.

## Independent corroboration that there is no signal

Offline analysis of the persisted ratios (`fba_eflux_bridge_2026-08-17_ratios.json`), which does not
depend on which run you pick:

- **Outside the lethal spike, the model is at chance.** Restricted to the 1,182 non-lethal cells,
  baseline AUROC is **0.4963** — indistinguishable from 0.5. Purity there is 0.304 against a 0.376 base
  rate, i.e. *below* base rate.
- **On the cells that actually moved** (311 of 1,441, 21.6 %), E-Flux is **worse**: AUROC 0.4542 →
  0.4313. Both below chance.
- **The lethal spike is what carries all the discrimination** (0.70 purity, 183/261 essential), and
  E-Flux barely touches it — 9 cells change lethal status, and excluding them changes the full-cohort
  delta not at all (+0.0212 → +0.0213), so the apparent gain was never coming from there either.

So the coherent picture is: **the model's essentiality discrimination lives entirely in the lethal spike,
which expression constraints do not move, and everywhere else it is at chance.** There is no graded
signal for E-Flux to improve.

## What now stands

- **The binary null stands and is exactly reproducible** — 0 of 1,441 cells change, in every run.
- **The bimodality mechanism stands** (`a880001`: 0 of 16,676 cells between 0.001 and 0.05).
- **The graded escape hatch is now CLOSED, not open.** The pre-registration named `NO_SIGNAL_EVEN_GRADED`
  as "the stronger, more useful result: it converts an open 'maybe if we change the metric' into a closed
  negative." That is what landed.
- **The metric-redesign recommendation is retired for this purpose.** A graded readout does not rescue
  expression constraints, because there is no sub-threshold signal to recover.

## The lesson, recorded

**A pre-registered bar cleared by a run whose run-to-run variance exceeds the effect is not a result.**
The pre-registration fixed the *threshold* but never specified *how many independent runs* had to clear
it. One run cleared +0.02 by +0.001. Two more runs put the estimate on both sides of zero.

Concretely, for this repo: **establish determinism BEFORE measuring an effect of this size, and require
at least two deterministic runs in agreement before publishing a direction.** The `processes=1` pin now
makes that cheap; it should have preceded the claim, not followed it.

This is the same failure class as the two earlier reversals in this epoch (the `min_abs_t` filter, the
verdict function comparing against nulls instead of the baseline arm) — a correctly-fired rule sitting on
an unexamined premise.
