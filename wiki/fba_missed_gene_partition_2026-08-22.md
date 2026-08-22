# The deficit is fully partitioned — and the condition-modelling hypothesis is empty

**Date:** 2026-08-22 · **Artifact:** `wiki/fba_missed_partition_2026-08-22.json`
**Script:** `scripts/fba_missed_gene_partition.py` · **Tests:** `tests/test_fba_missed_partition.py` (5)
Model + `feba.db`. 25 labelled carbon conditions × 217 two-sided genes.

## The result

The 131 conditionally-essential genes (11-condition panel, the denominator every prior artifact in this
arc uses) now partition **completely**, with no residue:

| class | genes | share | what would fix it |
|---|---|---|---|
| **PROVABLY UNCALLABLE** | **31** | 23.7 % | nothing — no medium, objective or constraint layer reaches them |
| **NEVER FIRES** | **63** | 48.1 % | the objective / network — the model never demands the product |
| **PARTIAL OVERLAP** | **37** | 28.2 % | already largely working |
| **MIS-CONDITIONED** | **0** | **0 %** | — |

The 31 reproduces the prior number **exactly** (23.7 %), which is a real consistency check across three
independently-written scripts.

## The finding: zero mis-conditioned genes

**Not one gene of 217 fires only in the wrong conditions.** The plausible, comfortable hypothesis — *"the
model basically knows these genes matter, our medium modelling just puts them in the wrong place"* — has
**zero support**. It was the cheapest remaining explanation and it is empty.

What the model does instead is **bimodal**. For genes it fires for at all, it is almost right:

> of the 749 true cells belonging to PARTIAL_OVERLAP genes, **721 are caught — 96.3 %**

and for the rest it is completely silent — 63 of 131 never predicted essential in **any** of the 25
labelled conditions, despite being structurally capable of it (they are not in the uncallable set).

So the model does not degrade gracefully. It either has a gene or it does not.

## Why that matters for what comes next

The four failed levers — gap-fill, threshold retune, pFBA, E-Flux — all reshape **flux under a fixed
objective**. This partition says that is aimed at the wrong layer for essentially the entire residue:

- **23.7 %** is unreachable by any of them (proved, validated over ~40k deletions).
- **48.1 %** never fires at all, so there is no flux to reshape — the biomass objective simply never
  demands whatever these genes make.
- **0 %** is a calibration problem, which is the only thing that class of lever could have fixed.
- **28.2 %** already works at 96.3 % per-cell.

This is the same shape as the earlier FLAT finding (for ~91 % of missed cells the deletion changes growth
by *nothing*) — but where that was a symptom, this names the layer: **objective incompleteness**, and it
is now quantified rather than inferred.

## Honest limits

1. **Denominator discipline.** This script's own panel is 25 conditions and yields **217** two-sided
   genes; the 131 comes from the 11-condition expression panel. Both are reported and the artifact
   carries an explicit `denominator_warning`. On the 217-gene panel the shares are 22.6 % / 44.2 % /
   33.2 % / 0 % — same story, different base.
2. **Sole-carbon-source conditions only.** A gene whose true dependence is on aeration, pH, or a nitrogen
   source cannot be caught by this panel and will read as NEVER_FIRES. That inflates NEVER_FIRES and is
   the main caveat on the 48.1 %.
3. **The "ever fires" test is deliberately generous** — it spans all 25 labelled conditions, not the 11
   with expression, so a gene lands in NEVER_FIRES only after failing everywhere we could look.
4. **Model reach, not biology.** This partitions the *model's* failure modes.
5. `MIS_CONDITIONED = 0` is an empirical result on this panel, not a theorem. The class is real and
   reachable — `tests/test_fba_missed_partition.py` pins it so it does not quietly become dead code.

## Where the arc stands

| | genes of 131 | status |
|---|---|---|
| provably uncallable | 31 | **closed** — proved + validated |
| …of which isozyme-masked | 8 | **closed** — expression-gating measured and rejected |
| never fires | 63 | **the frontier**, and now correctly labelled: objective, not constraints |
| partial overlap | 37 | working at 96.3 % per-cell |
| mis-conditioned | 0 | **empty** |

The remaining question is no longer "which constraint-based lever next?" — that family is exhausted and
its ceiling is measured. It is **what does biomass fail to demand?**, which is a different kind of
question about the objective function itself.
