# The deficit is fully partitioned — and the condition-modelling hypothesis is empty

**Date:** 2026-08-22 · **Artifact:** `wiki/fba_missed_partition_2026-08-22.json`
**Script:** `scripts/fba_missed_gene_partition.py` · **Tests:** `tests/test_fba_missed_partition.py` (5)
Model + `feba.db`. 25 labelled carbon conditions × 217 two-sided genes.

## The result

The 131 conditionally-essential genes (11-condition panel, the denominator every prior artifact in this
arc uses) now partition **completely**, with no residue:

| class | genes | share | what would fix it |
|---|---|---|---|
| **PROVABLY UNCALLABLE** | **31** | 23.7 % | nothing *within the standard-FBA setting* — see the scope correction below |
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

## SCOPE CORRECTION (2026-08-22, after adversarial review) — the uncallability claim was too broad

The row above originally read *"nothing — no medium, objective or constraint layer reaches them."*
**That is false as written**, and the same overbroad phrasing appeared in
`wiki/fba_expression_gated_gpr_result_2026-08-22.md`.

The proofs assume **boolean GPR, free enzyme capacity, fixed internal bounds, exchange-defined media, and
an objective linear over the retained reactions**. An enzyme-constrained reformulation (GECKO / ecFBA)
breaks the second assumption: isozymes carry different kcat against a shared proteome budget, so deleting
the *efficient* isozyme can reduce growth even when the boolean GPR is still satisfied.

**Corrected claim:** the 31 is a floor over **any constraint layer applied to iML1515's own reaction
set**, for any exchange-defined medium and any objective linear over the retained reactions. It does
**not** survive a *reaction-splitting* enzyme-constrained reformulation.

**How far that actually bites here — measured, not assumed:**

| fact | measured |
|---|---|
| iML1515 reaction annotations | `ec-code`, `sabiork` present — **no kcat values** |
| enzyme-pool / `prot_` / `draw_` reactions | **none** |
| isozyme representation | both complexes share **ONE** reaction via an `or` rule (`ACHBS` = `(b3670 and b3671) or (b0077 and b0078)`) |
| uncallable-essential genes sitting on an `or` rule | **22 of 23** |

So a plain GEM **structurally cannot express per-isozyme capacity at all** — the guarantee is lost only
under a model *transformation* (splitting each reaction per enzyme and assigning kcats), not under any
constraint added to iML1515 as shipped. This repo has no such variant, and the `ec-code`/`sabiork`
cross-references are the path to building one, not evidence that one exists.

**Robustness split:** the 6 `ALL_DISABLED_BLOCKED` genes survive even a reformulation (a zero-capacity
reaction stays zero under any tightening). The other 25 lose the **proof** — which is not the same as
becoming callable; an alternate enzyme with ample capacity, equal cost, or no kcat annotation would leave
them uncallable in practice. That remains untested.

## CORRECTION (same day) — "NEVER FIRES" is two classes, not one, and I labelled it too coarsely

The table above assigns all 63 NEVER_FIRES genes to "the objective / network". **That is half wrong.**
Measuring their worst-case deletion ratio across all 25 conditions splits them cleanly:

| | genes | worst ratio | what it means |
|---|---|---|---|
| **NO EFFECT** | **34** | ~1.000 everywhere | the model genuinely does not care → **objective incompleteness** |
| **SUB-THRESHOLD DEFECT** | **29** | 0.744 – 0.989 | the model **does** predict a real defect; the binary cutoff discards it |

The sub-threshold set is not marginal noise — it is coherent respiratory and central metabolism:

| genes | worst ratio | defect |
|---|---|---|
| `cyoA` `cyoB` `cyoC` (cytochrome *bo* oxidase) | 0.744 | **25.6 %** |
| `nuoA`–`nuoN` (NADH dehydrogenase I, 13 subunits) | 0.769 | **23.1 %** |
| `mdh` | 0.893 | 10.7 % |
| `sucC` `sucD` | 0.900 | 10.0 % |
| `rpe`, `aceE`, `galP`, `maeB`, `aceA` | 0.912–0.989 | 1–9 % |

So the corrected partition of the 131:

| class | genes | share | layer |
|---|---|---|---|
| provably uncallable | 31 | 23.7 % | unreachable |
| **sub-threshold defect** | **29** | **22.1 %** | **the READOUT** |
| **no effect anywhere** | **34** | **26.0 %** | **the OBJECTIVE** |
| partial overlap | 37 | 28.2 % | working |
| mis-conditioned | 0 | 0 % | empty |

**The class name is also misleading and is worth restating plainly:** `NEVER_FIRES` means *never crosses
the essentiality threshold*, **not** *never responds*. For 29 of 63 the model responds substantially and
the binary rule throws the response away.

**Restraint required on that 22.1 %.** This quantifies the earlier `nuo` observation — 13 genes at ratio
0.865 scored dispensable — as **29 of 131**. But it is *adjacent* to a lever already measured
quantitatively dead: the threshold sweep found ≤11 % of misses recoverable by moving the cutoff, and a
single global cutoff cannot separate a real 23 % defect from noise near 1.0 without flooding false
positives. The distinct, untried operation is transforming the prediction into the **label's own units**
(fitness ≈ n·(r−1)) rather than picking a prettier cutoff — and that is a **new endpoint**, so it must be
pre-registered before it is run. The expression-gating result earlier today is the reminder of why: a
primary endpoint can be met and still be worthless.

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
| sub-threshold defect | 29 | readout layer; adjacent to a dead lever — **pre-register before touching** |
| no effect anywhere | 34 | **the genuine frontier** — objective incompleteness |
| partial overlap | 37 | working at 96.3 % per-cell |
| mis-conditioned | 0 | **empty** |

The remaining question is no longer "which constraint-based lever next?" — that family is exhausted and
its ceiling is measured. It splits in two, and only one half is fresh territory:

- **34 genes (26.0 %) — what does biomass fail to demand?** LPS core biosynthesis, siderophore uptake,
  the stringent response. A question about the objective function itself, and genuinely untried.
- **29 genes (22.1 %) — the readout discards a real prediction.** Quantified, but neighbouring a lever
  already found dead; it earns a pre-registration, not a promise.
