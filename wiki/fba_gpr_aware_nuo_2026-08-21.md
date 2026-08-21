# The GPR fix lands, and the residue is one operon the binary cutoff throws away

**Date:** 2026-08-21 · **Artifact:** `wiki/fba_fva_requirement_2026-08-21.json`
**Script:** `scripts/fba_fva_requirement_class.py` · Model-only (D: still disconnected)

A review of the diagnostic program raised two defects and recommended a specific next experiment. Both
defects are now fixed, and fixing the second one exposed something concrete.

## Fix 1 — the zero-flux claim was overbroad, and my own data refuted it

`scripts/fba_flat_mechanism_partition.py` and `wiki/fba_innovate_2026-08-20.md` asserted *"No
constraint-based method can move such a cell."*

**That is false, and the counterexample is in this repo:** `wiki/fba_structural_blindspot_2026-08-21.md`
found that of 58 genes pinned at zero flux in glucose, **57 carry flux fine once other exchanges are
opened** — and opening an exchange *is* a constraint change.

Corrected claim, now in both files: *holding the objective and the medium FIXED, no intervention that
merely tightens or removes capacity on reactions carrying no flux can move the cell.* Changing medium,
objective, or adding a demand term are explicitly **outside** the claim. The numbers were unaffected.

## Fix 2 — the classifier attributed reactions a knockout does not disable

`fba_fva_requirement_class.py` built gene classes from raw `gene.reactions`. But cobrapy's
`Gene.knock_out()` zeros a reaction only `if not reaction.functional`, and `Reaction.functional`
evaluates the **full GPR** — so a reaction behind an isozyme `or` survives a single-gene deletion.

New helper `gpr_disabled_reactions()` simulates the knockout and takes only the reactions that actually
become non-functional. Effect:

| | raw `gene.reactions` | GPR-aware |
|---|---|---|
| REQUIRED | 55 | **39** |
| CAPABLE_BUT_IDLE | 18 | 17 |
| INACTIVE_IN_CONDITION | 58 | 58 |
| KO_DISABLES_NOTHING | — | **17** (new) |
| **FVA-vs-deletion disagreements** | **29** | **13** |
| …of which isozyme-explained | 16 | **0** |

**The fix removed exactly the disagreements it predicted.** 16 genes moved out of REQUIRED into
`KO_DISABLES_NOTHING` (fully isozyme-buffered — the KO disables no reaction at all, so the model can
never call them essential), and every isozyme-explained disagreement disappeared.

## What the residue turned out to be

The 13 remaining disagreements are **not** a mixed bag. They are **the entire nuo operon** —
`nuoA B C E F G H I J K L M N`, NADH dehydrogenase I — every subunit, identical growth ratio **0.865**.

This is the joint-feasibility limit made concrete: each nuo reaction can *individually* hit zero at the
optimum (so FVA says zero-attainable), but they cannot all be zero *simultaneously*, so the actual
deletion costs **13.5 % of growth**.

**And the binary rule throws that away.** The model scores essential at ratio ≤ 0.01; 0.865 is scored
**dispensable**. So for these 13 genes the model predicts a real, sizeable defect and the readout
discards it.

### How big is 13.5 % on the label's own scale?

RB-TnSeq fitness is roughly a log2 abundance change over competitive growth. To first order, for a mutant
growing at ratio *r* over *n* generations of wild-type growth:

> fitness ≈ n · (r − 1)

| growth ratio | ≈ fitness at n=25 | vs the `fit < −2` label |
|---|---|---|
| **0.865 (nuo)** | **−3.38** | **detected as essential** |
| 0.95 | −1.25 | not detected |
| 0.99 | −0.25 | not detected |

So a 13.5 % defect is comfortably inside RB-TnSeq's detection range, while a 5 % defect is not.

*(Correction: I first wrote this as n·log2(r), giving −5.23. That overstates the magnitude. The
first-order relation is n·(r−1). Both cross the −2 bar here, so the conclusion is unchanged, but the
formula matters and the linear one is right. **This is an approximate mapping** — the Fitness Browser's
exact normalisation is not verified here, so treat −3.38 as an order-of-magnitude argument, not a
calibrated prediction.)*

## Why this is not the failed threshold-retune lever again

Threshold retuning was already tested and found quantitatively dead (≤11 % of misses recoverable). That
swept a **cutoff on the ratio**, and a single global cutoff cannot separate a real 13.5 % defect from
noise near 1.0 without flooding false positives.

What the nuo case suggests is different: **transform the prediction into the label's own units** —
convert the growth ratio to an expected fitness via n·(r−1) and compare on that scale — rather than
picking a prettier cutoff. That is principled because it maps the prediction into the measurement's
scale, not because it improves the number.

**It is also not yet tested, and it must be pre-registered before it is.** Introducing a new endpoint
after the binary one disappoints is metric-shopping unless the labels, gene×condition universe,
transform, statistic and success criterion are frozen first — and this project has already retracted one
positive for exactly this class of discipline failure. Binary stays primary; any fitness-scale endpoint
is a declared secondary.

## Honest limits

1. **n=13, one operon, one condition.** A vivid instance, not a quantified share of the deficit.
2. **The fitness mapping is approximate and unverified** against the Fitness Browser's normalisation.
3. **`INACTIVE_IN_CONDITION` (58) is unchanged and still condition-confounded** — all genes evaluated in
   glucose though essential across an 11-condition panel. The GPR fix does not touch that.
4. **`KO_DISABLES_NOTHING` (17) is a new, useful class** — fully isozyme-buffered genes the model can
   *never* call essential — but it is measured in glucose only.

## Next

The recommended experiment is **bypass-closure**: for each missed gene, simulate the GPR-aware deletion,
then close the free inorganic/cofactor exchanges that bypass its pathway, and watch for a
**deletion-ratio flip**. Flip ⇒ medium realism; no flip ⇒ objective incompleteness. It is model-only and
needs no labels for the first pass, and it directly targets the free-iron observation
(`EX_fe2_e`/`EX_fe3_e`/`EX_zn2_e` open at `(-1000,1000)` with all 12 acquisition genes idle).
