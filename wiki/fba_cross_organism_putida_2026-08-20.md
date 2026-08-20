# All three findings replicate on a DIFFERENT ORGANISM and a DIFFERENT MODEL

**Date:** 2026-08-20 · **Verdict:** `ALL_THREE_REPLICATE` · **Determinism gate:** PASS
**Pre-registration:** `wiki/fba_cross_organism_prereg_2026-08-20.md` (frozen before any solve)
**Artifact:** `wiki/fba_cross_organism_putida_2026-08-20.json` · **Runner:** `scripts/fba_cross_organism.py`

## Result

*Pseudomonas putida* KT2440 on **`iJN1463`** (2,927 reactions / 1,462 genes) — a different organism and a
different genome-scale model from every prior result. Panel: 13 of 47 assay carbon sources. 1,148
complete-row genes → **168 conditionally essential**.

| prediction | E. coli carbon | E. coli nitrogen | **P. putida** | |
|---|---|---|---|---|
| **P1** threshold band near-empty | 0 / 16,676 | 0.0 % | **0.0 %** | REPLICATES |
| **P2** missed-essential cells FLAT | 76.9 % | 91.1 % | **70.4 %** (138/196) | REPLICATES |
| **P3** beats best-constant null | positive | +0.1444 | **0.794 vs 0.6932 (+0.1008)** | REPLICATES |

Determinism: **0 call flips**, max drift 1.4e-11, safety factor **7.18e8** (bar 1000), headline metric
**identical** across both passes (0.794 / 0.794).

Confusion totals: TP 474 · FP 254 · FN 196 · TN 1,260 — sensitivity 0.707, precision 0.651, on a
30.7 % essential-cell base rate.

## Why this is the strongest replication so far

Every previous result shared one organism and one model. This changes **both**, and the three findings
survive:

- **P1 is the load-bearing one**, because it underwrites the closed E-Flux negative. Exactly zero of
  2,184 cells land in the `[0.001, 0.05]` band on a *second* model — so bimodality is a property of the
  **readout** (a threshold on an FBA growth ratio), not of iML1515. The E-Flux closed negative stays
  closed, and now for a stated reason rather than an untested inference.
- **P2 replicates but is the weakest of the three** (70.4 %, vs 91.1 % on nitrogen and 76.9 % on carbon).
  The conditional deficit is still mostly deletions that changed nothing, so it is still a model problem —
  but the spread across three panels (70–91 %) is worth remembering before quoting any single figure.
- **P3's lift is +0.1008** over a high (0.6932) constant null. Real, comparable to nitrogen's +0.1444,
  and not an artefact of class imbalance: the base rate is stated above.

## Verification performed

**The generalized loader was checked against the audited one before any result was trusted.** The whole
P. putida number rests on `load_org_records`, which I wrote — `fitness_browser.load_records` hardcodes
`orgId='Keio'` and is heavily pinned, so rather than edit it I mirrored it. Run on E. coli, the two
produce **1,339 genes each, identical gene sets, 0 differing experimental calls**. Pinned as
`test_generalized_loader_equals_the_pinned_loader_on_ecoli` (skipped when feba.db is detached).

The other seven tests pin the contract that makes that equivalence possible — replicate averaging,
partial-row dropping, gene filtering, and the case that would be easiest to get wrong:
**`Keio` and `Putida` share experiment names** (`e1`), so a join that forgot `orgId` would silently import
another organism's fitness values. `test_another_organism_sharing_an_expName_cannot_bleed_in` fires on it.

## Named caveats

1. **n = 1 new organism.** A replication, not a survey. Two organisms agreeing is not "generalises across
   bacteria" — and the axis is *exhausted at two*, because of the 48 organisms in `feba.db` only E. coli
   and P. putida have a genome-scale model in the repo's BiGG registry.
2. **13 of 47 sources.** P. putida's assay is rich in compounds `iJN1463` does not carry (diols, alcohols,
   fatty acids), so the panel is biased toward central metabolism.
3. **Glucose is over-represented** in the assay (29 of 128 carbon experiments), so its condition is better
   determined than the others.
4. **`iJN1463` is less validated than `iML1515`.** That it performs comparably is reassuring, but a
   difference either way could reflect model quality rather than the method.

## Where the FBA epoch now stands

| axis | result |
|---|---|
| carbon (E. coli, 25 cond) | three findings established |
| nitrogen (E. coli, 13 cond) | all three replicate |
| **cross-organism (P. putida, 13 cond)** | **all three replicate** |
| stress (E. coli, 35 cond) | **structural NO-GO** — a stoichiometric model cannot represent a poison |
| E-Flux bridge | closed on both binary and graded readouts |

The findings hold across substrate **and** organism/model, and fail only at the perturbation-class
boundary — which is a property of the FBA formalism, not of the panels.

**The named axes are now exhausted.** Extending cross-organism further needs a genome-scale model for a
third Fitness-Browser organism, which is a model-acquisition task (BiGG has none for the other 46), not an
analysis task.
