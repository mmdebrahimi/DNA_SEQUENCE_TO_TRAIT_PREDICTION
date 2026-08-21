> # EXPLANATION CORRECTED 2026-08-21 (same day)
>
> The **INCONCLUSIVE verdict below stands** — closing `EX_fe2_e`/`EX_fe3_e` really does starve the model.
> But the **explanation given for it was wrong**, and repair-then-close showed why.
>
> This memo said: *every ferri-siderophore exchange is export-only, so no loaded siderophore can be
> imported.* True about the exchanges, **misleading as the cause** — the loaded siderophore never needs
> importing from the medium. `FEENTERexs` forms it EXTRACELLULARLY (`enter_e + fe3_e -> feenter_e`) from
> enterobactin the cell secretes, and `FEENTERtonex` (fepA/tonB/exbB/exbD) brings it in. The route is
> fully present, gene-associated, and traversable.
>
> The real reason closing those exchanges starves the model is that **`fe3_e` is the shared input to BOTH
> routes** — the siderophore path consumes extracellular ferric iron too, so closing the exchange kills
> the alternative as well as the shortcut.
>
> The actual cause of the idle machinery is **route redundancy plus orphan reactions**. See
> `wiki/fba_orphan_redundancy_2026-08-21.md`.

# Bypass closure returns INCONCLUSIVE — and diagnosing why beats the flip it was looking for

**Date:** 2026-08-21 · **Verdict:** `INCONCLUSIVE_WT_INFEASIBLE` (both bypasses)
**Artifact:** `wiki/fba_bypass_closure_2026-08-21.json` · **Script:** `scripts/fba_bypass_closure.py`
Model-only — D: still disconnected.

## The experiment

Two hypotheses both predict an idle gene, so idleness alone cannot separate them:

- **Medium realism** — the model gets a nutrient free, so the acquisition machinery idles.
- **Objective incompleteness** — biomass never demands the product, so the pathway idles.

The separating observable is a **deletion-ratio flip**: close the free shortcut, keep the wildtype
feasible, re-run the deletion. Flip ⇒ medium realism. No flip ⇒ objective incompleteness.

First case: iron/zinc, because `EX_fe2_e`, `EX_fe3_e`, `EX_zn2_e` sit in the default medium at
`(-1000, 1000)` while all 12 named Fe/Zn acquisition genes carry zero flux.

## Result: the pre-declared failure branch fired

| bypass | closed | wildtype open → closed | verdict |
|---|---|---|---|
| iron_zinc | `EX_fe2_e`, `EX_fe3_e`, `EX_zn2_e` | 0.876997 → **0.000000** | `INCONCLUSIVE_WT_INFEASIBLE` |
| iron_only | `EX_fe2_e`, `EX_fe3_e` | 0.876997 → **0.000000** | `INCONCLUSIVE_WT_INFEASIBLE` |

This outcome was **named in the script before running**: if closing the tap kills the wildtype, that is
starvation, not realism, and it is reported as inconclusive — never as a flip.

## Why it happened — the actual finding

iML1515 **does** carry a complete siderophore system. There are **26 reactions producing cytosolic
`fe2_c`**, covering enterobactin (`FEENTERR1/R2`), ferrichrome (`FECRMR2`), ferrioxamine (`FEOXAMR1`),
aerobactin (`ARBTNR3`), ferric citrate (`FE3DCITR5`) and dihydroxybenzoate (`FE3DHBZS3R`).

**But every ferri-siderophore exchange is export-only.** Measured bounds:

| exchange | bounds | in medium |
|---|---|---|
| `EX_feenter_e`, `EX_fecrm_e`, `EX_feoxam_e`, `EX_arbtn_fe3_e`, `EX_fe3dcit_e`, `EX_fe3hox_e`, `EX_fe3dhbzs_e` | **`(0.0, 1000.0)`** | no |
| `EX_fe2_e`, `EX_fe3_e` | `(-1000, 1000)` | **yes** |

A lower bound of `0.0` means **secretion only — no uptake**. So no loaded siderophore can ever enter the
cell, and with the direct tap closed there is **no route to iron at all**: closing `EX_fe2_e`/`EX_fe3_e`
while leaving *every other exchange open* still gives growth `-0.000000`.

**So the 12 acquisition genes are not idle because biomass ignores iron — biomass demands iron, and gets
it through a direct exchange that bypasses the machinery entirely. The machinery is structurally unable
to deliver iron in this reconstruction.**

That is a **third** category, distinct from both hypotheses the experiment was built to separate:

> The pathway is present, the product is demanded, and the pathway still cannot carry flux — because its
> substrate can never be imported.

## What this changes

- **For these 12 genes, the answer is neither "medium realism" nor "objective incompleteness."** It is an
  import-direction gap in the reconstruction. No demand term and no medium tweak within the current
  exchange bounds can fix it.
- **The bypass-closure design stands** — it just needs a nutrient whose alternative route is actually
  traversable. Iron in iML1515 is not such a case, and now we know why rather than guessing.
- **It sharpens the earlier free-iron observation.** Previously: "the model gets free iron, so the
  siderophore system idles." Now: "the model gets free iron *and could not use the siderophore system
  even if it didn't*."

## Honest limits

1. **One condition (glucose), one nutrient class.** Zinc was only tested jointly with iron; the
   `iron_only` arm isolates iron, and both fail the same way.
2. **Gene sets are curated by name**, not exhaustive — other genes may also sit on the same dead route.
3. **This is a model-structure claim, not a recovered true positive.** It explains why 12 genes can never
   be predicted essential; it does not by itself improve any score.
4. **Not yet tested:** whether opening a ferri-siderophore uptake (setting e.g. `EX_feenter_e` lower
   bound below zero) and *then* closing the direct tap produces the flip. That is the natural next
   version of this experiment and would make the discrimination actually run.

## Next

Two model-only follow-ups, in order:

1. **Repair-then-close:** open a ferri-siderophore uptake, verify the wildtype stays feasible, then close
   the direct tap and look for the flip. This is the version of the experiment that can return a verdict
   for iron.
2. **Pick a bypass with a traversable alternative** — a nutrient where the model already has an
   importable substrate for the alternate route — and run the original design there.
