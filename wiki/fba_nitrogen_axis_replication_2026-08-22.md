# The partition replicates on a second environmental axis — and the data ceiling is now measured

**Date:** 2026-08-22 · **Artifact:** `wiki/fba_missed_partition_nitrogen_2026-08-22.json`
**Script:** `scripts/fba_missed_gene_partition.py --axis nitrogen`
Resolves the **M2** issue left open by the adversarial review: *`MIS_CONDITIONED = 0` is a one-axis result
and cannot support "the condition-modelling hypothesis is empty" beyond the sampled axis.*

## 1 · The data ceiling, measured rather than asserted

The prior run banked on "the plateau is data-limited." That was a claim. It is now a measurement:

| | |
|---|---|
| Fitness Browser `Keio` carbon conditions | **28** (25 map to an iML1515 exchange) |
| PRECISE-1K distinct carbon values | **18** |
| **intersection** | **11** |
| of the 14 unmatched conditions, how many appear in PRECISE-1K under any name | **0** |

**The 14 are genuinely absent, not a name-normalisation bug** — I checked each against the compendium's
own carbon column. Both sides are fixed datasets, so no plumbing fix widens the intersection.

Worse, and previously unquantified: the matched 11 are extremely thin. Of 1,035 samples, **621 are
glucose** and 111 glycerol; the rest are `acetate 8 · fructose 8 · pyruvate 6 · xylose 5 · sorbitol 2 ·
ribose 2 · NAG 2 · galactose 2 · gluconate 2`. **Six of the eleven scored conditions rest on ≤5 samples,
four on exactly 2.** That is the concrete reason the composed selector recovered 1 of 8 — the per-condition
expression estimate for most conditions is a mean over two replicates.

## 2 · The second axis: `MIS_CONDITIONED = 0` replicates

The Fitness Browser `Keio` set also has **nitrogen source (32 experiments)** and **stress (55)**, and this
arc had only ever used carbon. 13 nitrogen conditions map to iML1515 exchanges. Running the identical
partition:

| class | carbon (217 genes / 25 cond) | **nitrogen (155 genes / 13 cond)** |
|---|---|---|
| NEVER_FIRES | 44.2 % | **41.9 %** |
| PARTIAL_OVERLAP | 33.2 % | **34.8 %** |
| PROVABLY_UNCALLABLE | 22.6 % | **23.2 %** |
| **MIS_CONDITIONED** | **0 %** | **0 %** |

> **Zero mis-conditioned genes on an independent axis, with a different gene set (155 vs 217) and
> different conditions.** The class shares land within ~2 points of the carbon partition.

So `MIS_CONDITIONED = 0` is **not** a one-axis artifact. The review's M2 is resolved in the direction that
strengthens the original finding: the model does not fire in the wrong place — it either fires roughly
right or stays silent.

## 3 · The 100 % catch rate is NOT a triumph — it is near-constitutive prediction

Nitrogen PARTIAL_OVERLAP genes catch **488 of 488** of their true cells (carbon: 721/749 = 96.3 %). That
number flatters, and inspecting it says why:

| axis | PARTIAL_OVERLAP genes | predicted essential in **every** condition | genuinely condition-specific (1 condition only) |
|---|---|---|---|
| nitrogen | 54 | **38 (70 %)**, and 53 of 54 in ≥12 of 13 | ~0 |
| carbon | 72 | 39 (54 %) | **25** |

On nitrogen the model is **almost entirely constitutive** among the genes it fires for, so catching every
true cell is close to automatic. Carbon has real condition-specificity for a third of them; nitrogen has
essentially none.

The likely reason is visible in the wildtype growths: **six of the thirteen nitrogen conditions give
*exactly* 0.92593** (D-alanine, D-serine, glycine, L-alanine, L-aspartate, L-serine). The model treats
those sources as metabolically interchangeable, so it has little room to be condition-specific in the
first place.

## 4 · Honest limits

1. **The nitrogen axis is not nitrogen-only.** `dna_decode/fba/nitrogen.py` documents this: several
   sources (alanine, serine, aspartate, glutamine) also supply carbon, matching the real assay (glucose
   minimal + test compound as sole N source). So it is a *partially overlapping* second axis, not an
   orthogonal one — the replication is weaker than a fully independent axis would be.
2. **Continuity view.** Restricted to genes shared with the carbon 131 (74 overlap): NEVER_FIRES 48.6 %,
   PARTIAL_OVERLAP 27.0 %, PROVABLY_UNCALLABLE 24.3 % — same shape.
3. ~~**Stress (55 experiments) is still untouched** and is the more genuinely independent axis.~~
   **CORRECTED, same day — this was wrong and would have sent a reader down a foreclosed path.** The
   stress axis was already probed on 2026-08-20 and is a **structural NO-GO**
   (`wiki/fba_stress_feasibility_2026-08-20.md`, verdict `STRESS_AXIS_NOT_REPRESENTABLE_IN_iML1515`):
   only 8 of 35 conditions have an exchange at all, **0 of those 8 reduce growth**, and acetate — a
   stressor in the assay — makes the model grow **28 % faster**. An exchange models *supplementation*, so
   the medium-swap contract has the **wrong sign** for stress. I recommended it here without checking.
4. **This is a partition of the model's failure modes**, not of biology.
5. The uncallable set was carried over from the carbon-derived artifacts; 36 of the 155 nitrogen genes
   fall in it. It was not re-derived on nitrogen.

## 5 · What it changes

- **M2 is resolved.** The condition-modelling hypothesis stays empty across two axes, not one.
- **The data ceiling is now a number, not an intuition:** 11 usable conditions, most on ≤5 expression
  samples, and 14 carbon sources with no expression data in principle. Expanding it needs a *different
  compendium*, which is an **external wall** — no amount of modelling closes it.
- **The environmental-axis expansion is EXHAUSTED.** Not a judgement — an enumeration of every
  `expGroup` the `Keio` organism has:

| axis | distinct conditions | status |
|---|---|---|
| carbon source | 28 (25 mappable) | **done** — the arc's primary panel |
| nitrogen source | 16 (13 mappable) | **done — replicates (this memo)** |
| stress | 35 | **structural NO-GO** (2026-08-20; wrong sign) |
| sulfur source | **2** (`L-Cysteine`, `Sodium sulfate`) | too few — a two-condition panel cannot support a meaningful two-sided partition |
| motility | 1 (`Agar`) | not a metabolic axis |
| lb | 0 | — |

There is no third usable axis. The aeration/pH caveat carried since the first partition memo therefore
**cannot be discharged with this dataset at all** — it is not a task waiting to be done, it is outside
what `feba.db` + iML1515 can express together.
