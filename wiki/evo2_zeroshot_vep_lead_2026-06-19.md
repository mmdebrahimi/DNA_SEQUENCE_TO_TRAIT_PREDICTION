# Discovery-tier lead — Evo 2 zero-shot variant-effect by likelihood (2026-06-19)

**Source:** Silvana Konermann (Arc Institute) TED talk "The Human Cell Is Wildly Complex. Can AI Decode
It?" (the Arc "virtual cell" pitch) → researched the underlying tech (Evo 2, STATE, virtual cell stack).
**Tier:** discovery / parked horizon. **NOT a pivot.** Captured so the one genuinely-untried angle isn't lost.

## What was rejected (and why) — the bulk of the talk's tech

| Arc tech | Useful here? | Why not |
|---|---|---|
| STATE / virtual cell / scBaseCount | No | single-cell RNA-seq, human, needs perturbational data we can't generate (no wet lab); non-commercial license |
| MULTI-evolve | No | protein design / directed evolution — out of scope |
| Evo 2 as a TRAINED embedding | No (CLOSED) | same class we falsified — 3 de-confounded embedding negatives across the kingdom boundary; our own rule: "do NOT scale embeddings on a bigger/paid GPU — a negative de-confounded metric is signal-vs-structure, not window-budget" |

## The ONE untried angle — Evo 2 zero-shot variant-effect-by-likelihood

Evo 2 scores a SNV's effect by the **reference-vs-variant likelihood delta** — **label-free, no trained
classifier**. We have NEVER tried this method; it is DIFFERENT from the mean-pool-embedding-XGBoost we
falsified. Falsifiable hypothesis:

> Does Evo 2's reference-vs-variant likelihood delta separate KNOWN cipro-resistance QRDR mutations
> (gyrA S83L / parC) from neutral SNVs on the existing N=147 cipro cohort?

**⚠️ Conceptual caveat (the reason expected value is LOW):** delta-likelihood scores *evolutionary
disruptiveness / conservation*. AMR resistance mutations are frequently **common under positive selection**
(low disruptiveness) — the OPPOSITE of a rare-deleterious variant. So Evo 2 may be measuring conservation,
not resistance. That is the test: a NEGATIVE result is informative (confirms conservation≠resistance for
AMR); a POSITIVE result would be a genuinely new label-free determinant signal. Either way it's cheap.

## Constraints
- **GPU-gated, not on the GTX 860M.** Evo 2 7B runs without fp8 on older devices but needs real GPU memory;
  use the hosted NVIDIA BioNeMo API or a Databricks / Precision-7780 GPU. Open weights + inference.
- Does NOT reopen the closed embedding arm (different method: zero-shot likelihood, not a trained embedding
  classifier). If run, it's a scoped probe with a pre-registered hypothesis, not a strategic pivot.

## The deeper signal — the talk VALIDATES our strategic position
Arc — $billion-funded, Evo 2 (40B), STATE (270M cells) — frames its own core unsolved problem as
"generalize to ENTIRELY NEW CONTEXTS without context-specific data," "models still fall short of real
experiments," "training-data quality is the critical unresolved factor." That is **exactly** our finding
(`wiki/negative_results_map_2026-06-13.md`: the bottleneck is LABELS not models). A frontier institute is
hitting the same wall we documented — external confirmation that the "ship the deterministic decoder + the
genome-map honesty tool, target embeddings only at sampling-independent lab labels" pivot was correct.

## Disposition
PARKED. If/when on a GPU host with idle cycles, run the scoped probe above (pre-register the hypothesis +
the conservation≠resistance caveat; report the negative honestly). Do NOT let "Evo 2" become a reason to
reopen the closed embedding bet. Sources: ted.com/talks/silvana_konermann_… · arcinstitute.org/tools/evo ·
nature.com/articles/s41586-026-10176-5 · docs.nvidia.com/bionemo-framework (Evo 2 zero-shot BRCA1 VEP).
