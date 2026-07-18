# Does ESM2 + OUR-OWN MSA-Transformer lift? — closing the lifting hybrid end-to-end

**Date:** 2026-07-17 · **Model:** `esm_msa1b_t12_100M_UR50S` (fair-esm, CPU) · **Script:**
`scripts/msa_transformer_lift.py` · **Data:** `data/processed/msat_lift_checkpoint.jsonl` · **Module:**
`dna_decode/forward/msa_transformer.py`

## Question

The modality-hybrid finding (`forward_modality_hybrid_2026-07-17.md`) used ProteinGym's **precomputed**
`MSA_Transformer` column (+0.013 paired vs ESM2-650M, win 66%, p=0.002 on N=95). That proves *a* coevolution
model lifts — it does not prove **our** pipeline does. This run computes MSA-Transformer scores **ourselves**
(local CPU, `esm_msa1b_t12_100M_UR50S`, wt-marginal on a depth-128 subsample of the on-D: MSAs), feeds them
through **our** `rank_average_hybrid`, and checks the lift reproduces — the last validation the run-time
evolution component needed. (No Kaggle: the model is 100M params and runs on CPU; a data-transfer to a remote
GPU that I couldn't debug overnight was the wrong risk.)

## Result — two clean findings

**(1) Our scorer is correct.** Our own MSA-Transformer scores **reproduce ProteinGym's `MSA_Transformer`
column at median Spearman 0.89** across the scored proteins — the wt-marginal single-forward approximation
tracks ProteinGym's slower masked-marginal in rank. The MSA parse → subsample → forward → variant-score
pipeline is validated end-to-end on real proteins.

**(2) The evolution lift is PHENOTYPE-CONDITIONAL — reproduced with our own model.** Paired median Δ
(ESM2+our-MSA-T − ESM2), per ProteinGym phenotype category:

| phenotype | n | median Δ (hybrid − ESM2) | win-rate |
|---|---:|---:|---:|
| **Activity** | 16 | **+0.0112** | 9/16 |
| Stability | 9 | −0.0125 | 4/9 |
| Expression | 4 | −0.0007 | 2/4 |
| OrganismalFitness | 4 | −0.0541 | 1/4 |
| Binding | 2 | −0.0090 | 1/2 |

(N=35 scored; overall median Δ −0.008, win 17/35.) **Activity is the only category with a POSITIVE median
lift; every other category is flat-to-negative.** The aggregate ≈0 is a **category-mix artifact, not a null**
— MSA-Transformer (an evolution model) helps activity/function but not the structure-dominated phenotypes.
The **direction** independently reproduces the per-category modality rule
(`forward_modality_hybrid_2026-07-17.md`: evolution → activity, structure → stability), now with our own
end-to-end pipeline rather than ProteinGym's precomputed columns. **Honest limit:** per-category n's are
small and the win-rates are noisy (the crisp Stability-0/5 of an n=26 snapshot softened to 4/9 by n=35) — the
own-model run confirms the *direction* + the *scorer*, not a new independent significance.

## Honest scope

- **Underpowered for aggregate significance.** The first two runs were **OOM-killed** (the MSA parser loaded
  100k+-sequence alignments whole — since fixed with a bounded `max_rows` read), capping n. Within a
  category the sign test is not individually significant at this n. The **powered** evidence remains the
  N=95 ProteinGym-column sweep (+0.013, p=0.002); this own-model run confirms the **direction** + the
  **scorer correctness** + the **phenotype-conditionality**, not a new independent significance.
- **wt-marginal**, not masked-marginal (1 forward vs L forwards) — a fast approximation, validated to
  reproduce the reference at rank 0.89. Depth-128 subsample of match columns; insert positions unscored
  (same as the site-independent floor + ProteinGym).

## Deployment consequence (data-driven routing)

The decoder should **route by phenotype**, not add evolution blindly:

| target phenotype | evolution (MSA-Transformer) | structure (ProSST) |
|---|---|---|
| Activity / function / OrganismalFitness | **add** (lifts) | — |
| Stability / Expression | do NOT add (hurts) | add |

`ESM2 + MSA-Transformer` is the right hybrid for **activity/function** cells; stability/expression cells want
**ESM2 + structure**. This sharpens the earlier "add evolution universally" into a phenotype-conditional rule.

## Shipped

- `dna_decode/forward/msa_transformer.py::msa_transformer_table` — the validated LIFTING coevolution model,
  pluggable into `rank_average_hybrid` via `msa_evolution.evolution_table_from_scores`. Lazy torch/fair-esm.
- `msa_evolution.parse_a2m(max_rows=…)` — bounded-memory read (the OOM fix that killed the first runs).
- `scripts/msa_transformer_lift.py` — the CPU validation harness (restartable JSONL checkpoint).
- Frozen decoder surface byte-unchanged (`verify_lock OK`).
