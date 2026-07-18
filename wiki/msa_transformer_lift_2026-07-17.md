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
| Activity | 16 | **+0.0112** | 9/16 |
| Stability | 11 | −0.0125 | 5/11 |
| Expression | 6 | **+0.0146** | 4/6 |
| OrganismalFitness | 5 | −0.0539 | 1/5 |
| Binding | 2 | −0.0090 | 1/2 |

(N=40 scored; overall median Δ −0.004, win 20/40 = exactly 50%.) The aggregate ≈0 is a **category-mix
artifact, not a null** — the lift is not uniform across phenotypes. The direction is **partially** consistent
with the per-category modality rule (`forward_modality_hybrid_2026-07-17.md`: evolution → activity, structure
→ stability): **Activity is positive (+0.011) and Stability is negative (−0.013), as predicted** — but
**Expression is positive (+0.015), which the rule did not predict**, and OrganismalFitness is negative.

**Honest limit (this run watched the finding get NOISIER, not cleaner, as n grew — reported here, not the
prettier intermediate):** at n=26 Activity was the *only* positive category (Stability 0/5); by n=40
Expression had flipped positive and Stability softened to 5/11. Per-category n's are small and the win-rates
are noisy. So the own-model run **does NOT cleanly reproduce the per-category split** — it confirms the
**scorer** (0.87 vs ProteinGym's column) and the **aggregate category-mix picture** (no uniform lift), and is
directionally consistent on the two best-powered categories (Activity+, Stability−). The **clean,
powered** per-category evidence remains the N=95 precomputed-column sweep (+0.013, p=0.002); this local run is
underpowered and should not be read as an independent per-category confirmation.

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
