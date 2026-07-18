# Does ESM2 + OUR-OWN MSA-Transformer lift? — a correct scorer that does NOT realize the lift

**Date:** 2026-07-18 (final N=66) · **Model:** `esm_msa1b_t12_100M_UR50S` (fair-esm, CPU, wt-marginal,
depth-128) · **Script:** `scripts/msa_transformer_lift.py` · **Data:** `data/processed/msat_lift_checkpoint.jsonl`
· **Module:** `dna_decode/forward/msa_transformer.py`

## Question

The modality-hybrid finding (`forward_modality_hybrid_2026-07-17.md`) used ProteinGym's **precomputed**
`MSA_Transformer_ensemble` column (+0.013 paired vs ESM2-650M, win 66%, p=0.002 on N=95). That proves *a*
coevolution model lifts. This run computes MSA-Transformer scores **ourselves** on CPU — the **fast**
zero-shot form (wt-marginal, one forward pass; depth-128 subsample) — and asks whether **our** pipeline
(fetch → MSA-T → `rank_average_hybrid`) reproduces the lift.

## Result — the scorer is correct, but the FAST approximation does not lift

**(1) Our scorer is correct.** Our own MSA-Transformer scores **reproduce ProteinGym's `MSA_Transformer`
column at median Spearman 0.84** (N=66). The MSA parse → subsample → forward → variant-score pipeline is a
faithful implementation.

**(2) But the ESM2 (+) our-MSA-T hybrid does NOT lift at scale.** Over N=66 proteins the paired median Δ
(hybrid − ESM2) is **+0.0008, win-rate 34/66 = 52% — indistinguishable from no lift.** Per phenotype
category, the picture is weak and has **no clean structure**:

| phenotype | n | median Δ (hybrid − ESM2) | win-rate |
|---|---:|---:|---:|
| Stability | 22 | +0.0044 | 13/22 |
| Activity | 17 | +0.0111 | 9/17 |
| OrganismalFitness | 13 | −0.0336 | 3/13 |
| Expression | 10 | +0.0146 | 7/10 |
| Binding | 4 | −0.0021 | 2/4 |

The only category with a clear signal is **OrganismalFitness (−0.034)** — evolution-coupling is *negatively*
useful for the most polygenic phenotype, which is sensible. Everything else hovers near +0.01 (weak). The
"evolution → activity, structure → stability" split that looked clean at small n **dissolved**: Stability,
the clearest "evolution doesn't help" case at n=5 (0/5), is **positive** (13/22) at n=22.

## Why the fast approximation loses the lift

Our MSA-T rank-correlates 0.84 with ProteinGym's column but does **not** inherit its +0.013 hybrid lift. The
difference is model quality: ProteinGym's `MSA_Transformer_ensemble` is the **masked-marginal** (one forward
per position) over a **depth-384** subsample, ensembled; ours is the **wt-marginal** (single forward,
query-row logits) over **depth-128**. The 16% rank gap sits exactly in the hard variants where the ESM2
hybrid lift comes from — so the fast form reproduces the *ranking* but not the *complementarity*. This is the
**same shape as the site-independent floor** one level up: a cheaper version of the evolution modality is a
correct-but-non-lifting scorer. Realizing the +0.013 needs ProteinGym-grade masked-marginal MSA-T.

## Honest process note (the lesson, applied to me three times)

This run watched the finding get **weaker as n grew**, and I twice wrote a premature conclusion:
- **n=26:** Activity-only-positive, Stability 0/5 → looked like a clean phenotype-conditional lift.
- **n=40:** Expression flipped positive → "clean split" already false; I corrected to "noisy".
- **n=66 (true final, run actually completed here):** aggregate ≈0, Stability positive → **no lift, no clean
  structure.** The N=40 I committed was *another* mistaken-for-final snapshot (I misread a between-protein
  moment as the run finishing).

The rule — *report the final n, never the prettier intermediate* — is exactly what this cell measures, and it
bit the measurement itself. The **committed conclusion is this N=66 one.**

## What stands

- **The scorer + the pipeline are validated correct** (0.84 reproduction; `msa_transformer_table` +
  `fetch_msa` + `rank_average_hybrid` compose end-to-end on a novel sequence). The infrastructure is real and
  reusable.
- **The deployable lift is NOT delivered by the fast local MSA-T.** The +0.013 evolution lift remains a
  property of ProteinGym-grade masked-marginal MSA-T (the N=95 precomputed-column sweep) — not reproduced by
  the cheap wt-marginal form here.
- **Robustness fix:** `parse_a2m(max_rows=…)` bounds memory (a 100k-row MSA loaded whole OOM-killed two runs).
- Frozen decoder surface byte-unchanged (`verify_lock OK`).

## Consequence

The run-time evolution pipeline is **built and correct**, but to actually LIFT it needs the **masked-marginal
ensemble** MSA-T (heavier: L forwards/protein × ensemble), or GEMME. The cheap single-forward form is a
correct scorer that does not clear the bar — a genuine negative that saves the effort of shipping it as "the
lift". `msa_transformer_table` is left as the honest fast baseline + the interface for the heavier scorer.
