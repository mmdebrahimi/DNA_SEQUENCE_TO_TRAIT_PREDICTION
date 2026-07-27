# GEUVADIS Stage-2: the DNA-encoder arm (Enformer) LOSES to data-fit linear (2026-07-27)

**The empirical close of the organism-multimodal loop, on our own cohort, on a Kaggle T4 (free).**

## Result — matched comparison, same 12 genes, within-population

Enformer zero-shot cis-variant effects (ref→alt Δ at the gene TSS, 24 CAGE lymphoblastoid/GM12878/blood
tracks) aggregated by per-individual dosage → predicted expression, vs a data-fit linear predictor on the
IDENTICAL 12 top cis-eQTL genes × 50 TSS-proximal variants × 462 individuals:

| Predictor | mean within-pop \|Spearman ρ\| | beats linear |
|---|---|---|
| **Enformer DNA-encoder arm** (zero-shot) | **0.273** | — |
| **Elastic-net data-fit linear** (5-fold CV, out-of-sample) | **0.483** | Enformer 0/12 |
| Best-SNP (in-sample) | 0.479 | — |

**Enformer − linear = −0.21; Enformer wins on 0 of 12 genes.**

## Verdict

The DNA-encoder arm carries **real** cross-individual signal (0.273 ≫ permutation null ≈0.08) but
**decisively LOSES to a simple data-fit linear model** on the same genes, within-population. This is the
empirical confirmation — on the project's own GEUVADIS cohort, with a real Enformer GPU run — of:
- **Nat Genet 2023 / Variformer 2026**: sequence-to-expression models tie/lose to linear cross-individual.
- **Row 572**: the organism-multimodal DNA arm adds nothing over linear; acquiring dbGaP/UKB would not change
  this (the wall is the regime, not data access).

**The organism-multimodal question is now closed with a project-owned, GPU-validated number, not a citation.**

## Why this is the fair test (verify-in-batch)

The raw Enformer within-pop 0.273 is meaningless alone — the 12 genes are the STRONGEST cis-eQTLs (top |r|),
so any predictor scores high on them. The decisive number is the MATCHED comparison: same genes, same
variants, same individuals, same within-population evaluator. Against that, Enformer's zero-shot 0.273 vs
data-fit linear's 0.483 is unambiguous — the DNA arm underperforms. (Earlier Stage-1 ceilings — single-SNP
0.29 in-sample / elastic-net 0.19 out-of-sample — were on the ALL-genes set and are NOT comparable to these
top-12; the matched linear here is 0.483.)

## Honesty / scope

- Enformer is **zero-shot** (weights not fit to this data); elastic-net is **fit** (out-of-sample CV). That IS
  the multimodal question — does the sequence model's learned knowledge beat a data-fit linear model? No.
- 12 genes, chr1, TSS-centered 196 kb windows, CAGE track-mean readout. A larger gene panel would tighten the
  estimate but not change the direction (0/12, −0.21 is not marginal).
- N Kaggle T4 run: kernel `emanueleebrahimi/stage2-push-test` v4, ~15 min GPU, $0.

## Reproduce
- Bundle: `scripts/stage2_prep_kaggle_bundle.py`; kernel: `scripts/kaggle/stage2_enformer_kernel.py`.
- Matched linear baseline computed locally from the same bundle.
