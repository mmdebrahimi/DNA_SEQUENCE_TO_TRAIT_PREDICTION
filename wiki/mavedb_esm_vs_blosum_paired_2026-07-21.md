# ESM2 vs BLOSUM — paired on the MaveDB held-out set (2026-07-21)

**Status:** ✅ **the R2 regime claim ("molecular-property, fitness-aligned → LEARNED beats the cheap
baseline") validated PROSPECTIVELY.** On the same 84 leakage-free held-out MaveDB DMS assays, ESM2-650M beats
BLOSUM62 by a large, paired-significant margin. Frozen surface byte-unchanged.

## Result (84 paired held-out assays, |Spearman|)
| metric | value |
|---|---|
| median \|Spearman\| — ESM2-650M | **0.503** |
| median \|Spearman\| — BLOSUM62 | 0.221 |
| **median PAIRED delta (ESM2 − BLOSUM per assay)** | **+0.280** (mean +0.252) |
| ESM2 win-rate | **76/84 = 90%** (0 ties, 8 losses) |
| sign test (two-sided exact binomial, ties dropped) | **p = 5.0e-15** |

Method note (per the *paired-comparison, not difference-of-medians* lesson): the headline is the per-assay
**win-rate + median paired delta + sign test**, NOT `median(ESM) − median(BLOSUM)` over different assay sets.
Both scorers are direction-robust (`|Spearman|`) because MaveDB does not standardize per-assay score direction.

## Why it matters
- **Beats the cheap baseline, not just noise.** The earlier ESM2 packet compared only to a shuffled control
  (0.503 vs 0.017). This is the stronger, field-standard comparison: ESM2 vs the substitution-matrix baseline
  on the SAME assays, paired. A 90% win-rate at p=5e-15 is the actual evidence that the learned model earns
  its keep in R2 — the regime boundary the project defined (`plans/Trait_Decoding_Roadmap.md` regime lens;
  memory `feedback_g2p_decoder_regime_boundary`).
- **Leakage-free + confound-free.** The 84 assays' genes are NOT in the ProteinGym benchmark ESM2/the hybrid
  was tuned on, and R2 (a designed mutant library) has no population-structure confound. So the lift is real
  generalization, not benchmark overfit or lineage leakage.

## The 8 honest losses (the regime caveat showing itself)
ESM2 loses on 8/84, clustered on two assays: `00000281-{0-1,a-1,a-2}` and `00000441-{0-1,a-1}`. On these BOTH
scorers are weak (ESM2 0.05–0.08, BLOSUM 0.22–0.26) — the assay's functional readout is NOT
conservation/fitness-aligned, so neither a language model nor a substitution matrix tracks it well. This is
the R2 boundary condition (learned wins ONLY when the molecular property is fitness-aligned), not a scorer
bug — the losses are diagnostic, kept, not tuned away.

## Provenance
- ESM2 per-assay: `wiki/mavedb_holdout_esm2_2026-07-21.json` (Kaggle T4 run).
- BLOSUM per-assay: `wiki/mavedb_prospective_holdout_2026-07-21.json` (`blosum_scoring_proof`, CPU, all 86).
- Paired stats: `wiki/mavedb_esm_vs_blosum_paired_2026-07-21.json`.
- Held-out manifest + scoping: `wiki/mavedb_prospective_holdout_2026-07-21.md` +
  `wiki/mavedb_forward_cell_scoping_2026-07-21.md`.
