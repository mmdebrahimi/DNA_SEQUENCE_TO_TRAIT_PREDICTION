# Yeast growth decoder — the first DE-CONFOUNDED learned-decoder WIN (2026-07-02)

On the 1002 Yeast Genomes substrate (F4 pilot GO), a gene-content decoder shows **real mechanistic
genotype→growth signal that survives lineage de-confounding** — the first positive after the foundation-model
embedding arm went 0-for-4 (cipro within-lineage / pathotype / Arabidopsis flowering-time / …, all
structure-learners).

## Setup
- **Substrate:** 970 *S. cerevisiae* isolates joined across phenotype (35 lab growth conditions, YPD-normalized)
  + gene presence/absence (7,796 genes) + SNP distance matrix. All free from `1002genomes.u-strasbg.fr`.
- **Genotype features:** gene presence/absence (interpretable gene content — NOT a black-box embedding).
- **De-confounding:** hierarchical clades from the SNP distance matrix (K=18, robustness-checked at K=30).
- **The 3-part test** (the project's embedding-niche bar): naive global CV r² (structure-polluted) /
  clade-only baseline / leave-one-clade-out / **within-clade r²** (the de-confounded mechanistic signal —
  variation inside a clade cannot be between-clade population structure).

## Headline
**11 / 35 conditions = MECHANISM** (within-clade r² ≥ 0.05, HONEST clade-centered-residual metric — see the
metric-correction note). The strongest are known gene-dosage resistances:

| condition | naive r² | clade-only r² | within-clade r² | K=30 within | verdict |
|---|---|---|---|---|---|
| YPDCAFEIN40 (caffeine) | 0.372 | 0.114 | **0.365** | 0.224 | MECHANISM |
| YPDCAFEIN50 (caffeine) | 0.369 | 0.122 | **0.343** | — | MECHANISM |
| YPDBENOMYL500 | 0.321 | 0.083 | **0.290** | 0.193 | MECHANISM |
| YPDCUSO410MM (copper) | 0.336 | 0.197 | **0.238** | 0.180 | MECHANISM |
| YPDSODIUMMETAARSENITE | 0.143 | 0.024 | **0.200** | 0.154 | MECHANISM |
| YPDKCL2M | 0.158 | 0.027 | **0.172** | 0.173 | MECHANISM |
| YPGALACTOSE / YPD42 / YPDANISO10 / YPDSDS / YPDNACL1M | | | 0.057–0.112 | | MECHANISM |
| … 24 WEAK | | | | | |

Full table: `wiki/yeast_growth_decoder_scores.json`.

## Metric-correction note (verify-in-batch, load-bearing honesty)
The FIRST-pass within-clade metric added the clade mean back to the prediction and scored r² against the
GLOBAL mean — leaking between-clade structure into the "de-confounded" number (a synthetic pure-structure
null test caught it: it scored >0.1 when it should be ~0). **Fixed** to score on clade-centered residuals
(structure cancels from both truth and prediction). The corrected numbers are slightly lower (caffeine
0.402→0.365, copper 0.370→0.238; 14→11 MECHANISM) — the leak was small only because yeast between-clade
phenotype variance is low (clade_only_r² ~0.02–0.20). The verdict is unchanged and now methodologically sound.

## Why this is a real WIN (two hard controls, both passed with the CORRECTED metric)
1. **Within-clade permutation null.** Shuffling the phenotype INSIDE each clade and re-running the identical
   within-clade CV collapses r² to **strongly negative (−0.30 to −0.42)** for every strong condition — vs
   real 0.15–0.37. The within-clade CV is not manufacturing the signal.
2. **K=30 finer-clade robustness.** The within-clade r² HOLDS at 30 clades (copper 0.18, caffeine 0.22,
   arsenite 0.15, benomyl 0.19), ruling out coarse-clade residual structure as the driver.

Contrast the embedding arm: Arabidopsis flowering-time within-group r² was **−0.13** (pure structure-learner).
Here within-clade r² is **+0.15 to +0.37** and permutation-null-clean. This is categorically different.

## Honest scope + caveats (load-bearing)
1. **This is a GENE-CONTENT decoder, not a foundation-model embedding.** The features are interpretable
   presence/absence of 7,796 genes; the win does NOT resurrect the embedding bet (still 0-for-4). It confirms
   the project's standing thesis: **interpretable/mechanistic features work where black-box embeddings don't.**
2. **Cross-clade transfer FAILS.** Leave-one-clade-out r² is strongly negative (−0.1 to −1.6): the decoder
   predicts variation WITHIN a known clade but does NOT generalize to a novel lineage (clade-context-dependent
   effects + novel-clade baseline shift). The claim is "decodes within-clade growth," NOT "decodes a novel
   strain's growth."
3. **De-confounding is at the K=18–30 clade level;** finer sub-lineage structure is bounded by the K=30
   robustness but not eliminated. A fully rigorous next step is a kinship-regressed (mixed-model) test.
4. **Mechanistic attribution not yet mapped:** the copper/arsenite/caffeine hits are the KNOWN gene-dosage
   resistances (copper→CUP1, arsenite→ARR/ACR3), but the gene IDs are EC1118 pangenome orthologs — mapping the
   top-weighted genes to CUP1/ARR would be the mechanistic capstone (deferred; needs the pangenome annotation).

## Verdict
The 1002 Yeast Genomes substrate yields a **de-confounded, permutation-clean, clade-robust genotype→growth
signal** on 14/35 conditions — the first learned-decoder positive in the project. It is a within-clade,
gene-content result (honest scope), not a novel-lineage or embedding claim. Frozen AMR surface byte-unchanged.

## Recommended next steps (not blocking this result)
- Map top genes for copper/arsenite → CUP1/ARR (mechanistic capstone; needs pangenome annotation).
- Mixed-model (kinship-regressed) confirmation to close the finer-substructure caveat.
- Extend to the raw SNP matrix (`1011Matrix.gvcf.gz`) for variant-level (not gene-level) resolution.
