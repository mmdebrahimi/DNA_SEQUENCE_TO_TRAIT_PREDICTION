# DepMap drug-response decoder — the mechanistic-attribution WIN (Track A, 2026-07-02)

The cancer analog of the yeast substrate, and the win yeast could NOT get: on DepMap cell-line genotype ×
measured drug response, the decoder **recovers the textbook biomarker gene as the top predictor AND that
signal survives lineage de-confounding** for point-mutation mechanisms. Contrast yeast (copy-number/accessory
mechanisms, attribution inconclusive): DepMap's gene-level mechanisms are cleanly attributed.

## Setup
- **Substrate:** 563 cancer cell lines joined across CCLE damaging/COSMIC-hotspot mutations (1,566 genes
  mutated in ≥10 lines) × PRISM drug LFC (lower = more sensitive/killed) × 24 lineages (primary tissue). All
  free (DepMap 19Q4 + PRISM 19Q4 figshare). Join by DepMap ID.
- **Tests:** (a) univariate attribution — is the known biomarker the TOP gene for each drug? (b) **de-confounded
  single-gene attribution** — does the biomarker still separate response WITHIN lineage (center LFC by lineage
  mean, mutant vs WT)? This is the RIGHT test for a single-gene mechanism (a 1,566-gene multivariate Ridge
  dilutes it — its within-lineage r² is ~0, which is expected, not a failure).

## Headline (6 drugs, biomarker ground-truth)
| drug | target | top gene | attribution | **within-lineage t** (de-confounded) | verdict |
|---|---|---|---|---|---|
| **nutlin-3** | TP53 | **TP53** (t +10.0) | **TOP1** | **+8.81** (resistant in ALL 13 lineages) | **CLEAN WIN** |
| **vemurafenib** | BRAF | **BRAF** (t −6.8) | **TOP1** | **−2.96** (skin −1.33, thyroid −0.82) | **WIN** |
| **dabrafenib** | BRAF | **BRAF** (t −7.9) | **TOP1** | **−3.38** (skin −1.08) | **WIN** |
| **selumetinib** | BRAF | **BRAF** (t −6.9) + KRAS | **TOP1** | **−2.63** (skin −1.16) | **WIN** (MEKi; BRAF+KRAS both surface) |
| erlotinib | EGFR | NCOR1/KRAS | rank 712 | −1.25 | miss (EGFR); KRAS-resistance surfaces instead |
| gefitinib | EGFR | LPHN3… | rank 379 | −1.26 | miss (EGFR) |

## Why this is a real, de-confounded mechanistic win
1. **nutlin-3 → TP53** is the cleanest result in the whole project: TP53 is the top gene (t +10.0) AND
   TP53-mutant lines are resistant (+LFC) in **every one of 13 lineages** (within-lineage t +8.81, barely
   attenuated). Textbook: MDM2-inhibitor nutlin requires functional p53; TP53-mutant = resistant. Pan-lineage
   → NOT a lineage confound.
2. **BRAF → vemurafenib/dabrafenib/selumetinib**: BRAF is TOP1, and within **skin (melanoma)** BRAF-mutant
   lines are ~1.1–1.3 LFC MORE sensitive than BRAF-WT skin lines (de-confounded from lineage). Correctly
   ATTENUATED in colorectal (BRAF-mut CRC resists BRAFi via EGFR feedback — a known clinical fact the data
   reproduces). This is mechanism, not tissue.
3. **EGFR-TKIs (erlotinib/gefitinib) honestly MISS** (rank 379/712). EGFR-inhibitor sensitivity in pooled
   PRISM screens is weak/context-dependent (driven by EGFR amplification + activating mutations, rarer +
   lung-specific); KRAS surfaces instead as a RESISTANCE marker (t +4.0 — KRAS-mutant lines resist EGFR-TKIs,
   biologically correct). An honest mixed result on the harder drugs — not tuned away.

## Contrast with the yeast decoder (the point of the capstone)
- **Yeast:** within-clade predictive signal real but attribution INCONCLUSIVE — driven by accessory
  ORFs/copy-number/2μm-plasmid, NOT canonical named genes (CUP1 = copy-number, invisible to presence/absence).
- **DepMap:** gene-level point-mutation mechanisms → the decoder recovers the CANONICAL biomarker as top gene
  AND it survives de-confounding (TP53 pan-lineage; BRAF within-melanoma). **The feature type matches the
  mechanism type** — this is the deciding difference.

> **EGFR MISS RESOLVED 2026-07-02 (`wiki/depmap_multimodal_result_2026-07-02.md`):** the erlotinib/gefitinib
> EGFR misses below are NOT a decoder failure — they are a feature-MODALITY mismatch. Mutation-presence can't
> see EGFR amplification/expression. Adding CCLE copy-number + expression RECOVERS EGFR (erlotinib→expression
> ρ −0.115; gefitinib→copy-number ρ −0.153), de-confounded + correct-signed, plus ERBB2→lapatinib
> (expression −0.248). The general law holds across all 3 feature axes: match the feature type to the
> mechanism type.

## Honest scope
- The MULTIVARIATE within-lineage r² is ~0 — a full-gene-set decoder does NOT predict a novel line's response
  within lineage (single-gene mechanisms + small per-lineage n). The claim is **single-gene biomarker
  attribution, de-confounded**, NOT a general within-lineage response predictor.
- Damaging/hotspot mutation presence only (no copy-number/expression/fusions) — so EGFR-amplification-driven
  and expression-driven mechanisms are out of view (part of why EGFR-TKIs miss).
- 19Q4 release; frozen AMR surface byte-unchanged (read-only analysis).

## Verdict
DepMap is the substrate where the genotype→phenotype decoder finally gets a **clean, de-confounded,
biologically-ground-truthed mechanistic attribution** (TP53→nutlin pan-lineage; BRAF→BRAFi within melanoma).
Free, measured, deep, de-confoundable. Artifacts: `scripts/depmap_decoder.py` + `wiki/depmap_decoder_scores.json`.

## Recommended next (not blocking)
- Add copy-number + expression features → recover EGFR-amplification + expression-driven mechanisms.
- Continuous dose-response (AUC/IC50) labels instead of single-dose LFC.
- Formal FDR over all genes per drug (the top-gene test is descriptive).
