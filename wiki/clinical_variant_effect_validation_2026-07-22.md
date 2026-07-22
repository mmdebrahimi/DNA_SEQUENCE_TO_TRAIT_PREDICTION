# Clinical-significance validation of the R2 `forward` cell — actionable human genes (2026-07-22)

**Status:** ✅ first ClinVar (clinical) validation of the R2 molecular variant-effect decoder on the human
genes clinicians actually interpret. The R2 cell was previously validated only on DMS *fitness* (Spearman);
this asks the **VUS-resolution** question: does the decoder's variant-effect score separate ClinVar
**pathogenic vs benign** missense? Frozen AMR surface byte-unchanged (READ-only).

Epoch context: user-ratified strategic direction "deepen R2 into clinical human proteins" (2026-07-22). R2
(molecular property, fitness-aligned) is the one regime where the learned decoder wins, is species-agnostic
(so "human" is already true), free, and confound-free by construction — see
`feedback_g2p_decoder_regime_boundary` + `plans/Trait_Decoding_Roadmap.md` (Regime-first lens).

## The question + the three numbers

For each actionable gene with a free MaveDB functional DMS, join **MaveDB-DMS ⋈ ClinVar path/benign** on the
protein change and measure how well three predictors separate the clinical labels:

1. **DMS-functional AUROC** — the **fitness-alignment CEILING**. Does the wet-lab molecular assay *itself*
   separate ClinVar path/benign? High → the phenotype is fitness-aligned → R2 is the right regime.
2. **BLOSUM62 AUROC** — the deterministic, no-GPU, no-network decoder **FLOOR** (the forward cell's offline
   baseline).
3. *[AlphaMissense (full-proteome) + ESM2+ProSST hybrid AUROC]* — the deployable learned decoder; a
   Kaggle/full-AM **follow-up** (the committed `am_filtered.tsv` covers only ProteinGym-overlap variants).
   Named, not faked.

## Result (`scripts/clinical_variant_effect_validate.py`, real MaveDB + ClinVar E-utilities)

| Gene | joined (path/benign) | **DMS-AUROC (ceiling)** | **BLOSUM-AUROC (floor)** | in ProteinGym |
|---|---|---|---|---|
| **TP53** | 261 (171/90) | **0.996** | 0.707 | (name overlap) |
| **MSH2** | 297 (64/233) | **0.955** | 0.832 | yes |
| PTEN | 124 / 1 | — (single-class) | — | — |
| BRCA1 | 0 / 1 | — (single-class) | — | — |

**Two clean, honest findings:**
1. **The fitness-alignment ceiling is high (0.955–0.996)** for TP53 + MSH2 → these clinical phenotypes ARE
   fitness-aligned → R2 is the correct regime → a good molecular decoder *should* recover clinical
   pathogenicity. This is the regime question answered **per gene**, on clinical labels.
2. **The gap ceiling→floor (0.996→0.707, 0.955→0.832) is the headroom the LEARNED decoder captures** — the
   R2 thesis exactly (learned > deterministic in this regime). The BLOSUM floor is decent but well below the
   assay ceiling; AlphaMissense/hybrid is expected to close most of that gap (the deployable follow-up).

## Honest rails (load-bearing)

- **Class balance is gene-specific and reported, never faked.** BRCA1 pathogenicity is truncating-dominated
  (missense-path rare) and the Findlay SGE assay covers only the RING domain → 0 path joined; PTEN
  missense-path dominate ClinVar (124/1). Single-class genes are **AUROC-INAPPLICABLE** (a finding about the
  gene's clinical-missense distribution), not a decoder failure. TP53 + MSH2 carry both classes.
- **Tier = in-distribution clinical, NOT held-out.** These are canonical ProteinGym proteins, so the
  *hybrid* is not leakage-free on their DMS. BUT ClinVar path/benign labels are **independent of the
  DMS-fitness tuning** (the decoder was never fit to ClinVar), so the ClinVar-AUROC is a partly-independent
  clinical readout. The leakage-free tier remains the gene-agnostic MaveDB prospective holdout
  (`wiki/mavedb_prospective_holdout_full_2026-07-22.md`).
- **Partial circularity, flagged:** ClinGen may use DMS as PS3/BS3 functional evidence for a minority of
  recent ClinVar calls; most path/benign rest on segregation/population/clinical evidence, so DMS-AUROC is
  largely (not perfectly) independent. The near-ceiling TP53 0.996 should be read with this caveat.
- **Orientation is LABEL-FREE.** MaveDB does not standardize functional-score direction, so each DMS is
  oriented by its rank-agreement with BLOSUM62 (both → higher = preserved) — an orientation that never
  consults the clinical labels and so cannot inflate the AUROC toward them.

## Reproduce + extend

```
uv run python scripts/clinical_variant_effect_validate.py            # all genes (real network)
uv run python scripts/clinical_variant_effect_validate.py --gene TP53
```

ClinVar E-utilities results cache to `D:/dna_decode_cache/clinvar/eutils/<gene>.json` (reruns offline).
The gene set (`CLINICAL_GENES`) is a dict → adding a gene = one row (uniprot + MaveDB URN). Extension
candidates with likely ClinVar balance: MLH1 (Lynch, MaveDB 00001218), VHL (00000675), CALM1, SCN5A.

**Follow-up (named):** run the AlphaMissense (full human proteome) + ESM2+ProSST hybrid AUROC on the same
joined sets to place the deployable learned decoder between the floor and ceiling — the Kaggle/full-AM
pattern already established for the forward cell. Tests: `tests/test_clinical_variant_effect.py` (9, offline).
Artifact JSON: `wiki/clinical_variant_effect_validation_2026-07-22.json`.
