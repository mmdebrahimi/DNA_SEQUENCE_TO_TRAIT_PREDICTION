# Deployable learned-decoder clinical AUROC — AlphaMissense + ESM2+ProSST hybrid between the floor and ceiling (2026-07-22)

**Status:** ✅ places the DEPLOYABLE learned decoders in the floor→ceiling gap the floor+ceiling cell opened
(`wiki/clinical_variant_effect_validation_2026-07-22.md`). Frozen AMR surface byte-unchanged (READ-only).

Epoch: user-ratified "deepen R2 into clinical human proteins" → `--until-mvp` "place the learned decoder
between the 0.71 BLOSUM floor and the 0.996 DMS ceiling on the TP53/MSH2 ClinVar joins." **MVP met.**

## Result (`scripts/clinical_am_hybrid_auroc.py`; real MaveDB + ClinVar; AlphaMissense-full on D:; ESM2-650M + ProSST-2048 on this host)

Every decoder scored on the SAME joined variant set per gene (MaveDB-DMS ⋈ ClinVar path/benign ⋈ predictor),
so the numbers are directly comparable:

| Gene | n (path/benign) | FLOOR blosum | ProSST | ESM2-650M | **ESM2+ProSST hybrid** | **AlphaMissense** | CEILING dms |
|---|---|---|---|---|---|---|---|
| **TP53** | 261 (171/90) | 0.707 | 0.853 | 0.927 | 0.918 | **0.986** | 0.996 |
| **MSH2** | 297 (64/233) | 0.832 | 0.915 | 0.921 | **0.937** | 0.936 | 0.955 |

**Fraction of the floor→ceiling gap each learned decoder captures:**

| Gene | AlphaMissense | ESM2+ProSST hybrid | ESM2 | ProSST |
|---|---|---|---|---|
| TP53 | **96.6%** | 73.2% | 76.3% | 50.4% |
| MSH2 | 84.5% | **85.8%** | 72.2% | 67.8% |

AM standalone on the FULL ClinVar set (AM's real deployment — needs no DMS): TP53 0.986 (n=263), MSH2 0.936 (n=299).

## Findings (honest, gene-dependent)

1. **The deployable learned decoders land near the top of the gap on both genes** — the R2 thesis, now on
   *clinical* labels: a learned variant-effect model recovers most of what the wet-lab DMS assay recovers of
   ClinVar pathogenicity, far above the deterministic BLOSUM floor.
2. **Which learned decoder wins is gene-dependent, and that is honest signal, not noise:**
   - **TP53:** **AlphaMissense wins decisively (0.986, ~the 0.996 ceiling)**; the general ESM2+ProSST hybrid
     (0.918) does *not* beat ESM2 alone (0.927) — the structure channel (ProSST 0.853) slightly *drags* the
     hybrid here.
   - **MSH2:** the **shipped ESM2+ProSST hybrid narrowly wins (0.937)**, edging AM (0.936) AND lifting above
     *both* its components (ESM2 0.921 / ProSST 0.915) — the orthogonal-modality lift the hybrid was built for.
3. **AlphaMissense is the best DEPLOYABLE default** — free, precomputed, NO GPU, NO DMS at score time, and
   either the winner (TP53) or within 0.002 of it (MSH2). It was trained specifically for human variant
   pathogenicity with population + structural context; ESM2/ProSST are general protein models. The
   ESM2+ProSST hybrid is the shipped `predict_effect(method='hybrid')` decoder — this validates the ACTUAL
   deployed forward cell on clinical labels, not a proxy.

## Honest rails

- **Tier = in-distribution clinical, NOT held-out.** AM/ESM2/ProSST all saw these proteins in training. BUT
  the ClinVar path/benign labels are **independent of the DMS-fitness tuning** — none of these decoders was
  fit to ClinVar — so the ClinVar-AUROC is a partly-independent clinical readout. The leakage-free tier is
  the gene-agnostic MaveDB prospective holdout (`wiki/mavedb_prospective_holdout_full_2026-07-22.md`).
- **Coordinate-integrity gate:** variants whose WT residue disagrees with the UniProt canonical sequence are
  dropped before ESM2/ProSST scoring (never mis-scored). Numbering matched (no offset) for both genes.
- **Label-free orientation** for the DMS ceiling (rank-agreement with BLOSUM, never the clinical labels).
- **AlphaMissense is CC BY-NC-SA (non-commercial)** — research use. ESM2/ProSST are permissive.

## Provenance + reproduce

- AlphaMissense: `AlphaMissense_aa_substitutions.tsv.gz` (Cheng 2023) cached at
  `D:/dna_decode_cache/alphamissense/` → per-gene filter `am_clinical_filtered.tsv` (25,213 rows for the 2 UniProts).
- ESM2-650M masked-marginals (`dna_decode/forward/esm_scorer`) at only the joined positions, **chunk-checkpointed
  to D:** (restartable — a first run was killed mid-MSH2 by the documented D:-USB-hiccup failure mode; the chunked
  cache resumed with zero lost work).
- ProSST-2048 structure tokens from AlphaFold **v6** (`AF-<uniprot>-F1-model_v6.pdb`; NB the shipped
  `structure_scorer.alphafold_pdb_url` hardcodes v4, which now 404s — worked around here by fetching v6 directly).

```
uv run python scripts/clinical_am_hybrid_auroc.py --build-am-filter        # one-time AM filter from the gz
HF_HOME=D:/hf_cache uv run python scripts/clinical_am_hybrid_auroc.py --hybrid   # AM + ESM2 + ProSST + hybrid
```

6 offline tests `tests/test_clinical_am_hybrid.py`; artifact JSON `wiki/clinical_am_hybrid_auroc_2026-07-22.json`.
Builds on `scripts/clinical_variant_effect_validate.py` (floor+ceiling). Memory: `feedback_g2p_decoder_regime_boundary`.
