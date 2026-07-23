# Full leakage-free ESM2+ProSST hybrid at scale — held-out MaveDB (2026-07-23)

**Status:** ✅ the Kaggle-GPU piece the data hunt named, LANDED after a 4-run debugging arc. First at-scale,
**leakage-free** run of the shipped ESM2-650M + ProSST-2048 modality hybrid on held-out MaveDB DMS assays
(genes NOT in the ProteinGym benchmark the hybrid was tuned on). Frozen AMR surface byte-unchanged (READ-only).

## Result (Kaggle T4, N=38 held-out human assays)

| Decoder | median \|Spearman\| |
|---|---|
| ESM2-650M (sequence) | 0.519 |
| **ProSST-2048 (structure)** | **0.601** |
| **ESM2+ProSST hybrid** | 0.586 |

Comparators: ESM2 full holdout **0.478** · AlphaMissense held-out **0.502** (`wiki/mavedb_am_holdout_2026-07-23`).

**PAIRED per-assay (the correct statistic — medians are over different assays, not a lift):**

| Comparison | paired wins | median paired delta |
|---|---|---|
| **hybrid > ESM2** | **34/38 (89%)** | **+0.060** |
| **hybrid > ProSST** | **26/38 (68%)** | **+0.006** |
| ProSST > ESM2 | 30/38 (79%) | +0.052 |

## Findings

1. **The hybrid beats BOTH components paired** — 34/38 vs ESM2 (+0.060) and 26/38 vs ProSST (+0.006). This
   CONFIRMS the ProteinGym modality-hybrid finding on a **leakage-free** held-out set at scale: a naive
   rank-average of orthogonal modalities (sequence ⊕ structure) beats either alone.
2. **Do not read the medians as a lift.** ProSST's median (0.601) is *above* the hybrid's (0.586), yet the
   hybrid **wins paired** — the medians are computed over different assays. This is exactly the
   difference-of-medians trap the project documented before (`feedback_paired_comparison_not_difference_of_medians`);
   the paired delta is authoritative. The hybrid's margin over ProSST is real but **small** (+0.006).
3. **Structure is the stronger single modality here** — ProSST 0.601 beats ESM2 0.519 (30/38 paired). It also
   beats every prior held-out number in the project (ESM2 0.478, AM 0.502), making ProSST the strongest
   single deployable predictor measured so far on held-out DMS fitness.
4. **Structure conditioning is verified, not assumed** — the in-run self-check scored one assay with REAL vs
   SHUFFLED structure tokens: correlation **0.4147** (≪ 1.0), proving `ss_input_ids` genuinely reaches the
   model and changes scores.

## The 4-run debugging arc (why the first three numbers were never published)

| Run | fix applied | ProSST median | verdict |
|---|---|---|---|
| 1 | — | 0.041 | decoder weight loaded RANDOM (`MISSING \| newly initialized`) |
| 2 | force-tie (verified True) | 0.053 | tie was NOT the cause |
| 3 | + pinned ProSST remote-code revision | 0.082 | revision was NOT the cause |
| **4** | **+ pinned transformers==5.13.0** | **0.601** | ✅ **the library version WAS the cause** |

- A `prosst_degenerate` guard (median < 0.10 → stamp artifact + print "hybrid number is INVALID, do NOT
  publish") blocked runs 1–3 from ever being reported as a result. Without it, 0.34 would have shipped twice.
- The decisive diagnostic was **local**: the SAME gene + structure tokens (MTHFR) scored Spearman **0.4025**
  locally, proving ProSST worked and isolating the fault to the Kaggle path.
- Root cause: `trust_remote_code` runs the model repo's own `modeling_prosst.py` against the INSTALLED
  transformers; Kaggle's newer transformers broke that interaction (structure conditioning not applied).
  Pinning transformers to the local working version (5.13.0) fixed it. Note torch still differs
  (2.10.0+cu128 Kaggle vs 2.12.1 local) and is evidently not load-bearing.

## Honest scope

- **N=38** held-out human assays: gene NOT in ProteinGym, UniProt + AlphaFold available, seq ≤ 1022 (ESM2
  context), non-giant structure, and **verified sequence identity ≥0.95** between the assay sequence and the
  AlphaFold structure at the MaveDB offset (7 genuinely misaligned assays were dropped — PSD95 0.06, OSTF1
  0.05, AβB42 0.12, HNRNPUL1 0.09, KRAS 0.89, HECTD1 0.92, NKX3-1 0.93).
- Structure tokens were quantized LOCALLY (the `torch_geometric` path) and sliced to the assay region
  `tokens[offset:offset+L]`; only pre-quantized tokens shipped to Kaggle.
- |Spearman| is direction-robust (MaveDB does not standardize per-assay score direction).
- The hybrid is an equal-weight rank-average — no fitting, no calibration (the deployability class of the
  shipped `rank_average_hybrid`).

Reproduce: `scripts/build_holdout_hybrid_manifest.py` (manifest + local quantize) →
`notebooks/mavedb_holdout_hybrid_kaggle.py` (inject manifest, push via `scripts/kaggle_push_poll.py --gpu`).
Artifacts: `wiki/mavedb_holdout_hybrid_2026-07-23.json` (with computed paired stats) + the three failed-run
artifacts kept as evidence (`*_run1_broken_prosst`, `*_run2_tied_still_degenerate`, `*_run3_pinned_still_degenerate`).
