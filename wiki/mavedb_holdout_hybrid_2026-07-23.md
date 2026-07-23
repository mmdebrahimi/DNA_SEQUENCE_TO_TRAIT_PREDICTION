# Full leakage-free ESM2+ProSST hybrid at scale — held-out MaveDB (2026-07-23)

**Status:** ✅ the Kaggle-GPU piece the data hunt named, LANDED and then WIDENED. First at-scale,
**leakage-free** run of the shipped ESM2-650M + ProSST-2048 modality hybrid on held-out MaveDB DMS assays
(genes NOT in the ProteinGym benchmark the hybrid was tuned on). N grown 38 → **76** to test whether the
hybrid's small margin over ProSST was real. Frozen AMR surface byte-unchanged (READ-only).

## Result (Kaggle T4, N=76 held-out human assays)

| Decoder | median \|Spearman\| |
|---|---|
| ESM2-650M (sequence) | 0.538 |
| ProSST-2048 (structure) | 0.596 |
| **ESM2+ProSST hybrid** | **0.602** |

Comparators: ESM2 full holdout **0.478** · AlphaMissense held-out **0.502** (`wiki/mavedb_am_holdout_2026-07-23`).

**PAIRED per-assay (the correct statistic):**

| Comparison | paired wins | median paired delta | sign-test p |
|---|---|---|---|
| **hybrid > ESM2** | **70/76 (92%)** | **+0.063** | — |
| **hybrid > ProSST** | **52/76 (68%)** | **+0.011** | **0.0009** |
| ProSST > ESM2 | 58/76 (76%) | +0.052 | — |

## Findings

1. **The hybrid beats BOTH components paired, and beating ProSST is now statistically significant** — 52/76
   wins over the structure model, one-sided sign-test **p = 0.0009**. At the first N=38 run this margin was
   +0.006 (26/38, borderline); doubling N confirmed it is **real, not noise**. This is the decisive answer
   to "is the naive rank-average hybrid genuinely better than its strongest single component, or a wash?" —
   genuinely better, on a leakage-free held-out set.
2. **The difference-of-medians tension resolved with more N.** At N=38, ProSST's *median* (0.601) exceeded
   the hybrid's (0.586) even though the hybrid won paired — the documented difference-of-medians trap. At
   N=76 the hybrid median (0.602) now *also* edges ProSST (0.596), so median and paired agree. The paired
   statistic was right at N=38; the median caught up. (Lesson reaffirmed:
   `feedback_paired_comparison_not_difference_of_medians`.)
3. **Ranking on held-out DMS fitness: hybrid > ProSST > ESM2 > AlphaMissense > ESM2-full-holdout > BLOSUM.**
   Structure (ProSST 0.596) is the strongest single modality; the sequence⊕structure hybrid is the best
   overall. This CONFIRMS the ProteinGym modality-hybrid finding on data the methods never saw.
4. **Structure conditioning verified, not assumed** — the in-run self-check scores one assay with REAL vs
   SHUFFLED structure tokens: correlation **0.4147** (≪ 1.0), proving `ss_input_ids` reaches the model.

## Honest scope

- **N=76** held-out human assays: gene NOT in ProteinGym, UniProt + AlphaFold, seq ≤ 1022, non-giant
  structure, and **verified sequence identity ≥ 0.95** between the assay sequence and the AlphaFold structure
  at the MaveDB offset (all 76 ≥ 0.958). The identity gate now runs at manifest BUILD time
  (`build_holdout_hybrid_manifest.struct_identity`), before quantize — a length fit is not an alignment.
- Structure tokens quantized LOCALLY and sliced to the assay region `tokens[offset:offset+L]`; only
  pre-quantized tokens shipped to Kaggle (no `torch_geometric` on Kaggle).
- The hybrid is an equal-weight rank-average — no fitting, no calibration (the deployability class of the
  shipped `rank_average_hybrid`).
- |Spearman| is direction-robust (MaveDB does not standardize per-assay score direction).

## The debugging arc that produced the first number (kept for the record)

The first run of this pipeline took 4 attempts because ProSST scored ~0.04–0.08 (chance) on Kaggle while
scoring 0.60 locally. Root cause: `trust_remote_code` runs the model repo's `modeling_prosst.py` against the
**installed** transformers, and Kaggle's newer transformers broke structure conditioning; pinning
`transformers==5.13.0` fixed it (the decoder-tie and remote-code-revision fixes were verified-correct and NOT
the cause). A `prosst_degenerate` guard (median < 0.10 → "INVALID, do NOT publish") blocked the three broken
runs from ever being reported. Full arc + the three failed-run artifacts (`*_run1/2/3`) retained as evidence.
Reusable lesson: `feedback_trust_remote_code_library_version_drift`.

Reproduce: `scripts/build_holdout_hybrid_manifest.py` (manifest + local quantize + identity gate) →
`notebooks/mavedb_holdout_hybrid_kaggle.py` (inject manifest, push via `scripts/kaggle_push_poll.py --gpu`).
Artifact: `wiki/mavedb_holdout_hybrid_2026-07-23.json` (N=76, with computed paired stats + sign test).
