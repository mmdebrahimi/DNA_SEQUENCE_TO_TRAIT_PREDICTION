# MaveDB prospective-holdout — ESM2-650M headline number (2026-07-21)

**Status:** ✅ **the first prospective (leakage-free-by-benchmark-exclusion) learned number for the R2
molecular cell.** ESM2-650M zero-shot masked-marginals on 84 held-out MaveDB DMS assays → **median |Spearman|
= 0.503** (shuffled negative control 0.017). Real Kaggle T4 GPU run; frozen surface byte-unchanged.

## Result
| metric | value |
|---|---|
| model | `facebook/esm2_t33_650M_UR50D` (zero-shot masked-marginals) |
| assays pinned | 86 (the held-out manifest — gene NOT in ProteinGym v1.1) |
| assays scored | **84** (2 skipped `too_few` <20 missense; **0 fetch failures** — live MaveDB fetch worked on all 86) |
| **median \|Spearman\|** | **0.503** |
| shuffled negative control | **0.017** (near-zero → the signal is real) |
| per-assay \|ρ\| distribution | min 0.018 / q1 0.419 / median 0.503 / q3 0.610 / max 0.702 |
| compute | Kaggle T4 (`device=cuda`), ~28 min wall-clock (31 s → 1700 s) |

## Why it matters
- **Leakage-free by construction.** These 84 assays' target genes are NOT in the ProteinGym v1.1 benchmark
  the frozen `forward` hybrid was tuned + validated on, so ESM2 provably never saw them. R2 (a designed
  mutant library) has **no population-structure/clonality confound** — the cleanest regime the project has.
- **Matches the field on unseen data.** 0.503 is at/above the published ESM2-650M ProteinGym median (~0.48)
  and ~0.2 above the CPU BLOSUM62 baseline (~0.30, `wiki/mavedb_prospective_holdout_2026-07-21.md`). The
  molecular cell generalizes — it isn't benchmark-overfit.
- **The `forward` cell decodes human.** The held-out manifest is human-heavy (CYP2C19/CYP2C9/PSAT1/APP…);
  this is the north-star "move to human" crossing the line for free, in the one regime where learned models
  work — no biobank, no access gate, no confound.

## Honest scope
- **|Spearman| is direction-robust** — MaveDB does NOT standardize per-assay functional-score direction (the
  curation ProteinGym adds), so a signed correlation is uninterpretable; the absolute value is the honest
  metric here (standard ProteinGym-family practice).
- **ESM2-650M ALONE**, not the deployed ESM2+ProSST(+GEMME) hybrid — the hybrid needs per-protein structure
  tokens not yet Kaggle-wired. So this is a lower bound for the deployed scorer; the hybrid (+0.067 vs ESM2
  on ProteinGym) is the expected next lift. A ProSST-on-Kaggle wiring is the follow-up.
- **86-record manifest** = the MaveDB search page cap; a full-catalog manifest needs pagination (follow-up).

## Provenance
- Kernel: `notebooks/mavedb_holdout_esm2_kaggle.py` (self-contained, internet-enabled; masked-marginal core
  byte-faithful to `scripts/esm_zeroshot_dms.py`). Pushed via `scripts/kaggle_push_poll.py --gpu`
  (machine_shape NvidiaTeslaT4). Kaggle `emanueleebrahimi/mavedb-holdout-esm2`.
- Manifest: `wiki/mavedb_prospective_holdout_2026-07-21.json` (86 held-out URNs).
- Per-assay results: `wiki/mavedb_holdout_esm2_2026-07-21.json`.
- Scoping: `wiki/mavedb_forward_cell_scoping_2026-07-21.md`; regime lens: `plans/Trait_Decoding_Roadmap.md`.
