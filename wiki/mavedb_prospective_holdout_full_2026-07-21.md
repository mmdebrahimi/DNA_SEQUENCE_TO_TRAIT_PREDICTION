# MaveDB held-out manifest — FULL catalog (paginated) (2026-07-21)

**Status:** ✅ MVP reached. Paginated the whole MaveDB catalog (28 pages / 2797 score sets) → the leakage-free
held-out manifest grows **86 → 2569 assays** (30×). More than confirmation: it tightens any future estimate,
gives a big Kaggle set, and surfaces the R1 pharmacogene seed. Frozen surface byte-unchanged.

## What pagination delivered
| metric | single-page (before) | FULL catalog (now) |
|---|---|---|
| MaveDB score sets fetched | 100 (page cap) | **2797** (28 pages via `offset`) |
| held-out assays (protein_coding, NOT in ProteinGym) | 86 | **2569** |
| — human | 58 | **1072** |
| — published ≥2024 | 37 | 819 |
| distinct genes | ~52 | **1036** (603 human) |
| ProteinGym leakage (held-out gene in the benchmark) | 0 | **0** (dedup clean at scale) |

Pagination mechanics (R2 pre-bar check): the search endpoint caps a page at 100 (`limit`>100 → HTTP 422); the
only lever is `offset`, which returns disjoint pages (offset:100 → 0 overlap with base). `_fetch_score_sets(
all_pages=True)` loops offset until a short page. The full run is namespace-separated
(`mavedb_prospective_holdout_full_*.json`) so it does NOT clobber the 86-assay manifest that carries the
committed `blosum_scoring_proof` the paired-comparison artifact reads (shared-key-overwrite trap avoided).

## Why it matters (feeds BOTH lanes)
- **R2 (molecular cell):** a 2569-assay leakage-free held-out corpus (1072 human). The ESM2 0.503 headline was
  on 84; the enlarged set is the substrate for a tighter, larger prospective number (a full-manifest Kaggle
  ESM run — the kernel embeds URNs, so it regenerates with 2569; ~1h+ T4).
- **R1 (pharmacogenomics catalog cell):** the held-out human set contains the exact seed genes —
  **CYP2C19, CYP2C9, NUDT15, G6PD, VKOR** (+ clinical BRCA1, TP53, MLH1, KRAS, PMS1). These are curated-catalog
  drug-response genes with real DMS scores — the validation substrate for the R1 human cell.

## Honest scope
- The deliverable is the ENLARGED MANIFEST (the leakage-free held-out set), not scores. Full scoring is the
  follow-up: BLOSUM over 2569 assays is a multi-hour local API job; the ESM headline needs a full-manifest
  Kaggle run. The prior 84-assay numbers (ESM2 0.503, 90% BLOSUM beat) stand as the current prospective result.
- Manifest is `wiki/mavedb_prospective_holdout_full_2026-07-21.json`; build/rebuild:
  `uv run python scripts/mavedb_prospective_holdout.py --all-pages`.
