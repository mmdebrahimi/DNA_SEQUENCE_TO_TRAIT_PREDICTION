# Staph + Pneumo external validation on NCBI-PD — 4 more cells CONFIRMED (2026-07-21)

**Status:** ✅ **all 4 SCORED_ENDORSED + lineage-resolved** · extends the reusable no-compute NCBI-PD
substrate to organisms 3 + 4 (`--advance` follow-up (b)) · NON-FROZEN cells · **Frozen surface:** byte-unchanged.

## Result (raw + lineage-collapsed, both no-compute)

| organism | drug | n | R/S | RAW sens/spec | **LINEAGE sens/spec** | verdict |
|---|---|---|---|---|---|---|
| S. aureus | ciprofloxacin (gyrA/grlA QRDR) | 82 | 22R/60S | 1.00/0.967 | **1.00/0.976** | SCORED_ENDORSED |
| S. aureus | rifampin (rpoB RRDR) | 109 | 8R/101S | 0.875/0.980 | **0.857/0.986** | SCORED_ENDORSED |
| S. pneumoniae | erythromycin (erm/mef/msr) | 113 | 43R/70S | 0.977/0.971 | **1.00/0.962** | SCORED_ENDORSED |
| S. pneumoniae | tetracycline (tet(M)/tet(O)) | 113 | 10R/103S | 1.00/1.00 | **1.00/1.00** | SCORED_ENDORSED |

All 4 hold at the lineage-collapsed level (not clonally inflated). Pneumo erythromycin was previously
validated on the GPS-Poland cohort (acc 0.961) — this is a **second, larger, independent** confirmation.

## The reusable substrate now spans 4 organisms / 9 cells

`scripts/score_ncbipd_extval.py --organism {gono,campylobacter,staph,pneumo}` — a pure metadata join + score
(no Docker/Kaggle/assembly) with no-compute lineage-collapse via NCBI-PD's own SNP clusters:
- **gono** (3): ciprofloxacin, cefixime v0.1, penicillin v0.2
- **campylobacter** (2): tetracycline, gentamicin
- **staph** (2): ciprofloxacin, rifampin
- **pneumo** (2): erythromycin, tetracycline

## Determinant-format catch (verify-in-batch)

The pneumo cell's gene-presence rule keys (`ermb`/`mefa`/`tetm`, no parens) do NOT substring-match NCBI-PD's
`erm(B)`/`tet(M)` symbols → the raw call would false-negative every isolate (returns S). Fixed with a
**normalization adapter** (strip non-alphanumerics: `erm(B)` → `ermb`) in the scorer's pneumo config. Verified:
raw `erm(B)` → S (wrong); normalized `ermb` → R (correct). The staph cells (point-mutation dicts) needed no
adapter.

## Honest caveats

- **Provenance-disjoint** (measured-AST NCBI-PD isolates with assemblies) but **NOT methodology-independent**
  (same AMRFinderPlus + same cells).
- Both S. aureus + S. pneumoniae measured AST is partly **surveillance-linked** — the provenance-disjointness
  is cell-level (these organisms are not in the frozen surface's provdisjoint cohorts), not
  surveillance-ecosystem exclusion.
- NON-FROZEN cells; frozen decoder surface byte-unchanged (`verify_lock` OK).
