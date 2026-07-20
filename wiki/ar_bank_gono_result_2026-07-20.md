# AR Bank N. gonorrhoeae scoring — first gono numbers (1 clean win, 4 diagnostic rule failures)

**Date:** 2026-07-20 · **Status:** SCORED (curated NON-FROZEN gono rule vs CDC measured S/I/R) · **Cohort:**
52 gono isolates (CORRECTED from the earlier "94" genus miscount), 49 FREE, **zero leaks** (BioSample +
throttle-proof assembly-base) · **Rule:** `organism_rules/neisseria_amr.call_ng_amr` (CURATED_NONFROZEN,
sourced from Pathogenwatch 485.toml) · **Frozen decoder surface:** byte-unchanged.

## Results (strict, per drug)

| Drug | scored n | TP/FN/TN/FP | sens | spec | verdict |
|---|---|---|---|---|---|
| **ciprofloxacin** | 25 | 22/0/3/0 | **1.00** | **1.00** | ✅ **ENDORSED** (two-sided, perfect) |
| **azithromycin** | 20 | 0/0/20/0 | — | **1.00** | ✅ **ENDORSED** (specificity; no R in scored set) |
| cefixime | 23 | 20/0/0/3 | 1.00 | 0.00 | ⚠ over-call (sens perfect, spec 0 — 3 FP) |
| ceftriaxone | 25 | 0/0/0/25 | — | 0.00 | ❌ over-call (all 25 S → FP) |
| penicillin | 24 | 0/24/0/0 | 0.00 | — | ❌ under-call (all 24 R → FN) |
| tetracycline | 21 | 0/21/0/0 | 0.00 | — | ❌ under-call (all 21 R → FN) |

(~49% indeterminate per drug = the same NCBI GCA-unavailability seen on E. coli/Klebsiella; ~25 of 49 FREE
scored. Endorsement gate = spec ≥ 0.85 on the powered provenance-disjoint set, per the AMR-Portal cell.)

## Honest interpretation — the validation DISCRIMINATES

**Two rules validated:**
- **Ciprofloxacin — perfect, two-sided (sens 1.0 / spec 1.0).** The `gyrA` QRDR (Ser91/Asp95) rule is the
  clean single-marker case (Eyre 32/32) and it holds on independent CDC gono MICs: 22/22 R called R, 3/3 S
  called S. This is the strongest single-cell result in the whole AR-Bank arm — genuinely two-sided.
- **Azithromycin — spec 1.0.** The 23S rRNA rule correctly calls S when no macrolide-target mutation is
  present (20/20). No azithromycin-R in the scored subset, so sensitivity is untested here.

**Four rules exposed as WRONG — each with a clear mechanism (the value of an external test):**
- **Cefixime + ceftriaxone — penA OVER-call.** "Any curated penA point → R" is too aggressive for the ESCs:
  penA mosaics are *common* in gono and confer *reduced susceptibility*, not reliably full R — and
  ceftriaxone-R is globally rare. So near-every isolate carries a penA point → called R → **specificity
  collapses** (cefixime 0/3, ceftriaxone 0/25). The rule needs the *specific high-level mosaic combinations*
  (e.g. penA-60 A501/G545S + porB + mtrR context), NOT any penA point. The R3 real-symbol fix (match any
  penA point) that helped SENSITIVITY is exactly what sinks SPECIFICITY here — the honest trade the MIC
  data reveals.
- **Penicillin + tetracycline — chromosomal UNDER-call.** "Plasmid gene only" (blaTEM / tet(M)) MISSES
  gono's dominant CHROMOSOMAL resistance: penicillin-R is largely penA/mtrR/ponA-mediated (not blaTEM), and
  tet-R is largely rpsJ V57M / mtrR efflux-mediated (not tet(M)). My "primary determinant only, demote
  chromosomal to accessory" discipline — which was correct for E. coli (I reused the frozen validated rule)
  — is BACKWARDS for gono pen/tet, where the chromosomal path IS primary. Every R isolate here is R via the
  demoted accessory → 24/24 + 21/21 FN.

## Why this is a good result, not a bad one

This is the "suspect the rule when the metric is degenerate" pattern working as designed. A curated rule
built from a catalogue's *primary-determinant* framing met an independent measured-MIC cohort, and the cohort
**pinpointed exactly where the framing is wrong** — over-calling ESCs (penA presence ≠ resistance) and
under-calling chromosomal pen/tet. Two cells endorse cleanly (cipro perfect, azithromycin spec); four have a
concrete, mechanistic v0.1 fix. NONE of this touches the frozen decoder surface (CURATED_NONFROZEN).

## v0.1 fixes (named, deferred)

1. **Cef/ceftriaxone:** replace "any penA point → R" with the specific high-level mosaic *combinations*
   (485.toml's multi-member penA+porB+mtrR rules) — accept the sens↓ for the spec↑ the MIC data demands.
2. **Penicillin:** promote chromosomal penA/mtrR from accessory → primary (blaTEM stays high-level).
3. **Tetracycline:** promote rpsJ V57M / mtrR from accessory → primary (tet(M) stays high-level).
Each re-scored against this same cohort (spec≥0.85 / sens≥0.85 endorsement).

## Honest scope

- Curated NON-FROZEN rule (not the deployed surface); genotype = AMRFinder `-O Neisseria_gonorrhoeae`;
  phenotype = CDC's own gonococcal S/I/R (correct species breakpoints). Provenance-disjoint (0 leaks,
  BioSample + assembly-base). NOT methodology-independent (rule curated from 485.toml). One-sided for 4 drugs
  (resistance-enriched bank); cipro + cefixime are two-sided.
- Artifacts: `wiki/ar_bank_gono_validation_<drug>_2026-07-20_*.json`. Frozen surface `verify_lock` unchanged.
