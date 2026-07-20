# AR Bank S. aureus Levofloxacin — FQ-QRDR rule confirmed 10/10 (perfect, underpowered)

**Date:** 2026-07-20 · **Status:** SCORED, perfect but UNDERPOWERED · **Cell:** `organism_rules/
staphylococcus_amr.call_sa_ciprofloxacin` (gyrA/grlA QRDR; NON-FROZEN) via the generalized registry ·
**Frozen surface:** byte-unchanged (verify_lock OK).

## Result

| metric | value |
|---|---|
| scored n | 10 (6 R / 4 S) |
| sensitivity | **1.00** (6/6 R) |
| specificity | **1.00** (4/4 S) |
| accuracy | **1.00** (0 errors) |
| powering | **UNDERPOWERED** (4 S < min-per-class 5) |

The Staph fluoroquinolone rule (gyrA / grlA QRDR point mutation → R) matched CDC's measured **Levofloxacin**
S/I/R on every disjoint isolate it could score — a clean two-sided independent confirmation. (Levofloxacin
is scored against the cell's ciprofloxacin rule: both are fluoroquinolones driven by the same gyrA/grlA QRDR
mechanism — a documented same-mechanism cross-drug application.)

## Why UNDERPOWERED + only 1 drug (honest, both are data-shape ceilings)

- **Only Levofloxacin is scorable.** The AR Bank S. aureus panel has NO Rifampin (so `call_sa_rifampicin` has
  no label) and NO ciprofloxacin — Levofloxacin (FQ) is the only cell-drug with a matching label. The other
  panel drugs (cefoxitin/oxacillin) hit the mecA/oxacillin label-confounding wall (already flagged
  `LABEL_CONFOUNDED` in the report card).
- **128 → 10 scored.** Of 128 S. aureus, only 28 have a downloadable assembly; the hardened preflight
  excluded **15 tuning leaks** (some AR Bank S. aureus ARE in the decoder's tuning set — the resolution-free
  assembly-base check caught them); of the ~13 disjoint-FREE, 3 failed download/smoke → 10 scored. The S side
  (4) lands just below the 5-per-class powering floor.

## Context — the AR-Bank Gram-positive/fungal assembly ceiling

This is the third small-N-but-perfect result in this arc (C. auris fungal 8/8, S. aureus 10/10), all capped
by the same cause: the CDC AR Bank deposits downloadable assemblies for the priority Gram-negatives
(E. coli/Klebsiella/gono ~half FREE, all POWERED) but SRA-reads-only for Enterococcus (0 FREE, unscorable)
and much of C. auris/S. aureus. The curated cells score PERFECTLY on the disjoint isolates they can assemble;
a powered Gram-positive/fungal validation needs SRA-read-mapping (`assemble_sra_cohort --method map`), not the
assembly-download path.

## Scope

- CURATED_NONFROZEN (Staph FQ-QRDR rule) via the generalized `ar_bank_registry` harness; genotype = AMRFinder
  `-O Staphylococcus_aureus`; phenotype = CDC Levofloxacin S/I/R; provenance-disjoint (15 leaks excluded).
- Artifact: `wiki/ar_bank_staphylococcus_aureus_validation_levofloxacin_staph_2026-07-20_*.json`.
  Frozen surface `verify_lock` OK.
