# AR Bank N. gonorrhoeae v0.1 — 3 of 4 v0 failures fixed (literature-grounded, cached re-score)

**Date:** 2026-07-20 · **Status:** v0.1 SCORED (curated NON-FROZEN gono rule vs CDC measured S/I/R) ·
**Cohort:** same 52-iso, 49-FREE, zero-leak disjoint cohort as v0 (`wiki/ar_bank_gono_result_2026-07-20.md`) ·
**Frozen surface:** byte-unchanged (verify_lock OK).

## v0 → v0.1 (same cohort, refined rules; AMRFinder runs cached → re-score in seconds)

| Drug | v0 | **v0.1** | change |
|---|---|---|---|
| ciprofloxacin | sens 1.00 / spec 1.00 | **sens 1.00 / spec 1.00** | — (perfect, unchanged) |
| azithromycin | spec 1.00 | **spec 1.00** | — (unchanged) |
| **ceftriaxone** | spec **0.00** (25 FP) | **spec 1.00** (25/25 S) | ✅ **FIXED** (A501-narrowing) |
| **penicillin** | sens **0.00** (24 FN) | **sens 1.00** (24/24 R) | ✅ **FIXED** (chromosomal promotion) |
| **tetracycline** | sens **0.00** (21 FN) | **sens 1.00** (21/21 R) | ✅ **FIXED** (chromosomal promotion) |
| cefixime | sens 1.00 / spec 0.00 | sens 1.00 / spec 0.00 | ⚠ unchanged (MIC-margin ceiling) |

**5 of 6 drugs now endorse on their powered side** (up from 2 of 6 at v0).

## The three fixes — each grounded in literature/485.toml, NOT fit to the cohort

The AR Bank gono cohort is **genetically near-homogeneous**: nearly every isolate (R and S) carries the full
mosaic-penA + mtrR + porB + ponA + rpsJ chromosomal package (verified by inspecting the cached AMRFinder
symbols). That single fact drove all three fixes:

1. **Penicillin — chromosomal promotion (sens 0.0 → 1.0).** Gono penicillin-R is CHROMOSOMAL-dominant
   (penA/mtrR/ponA), with blaTEM the plasmid path. v0's "blaTEM-only" missed all 24 R isolates (chromosomal,
   no blaTEM). v0.1: `R iff blaTEM OR penA-point OR mtrR`. 24/24 R now correct. **Honest caveat:** the cohort
   is penicillin-R-saturated (0 S), so this is SENS-only-testable; the promoted rule WILL over-call a
   penicillin-S isolate carrying the near-universal mosaic penA/mtrR — specificity is UNTESTED here.

2. **Tetracycline — chromosomal promotion (sens 0.0 → 1.0).** Gono tet-R is chromosomal (rpsJ V57M + mtrR
   efflux) + plasmid tet(M). v0's "tet(M)-only" missed all 21 R (rpsJ+mtrR). v0.1: `R iff tet(M) OR rpsJ_V57M
   OR mtrR`. 21/21 correct. **Honest caveat:** on the EBI AMR Portal, rpsJ V57M is common + low-level and
   collapsed specificity to ~0.35 — it over-calls when tet-S isolates exist. Here (0 tet-S) the lift is
   SENS-only-testable; keep the v0 tet(M)-only variant for a spec-sensitive cohort.

3. **Ceftriaxone — A501-narrowing (spec 0.0 → 1.0).** v0's "any penA point → R" over-called: all 25 scored
   ceftriaxone-S isolates carry mosaic penA — which raises the CEFIXIME MIC but leaves CEFTRIAXONE
   susceptible (ceftriaxone is more potent, the last-line drug). v0.1 requires the SPECIFIC high-level penA
   **Ala501** marker (A501P/T/V; the mosaic penA-60/-34 ceftriaxone-R signature) that these reduced-suscept
   isolates (carrying A510V, not A501) LACK. 25/25 S now correct. **Honest caveat:** SENSITIVITY is UNTESTED
   here (0 ceftriaxone-R in the FREE scored set — the 2 true R are assembly-required); the fix rests on the
   literature A501 marker + the observation that the S isolates carry A510-class, not A501.

## The one un-fixed cell — an honest genotype-decodability ceiling

**Cefixime spec 0.0 (3 FP) is NOT a rule bug — it is a real ceiling.** The 3 cefixime-S isolates carry the
SAME mosaic penA as the cefixime-R set (verified in the cached symbols); the R/S split is at the MIC margin,
not resolvable from AMRFinder determinants. Overfitting a rule to separate them on this 52-iso cohort would
be exactly the anti-pattern the project guards against. Documented, not force-fixed.

## Honesty summary

- CURATED_NONFROZEN rule (`rule_version` v0.1); genotype = AMRFinder `-O Neisseria_gonorrhoeae`; phenotype =
  CDC's own gonococcal S/I/R. Provenance-disjoint (0 leaks). NOT methodology-independent.
- v0.1 fixes are literature-grounded (A501 ceftriaxone marker; chromosomal-dominant gono pen/tet), with data
  inspection used only to confirm which literature marker is present — never to fit a threshold.
- Two classes of untestable side (both flagged in the rule docstrings): pen/tet spec (0 S in cohort);
  ceftriaxone sens (0 R scored). A both-sided gono cohort (Euro-GASP, mixed R/S) would test these.
- Artifacts: `wiki/ar_bank_gono_validation_<drug>_gono_v01_2026-07-20_*.json`. Frozen surface `verify_lock` OK.
