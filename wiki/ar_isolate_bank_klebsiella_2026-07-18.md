# CDC AR Isolate Bank → frozen-decoder external re-validation, 2nd organism: K. pneumoniae

**Date:** 2026-07-18 · **Status:** labels built + scorer organism-parameterized; disjointness preflight
running, then scoring · **Chosen path:** phase 2 of "both, strengthen first" — reuse the AR Bank
ingester + the (now organism-parameterized) external-cohort arm for a 2nd bacterial organism ·
**Frozen surface:** byte-unchanged.

## Why this was low-friction

The E. coli re-validation built the whole pipeline; Klebsiella needed only:
1. **Scorer parameterization** (commit c67d886): `external_cohort_revalidate --amrfinder-organism
   Klebsiella_pneumoniae --registry-organism Klebsiella` (the VERBATIM triple from the frozen
   `provenance_disjoint_validation_klebsiella_*` cells), E. coli defaults preserved.
2. **The ingester already generalizes** — `build_ar_bank_labels --organism Klebsiella`.
3. **Breakpoints transfer for free** — `mic_tiers.breakpoints_for` is drug-keyed and the cipro/cef/gent
   values are CLSI M100 **Enterobacterales** breakpoints (not species-specific), valid for K. pneumoniae.

## Cohort (before leak exclusion)

157 Klebsiella isolates enumerated; 153 with BioSample + MIC. Strict-tier labels — even more
resistance-enriched than E. coli:

| Drug | strict N | R | S | Powered side |
|---|---|---|---|---|
| ceftriaxone | 143 | 143 | 0 | **sensitivity** |
| ciprofloxacin | 114 | 114 | 0 | **sensitivity** |
| gentamicin | 49 | 0 | 49 | **specificity** |

(The AR Bank Klebsiella panels are ESBL/carbapenemase-dominated → cef+cipro all-R; gent all-S because
those mechanisms don't confer gent resistance. Pure one-sided tests, large N.)

## Disjointness + scoring

<!-- KLEB_PREFLIGHT_AND_SCORING -->
Preflight (BioSample + the new resolution-free assembly-base check) running; leaks excluded like the
E. coli run; then score each powered drug with the Klebsiella organism triple. Results filled on
completion.

## Honest scope

- Same as the E. coli cell: free public data; provenance-separable at BioSample + assembly level
  (preflight-enforced); one-sided by design (resistance-enriched); reference-BMD-MIC labels (real G1);
  NOT methodology-independent (same AMRFinder `-O` + frozen `call_resistance`), NOT lineage-corrected.
- Frozen decoder surface byte-unchanged; READ-only external-validation adapter.
