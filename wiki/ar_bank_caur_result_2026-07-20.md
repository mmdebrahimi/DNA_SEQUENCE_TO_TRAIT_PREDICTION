# AR Bank C. auris fluconazole — fungal ERG11 cell independently confirmed (perfect, underpowered)

**Date:** 2026-07-20 (autonomous run) · **Status:** SCORED, perfect but UNDERPOWERED · **Cell:** the
NON-FROZEN fungal ERG11 cell (`dna_decode/data/fungal_amr` + `scripts/fungal_erg11_caller`), G1-validated
2026-06-08 · **Frozen surface:** byte-unchanged.

## Result

| metric | value |
|---|---|
| scored n | 8 (5 R / 3 S) |
| sensitivity | **1.00** (5/5 R) — Wilson95 [0.57, 1.0] |
| specificity | **1.00** (3/3 S) — Wilson95 [0.44, 1.0] |
| accuracy | **1.00** (0 errors) |
| powering | **UNDERPOWERED** (3 S < min-per-class 5) |

**Every call is mechanistically exact** — the fungal cell's ERG11 target-site scan (BLAST) recovered the
canonical C. auris fluconazole-R substitutions and matched the CDC MIC every time:

| isolate | ERG11 determinants | MIC-label | call |
|---|---|---|---|
| SAMN05379621 / SAMN05379619 | Y132F (+E343D,K177R,N335S) | R | ✅ R |
| SAMN05379609 | F126L | R | ✅ R |
| SAMN05379594 | K143R | R | ✅ R |
| SAMN11570381 | E343D (+K177R,N335S) | R | ✅ R |
| SAMN05379608 / SAMN05379624 / SAMN13294127 | none (WT ERG11) | S | ✅ S |

## What this is (honest)

- **An independent, provenance-disjoint confirmation** of the G1-validated fungal cell on a SEPARATE
  cohort: the CDC AR Isolate Bank C. auris, **0 overlap** with the fungal G1 tuning cohort (90 ids vs the
  32 AR Bank BioSamples). Genotype = BLAST ERG11 target-site; phenotype = CDC measured fluconazole MIC via
  the CDC tentative breakpoint (>=32 -> R) — the phenotype side is independent of the ERG11 genotype call,
  so the test is NOT circular.
- **A kingdom-jump external result:** the determinant-scan decoder, validated on bacteria, holds on a fungal
  pathogen against independent CDC MICs — the fungal analogue of the bacterial AR-Bank arm.

## Why UNDERPOWERED (a data-availability ceiling, not a cell failure)

Of 32 AR Bank C. auris, only **8 have a downloadable assembly** (24 are assembly-required / SRA-reads-only;
2 versionless-accession download failures were recovered via a `.1/.2` fallback). The fluconazole-**S**
isolates in particular mostly lack assemblies, capping the scored S at 3 (< the 5-per-class powering floor).
So the cell scores perfectly on what can be assembled, but the AR Bank C. auris cannot supply a powered
both-class set. A both-sided powered fungal validation needs a cohort with assembled S isolates (e.g. the
G1 SRA-mapping cohort, or a broader C. auris genome set).

## Scope

- CURATED_NONFROZEN fungal cell; voriconazole NOT scored (no CDC tentative breakpoint configured — cannot
  label). Only fluconazole is labelable + scored. Echinocandins (FKS1) deferred (need the FKS1 ref path).
- Labels: fluconazole MIC (AR Bank page; INT blank for antifungals) -> `fungal_amr.mic_to_phenotype`.
- Artifact: `wiki/ar_bank_caur_validation_fluconazole_caur_2026-07-20_*.json`. `scripts/ar_bank_caur_validate.py`.
  Frozen surface `verify_lock` OK.
