# Campylobacter external validation on NCBI-PD — tet + gent CONFIRMED (2026-07-21)

**Status:** ✅ **tetracycline + gentamicin both SCORED_ENDORSED** — the non-frozen `campylobacter_amr`
cells externally validated · demonstrates the **generic no-compute NCBI-PD substrate** (follow-up #2 of the
gonococcus external validation) · **Frozen surface:** byte-unchanged.

## What this demonstrates

The gonococcus run discovered that NCBI Pathogen Detection publishes its OWN AMRFinderPlus determinant calls
(the `AMR_genotypes` field) + isolates have downloadable assemblies → external validation of any
determinant-based cell is a pure **metadata join + score, offline, no compute**. This applies that substrate
to a *second* organism via a **generic scorer** (`scripts/score_ncbipd_extval.py --organism {gono,campylobacter}`).

## Result (121-isolate provenance-disjoint cohort; frozen cipro accessions excluded)

| drug | n | R/S | sens | spec | acc | verdict |
|---|---|---|---|---|---|---|
| **tetracycline** (tet(O)) | 121 | 66R/55S | 1.00 | 0.945 | 0.975 | **SCORED_ENDORSED** ✅ |
| **gentamicin** (aph/aac enzymes) | 121 | 31R/90S | 0.968 | 1.00 | 0.992 | **SCORED_ENDORSED** ✅ |

Both `campylobacter_amr` cells were NON-FROZEN additions that had never been independently scored (only the
frozen *ciprofloxacin* cell was validated). This is their first external validation:
- **tet(O)** ribosomal-protection gene presence → tet R: sens 1.0, spec 0.945 (3 FP — likely non-tet(O)
  tet-R or a borderline MIC). Clean gene-presence determinant.
- **gentamicin** true-enzyme rule (aph(2'')/aac(3)/aac(6')-aph, with the aad9/spw non-gentamicin exclusion):
  sens 0.968, spec 1.00 — the intrinsic-gene-exclusion discipline holds on 90 independent gent-S isolates.

## Honest caveats

- **Provenance-disjoint** (128 frozen campylobacter identifiers excluded at cohort-build; and the frozen campy
  cell is *ciprofloxacin* — a different drug — so the tet/gent rules never saw these isolates regardless) but
  **NOT methodology-independent** (same AMRFinderPlus + same cell).
- **Lineage-collapse (clonality-corrected) — no compute:** via NCBI-PD's own SNP clusters (`PDS_acc`) +
  `clonality.cluster_weighted_confusion`. **Both cells HOLD:** tetracycline lineage sens 1.0 / spec 0.933
  (1 discordant), gentamicin lineage sens 1.0 / spec 1.0 (1 discordant). Not clonally inflated — the rules
  decode mechanism, not clonal structure.
- Campylobacter measured AST is **surveillance-dominated** (NARMS/food); this cohort is a subset of that
  ecosystem (the isolates carry measured AST because they were surveilled) — the provenance-disjointness is
  from the *frozen cohort*, not from the surveillance ecosystem.

## The reusable takeaway

`score_ncbipd_extval.py` is now a **generic no-compute external-validation harness** — point it at any
NCBI-PD-covered organism's cohort + determinants (built from `<PDG>.amr.metadata.tsv`) + its non-frozen cell.
Demonstrated on gonococcus (cipro/cefixime/penicillin) + campylobacter (tet/gent). The DEGENERATE guard (a
cell predicting all-one-class is never endorsed) is shared. NON-FROZEN cells; frozen surface byte-unchanged.
