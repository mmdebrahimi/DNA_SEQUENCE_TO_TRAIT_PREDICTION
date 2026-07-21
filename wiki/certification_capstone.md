# Decoder certification capstone (2026-07-17)

> _A thin presentation layer over existing milestone evidence — NOT a clinical tool, NOT a new gate, NOT an aggregate verdict. Each cell + card keeps its own honest tier._

**No aggregate verdict.** Cells span 6 honest tiers; an aggregate pass/fail would certify the weakest cell as strongly as the strongest. No such field exists by design.

## Registry census (the per-cell spine)

- **Total cells:** 84
- **By track:** amr=25, finder=6, hla=1, mendelian=1, pgx=14, typing=8, viral=29
- **By honest evidence tier:** faithful_to_tool=13, independent_measured=27, knowledge_baseline=11, near_independent=21, no_free_source=11, not_censused=1
- **Routable-but-NOT_CENSUSED (1):** viral:HIV-1:delavirdine

## Per-domain validated headlines (verbatim — NOT averaged)

| card | status | headline |
|---|---|---|
| decoder_validation_report_card.json | present | {"honest_tier": "isolate-level provenance-disjoint stress test (different submitter/lab/country); R classes clonally dominated \u2014 lineage-effective N + cluster-weighted metrics (with Wilson CI) disclosed in the lineage table; NOT methodology-independent (most submitters use CLSI broth microdilution) and NOT lineage-independent external clinical validation", "no_aggregate_headline": true, "state_counts": {"ABSTAINS_BY_DESIGN": 2, "SCORED": 10, "NO_FREE_PHENOTYPE_SOURCE": 11, "UNDERPOWERED": 3, "LABEL_CONFOUNDED": 1}} |
| amr_portal_independent_report_card.json | present | {"n_cells": 27, "n_scored_independent": 23, "n_underpowered": 4, "status_field": "PROVENANCE_DISJOINT_INDEPENDENT_ACCESSION_LEVEL"} |
| external_validation_report_card.json | present | {"note": "external clinical re-validation of the frozen decoder; strict-tier is the primary metric, relaxed secondary; raw sens/spec is clonality-inflated \u2014 see the cluster-weighted block. Separate from the frozen decoder report card."} |
| tb_report_card.json | present | {"headline_rule": "RAW per-isolate is the TB-AMR headline; lineage is a clonality DISCLOSURE (resistance is homoplasic). NAMESPACE-SEPARATE from the frozen bacterial cards.", "n_independent_drugs": 2, "n_in_distribution_drugs": 2} |
| hiv_decoder_report_card.json | present | {"label_independence": "PhenoSense fold-change is NOT HIVDB's own Sierra interpretation (non-circular)", "n_cells": 25, "honest_caveats": ["in-distribution (HIVDB), NOT provenance-disjoint -> a lower external-rigour bar than the bacterial card", "NNRTI = mutant-specific (excellent on 1st-gen EFV/NVP); NRTI v0 = position-based (over-calls, fixed by the deconfounded mutant-specific v0.1 for 5/6 drugs; ddI keeps position-based)", "non-B subtype transfer is under-powered (data ~96% subtype B)", "PI/INSTI = position-based v0 (PI AUC 0.78-0.96; INSTI 0.74-1.0, 2nd-gen DTG/BIC lower as the class-level over-call predicts); CAI/lenacapavir = mutant-level (AUC 0.91) on a small resistance-enriched dataset (n=140, 11 S)", "OLS underlying-tool baseline now run for PI/INSTI/CAI (uniform illustrative fold>=3 cutoff, delta is the wrapper-vs-tool signal): PI catalog is high-sens/low-spec so OLS recovers +0.06..+0.13 balacc (real v0.1 mutant-specific headroom, like NRTI); INSTI catalog is competitive (+-0.07, ties/beats OLS on EVG/BIC); CAI catalog BEATS OLS +0.112 (OLS overfits the tiny resistance-enriched set)", "WITHIN-SUBTYPE de-confounding now cleared for ALL 4 classes (2026-07-03): the catalogs decode MECHANISM not subtype structure \u2014 median within-B AUC NNRTI 0.795 / PI 0.921 / INSTI 0.898 (NRTI held earlier); pooled-minus-within-B ~0 for every class -> the class-mixed numbers were NOT subtype-inflated. Subtype-transfer column now populated for all classes (within-B AUC).", "PI v0.1 (2026-06-23) + INSTI v0.1 (2026-06-27) deconfounded mutant-specific catalogs SHIPPED (same OLS-coef + 5-fold-CV arc as NRTI): PI 8/8 improve-or-hold (mean +0.056), INSTI 5/5 improve-or-hold (mean +0.087). The HIV class-deconfounding arc NRTI->PI->INSTI is COMPLETE.", "v0.2 ABSOLUTE-CUTOFF calibration DONE (2026-07-03) by SOURCING per-drug cutoffs from Stanford DRMcv.R (not fabricated): PI calibrated 8/8 (real per-drug cutoffs LPV 9/TPV 2/DRV 10; the position-based over-call now quantified as low spec at the real cutoff); NNRTI confirmed 4/5 (all DRMcv cutoffs = the prior illustrative 3; DOR postdates DRMcv -> walled); INSTI 0/5 CUTOFF_UNAVAILABLE (integrase inhibitors postdate DRMcv.R -> external wall, reported not guessed). wiki/hiv_{nnrti,pi,insti}_absolute_cutoff_2026-07-03."]} |
| pgx_report_card.json | present | {"note": "Standing PGx trust surface -- a roll-up, NOT a gate (exit 0 always). No aggregate headline; each cell's honest tier stands alone. CALLING is independently validatable vs GeT-RM (free consensus panel); PHENOTYPE is faithful-to-CPIC (assigned, not measured).", "sources": {"getrm_cyp2c19": true, "getrm_cyp2c9": true, "getrm_cyp2c8": true, "getrm_cyp3a5": true, "getrm_tpmt": true, "getrm_cyp2b6": true, "getrm_cyp2d6": true, "pharmcat_cyp2c19": true, "functional_evidence": true, "trio_mendelian": true}} |

## Boundaries (what is closed + the label wall)

- **reproducibility_freeze_2026-06-13.md** (present): Reproducibility freeze — deterministic AMR decoder (2026-06-13)
- **negative_results_map_2026-06-13.md** (present): Operational negative-results map — why each G2P expansion was rejected (2026-06-13)

_A presentation of existing evidence; it certifies nothing the underlying cards do not. NOT a clinical tool._
