# First genuinely-independent numbers for the frozen bacterial cells (AMR Portal, 2026-06-23)

The recommendation from `wiki/amr_portal_feasibility_result_2026-06-23.md` — score the frozen rule on the
AMR Portal's provenance-disjoint, measured-phenotype isolates — RAN, with no Docker and no genome fetch. The
AMR Portal ships a per-isolate **AMRFinderPlus-format genotype table** (incl. 99 k POINT mutations: gyrA/parC
QRDR etc.), so `scripts/amr_portal_score_independent.py` reconstructs a faithful `main.tsv` per isolate and
calls the **FROZEN `amr_rules.call_resistance` unchanged** → confusion vs the measured `resistance_phenotype`.

## Result — E. coli, provenance-disjoint (not in CRyPTIC / our cohorts), measured AST
| Cell | nR / nS | sens (95% CI) | spec (95% CI) | accuracy | tuning-cohort acc (for comparison) |
|---|---|---|---|---|---|
| **gentamicin** | 1,944 / 13,762 | 0.927 [0.915, 0.938] | 0.995 [0.994, 0.996] | **0.987** | 0.945 (N=128) |
| **tetracycline** | 3,310 / 5,113 | 0.970 [0.963, 0.975] | 0.987 [0.984, 0.990] | **0.980** | 0.833 (N=12, noisy) |
| **ciprofloxacin** | 3,068 / 12,951 | 0.837 [0.823, 0.849] | 0.977 [0.974, 0.979] | **0.950** | 0.925 (N=147) |
| **ceftriaxone** | 1,544 / 7,234 | 0.921 [0.906, 0.933] | 0.918 [0.912, 0.924] | **0.919** | 0.933 (N=60) |
| **meropenem** | 98 / 12,670 | 0.776 [0.683, 0.847] | 0.988 [0.986, 0.989] | **0.986** | — |

**The frozen deterministic rule holds up on independent data** — 0.92–0.99 accuracy across 5 drugs on
8 k–16 k disjoint isolates each, every cell powered, with Wilson CIs. This is the external-cohort validation
HIV has, now delivered for the deployed bacterial core — for free. The tet number especially: the deployed
config was tuned on a noisy N=12 (acc 0.833); the independent N=8,423 confirms it at 0.980.

## Honest reading (the rule's known shape shows through, correctly)
- **cipro sens 0.837 < spec 0.977** — the QRDR-threshold-2 rule UNDER-calls some R (single-QRDR + plasmid-
  mediated quinolone resistance qnr/aac(6')-Ib-cr are below the 2-point threshold by design). A documented
  blind spot, not a regression — and the validation quantifies it at scale for the first time.
- **meropenem sens 0.776 on only 98 R** — carbapenem-R is rare in E. coli; the sens CI is wide
  [0.683, 0.847]. Powered at ≥10 but R-starved; report as such.
- **cef sens/spec ≈ 0.92 balanced** — the extended-spectrum subclass refinement generalizes (the broad-class
  version would over-call; this confirms the refinement was the right call).

## Honesty rails (what this is, and is NOT)
1. **Independent at the ACCESSION level (upper bound).** Leakage is BioSample/ERS/GCA-disjoint vs CRyPTIC +
   our cohorts. BioSample cross-archive resolution (`biosample_resolver`, the external-arm Gate-0) would only
   TIGHTEN this — a refinement, not a correction. The overlap with our tiny tuning cohorts was ~400/16 k, so
   the disjoint set is overwhelmingly the AMR Portal's own (largely NARMS / literature / PATRIC).
2. **Genotype = the AMR Portal's OWN AMRFinder run** (a different operator, possibly a different AMRFinder
   version than our pinned image). That makes it MORE independent (different pipeline), but the AMRFinder
   version is a named caveat — point-mutation calling can differ slightly across versions.
3. **The rule is FROZEN + applied UNCHANGED** (`organism=None` → the validated `DRUG_RULE` path, exactly as
   the E. coli cells were validated). The frozen surface (`amr_rules.py` + `calibrated_amr_rules.json`) is
   byte-unchanged (leak-guard green).
4. **Measured phenotype = non-circular** (wet-lab MIC/disk, not a prediction).

## Next (documented, not done)
- **TB independent number.** The same data has 26 k disjoint M. tuberculosis isolates with measured DST +
  AMRFinder rpoB/katG POINT calls. TB uses a DIFFERENT rule (the WHO catalogue on VCF, `organism_rules/tb_amr`),
  so scoring it needs an AMR-Portal-POINT → WHO-catalogue adapter (rpoB_S450L → the catalogue determinant).
  That is the clear next step — it would finally deliver the independent TB number the gold-set saga blocked.
- **Tighten with BioSample resolution** (Gate-0) for the headline-publishable number.
- **Candidate NEW cells.** Salmonella, Shigella, N. gonorrhoeae, Klebsiella, Pseudomonas, Acinetobacter each
  have a powered disjoint measured-AST set here (the 74 cells) — each a free new cell.

## Provenance
`scripts/amr_portal_score_independent.py` (pure logic + a frozen-rule faithfulness test in
`tests/test_amr_portal_score_independent.py`). Scores `wiki/amr_portal_independent_scores.json`. Inputs: AMR
Portal 2025-12 phenotype + genotype parquet (on D:); leak set = `tb_goldset.cryptic_accessions` ∪
`cohort_manifest.prior_accessions`. Feasibility `wiki/amr_portal_feasibility_result_2026-06-23.md`.
