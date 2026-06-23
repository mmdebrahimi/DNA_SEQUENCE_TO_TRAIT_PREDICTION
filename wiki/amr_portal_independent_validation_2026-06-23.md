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

## Cross-organism independent validation (added — the rules + calibrated registry GENERALIZE)
The same frozen scorer, routed per organism (E. coli/Shigella → DRUG_RULE default; Klebsiella/Salmonella →
the OPT-IN calibrated registry for cipro), on the AMR Portal disjoint sets. **This is the first INDEPENDENT
validation of the calibrated organism registry** — whose own provenance says "CALIBRATED configs need an
INDEPENDENT cohort before becoming a default" (it was IN-SAMPLE N~30).

| Cell | N (R/S) | sens | spec | accuracy | note |
|---|---|---|---|---|---|
| **Salmonella ciprofloxacin** | 2,434 / 22,538 | 0.907 | 0.964 | **0.959** | calibrated (broad) — VALIDATED at N=24,972 |
| Salmonella ceftriaxone | 1,629 / 21,259 | 0.952 | 0.998 | **0.995** | default rule generalizes |
| Salmonella gentamicin | 1,122 / 24,608 | 0.912 | 0.991 | **0.987** | |
| Salmonella tetracycline | 6,941 / 19,411 | 0.958 | 0.992 | **0.983** | |
| **Klebsiella ciprofloxacin** | 2,970 / 1,415 | 0.755 | **0.994** | 0.832 | calibrated (qrdr_point + oqxAB-exclusion) — spec 0.994 CONFIRMS the efflux-exclusion design; sens 0.755 = the non-QRDR (plasmid qnr / porin) blind spot, now quantified |
| Klebsiella gentamicin | 2,026 / 2,904 | 0.928 | 0.963 | 0.949 | default generalizes |
| Klebsiella ceftriaxone | 2,496 / 874 | 0.964 | 0.905 | 0.948 | |
| Klebsiella meropenem | 1,724 / 2,898 | 0.948 | 0.859 | 0.892 | |
| Klebsiella tetracycline | 1,059 / 936 | 0.731 | 0.980 | 0.848 | sens-limited |
| Shigella sonnei ceftriaxone | 426 / 910 | 0.988 | 0.995 | **0.993** | Shigella shares E. coli rules |
| Shigella sonnei tetracycline | 951 / 475 | 0.950 | 0.975 | 0.958 | |
| Shigella sonnei ciprofloxacin | 625 / 961 | 0.730 | 0.997 | 0.892 | sens-limited (plasmid/single-QRDR) |
| Shigella sonnei gentamicin | 513 / 967 | 0.472 | 0.999 | 0.816 | low sens — aminoglycoside mechanism the gene rule misses |

**Headline:** the deterministic decoder generalizes across the Enterobacterales + Salmonella on independent
data — Salmonella 0.96–0.995, E. coli 0.92–0.99, Klebsiella 0.83–0.95, Shigella 0.82–0.99. The calibrated
Salmonella/Klebsiella cipro configs are now independently validated (Salmonella excellent; Klebsiella highly
SPECIFIC with the intended oqxAB-exclusion, but sens-limited by non-QRDR mechanisms — a real, now-measured
ceiling, not a regression). Salmonella + Shigella meropenem have ZERO resistant isolates (not powered;
carbapenem-R absent in these species' disjoint sets) — reported, not hidden.

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

## Next (documented, not done) — and a TB honesty caveat
- **Why the bacterial cells are a CLEAN test here:** our deployed bacterial rule IS an AMRFinder-determinant
  rule, so the AMR Portal's AMRFinder genotype table is the EXACT right input — scoring it tests OUR rule.
- **TB is NOT clean via this genotype table (a real confound, not just "an adapter").** Our TB rule is the
  WHO-2023-catalogue on RAW per-isolate VCF (`organism_rules/tb_amr`). The AMR Portal genotype is AMRFinder's
  own (narrower, catalogue-derived) TB POINT calls — scoring our WHO rule against AMRFinder's pre-filtered TB
  calls would test AMRFinder's TB caller, not the WHO catalogue, AND mismatch the rule's raw-variant input.
  So the 26 k disjoint M. tuberculosis isolates here give a free independent **phenotype** + **accession set**,
  but the TB independent number needs the RAW variants (fetch the assemblies the AMR Portal links → call
  variants → score the WHO rule), i.e. the heavier genome-fetch path — not this table. This is the honest TB
  next step (and it's now unblocked on the LABEL, which was the binding constraint).
- **Tighten with BioSample resolution** (Gate-0) for the headline-publishable number.
- **Candidate NEW cells.** Salmonella, Shigella, N. gonorrhoeae, Klebsiella, Pseudomonas, Acinetobacter each
  have a powered disjoint measured-AST set here (the 74 cells) — each a free new cell.

## Provenance
`scripts/amr_portal_score_independent.py` (pure logic + a frozen-rule faithfulness test in
`tests/test_amr_portal_score_independent.py`). Scores `wiki/amr_portal_independent_scores.json`. Inputs: AMR
Portal 2025-12 phenotype + genotype parquet (on D:); leak set = `tb_goldset.cryptic_accessions` ∪
`cohort_manifest.prior_accessions`. Feasibility `wiki/amr_portal_feasibility_result_2026-06-23.md`.
