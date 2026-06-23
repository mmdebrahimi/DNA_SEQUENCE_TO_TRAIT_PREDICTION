# AMR Portal feasibility + leakage — the label wall has a free door (2026-06-23)

`scripts/amr_portal_feasibility.py` measured, per (organism, drug), the PROVENANCE-DISJOINT subset of the EBI
AMR Portal (2025-12, 1.71 M measured-AST rows) — isolates whose BioSample / ERS / GCA accession is NOT in
CRyPTIC and NOT in any of our existing cohorts (`cohort_manifest`, 1,039 accessions; total leak set 37,900).
Alias-aware (the `assert_independent_aliased` discipline). **Result: a free independent number is in reach for
the headline cells AND 74 cells total — the binding label constraint of the whole project is broken.**

## Headline cells (all POWERED + provenance-disjoint, measured AST)
| Cell | disjoint (R / S) | leaked (in CRyPTIC/our cohorts) |
|---|---|---|
| **M. tuberculosis rifampicin** | **26,677** (7,926 R / 18,751 S) | 12,097 |
| **M. tuberculosis isoniazid** | **25,773** (8,757 R / 17,016 S) | 12,051 |
| **E. coli ciprofloxacin** | **16,022** (3,068 R / 12,954 S) | 420 |
| E. coli gentamicin | 15,706 (1,944 / 13,762) | 418 |
| E. coli ceftriaxone | 8,778 (1,544 / 7,234) | 396 |
| E. coli tetracycline | 8,427 (3,312 / 5,115) | 306 |
| E. coli meropenem | 12,768 (98 / 12,670) | 391 |

**74 provenance-disjoint POWERED cells (≥10 R and ≥10 S)** across M. tuberculosis, E. coli, Klebsiella,
Salmonella, Shigella, Pseudomonas, Acinetobacter, S. aureus, S. pneumoniae, N. gonorrhoeae, Enterococcus,
Enterobacter, Serratia, C. difficile, H. influenzae … (full table in `wiki/amr_portal_feasibility.json`).

## What this means (and the two things it ends)
1. **The independent TB number is FREE and in reach — the author emails are moot.** 26,677 M. tuberculosis
   rifampicin isolates with MEASURED DST that are NOT in CRyPTIC (CRyPTIC swept ~12 k of them; the rest come
   from NARMS / literature / other sources). The TB gold-set saga (5 sources, all author-request / DUA /
   circular) is resolved by a free FTP download. The Thorpe / India / Ethiopia emails are no longer the path.
2. **Every frozen bacterial cell gets a 2nd, larger external validation set** — E. coli cipro 16 k disjoint
   (only 420 overlap our cohorts), + tet / cef / gent / mero — and 60+ candidate NEW cells (Salmonella,
   Shigella, N. gonorrhoeae, Pseudomonas, Acinetobacter …) each with a free measured label.

## Honest caveats (don't overclaim — the lessons)
- **Accession-string match is an UPPER BOUND on independence.** Leakage here is BioSample/ERS/GCA string
  overlap. CRyPTIC's reuse table is keyed by ERR/ERS/UNIQUEID; AMR Portal TB rows by SAMEA/ERS/GCA. The ERS↔ERS
  match works, but a CRyPTIC isolate present under an accession alias the AMR Portal row doesn't carry could
  slip into "disjoint". The final independent number needs the BioSample-level cross-archive resolution we
  already built (`biosample_resolver` / the external-cohort Gate-0 preflight) — exactly the next step. Even
  so, 26 k disjoint almost certainly contains a large genuinely-independent core.
- **Deterministic-decoder validation only.** The AMR Portal aggregates published surveillance/literature, so
  the negative-results-map's study==class / sampling caveats apply to any EMBEDDING use. For scoring our
  FROZEN deterministic rules against measured AST, it is exactly the right substrate (non-circular by
  construction — the phenotype is wet-lab MIC/disk, not a prediction).
- **Powering is per-class.** A few cells are R- or S-starved (E. coli meropenem 98 R) — powered at ≥10 but
  report the imbalance + Wilson CI when scored.

## Recommended next move (the actual scoring)
Wire the AMR Portal disjoint subsets into the EXISTING external-cohort revalidation arm
(`scripts/external_cohort_revalidate.py` + the Gate-0 BioSample preflight): for each headline cell, take the
disjoint isolates, BioSample-resolve to confirm true independence, fetch the genomes, and SCORE the frozen
rule → the first genuinely-independent numbers for TB + the bacterial cells, free. This is the highest-VOI
move in the project: it converts "74 disjoint powered cells" into actual validated independent numbers.

## Provenance
`scripts/amr_portal_feasibility.py` (pure logic unit-tested `tests/test_amr_portal_feasibility.py`). Leak set:
`tb_goldset.cryptic_accessions` (CRyPTIC reuse table) ∪ `cohort_manifest.prior_accessions`. AMR Portal 2025-12
`phenotype.parquet` (8.7 MB, 1.71 M rows; on D:). Finding `wiki/ebi_amr_portal_finding_2026-06-23.md`. Frozen
AMR surface byte-unchanged.
