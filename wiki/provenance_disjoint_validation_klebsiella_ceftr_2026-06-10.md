# Provenance-disjoint validation — Klebsiella x ceftriaxone — 2026-06-10

Deployed decoder `call_resistance(organism=Klebsiella, drug=ceftriaxone)` scored on a PROVENANCE-DISJOINT NCBI-PD subset (submitters OUTSIDE NARMS/CDC/FDA/GenomeTrakr/PulseNet/USDA).

## Result

| metric | value |
|---|---|
| n scored | 60 (TP 30 FP 3 TN 27 FN 0; abstain 0) |
| accuracy | 0.95 |
| sensitivity (R) | 1.0 |
| specificity (S) | 0.9 |

## Independence tier (DO NOT inflate)
This is **provenance-disjoint** (different submitter / lab / country than the BV-BRC/NCBI-PD records the decoder was tuned + cross-source-validated on) — a stress-test against provenance leakage. It is NOT methodology-independent (most NCBI submitters use CLSI broth microdilution) and NOT external clinical validation. Headline accordingly. Excluded ecosystem submitters: ['narms', 'genometrakr', 'pulsenet', 'cdc', 'centers for disease', 'fda', 'food and drug', 'usda', 'national antimicrobial'].

## Leakage control
All 60 accessions are FRESH — excluded 272 accessions used in ANY prior klebsiella_* cohort (registry calibration + prior validation), so the score is on strains never seen in tuning or earlier validation.
