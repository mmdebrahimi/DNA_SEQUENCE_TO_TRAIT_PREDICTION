# Provenance-disjoint validation — Escherichia_coli_Shigella x gentamicin — 2026-06-12

Deployed decoder `call_resistance(organism=Escherichia_coli_Shigella, drug=gentamicin)` scored on a PROVENANCE-DISJOINT NCBI-PD subset (submitters OUTSIDE NARMS/CDC/FDA/GenomeTrakr/PulseNet/USDA).

## Result

| metric | value |
|---|---|
| n scored | 60 (TP 27 FP 0 TN 30 FN 3; abstain 0) |
| accuracy | 0.95 |
| sensitivity (R) | 0.9 |
| specificity (S) | 1.0 |

## Independence tier (DO NOT inflate)
This is **provenance-disjoint** (different submitter / lab / country than the BV-BRC/NCBI-PD records the decoder was tuned + cross-source-validated on) — a stress-test against provenance leakage. It is NOT methodology-independent (most NCBI submitters use CLSI broth microdilution) and NOT external clinical validation. Headline accordingly. Excluded ecosystem submitters: ['narms', 'genometrakr', 'pulsenet', 'cdc', 'centers for disease', 'fda', 'food and drug', 'usda', 'national antimicrobial'].

## Leakage control
All 60 accessions are FRESH — excluded 919 accessions used in ANY prior escherichia_coli_shigella_* cohort (registry calibration + prior validation), so the score is on strains never seen in tuning or earlier validation.
