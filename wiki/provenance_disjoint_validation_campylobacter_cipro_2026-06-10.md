# Provenance-disjoint validation — Campylobacter x ciprofloxacin — 2026-06-10

Deployed decoder `call_resistance(organism=Campylobacter, drug=ciprofloxacin)` scored on a PROVENANCE-DISJOINT NCBI-PD subset (submitters OUTSIDE NARMS/CDC/FDA/GenomeTrakr/PulseNet/USDA).

## Result

| metric | value |
|---|---|
| n scored | 40 (TP 20 FP 0 TN 20 FN 0; abstain 0) |
| accuracy | 1.0 |
| sensitivity (R) | 1.0 |
| specificity (S) | 1.0 |

## Independence tier (DO NOT inflate)
This is **provenance-disjoint** (different submitter / lab / country than the BV-BRC/NCBI-PD records the decoder was tuned + cross-source-validated on) — a stress-test against provenance leakage. It is NOT methodology-independent (most NCBI submitters use CLSI broth microdilution) and NOT external clinical validation. Headline accordingly. Excluded ecosystem submitters: ['narms', 'genometrakr', 'pulsenet', 'cdc', 'centers for disease', 'fda', 'food and drug', 'usda', 'national antimicrobial'].

## Leakage control
All 40 accessions are FRESH — excluded 60 accessions used in ANY prior campylobacter_* cohort (registry calibration + prior validation), so the score is on strains never seen in tuning or earlier validation.
