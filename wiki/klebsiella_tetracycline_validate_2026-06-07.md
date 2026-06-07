# Klebsiella pneumoniae tetracycline — 2nd-organism drug-matrix validation — 2026-06-07

> Phase 3 matrix. Deployed dna-amr tetracycline rule applied UNCHANGED from E. coli.
- Source: NCBI Pathogen Detection `PDG000000012.2431`; cohort 30 K. pneumoniae (15R/15S), 30 with runs; `-O Klebsiella_pneumoniae`; AMRFinder `ncbi/amr:4.2.7-2026-03-24.1`

## VERDICT: FAILS_BAR

| caller | N | acc | sens | spec |
|---|---:|---:|---:|---:|
| **dna-amr (tetracycline rule, unchanged)** | 30 | **0.8** | 0.6 | 1.0 |
| naive AMRFinder (any drug-class determinant) | 30 | 0.6 | 1.0 | 0.2 |

## Discordance
- FN (R missed): 6
- FP (called R, susceptible): 0

## Honest scope
1 organism, 1 drug, N=30, NCBI labels (different source/curation, not a different-lab study).