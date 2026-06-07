# Klebsiella pneumoniae gentamicin — 2nd-organism drug-matrix validation — 2026-06-07

> Phase 3 matrix. Deployed dna-amr gentamicin rule applied UNCHANGED from E. coli.
- Source: NCBI Pathogen Detection `PDG000000012.2431`; cohort 30 K. pneumoniae (15R/15S), 30 with runs; `-O Klebsiella_pneumoniae`; AMRFinder `ncbi/amr:4.2.7-2026-03-24.1`

## VERDICT: VALIDATED

| caller | N | acc | sens | spec |
|---|---:|---:|---:|---:|
| **dna-amr (gentamicin rule, unchanged)** | 30 | **0.867** | 0.867 | 0.867 |
| naive AMRFinder (any drug-class determinant) | 30 | 0.667 | 1.0 | 0.333 |

## Discordance
- FN (R missed): 2
- FP (called R, susceptible): 2

## Honest scope
1 organism, 1 drug, N=30, NCBI labels (different source/curation, not a different-lab study).