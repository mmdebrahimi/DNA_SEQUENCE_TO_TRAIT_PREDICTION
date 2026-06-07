# Pseudomonas_aeruginosa ciprofloxacin — cross-organism validation — 2026-06-07

> Deployed dna-amr ciprofloxacin rule applied UNCHANGED. AMRFinder `-O Pseudomonas_aeruginosa`.
- NCBI group `Pseudomonas_aeruginosa`; cohort 30 (15R/15S), 30 runs; `ncbi/amr:4.2.7-2026-03-24.1`

## VERDICT: VALIDATED

| caller | N | acc | sens | spec |
|---|---:|---:|---:|---:|
| **dna-amr (ciprofloxacin)** | 30 | **0.867** | 0.8 | 0.933 |
| naive AMRFinder | 30 | 0.767 | 1.0 | 0.533 |

## Discordance
- FN (R missed): 3
- FP (called R, susceptible): 1

## Honest scope
1 organism, 1 drug, N=30, NCBI labels (different source/curation, not a different-lab study).