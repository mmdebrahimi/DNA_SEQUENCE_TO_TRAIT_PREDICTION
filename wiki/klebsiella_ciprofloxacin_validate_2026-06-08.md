# Klebsiella ciprofloxacin — cross-organism validation — 2026-06-08

> Deployed dna-amr ciprofloxacin rule applied UNCHANGED. AMRFinder `-O Klebsiella_pneumoniae`.
- NCBI group `Klebsiella`; cohort 30 (15R/15S), 30 runs; `ncbi/amr:4.2.7-2026-03-24.1`

## VERDICT: VALIDATED

| caller | N | acc | sens | spec |
|---|---:|---:|---:|---:|
| **dna-amr (ciprofloxacin)** | 30 | **1.0** | 1.0 | 1.0 |
| naive AMRFinder | 30 | 0.5 | 1.0 | 0.0 |

## Discordance
- FN (R missed): 0
- FP (called R, susceptible): 0

## Honest scope
1 organism, 1 drug, N=30, NCBI labels (different source/curation, not a different-lab study).