# Klebsiella pneumoniae meropenem — 2nd-organism, carbapenem mechanism — 2026-06-07

> Phase 3 slice 2. Carbapenem (KPC/NDM/OXA-48) — a mechanism class E. coli AMR never covered.
- Source: NCBI Pathogen Detection `PDG000000012.2431`; cohort 30 K. pneumoniae (15R/15S), 30 with runs
- AMRFinder `ncbi/amr:4.2.7-2026-03-24.1` `-O Klebsiella_pneumoniae`; rule: carbapenemase (CARBAPENEM-subclass) >=1

## VERDICT: VALIDATED

| caller | N | acc | sens | spec |
|---|---:|---:|---:|---:|
| **dna-amr (CARBAPENEM-subclass)** | 30 | **0.867** | 1.0 | 0.733 |
| naive AMRFinder (any beta-lactam/carbapenem determinant) | 30 | 0.533 | 1.0 | 0.067 |

## Discordance
- FN (R missed — porin-loss/ESBL+impermeability/low-level): 0
- FP (called R, susceptible): 4

## Honest scope
1 organism, 1 drug, NCBI labels (different source/curation, not a different-lab study). The rule is
blind to porin-loss-mediated carbapenem resistance (no carbapenemase gene) — expected FN mode.