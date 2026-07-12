# non-plasmid → plasmid cross-axis — lineage de-confound (leave-one-clade-out) (2026-07-12)

**Verdict: CROSS_AXIS_IS_LINEAGE_MEDIATED** — median clade-grouped AUC = 0.603 (vs naive 0.819); 9/23 plasmid replicons still predicted at AUC >= 0.7 when their clade is held out.

Mash: 232 escherichia_coli_shigella genomes → 80 clades at threshold 0.005 (largest clade 0.353).

Does the E. coli non-plasmid -> plasmid cross-axis signal survive leave-one-clade-out CV (real beyond lineage), or collapse (lineage-mediated)?

| plasmid replicon | n | AUC naive | **AUC clade-grouped** | drop | generalizes |
|---|---|---|---|---|---|
| IncX1 | 18 | 0.924 | **0.859** | 0.064 | YES |
| IncFIC(FII) | 23 | 0.935 | **0.854** | 0.081 | YES |
| IncQ1 | 17 | 0.907 | **0.814** | 0.093 | YES |
| Col156 | 66 | 0.869 | **0.801** | 0.068 | YES |
| IncFII(pSE11) | 27 | 0.895 | **0.796** | 0.099 | YES |
| IncN | 17 | 0.715 | **0.763** | -0.048 | YES |
| IncFIA | 133 | 0.869 | **0.749** | 0.12 | YES |
| IncB/O/K/Z | 13 | 0.819 | **0.745** | 0.074 | YES |
| IncFIB(AP001918) | 166 | 0.816 | **0.715** | 0.101 | YES |
| IncY | 24 | 0.654 | **0.63** | 0.024 | no (lineage) |
| IncI1-I(Alpha) | 44 | 0.716 | **0.612** | 0.103 | no (lineage) |
| IncI(Gamma) | 43 | 0.702 | **0.603** | 0.098 | no (lineage) |
| IncFII(29) | 18 | 0.709 | **0.602** | 0.107 | no (lineage) |
| ColpEC648 | 18 | 0.478 | **0.54** | -0.062 | no (lineage) |
| IncFII | 94 | 0.788 | **0.502** | 0.286 | no (lineage) |
| p0111 | 13 | 0.682 | **0.447** | 0.235 | no (lineage) |
| IncX4 | 11 | 0.426 | **0.426** | -0.0 | no (lineage) |
| IncFII(pAMA1167-NDM-5) | 66 | 0.905 | **0.422** | 0.483 | no (lineage) |
| Col(BS512) | 43 | 0.811 | **0.401** | 0.41 | no (lineage) |
| IncFII(pRSB107) | 68 | 0.872 | **0.381** | 0.491 | no (lineage) |
| Col(MG828) | 11 | 0.836 | **0.38** | 0.456 | no (lineage) |
| IncFII(pHN7A8) | 67 | 0.846 | **0.359** | 0.487 | no (lineage) |
| IncFII(pCoo) | 12 | 0.826 | **0.205** | 0.62 | no (lineage) |

## Literature / mechanism
- PREDICTION: IncF/IncQ/IncN plasmid replicons are MOBILE (conjugative), so AMR<->plasmid co-occurrence should GENERALIZE across held-out clades MORE than the chromosomal-PAI virulence axis did — mobile elements transfer horizontally between lineages, chromosomal islands do not. A clade-surviving AUC here is the horizontal-co-transfer signal; a collapse would mean the replicon is clade-fixed (vertically inherited on a stable backbone, e.g. the ST131 IncFII resistance plasmid).

## Honest caveats
- Mash-clade collapse is the proper lineage control; a within-lineage-surviving AUC means the cross-link is not PURELY clonal co-inheritance (but sub-clade structure below Mash resolution remains).
- Organism = escherichia_coli_shigella (the plasmid axis). Cohorts drug-R/S-selected. Associational.
- Generalization tracks PREVALENCE of the target feature (a common feature present across many clades has cross-clade signal to learn; a clade-restricted accessory feature does not).