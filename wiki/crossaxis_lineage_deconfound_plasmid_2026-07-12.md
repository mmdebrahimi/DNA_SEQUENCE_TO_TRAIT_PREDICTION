# non-plasmid → plasmid cross-axis — lineage de-confound (leave-one-clade-out) (2026-07-12)

**Verdict: CROSS_AXIS_IS_LINEAGE_MEDIATED** — median clade-grouped AUC = 0.615 (vs naive 0.803); 8/23 plasmid replicons still predicted at AUC >= 0.7 when their clade is held out.

Mash: 240 E. coli genomes → 85 clades at threshold 0.005 (largest clade 0.342).

Does the E. coli non-plasmid -> plasmid cross-axis signal survive leave-one-clade-out CV (real beyond lineage), or collapse (lineage-mediated)?

| plasmid replicon | n | AUC naive | **AUC clade-grouped** | drop | generalizes |
|---|---|---|---|---|---|
| IncQ1 | 17 | 0.938 | **0.829** | 0.108 | YES |
| Col156 | 66 | 0.879 | **0.807** | 0.072 | YES |
| IncFIC(FII) | 23 | 0.915 | **0.802** | 0.113 | YES |
| IncN | 17 | 0.773 | **0.767** | 0.006 | YES |
| IncX1 | 19 | 0.936 | **0.764** | 0.172 | YES |
| IncFIB(AP001918) | 167 | 0.826 | **0.762** | 0.064 | YES |
| IncFII(pSE11) | 27 | 0.885 | **0.761** | 0.124 | YES |
| IncFIA | 134 | 0.89 | **0.757** | 0.133 | YES |
| IncB/O/K/Z | 13 | 0.787 | **0.648** | 0.139 | no (lineage) |
| IncY | 24 | 0.619 | **0.621** | -0.002 | no (lineage) |
| IncI1-I(Alpha) | 45 | 0.672 | **0.618** | 0.054 | no (lineage) |
| IncI(Gamma) | 44 | 0.685 | **0.615** | 0.07 | no (lineage) |
| IncFII(29) | 18 | 0.683 | **0.519** | 0.164 | no (lineage) |
| ColpEC648 | 18 | 0.529 | **0.513** | 0.016 | no (lineage) |
| IncFII | 95 | 0.796 | **0.508** | 0.288 | no (lineage) |
| p0111 | 13 | 0.651 | **0.431** | 0.22 | no (lineage) |
| Col(BS512) | 43 | 0.803 | **0.408** | 0.395 | no (lineage) |
| IncFII(pAMA1167-NDM-5) | 67 | 0.897 | **0.404** | 0.493 | no (lineage) |
| IncFII(pRSB107) | 69 | 0.867 | **0.369** | 0.498 | no (lineage) |
| IncFII(pHN7A8) | 69 | 0.843 | **0.353** | 0.49 | no (lineage) |
| Col(MG828) | 11 | 0.847 | **0.343** | 0.504 | no (lineage) |
| IncX4 | 11 | 0.369 | **0.329** | 0.04 | no (lineage) |
| IncFII(pCoo) | 14 | 0.793 | **0.24** | 0.553 | no (lineage) |

## Literature / mechanism
- PREDICTION: IncF/IncQ/IncN plasmid replicons are MOBILE (conjugative), so AMR<->plasmid co-occurrence should GENERALIZE across held-out clades MORE than the chromosomal-PAI virulence axis did — mobile elements transfer horizontally between lineages, chromosomal islands do not. A clade-surviving AUC here is the horizontal-co-transfer signal; a collapse would mean the replicon is clade-fixed (vertically inherited on a stable backbone, e.g. the ST131 IncFII resistance plasmid).

## Honest caveats
- Mash-clade collapse is the proper lineage control; a within-lineage-surviving AUC means the cross-link is not PURELY clonal co-inheritance (but sub-clade structure below Mash resolution remains).
- E. coli/Shigella only (the plasmid axis). Cohorts drug-R/S-selected. Associational.
- Generalization tracks PREVALENCE of the target feature (a common feature present across many clades has cross-clade signal to learn; a clade-restricted accessory feature does not).