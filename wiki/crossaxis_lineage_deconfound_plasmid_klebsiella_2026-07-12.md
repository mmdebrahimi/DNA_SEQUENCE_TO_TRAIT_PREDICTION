# non-plasmid → plasmid cross-axis — lineage de-confound (leave-one-clade-out) (2026-07-12)

**Verdict: CROSS_AXIS_IS_LINEAGE_MEDIATED** — median clade-grouped AUC = 0.668 (vs naive 0.84); 9/21 plasmid replicons still predicted at AUC >= 0.7 when their clade is held out.

Mash: 307 klebsiella genomes → 118 clades at threshold 0.005 (largest clade 0.254).

Does the E. coli non-plasmid -> plasmid cross-axis signal survive leave-one-clade-out CV (real beyond lineage), or collapse (lineage-mediated)?

| plasmid replicon | n | AUC naive | **AUC clade-grouped** | drop | generalizes |
|---|---|---|---|---|---|
| IncFIB(pKPHS1) | 22 | 0.995 | **0.991** | 0.004 | YES |
| ColRNAI | 106 | 0.938 | **0.912** | 0.026 | YES |
| IncFIB(pQil) | 53 | 0.947 | **0.817** | 0.13 | YES |
| RepB | 25 | 0.84 | **0.764** | 0.076 | YES |
| IncFII(pKP91) | 23 | 0.812 | **0.759** | 0.053 | YES |
| IncFIB(K)(pCAV1099-114) | 38 | 0.82 | **0.749** | 0.071 | YES |
| repB | 38 | 0.82 | **0.749** | 0.071 | YES |
| IncC | 11 | 0.912 | **0.739** | 0.173 | YES |
| IncFII(K) | 135 | 0.858 | **0.712** | 0.146 | YES |
| IncFII(Yp) | 29 | 0.951 | **0.684** | 0.268 | no (lineage) |
| IncHI1B(pNDM-MAR) | 18 | 0.838 | **0.652** | 0.186 | no (lineage) |
| IncFIB(K) | 151 | 0.816 | **0.646** | 0.17 | no (lineage) |
| IncFIA(HI1) | 20 | 0.717 | **0.64** | 0.076 | no (lineage) |
| repB(R1701) | 61 | 0.894 | **0.638** | 0.257 | no (lineage) |
| Col440I | 16 | 0.726 | **0.631** | 0.095 | no (lineage) |
| IncFIA(pBK30683) | 26 | 0.887 | **0.627** | 0.26 | no (lineage) |
| IncR | 84 | 0.836 | **0.621** | 0.214 | no (lineage) |
| IncX3 | 52 | 0.958 | **0.467** | 0.491 | no (lineage) |
| IncN | 20 | 0.784 | **0.348** | 0.436 | no (lineage) |
| Col440II | 14 | 0.752 | **0.305** | 0.446 | no (lineage) |
| IncI2(Delta) | 24 | 0.963 | **clade-concentrated** | None | no (clade-concentrated) |

## Literature / mechanism
- PREDICTION: IncF/IncQ/IncN plasmid replicons are MOBILE (conjugative), so AMR<->plasmid co-occurrence should GENERALIZE across held-out clades MORE than the chromosomal-PAI virulence axis did — mobile elements transfer horizontally between lineages, chromosomal islands do not. A clade-surviving AUC here is the horizontal-co-transfer signal; a collapse would mean the replicon is clade-fixed (vertically inherited on a stable backbone, e.g. the ST131 IncFII resistance plasmid).

## Honest caveats
- Mash-clade collapse is the proper lineage control; a within-lineage-surviving AUC means the cross-link is not PURELY clonal co-inheritance (but sub-clade structure below Mash resolution remains).
- Organism = klebsiella (the plasmid axis). Cohorts drug-R/S-selected. Associational.
- Generalization tracks PREVALENCE of the target feature (a common feature present across many clades has cross-clade signal to learn; a clade-restricted accessory feature does not).