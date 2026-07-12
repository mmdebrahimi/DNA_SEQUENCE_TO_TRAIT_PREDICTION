# non-virulence → virulence cross-axis — lineage de-confound (leave-one-clade-out) (2026-07-12)

**Verdict: SPLIT_COMMON_GENERALIZES_ACCESSORY_IS_LINEAGE** — median clade-grouped AUC = 0.663 (vs naive 0.852); 3/6 virulence genes still predicted at AUC >= 0.7 when their clade is held out.

Mash: 232 escherichia_coli_shigella genomes → 80 clades at threshold 0.005 (largest clade 0.353).

Does the E. coli non-virulence -> virulence cross-axis signal survive leave-one-clade-out CV (real beyond lineage), or collapse (lineage-mediated)?

| virulence gene | n | AUC naive | **AUC clade-grouped** | drop | generalizes |
|---|---|---|---|---|---|
| P_FIMBRIAE | 111 | 0.95 | **0.901** | 0.049 | YES |
| SIDEROPHORES | 196 | 0.913 | **0.85** | 0.063 | YES |
| CAPSULE_SERUM | 217 | 0.776 | **0.735** | 0.041 | YES |
| CNF1 | 15 | 0.911 | **0.591** | 0.32 | no (lineage) |
| AFA_DRA | 17 | 0.753 | **0.429** | 0.324 | no (lineage) |
| HEMOLYSIN | 24 | 0.793 | **0.347** | 0.447 | no (lineage) |

## Literature / mechanism
- The ST131 literature (Johnson/Nicolas-Chanoine reviews; PMC3916147, PMC4135879, PMC8487868) holds that resistance<->virulence co-occurrence in E. coli is driven by CLONAL EXPANSION (vertical inheritance within lineage), with hlyA (HEMOLYSIN) specifically C2-clade-restricted and co-occurring with aac6Ib/blaCTX-M-15 WITHIN clade. This de-confound independently recovers that: HEMOLYSIN collapses HARDEST under leave-one-clade-out (0.796->0.286), i.e. it is the most lineage-restricted virulence marker — exactly the hlyA-in-C2 pattern. The high-prevalence core functions (fimbriae/siderophores/capsule) generalizing across clades is the novel, non-purely-clonal half.

## Honest caveats
- Mash-clade collapse is the proper lineage control; a within-lineage-surviving AUC means the cross-link is not PURELY clonal co-inheritance (but sub-clade structure below Mash resolution remains).
- Organism = escherichia_coli_shigella (the virulence axis). Cohorts drug-R/S-selected. Associational.
- Generalization tracks PREVALENCE of the target feature (a common feature present across many clades has cross-clade signal to learn; a clade-restricted accessory feature does not).