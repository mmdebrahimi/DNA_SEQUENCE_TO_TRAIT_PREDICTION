# Resistance→virulence cross-axis — lineage de-confound (leave-one-clade-out) (2026-07-12)

**Verdict: SPLIT_COMMON_GENERALIZES_ACCESSORY_IS_LINEAGE** — median clade-grouped AUC = 0.676 (vs naive 0.848); 3/6 virulence genes still predicted from AMR+plasmid at AUC >= 0.7 when their clade is held out.

Mash: 240 E. coli genomes → 85 clades at threshold 0.005 (largest clade 0.342).

Does the E. coli AMR+plasmid -> virulence cross-axis signal survive leave-one-clade-out CV (real beyond lineage), or collapse (lineage-mediated)?

| virulence gene | n | AUC naive | **AUC clade-grouped** | drop | generalizes |
|---|---|---|---|---|---|
| P_FIMBRIAE | 111 | 0.943 | **0.901** | 0.041 | YES |
| SIDEROPHORES | 201 | 0.9 | **0.867** | 0.033 | YES |
| CAPSULE_SERUM | 224 | 0.771 | **0.763** | 0.008 | YES |
| CNF1 | 15 | 0.917 | **0.59** | 0.327 | no (lineage) |
| AFA_DRA | 17 | 0.639 | **0.459** | 0.18 | no (lineage) |
| HEMOLYSIN | 24 | 0.796 | **0.286** | 0.51 | no (lineage) |

## Honest caveats
- Mash-clade collapse is the proper lineage control; a within-lineage-surviving AUC means the cross-link is not PURELY clonal co-inheritance (but sub-clade structure below Mash resolution remains).
- E. coli/Shigella only (the virulence axis). Cohorts drug-R/S-selected. Associational.
- Generalization tracks PREVALENCE: the 3 generalizing functions are all high-prevalence (n>=111); the 3 lineage-only genes are all low-prevalence accessory markers (n<=24). A common function present across many clades has cross-clade signal to learn; a clade-restricted accessory gene does not.