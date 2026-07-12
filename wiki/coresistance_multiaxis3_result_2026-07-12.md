# Multi-axis co-resistance — AMR determinant × plasmid Inc-type × virulence (2026-07-12)

**Verdict: PASS_MULTIAXIS_LINKAGE** — 166/169 within-organism features (determinant / Inc-type / virulence gene) impute from the joint set (fraction 0.982; bar 0.5). (raw)

Axes: amr_determinant, plasmid_inc_type, virulence_gene. Do AMR determinants, plasmid Inc-types, and (E. coli) virulence genes co-occur / impute each other, within organism? (which resistance rides which plasmid backbone / virulence context) Substrate: AMR determinants (AMRFinder) + plasmid Inc-replicons (PlasmidFinder) + virulence genes (VirulenceFinder, E. coli) blastn local sweeps over Enterobacterales cohort genomes.

| organism | genomes | Inc-types | vir genes | testable | **imputable** | frac | note |
|---|---|---|---|---|---|---|---|
| klebsiella | 307 | 62 | 0 | 81 | **81** | 1.0 |  |
| escherichia_coli_shigella | 240 | 55 | 13 | 74 | **71** | 0.959 |  |
| salmonella | 60 | 25 | 0 | 14 | **14** | 1.0 |  |

## Which resistance rides which plasmid backbone / virulence context (determinant → cross-axis, by lift)
### klebsiella
- **parC_S80I** → IncI2(Delta)(lift 2.72)
- **gyrA_S83I** → IncI2(Delta)(lift 3.16)
- **blaSHV-11** → IncN(lift 2.19), IncHI1B(pNDM-MAR)(lift 2.06), IncFIB(K)(pCAV1099-114)(lift 1.78), repB(lift 1.78)
### escherichia_coli_shigella
- **uhpT_E350Q** → Col156(lift 1.9), VIR/CNF1(lift 1.87)
- **tet(A)** → IncN(lift 1.68)
- **sul2** → IncQ1(lift 2.35)
- **aph(3'')-Ib** → IncQ1(lift 2.4)
### salmonella
- **tet(A)** → IncFIB(pN55391)(lift 2.17)
- **sul1** → IncFIB(pN55391)(lift 2.82)
- **gyrA_D87Y** → IncFIB(pN55391)(lift 3.16)
- **aadA1** → IncFIB(pN55391)(lift 3.32)
- **aph(3')-Ia** → IncFIB(pN55391)(lift 2.99)
- **aac(3)-IVa** → IncFIB(pN55391)(lift 2.94)
- **aph(4)-Ia** → IncFIB(pN55391)(lift 2.94)

## Cross-axis-only: does AMR+plasmid predict VIRULENCE? (excludes virulence-virulence PAI clustering)
- **escherichia_coli_shigella** (AMR+plasmid → virulence-gene AUC): P_FIMBRIAE 0.947, CNF1 0.901, SIDEROPHORES 0.871, HEMOLYSIN 0.821, CAPSULE_SERUM 0.796, AFA_DRA 0.788

## Honest caveats
- Plasmid axis is Enterobacterales-only; virulence axis is E. coli/Shigella-only (VirulenceFinder DB).
- PlasmidFinder / VirulenceFinder PRESENCE is a genotype axis (blastn), not a wet content readout.
- Cohorts are drug-R/S-selected; within-organism de-confound + dedup clonality proxy; associational.
- The cross-axis AMR+plasmid->virulence signal (AUC ~0.79-0.95 in E. coli) is likely LINEAGE-mediated (specific STs, e.g. ExPEC ST131, carry both particular CTX-M plasmids AND UPEC virulence PAIs); within-organism controls species but not sub-lineage. Mash-clade collapse is the proper control.