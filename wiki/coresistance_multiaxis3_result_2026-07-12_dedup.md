# Multi-axis co-resistance — AMR determinant × plasmid Inc-type × virulence (2026-07-12)

**Verdict: PASS_MULTIAXIS_LINKAGE** — 149/155 within-organism features (determinant / Inc-type / virulence gene) impute from the joint set (fraction 0.961; bar 0.5). (deduped)

Axes: amr_determinant, plasmid_inc_type, virulence_gene. Do AMR determinants, plasmid Inc-types, and (E. coli) virulence genes co-occur / impute each other, within organism? (which resistance rides which plasmid backbone / virulence context) Substrate: AMR determinants (AMRFinder) + plasmid Inc-replicons (PlasmidFinder) + virulence genes (VirulenceFinder, E. coli) blastn local sweeps over Enterobacterales cohort genomes.

| organism | genomes | Inc-types | vir genes | testable | **imputable** | frac | note |
|---|---|---|---|---|---|---|---|
| klebsiella | 259 | 62 | 0 | 78 | **77** | 0.987 | profile-deduped 307->259 |
| escherichia_coli_shigella | 197 | 55 | 13 | 71 | **66** | 0.93 | profile-deduped 240->197 |
| salmonella | 35 | 25 | 0 | 6 | **6** | 1.0 | profile-deduped 60->35 |

## Which resistance rides which plasmid backbone / virulence context (determinant → cross-axis, by lift)
### klebsiella
- **parC_S80I** → IncI2(Delta)(lift 2.64)
- **blaSHV-11** → IncN(lift 2.0), IncFII(pKP91)(lift 1.46)
- **gyrA_S83I** → IncI2(Delta)(lift 3.12)
### escherichia_coli_shigella
- **uhpT_E350Q** → Col156(lift 1.97), VIR/CNF1(lift 1.96)
- **sul1** → VIR/HEMOLYSIN(lift 1.64)
- **tet(A)** → IncN(lift 1.87)
- **blaCTX-M-15** → VIR/CNF1(lift 2.14), VIR/HEMOLYSIN(lift 1.92)
- **sul2** → IncQ1(lift 2.56)
### salmonella
- **tet(A)** → IncFIB(pN55391)(lift 1.85)
- **gyrA_D87Y** → IncFIB(pN55391)(lift 2.92)
- **sul1** → IncFIB(pN55391)(lift 2.62)
- **aadA1** → IncFIB(pN55391)(lift 3.15)

## Cross-axis-only: does AMR+plasmid predict VIRULENCE? (excludes virulence-virulence PAI clustering)
- **escherichia_coli_shigella** (AMR+plasmid → virulence-gene AUC): P_FIMBRIAE 0.922, CNF1 0.869, SIDEROPHORES 0.859, CAPSULE_SERUM 0.75, HEMOLYSIN 0.74, AFA_DRA 0.669

## Honest caveats
- Plasmid axis is Enterobacterales-only; virulence axis is E. coli/Shigella-only (VirulenceFinder DB).
- PlasmidFinder / VirulenceFinder PRESENCE is a genotype axis (blastn), not a wet content readout.
- Cohorts are drug-R/S-selected; within-organism de-confound + dedup clonality proxy; associational.
- The cross-axis AMR+plasmid->virulence signal (AUC ~0.79-0.95 in E. coli) is likely LINEAGE-mediated (specific STs, e.g. ExPEC ST131, carry both particular CTX-M plasmids AND UPEC virulence PAIs); within-organism controls species but not sub-lineage. Mash-clade collapse is the proper control.