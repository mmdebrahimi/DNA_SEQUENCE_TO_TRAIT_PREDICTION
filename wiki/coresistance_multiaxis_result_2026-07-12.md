# Multi-axis co-resistance — AMR determinant × plasmid Inc-type (2026-07-12)

**Verdict: PASS_MULTIAXIS_LINKAGE** — 160/163 within-organism features (determinant OR Inc-type) impute from the joint set (fraction 0.982; bar 0.5). (raw)

Do AMR determinants and plasmid Inc-types co-occur / impute each other, within organism? (which resistance rides which plasmid backbone) Substrate: AMR determinants (AMRFinder) + plasmid Inc-replicons (PlasmidFinder blastn, local sweep) over Enterobacterales cohort genomes.

| organism | genomes | Inc-types | testable | **imputable** | frac | note |
|---|---|---|---|---|---|---|
| klebsiella | 307 | 62 | 81 | **81** | 1.0 |  |
| escherichia_coli_shigella | 240 | 55 | 68 | **65** | 0.956 |  |
| salmonella | 60 | 25 | 14 | **14** | 1.0 |  |

## Which resistance rides which plasmid backbone (determinant → Inc-type, by lift)
### klebsiella
- **parC_S80I** → IncI2(Delta)(lift 2.72)
- **gyrA_S83I** → IncI2(Delta)(lift 3.16)
- **blaSHV-11** → IncN(lift 2.19), IncHI1B(pNDM-MAR)(lift 2.06), IncFIB(K)(pCAV1099-114)(lift 1.78)
### escherichia_coli_shigella
- **uhpT_E350Q** → Col156(lift 1.9)
- **tet(A)** → IncN(lift 1.68)
- **sul2** → IncQ1(lift 2.35)
- **aph(3'')-Ib** → IncQ1(lift 2.4)
### salmonella
- **tet(A)** → IncFIB(pN55391)(lift 2.17)
- **sul1** → IncFIB(pN55391)(lift 2.82)
- **gyrA_D87Y** → IncFIB(pN55391)(lift 3.16)
- **aadA1** → IncFIB(pN55391)(lift 3.32)
- **aph(3')-Ia** → IncFIB(pN55391)(lift 2.99)

## Honest caveats
- Plasmid axis is Enterobacterales-only (enterobacteriales.fsa); other organisms excluded.
- PlasmidFinder replicon PRESENCE is a genotype axis (blastn), not a wet plasmid-content readout.
- Cohorts are drug-R/S-selected; within-organism de-confound + dedup clonality proxy; associational.