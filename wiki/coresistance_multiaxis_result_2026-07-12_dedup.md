# Multi-axis co-resistance — AMR determinant × plasmid Inc-type (2026-07-12)

**Verdict: PASS_MULTIAXIS_LINKAGE** — 142/149 within-organism features (determinant OR Inc-type) impute from the joint set (fraction 0.953; bar 0.5). (deduped)

Do AMR determinants and plasmid Inc-types co-occur / impute each other, within organism? (which resistance rides which plasmid backbone) Substrate: AMR determinants (AMRFinder) + plasmid Inc-replicons (PlasmidFinder blastn, local sweep) over Enterobacterales cohort genomes.

| organism | genomes | Inc-types | testable | **imputable** | frac | note |
|---|---|---|---|---|---|---|
| klebsiella | 259 | 62 | 78 | **77** | 0.987 | profile-deduped 307->259 |
| escherichia_coli_shigella | 193 | 55 | 65 | **59** | 0.908 | profile-deduped 240->193 |
| salmonella | 35 | 25 | 6 | **6** | 1.0 | profile-deduped 60->35 |

## Which resistance rides which plasmid backbone (determinant → Inc-type, by lift)
### klebsiella
- **parC_S80I** → IncI2(Delta)(lift 2.64)
- **blaSHV-11** → IncN(lift 2.0), IncFII(pKP91)(lift 1.46)
- **gyrA_S83I** → IncI2(Delta)(lift 3.12)
### escherichia_coli_shigella
- **uhpT_E350Q** → Col156(lift 2.0)
- **tet(A)** → IncN(lift 1.86)
- **sul2** → IncQ1(lift 2.54)
### salmonella
- **tet(A)** → IncFIB(pN55391)(lift 1.85)
- **gyrA_D87Y** → IncFIB(pN55391)(lift 2.92)
- **sul1** → IncFIB(pN55391)(lift 2.62)
- **aadA1** → IncFIB(pN55391)(lift 3.15)

## Honest caveats
- Plasmid axis is Enterobacterales-only (enterobacteriales.fsa); other organisms excluded.
- PlasmidFinder replicon PRESENCE is a genotype axis (blastn), not a wet plasmid-content readout.
- Cohorts are drug-R/S-selected; within-organism de-confound + dedup clonality proxy; associational.