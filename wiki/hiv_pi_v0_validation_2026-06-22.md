# HIV PI v0 cell — validation vs Stanford HIVDB PhenoSense (2026-06-22)

Gene **PR**; position-based catalog at [30, 32, 33, 46, 47, 48, 50, 54, 76, 82, 84, 88, 90].
Label = Stanford HIVDB PhenoSense fold-change (independent wet-lab; NOT Sierra interpretation). Caller = dna_decode.data.hiv_amr v0 (PI position-based major-DRM catalog).
Dataset = 2171 isolates (`PI_DataSet.txt`).
Catalog source: Stanford HIVDB protease major drug-resistance positions; xref Rhee 2022 Pathogens; cite Rhee 2003 Nucleic Acids Res 31:298-303.

Primary metric (cutoff-free): **AUC = P(fold of a called-R isolate > fold of a called-S isolate)**.

| Drug | n (fold) | called R / S | **AUC** | median fold R | median fold S | sens@f3 | spec@f3 |
|---|---|---|---|---|---|---|---|
| fosamprenavir | 2052 | 1343/709 | **0.8916** | 6.5 | 0.7 | 0.998 | 0.605 |
| atazanavir | 1505 | 999/506 | **0.9571** | 17.0 | 0.9 | 0.999 | 0.712 |
| indinavir | 2098 | 1354/744 | **0.9288** | 12.7 | 0.8 | 0.996 | 0.706 |
| lopinavir | 1807 | 1191/616 | **0.933** | 27.0 | 0.8 | 0.998 | 0.709 |
| nelfinavir | 2133 | 1383/750 | **0.9477** | 24.0 | 1.1 | 0.987 | 0.805 |
| saquinavir | 2084 | 1355/729 | **0.8991** | 8.7 | 0.8 | 1.0 | 0.614 |
| tipranavir | 1226 | 793/433 | **0.7825** | 2.0 | 0.9 | 0.997 | 0.469 |
| darunavir | 993 | 575/418 | **0.8825** | 3.6 | 0.7 | 1.0 | 0.608 |

## Honest caveats
- v0 is POSITION-BASED -> over-calls non-resistant polymorphisms/revertants at a major position (the PI spec quantifies it); mutant-specific deconfounded v0.1 mirrors the NRTI arc
- class-level catalog -> per-drug differential resistance over-calls drugs a mutation spares
- fold>=3 sens/spec is illustrative, NOT a per-drug clinical breakpoint
- in-distribution vs HIVDB-PhenoSense (NOT provenance-disjoint external validation)

Citation: Rhee et al. 2003 Nucleic Acids Res 31:298-303; dataset public per HIVDB Terms of Use.