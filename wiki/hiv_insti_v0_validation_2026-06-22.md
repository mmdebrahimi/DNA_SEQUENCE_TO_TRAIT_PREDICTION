# HIV INSTI v0 cell — validation vs Stanford HIVDB PhenoSense (2026-06-22)

Gene **IN**; position-based catalog at [66, 92, 118, 138, 140, 143, 147, 148, 155, 263].
Label = Stanford HIVDB PhenoSense fold-change (independent wet-lab; NOT Sierra interpretation). Caller = dna_decode.data.hiv_amr v0 (INSTI position-based major-DRM catalog).
Dataset = 861 isolates (`INI_DataSet.txt`).
Catalog source: Stanford HIVDB integrase major drug-resistance positions; cite Rhee 2003 Nucleic Acids Res 31:298-303.

Primary metric (cutoff-free): **AUC = P(fold of a called-R isolate > fold of a called-S isolate)**.

| Drug | n (fold) | called R / S | **AUC** | median fold R | median fold S | sens@f3 | spec@f3 |
|---|---|---|---|---|---|---|---|
| raltegravir | 753 | 374/379 | **0.9053** | 18.85 | 1.0 | 0.939 | 0.814 |
| elvitegravir | 754 | 382/372 | **0.9134** | 34.5 | 1.4 | 0.87 | 0.857 |
| dolutegravir | 370 | 211/159 | **0.745** | 1.6 | 0.9 | 1.0 | 0.534 |
| bictegravir | 287 | 188/99 | **0.8456** | 1.9 | 0.8 | 1.0 | 0.419 |
| cabotegravir | 64 | 60/4 | **1.0** | 9.0 | 0.95 | 1.0 | 0.333 |

## Honest caveats
- v0 is POSITION-BASED -> over-calls non-resistant polymorphisms/revertants at a major position (the INSTI spec quantifies it); mutant-specific deconfounded v0.1 mirrors the NRTI arc
- class-level catalog -> per-drug differential resistance over-calls drugs a mutation spares
- fold>=3 sens/spec is illustrative, NOT a per-drug clinical breakpoint
- in-distribution vs HIVDB-PhenoSense (NOT provenance-disjoint external validation)

Citation: Rhee et al. 2003 Nucleic Acids Res 31:298-303; dataset public per HIVDB Terms of Use.