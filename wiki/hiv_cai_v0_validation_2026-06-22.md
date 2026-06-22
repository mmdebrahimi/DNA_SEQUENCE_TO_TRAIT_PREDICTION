# HIV CAI v0 cell — validation vs Stanford HIVDB PhenoSense (2026-06-22)

Gene **CA**; position-based catalog at [56, 66, 67, 70, 74, 105, 107].
Label = Stanford HIVDB PhenoSense fold-change (independent wet-lab; NOT Sierra interpretation). Caller = dna_decode.data.hiv_amr v0 (CAI position-based major-DRM catalog).
Dataset = 140 isolates (`CAI_DataSet.txt`).
Catalog source: lenacapavir CAPELLA treatment-emergent capsid substitutions (Margot et al. 2022 J Infect Dis 226:1985-1991).

Primary metric (cutoff-free): **AUC = P(fold of a called-R isolate > fold of a called-S isolate)**.

| Drug | n (fold) | called R / S | **AUC** | median fold R | median fold S | sens@f3 | spec@f3 |
|---|---|---|---|---|---|---|---|
| lenacapavir | 140 | 129/11 | **0.9098** | 20.0 | 1.0 | 0.982 | 0.31 |

## Honest caveats
- v0 is MUTANT-LEVEL (12 catalogued substitutions); the CAI dataset is resistance-enriched (treatment-selected isolates) so the susceptible contrast arm is small
- class-level catalog -> per-drug differential resistance over-calls drugs a mutation spares
- fold>=3 sens/spec is illustrative, NOT a per-drug clinical breakpoint
- in-distribution vs HIVDB-PhenoSense (NOT provenance-disjoint external validation)

Citation: Rhee et al. 2003 Nucleic Acids Res 31:298-303; dataset public per HIVDB Terms of Use.