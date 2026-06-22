# HIV NNRTI v0 cell — validation vs Stanford HIVDB PhenoSense (2026-06-21)

**The first validated viral cell.** Label = Stanford HIVDB PhenoSense fold-change (independent wet-lab; NOT Sierra interpretation).
Caller = dna_decode.data.hiv_amr v0 (class-level NNRTI major-DRM catalog). Dataset = 2272 isolates (`NNRTI_DataSet.txt`).

Primary metric (cutoff-free): **AUC = P(fold of a called-R isolate > fold of a called-S isolate)** — how well the genotypic call orders isolates by the independent lab phenotype.

| Drug | n (fold) | called R / S | **AUC** | median fold R | median fold S | sens@f3 | spec@f3 |
|---|---|---|---|---|---|---|---|
| efavirenz | 2168 | 1057/1111 | **0.9618** | 35.0 | 0.6 | 0.947 | 0.904 |
| nevirapine | 2052 | 963/1089 | **0.9853** | 100.0 | 0.7 | 0.906 | 0.991 |
| etravirine | 998 | 558/440 | **0.75** | 2.3 | 0.8 | 0.887 | 0.571 |
| rilpivirine | 311 | 185/126 | **0.7001** | 3.5 | 1.25 | 0.822 | 0.544 |
| doravirine | 130 | 93/37 | **0.5554** | 4.8 | 3.9 | 0.75 | 0.348 |

## Honest caveats
- v0 is CLASS-LEVEL -> per-drug differential resistance over-calls ETR/RPV (the AUC gap quantifies it)
- filtered dataset has no Subtype column -> per-subtype transfer check is v0.1 (needs the unfiltered set)
- Stanford R-script least-squares regression is the v0.1 'validate-vs-underlying-tool' baseline (not run here)
- fold>=3 sens/spec is illustrative, NOT a per-drug clinical breakpoint

Citation: Rhee et al. 2003 Nucleic Acids Res 31:298-303; dataset CC public per HIVDB Terms of Use.