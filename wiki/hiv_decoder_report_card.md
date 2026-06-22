# HIV viral-decoder report card (2026-06-21)

**Modality:** IN-DISTRIBUTION validation vs Stanford HIVDB PhenoSense fold-change (independent wet-lab IC50). DISTINCT from the bacterial NCBI-PD provenance-disjoint report card — NOT conflated (a different, more external rigour modality).
**Label independence:** PhenoSense fold-change is NOT HIVDB's own Sierra interpretation (non-circular).
**Underlying-tool baseline:** Stanford DRMcv.R OLS, reimplemented (sklearn).

| Class | Drug | n | AUC (call sep. fold) | catalog balacc | OLS balacc | Δ(OLS−cat) | v0.1 gain | subtype transfer |
|---|---|---|---|---|---|---|---|---|
| NNRTI | efavirenz | 2168 | 0.9618 | 0.926 | 0.941 | 0.015 | - | - |
| NNRTI | nevirapine | 2052 | 0.9853 | 0.948 | 0.925 | -0.023 | - | - |
| NNRTI | etravirine | 998 | 0.75 | 0.729 | 0.849 | 0.12 | - | - |
| NNRTI | rilpivirine | 311 | 0.7001 | 0.683 | 0.77 | 0.087 | - | - |
| NNRTI | doravirine | 130 | 0.5554 | 0.549 | 0.767 | 0.218 | - | - |
| NRTI | lamivudine | 1839 | 0.9754 | 0.818 | 0.938 | 0.12 | 0.06 | non-B balacc 0.923 (n=49) |
| NRTI | abacavir | 1731 | 0.9737 | 0.888 | 0.934 | 0.046 | 0.055 | non-B balacc 0.934 (n=42) |
| NRTI | zidovudine | 1853 | 0.7793 | 0.745 | 0.907 | 0.162 | 0.141 | non-B balacc 0.921 (n=49) |
| NRTI | stavudine | 1846 | 0.8166 | 0.741 | 0.868 | 0.127 | 0.122 | non-B balacc 0.877 (n=49) |
| NRTI | didanosine | 1849 | 0.9192 | 0.763 | 0.805 | 0.042 | -0.201 | non-B balacc 0.923 (n=49) |
| NRTI | tenofovir | 1548 | 0.6986 | 0.701 | 0.816 | 0.115 | 0.084 | non-B balacc 0.804 (n=31) |

## Honest caveats
- in-distribution (HIVDB), NOT provenance-disjoint -> a lower external-rigour bar than the bacterial card
- NNRTI = mutant-specific (excellent on 1st-gen EFV/NVP); NRTI v0 = position-based (over-calls, fixed by the deconfounded mutant-specific v0.1 for 5/6 drugs; ddI keeps position-based)
- non-B subtype transfer is under-powered (data ~96% subtype B)

Citation: Rhee 2003 Nucleic Acids Res 31:298-303; dataset CC public per HIVDB Terms of Use.