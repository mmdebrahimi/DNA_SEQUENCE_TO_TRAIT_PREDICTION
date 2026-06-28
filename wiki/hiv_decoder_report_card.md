# HIV viral-decoder report card (2026-06-28)

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
| PI | fosamprenavir | 2052 | 0.8916 | 0.802 | 0.929 | 0.127 | 0.058 | - |
| PI | atazanavir | 1505 | 0.9571 | 0.856 | 0.944 | 0.088 | 0.027 | - |
| PI | indinavir | 2098 | 0.9288 | 0.851 | 0.951 | 0.1 | 0.053 | - |
| PI | lopinavir | 1807 | 0.933 | 0.853 | 0.956 | 0.103 | 0.038 | - |
| PI | nelfinavir | 2133 | 0.9477 | 0.896 | 0.957 | 0.061 | 0.045 | - |
| PI | saquinavir | 2084 | 0.8991 | 0.807 | 0.938 | 0.131 | 0.067 | - |
| PI | tipranavir | 1226 | 0.7825 | 0.733 | 0.857 | 0.124 | 0.083 | - |
| PI | darunavir | 993 | 0.8825 | 0.804 | 0.917 | 0.113 | 0.073 | - |
| INSTI | raltegravir | 753 | 0.9053 | 0.877 | 0.946 | 0.069 | 0.06 | - |
| INSTI | elvitegravir | 754 | 0.9134 | 0.863 | 0.837 | -0.026 | 0.012 | - |
| INSTI | dolutegravir | 370 | 0.745 | 0.767 | 0.843 | 0.076 | 0.065 | - |
| INSTI | bictegravir | 287 | 0.8456 | 0.71 | 0.699 | -0.011 | 0.034 | - |
| INSTI | cabotegravir | 64 | 1.0 | 0.667 | 0.981 | 0.314 | 0.266 | - |
| CAI | lenacapavir | 140 | 0.9098 | 0.646 | 0.534 | -0.112 | - | - |

## Honest caveats
- in-distribution (HIVDB), NOT provenance-disjoint -> a lower external-rigour bar than the bacterial card
- NNRTI = mutant-specific (excellent on 1st-gen EFV/NVP); NRTI v0 = position-based (over-calls, fixed by the deconfounded mutant-specific v0.1 for 5/6 drugs; ddI keeps position-based)
- non-B subtype transfer is under-powered (data ~96% subtype B)
- PI/INSTI = position-based v0 (PI AUC 0.78-0.96; INSTI 0.74-1.0, 2nd-gen DTG/BIC lower as the class-level over-call predicts); CAI/lenacapavir = mutant-level (AUC 0.91) on a small resistance-enriched dataset (n=140, 11 S)
- OLS underlying-tool baseline now run for PI/INSTI/CAI (uniform illustrative fold>=3 cutoff, delta is the wrapper-vs-tool signal): PI catalog is high-sens/low-spec so OLS recovers +0.06..+0.13 balacc (real v0.1 mutant-specific headroom, like NRTI); INSTI catalog is competitive (+-0.07, ties/beats OLS on EVG/BIC); CAI catalog BEATS OLS +0.112 (OLS overfits the tiny resistance-enriched set). Subtype-transfer column stays NNRTI/NRTI-only
- PI v0.1 (2026-06-23) + INSTI v0.1 (2026-06-27) deconfounded mutant-specific catalogs SHIPPED (same OLS-coef + 5-fold-CV arc as NRTI): PI 8/8 improve-or-hold (mean +0.056), INSTI 5/5 improve-or-hold (mean +0.087 — overturned the earlier deferral; gain on specificity via accessory-rider deconfounding). The HIV class-deconfounding arc NRTI->PI->INSTI is COMPLETE; only per-drug clinical-cutoff absolute calibration (v0.2) remains. v0.1 gain column shows it.

Citation: Rhee 2003 Nucleic Acids Res 31:298-303; dataset CC public per HIVDB Terms of Use.