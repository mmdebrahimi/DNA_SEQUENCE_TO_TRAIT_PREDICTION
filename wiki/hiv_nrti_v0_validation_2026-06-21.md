# HIV NRTI v0 cell - validation vs HIVDB PhenoSense + OLS baseline (2026-06-21)

Caller = dna_decode.data.hiv_amr v0 (POSITION-BASED NRTI major-position catalog) (POSITION-BASED - deliberately over-calls).
Label = Stanford HIVDB PhenoSense fold-change (independent wet-lab; NOT Sierra). 1867 isolates; OLS features = 325.

| Drug | n | prev R | cutoff | catalog sens/spec/**balacc** (AUCsep) | OLS **balacc** (AUC) | d(OLS-cat) |
|---|---|---|---|---|---|---|
| lamivudine | 1839 | 0.572 | 5.0 | 0.998/0.638/**0.818** (0.9754) | **0.938** (0.9798) | 0.12 |
| abacavir | 1731 | 0.637 | 2.0 | 0.997/0.779/**0.888** (0.9737) | **0.934** (0.9829) | 0.046 |
| zidovudine | 1853 | 0.45 | 3.0 | 0.998/0.492/**0.745** (0.7793) | **0.907** (0.9688) | 0.162 |
| stavudine | 1846 | 0.454 | 1.5 | 0.992/0.49/**0.741** (0.8166) | **0.868** (0.9372) | 0.127 |
| didanosine | 1849 | 0.525 | 1.5 | 0.979/0.547/**0.763** (0.9192) | **0.805** (0.8927) | 0.042 |
| tenofovir | 1548 | 0.318 | 1.5 | 0.994/0.408/**0.701** (0.6986) | **0.816** (0.8877) | 0.115 |

## Honest caveats
- v0 NRTI is POSITION-BASED -> over-calls T215 revertants / V75 polymorphisms (spec hit, esp AZT/D4T)
- no Subtype column -> per-subtype transfer check is v0.1 (unfiltered set)
- fold>=lower-cutoff binarization uses DRMcv.R's per-drug clinical cutoffs
- v0.1: mutant-specific NRTI catalog (data-derived OLS coefficients / sourced SDRM list)

Citation: Rhee 2003 Nucleic Acids Res 31:298-303; cutoffs from Hedlin/Stanford DRMcv.R.