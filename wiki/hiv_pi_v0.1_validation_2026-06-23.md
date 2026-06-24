# HIV PI v0.1 - mutant-specific catalog (data-derived, CV held-out) (2026-06-23)

Method: data-derived resistant mutants = MULTIVARIATE OLS log10-fold coefficient >= log10(1.5) (independent >=1.5x effect, controlling for co-occurrence -> deconfounds accessory/polymorphic riders), >=5 carriers; 5-fold CROSS-VALIDATED held-out (derive on train, eval on test) -> NOT in-sample.

**Cutoff (delta-honest):** DELTA-HONEST: both v0 (position-based) and v0.1 (mutant-specific) scored at the SAME UNIFORM illustrative fold>=3 (no per-drug PI clinical cutoff sourced in-repo); the headline is the v0->v0.1 balanced-accuracy GAIN at that fixed cutoff, NOT an absolute-calibration claim. Per-drug PI clinical cutoffs would upgrade to absolute calibration (a v0.2 item; not fabricated here).

Label = Stanford HIVDB PhenoSense fold-change (independent wet-lab; NOT Sierra).

Position-based v0 over-calls accessory/polymorphic riders at the 13 major protease positions; v0.1 keeps only fold-associated mutants. Both metrics are 5-fold CV held-out (out-of-sample). **8/8 PI drugs improve-or-hold; mean balacc gain 0.056.**

| Drug | n | prev R | v0 pos-based sens/spec/**balacc** | v0.1 mutant sens/spec/**balacc** | balacc gain |
|---|---|---|---|---|---|
| fosamprenavir | 2052 | 0.431 | 0.998/0.605/**0.802** | 0.995/0.724/**0.86** | 0.058 |
| atazanavir | 1505 | 0.529 | 0.999/0.712/**0.856** | 0.994/0.773/**0.883** | 0.027 |
| indinavir | 2098 | 0.5 | 0.996/0.706/**0.851** | 0.993/0.815/**0.904** | 0.053 |
| lopinavir | 1807 | 0.521 | 0.998/0.709/**0.853** | 0.997/0.785/**0.891** | 0.038 |
| nelfinavir | 2133 | 0.572 | 0.987/0.805/**0.896** | 0.98/0.901/**0.941** | 0.045 |
| saquinavir | 2084 | 0.43 | 1.0/0.614/**0.807** | 0.998/0.75/**0.874** | 0.067 |
| tipranavir | 1226 | 0.249 | 0.997/0.469/**0.733** | 0.948/0.684/**0.816** | 0.083 |
| darunavir | 993 | 0.307 | 1.0/0.608/**0.804** | 0.99/0.763/**0.877** | 0.073 |

## Deliverable catalog (derived on all data)
- **fosamprenavir** (25): G48M, I47A, I47V, I50V, I54A, I54L, I54M, I54S, I54T, I54V, I84A, I84C, I84V, L33F, L33M, L76V, L90M, M46I, M46L, V32I, V82C, V82F, V82L, V82M, V82S
- **atazanavir** (24): D30N, G48M, G48V, I50L, I54A, I54L, I54M, I54S, I54T, I54V, I84A, I84C, I84V, L90M, M46I, M46L, N88D, N88S, V32I, V82A, V82C, V82F, V82S, V82T
- **indinavir** (22): G48M, G48V, I47A, I54A, I54M, I54S, I54T, I54V, I84A, I84C, I84V, L76V, L90M, M46I, M46L, N88D, N88S, V32I, V82A, V82F, V82S, V82T
- **lopinavir** (23): I47A, I47V, I50V, I54A, I54L, I54M, I54S, I54T, I54V, I84A, I84C, I84V, L33F, L76V, L90M, M46I, M46L, N88D, V82A, V82C, V82F, V82S, V82T
- **nelfinavir** (22): D30N, G48M, I47A, I54A, I54L, I54M, I54S, I54T, I54V, I84A, I84C, I84V, L90M, M46I, M46L, M46V, N88S, V32I, V82A, V82C, V82F, V82S
- **saquinavir** (18): G48M, G48V, I50V, I54A, I54L, I54M, I54S, I54T, I54V, I84A, I84C, I84V, L33I, L90M, N88D, N88S, V82C, V82S
- **tipranavir** (15): I47A, I47V, I54A, I54M, I54S, I54T, I54V, I84V, L33F, L33I, V82C, V82F, V82L, V82S, V82T
- **darunavir** (17): I47A, I47V, I50V, I54A, I54L, I54M, I54S, I54T, I54V, I84V, L33F, L76V, M46L, V32I, V82C, V82F, V82M

## Honest caveats
- DELTA-HONEST at the uniform fold>=3 cutoff (no per-drug PI clinical breakpoint sourced); the GAIN is the signal, not the absolute balacc
- deconfounded (multivariate-OLS-coefficient) derivation, mirroring the shipped NRTI v0.1 arc; 5-fold CV held-out -> out-of-sample, not optimistic in-sample
- in-distribution vs HIVDB-PhenoSense (NOT provenance-disjoint external validation)
- INSTI deferred: its OLS-vs-catalog deltas are thin/mixed (EVG -0.026, BIC -0.011, only CAB +0.314 at n=64 unstable); PI is where the headroom is uniform across all 8 drugs

Citation: Rhee 2003 Nucleic Acids Res 31:298-303; method Stanford DRMcv.R.