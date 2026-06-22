# HIV NRTI v0.1 - mutant-specific catalog (data-derived, CV held-out) (2026-06-21)

Method: data-derived resistant mutants = MULTIVARIATE OLS log10-fold coefficient >= log10(1.5) (independent >=1.5x effect, controlling for co-occurrence -> deconfounds revertants), >=5 carriers; 5-fold CROSS-VALIDATED held-out (derive on train, eval on test) -> NOT in-sample.
Label = Stanford HIVDB PhenoSense fold-change (independent wet-lab; NOT Sierra).

Position-based v0 over-calls revertants; v0.1 keeps only fold-associated mutants. Both metrics are 5-fold CV held-out (out-of-sample).

| Drug | n | prev R | v0 pos-based sens/spec/**balacc** | v0.1 mutant sens/spec/**balacc** | balacc gain |
|---|---|---|---|---|---|
| lamivudine | 1839 | 0.572 | 0.998/0.638/**0.818** | 0.99/0.765/**0.878** | 0.06 |
| abacavir | 1731 | 0.637 | 0.997/0.779/**0.888** | 0.981/0.905/**0.943** | 0.055 |
| zidovudine | 1853 | 0.45 | 0.998/0.492/**0.745** | 0.994/0.778/**0.886** | 0.141 |
| stavudine | 1846 | 0.454 | 0.992/0.49/**0.741** | 0.912/0.813/**0.863** | 0.122 |
| didanosine | 1849 | 0.525 | 0.979/0.547/**0.763** | 0.144/0.98/**0.562** | -0.201 |
| tenofovir | 1548 | 0.318 | 0.994/0.408/**0.701** | 0.972/0.598/**0.785** | 0.084 |

## Deliverable catalog (derived on all data)
- **lamivudine** (12): K65R, K70E, K70G, M184I, M184V, Q151M, T215E, T215F, T215I, T215L, T215Y, V75T
- **abacavir** (8): K65R, L74V, M184I, M184V, Q151M, T215F, T215Y, V75T
- **zidovudine** (12): K70R, K70T, L210W, M41L, Q151M, T215F, T215I, T215V, T215Y, V75A, V75I, V75M
- **stavudine** (8): K65R, K70N, L210W, Q151M, T215F, T215Y, V75A, V75I
- **didanosine** (4): K65R, M184I, Q151M, V75T
- **tenofovir** (8): K65R, K70N, K70R, K70T, L210W, Q151M, T215F, T215Y

Citation: Rhee 2003 Nucleic Acids Res 31:298-303; cutoffs from DRMcv.R.