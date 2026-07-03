# HIV PI within-subtype de-confounding check (2026-07-03)

**Verdict: `HOLDS_WITHIN_SUBTYPE`** (median within-B AUC = 0.9209; median pooled−within-B = -0.0055; 8 powered drugs).

Catalog = frozen dna_decode.data.hiv_amr (position-level, consensus-B numbering). Label = Stanford HIVDB PhenoSense fold (independent wet-lab IC50; NOT Sierra). Metric = cutoff-free AUC = P(fold|called-R > fold|called-S), per subtype group. Filter = Method=PhenoSense AND Type=Clinical; N = 2646.
**Subtype mix:** B = 2492, non-B (pooled) = 154 -> non-B under-powered; within-B is the test arm.

| Drug | all AUC (n) | **B AUC (n)** | non-B AUC (n) | pooled−B |
|---|---|---|---|---|
| fosamprenavir | 0.9002 (2502) | **0.9052 (2380)** | 0.7751 (122) | -0.005 |
| atazanavir | 0.9537 (1870) | **0.9541 (1827)** | 0.9589 (43) | -0.0 |
| indinavir | 0.9339 (2544) | **0.9413 (2422)** | 0.7147 (122) | -0.007 |
| lopinavir | 0.9334 (2265) | **0.9367 (2159)** | 0.8189 (106) | -0.003 |
| nelfinavir | 0.9496 (2600) | **0.9554 (2477)** | 0.8047 (123) | -0.006 |
| saquinavir | 0.8974 (2560) | **0.9025 (2438)** | 0.738 (122) | -0.005 |
| tipranavir | 0.7911 (1561) | **0.7985 (1501)** | 0.579 (60) | -0.007 |
| darunavir | 0.8711 (1282) | **0.8787 (1231)** | 0.5956 (51) | -0.008 |

## What the verdict means
- **`HOLDS_WITHIN_SUBTYPE`** — the deterministic call orders isolates by the independent lab phenotype INSIDE subtype B (AUC materially > 0.5) and the pooled number is not subtype-inflated. The catalog is mechanism, not subtype structure — the same rail that NRTI cleared (2026-06-21).
- **`SUBTYPE_INFLATED`** — the pooled AUC exceeds the within-B AUC by more than 0.1: the class-mixed number was riding subtype structure.

## Honest caveats
- the free HIVDB gp data is ~96% subtype B -> the non-B arm is UNDER-POWERED (a free-data limit, reported not hidden); within-B is the well-powered arm and is the de-confounding test that matters
- cutoff-free AUC needs no clinical breakpoint; the fold>=3 sens/spec is illustrative only
- PI/INSTI are POSITION-BASED v0 (deliberate over-call at major positions) -> their AUC has a built-in ceiling below a mutant-specific catalog; the within-vs-pooled DELTA (not the level) is the de-confounding readout
- a within-B AUC ~ the pooled AUC means the class-mixed number was NOT subtype-inflated; it does NOT prove non-B generalisation at scale (non-B is small)

Citation: Rhee et al. 2003 Nucleic Acids Res 31:298-303; dataset CC public per HIVDB Terms of Use.