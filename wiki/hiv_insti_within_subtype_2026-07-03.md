# HIV INSTI within-subtype de-confounding check (2026-07-03)

**Verdict: `HOLDS_WITHIN_SUBTYPE`** (median within-B AUC = 0.8981; median pooled−within-B = 0.0434; 4 powered drugs).

Catalog = frozen dna_decode.data.hiv_amr (position-level, consensus-B numbering). Label = Stanford HIVDB PhenoSense fold (independent wet-lab IC50; NOT Sierra). Metric = cutoff-free AUC = P(fold|called-R > fold|called-S), per subtype group. Filter = Method=PhenoSense AND Type=Clinical; N = 782.
**Subtype mix:** B = 640, non-B (pooled) = 39 -> non-B under-powered; within-B is the test arm.

| Drug | all AUC (n) | **B AUC (n)** | non-B AUC (n) | pooled−B |
|---|---|---|---|---|
| raltegravir | 0.9317 (689) | **0.9111 (629)** | — (12; under-powered (<15)) | 0.021 |
| elvitegravir | 0.9428 (708) | **0.8937 (596)** | — (12; under-powered (<15)) | 0.049 |
| dolutegravir | 0.9133 (174) | **0.8509 (84)** | 0.9709 (39) | 0.062 |
| bictegravir | 0.9402 (144) | **0.9024 (29)** | 0.9205 (15) | 0.038 |
| cabotegravir | 1.0 (66) | **— (4; under-powered (<15))** | — (10; under-powered (<15)) | — |

## What the verdict means
- **`HOLDS_WITHIN_SUBTYPE`** — the deterministic call orders isolates by the independent lab phenotype INSIDE subtype B (AUC materially > 0.5) and the pooled number is not subtype-inflated. The catalog is mechanism, not subtype structure — the same rail that NRTI cleared (2026-06-21).
- **`SUBTYPE_INFLATED`** — the pooled AUC exceeds the within-B AUC by more than 0.1: the class-mixed number was riding subtype structure.

## Honest caveats
- the free HIVDB gp data is ~96% subtype B -> the non-B arm is UNDER-POWERED (a free-data limit, reported not hidden); within-B is the well-powered arm and is the de-confounding test that matters
- cutoff-free AUC needs no clinical breakpoint; the fold>=3 sens/spec is illustrative only
- PI/INSTI are POSITION-BASED v0 (deliberate over-call at major positions) -> their AUC has a built-in ceiling below a mutant-specific catalog; the within-vs-pooled DELTA (not the level) is the de-confounding readout
- a within-B AUC ~ the pooled AUC means the class-mixed number was NOT subtype-inflated; it does NOT prove non-B generalisation at scale (non-B is small)

Citation: Rhee et al. 2003 Nucleic Acids Res 31:298-303; dataset CC public per HIVDB Terms of Use.