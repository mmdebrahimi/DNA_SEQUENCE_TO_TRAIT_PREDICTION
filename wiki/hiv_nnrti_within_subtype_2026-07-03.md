# HIV NNRTI within-subtype de-confounding check (2026-07-03)

**Verdict: `HOLDS_WITHIN_SUBTYPE`** (median within-B AUC = 0.795; median pooled−within-B = -0.0005; 5 powered drugs).

Catalog = frozen dna_decode.data.hiv_amr (mutant-level, consensus-B numbering). Label = Stanford HIVDB PhenoSense fold (independent wet-lab IC50; NOT Sierra). Metric = cutoff-free AUC = P(fold|called-R > fold|called-S), per subtype group. Filter = Method=PhenoSense AND Type=Clinical; N = 2585.
**Subtype mix:** B = 2482, non-B (pooled) = 97 -> non-B under-powered; within-B is the test arm.

| Drug | all AUC (n) | **B AUC (n)** | non-B AUC (n) | pooled−B |
|---|---|---|---|---|
| efavirenz | 0.9586 (2529) | **0.9614 (2444)** | 0.8731 (85) | -0.003 |
| nevirapine | 0.9795 (2423) | **0.9798 (2343)** | 0.9765 (80) | -0.0 |
| etravirine | 0.7631 (1175) | **0.7636 (1152)** | 0.7311 (23) | -0.0 |
| rilpivirine | 0.7836 (343) | **0.795 (324)** | — (13; under-powered (<15)) | -0.011 |
| doravirine | 0.6043 (127) | **0.5815 (122)** | — (5; under-powered (<15)) | 0.023 |

## What the verdict means
- **`HOLDS_WITHIN_SUBTYPE`** — the deterministic call orders isolates by the independent lab phenotype INSIDE subtype B (AUC materially > 0.5) and the pooled number is not subtype-inflated. The catalog is mechanism, not subtype structure — the same rail that NRTI cleared (2026-06-21).
- **`SUBTYPE_INFLATED`** — the pooled AUC exceeds the within-B AUC by more than 0.1: the class-mixed number was riding subtype structure.

## Honest caveats
- the free HIVDB gp data is ~96% subtype B -> the non-B arm is UNDER-POWERED (a free-data limit, reported not hidden); within-B is the well-powered arm and is the de-confounding test that matters
- cutoff-free AUC needs no clinical breakpoint; the fold>=3 sens/spec is illustrative only
- PI/INSTI are POSITION-BASED v0 (deliberate over-call at major positions) -> their AUC has a built-in ceiling below a mutant-specific catalog; the within-vs-pooled DELTA (not the level) is the de-confounding readout
- a within-B AUC ~ the pooled AUC means the class-mixed number was NOT subtype-inflated; it does NOT prove non-B generalisation at scale (non-B is small)

Citation: Rhee et al. 2003 Nucleic Acids Res 31:298-303; dataset CC public per HIVDB Terms of Use.