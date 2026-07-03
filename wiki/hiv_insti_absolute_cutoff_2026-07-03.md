# HIV INSTI v0.2 absolute-cutoff calibration (2026-07-03)

**Cutoffs SOURCED from Stanford HIVDB DRMcv.R cutoffmat (per-drug clinical lower cutoff) — SOURCED, not fabricated.** 0/5 drugs calibrated; the rest CUTOFF_UNAVAILABLE (external wall).
Catalog = frozen dna_decode.data.hiv_amr (position-level). Label = Stanford HIVDB PhenoSense fold (independent wet-lab IC50; NOT Sierra). Filter = Method=PhenoSense AND Type=Clinical; N = 782.

| Drug | cutoff | all sens/spec/**balacc** (n) | within-B sens/spec/**balacc** (n) |
|---|---|---|---|
| raltegravir | — | _CUTOFF_UNAVAILABLE (external)_ | — |
| elvitegravir | — | _CUTOFF_UNAVAILABLE (external)_ | — |
| dolutegravir | — | _CUTOFF_UNAVAILABLE (external)_ | — |
| bictegravir | — | _CUTOFF_UNAVAILABLE (external)_ | — |
| cabotegravir | — | _CUTOFF_UNAVAILABLE (external)_ | — |

## Honest caveats
- cutoffs SOURCED from DRMcv.R (the script the NRTI cell used); DOR + all INSTI are ABSENT from it (postdate integrase inhibitors / doravirine) -> CUTOFF_UNAVAILABLE, reported as a wall not guessed
- PI/INSTI are POSITION-BASED v0 (deliberate over-call at major positions) -> LOW spec at the cutoff is EXPECTED; the mutant-specific v0.1 catalogs (hiv_pi_mutant_catalog) lift it
- for NNRTI every DRMcv.R lower cutoff is 3 == the prior illustrative fold>=3 -> the illustrative choice already matched the clinical cutoff (a confirmation, not a change)
- within-B is the powered de-confound arm (~96% B data)

Citation: Rhee 2003 Nucleic Acids Res 31:298-303; cutoffs from Stanford HIVDB DRMcv.R.