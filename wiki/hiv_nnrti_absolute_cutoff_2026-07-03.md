# HIV NNRTI v0.2 absolute-cutoff calibration (2026-07-03)

**Cutoffs SOURCED from Stanford HIVDB DRMcv.R cutoffmat (per-drug clinical lower cutoff) — SOURCED, not fabricated.** 4/5 drugs calibrated; the rest CUTOFF_UNAVAILABLE (external wall).
Catalog = frozen dna_decode.data.hiv_amr (mutant-level). Label = Stanford HIVDB PhenoSense fold (independent wet-lab IC50; NOT Sierra). Filter = Method=PhenoSense AND Type=Clinical; N = 2585.

| Drug | cutoff | all sens/spec/**balacc** (n) | within-B sens/spec/**balacc** (n) |
|---|---|---|---|
| efavirenz | 3.0 | 0.962/0.855/**0.908** (2529) | 0.963/0.861/**0.912** (2444) |
| nevirapine | 3.0 | 0.926/0.962/**0.944** (2423) | 0.926/0.964/**0.945** (2343) |
| etravirine | 3.0 | 0.916/0.54/**0.728** (1175) | 0.918/0.539/**0.728** (1152) |
| rilpivirine | 3.0 | 0.853/0.595/**0.724** (343) | 0.88/0.586/**0.733** (324) |
| doravirine | — | _CUTOFF_UNAVAILABLE (external)_ | — |

## Honest caveats
- cutoffs SOURCED from DRMcv.R (the script the NRTI cell used); DOR + all INSTI are ABSENT from it (postdate integrase inhibitors / doravirine) -> CUTOFF_UNAVAILABLE, reported as a wall not guessed
- PI/INSTI are POSITION-BASED v0 (deliberate over-call at major positions) -> LOW spec at the cutoff is EXPECTED; the mutant-specific v0.1 catalogs (hiv_pi_mutant_catalog) lift it
- for NNRTI every DRMcv.R lower cutoff is 3 == the prior illustrative fold>=3 -> the illustrative choice already matched the clinical cutoff (a confirmation, not a change)
- within-B is the powered de-confound arm (~96% B data)

Citation: Rhee 2003 Nucleic Acids Res 31:298-303; cutoffs from Stanford HIVDB DRMcv.R.