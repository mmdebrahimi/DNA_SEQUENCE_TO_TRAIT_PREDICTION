# HIV PI v0.2 absolute-cutoff calibration (2026-07-03)

**Cutoffs SOURCED from Stanford HIVDB DRMcv.R cutoffmat (per-drug clinical lower cutoff) — SOURCED, not fabricated.** 8/8 drugs calibrated; the rest CUTOFF_UNAVAILABLE (external wall).
Catalog = frozen dna_decode.data.hiv_amr (position-level). Label = Stanford HIVDB PhenoSense fold (independent wet-lab IC50; NOT Sierra). Filter = Method=PhenoSense AND Type=Clinical; N = 2646.

| Drug | cutoff | all sens/spec/**balacc** (n) | within-B sens/spec/**balacc** (n) |
|---|---|---|---|
| fosamprenavir | 3.0 | 0.998/0.573/**0.786** (2502) | 0.999/0.567/**0.783** (2380) |
| atazanavir | 3.0 | 0.999/0.682/**0.841** (1870) | 0.999/0.678/**0.838** (1827) |
| indinavir | 3.0 | 0.997/0.689/**0.843** (2544) | 0.998/0.691/**0.844** (2422) |
| lopinavir | 9.0 | 0.999/0.58/**0.789** (2265) | 1.0/0.57/**0.785** (2159) |
| nelfinavir | 3.0 | 0.986/0.773/**0.88** (2600) | 0.987/0.779/**0.883** (2477) |
| saquinavir | 3.0 | 1.0/0.578/**0.789** (2560) | 1.0/0.571/**0.785** (2438) |
| tipranavir | 2.0 | 0.987/0.5/**0.743** (1561) | 0.994/0.498/**0.746** (1501) |
| darunavir | 10.0 | 1.0/0.492/**0.746** (1282) | 1.0/0.487/**0.744** (1231) |

## Honest caveats
- cutoffs SOURCED from DRMcv.R (the script the NRTI cell used); DOR + all INSTI are ABSENT from it (postdate integrase inhibitors / doravirine) -> CUTOFF_UNAVAILABLE, reported as a wall not guessed
- PI/INSTI are POSITION-BASED v0 (deliberate over-call at major positions) -> LOW spec at the cutoff is EXPECTED; the mutant-specific v0.1 catalogs (hiv_pi_mutant_catalog) lift it
- for NNRTI every DRMcv.R lower cutoff is 3 == the prior illustrative fold>=3 -> the illustrative choice already matched the clinical cutoff (a confirmation, not a change)
- within-B is the powered de-confound arm (~96% B data)

Citation: Rhee 2003 Nucleic Acids Res 31:298-303; cutoffs from Stanford HIVDB DRMcv.R.