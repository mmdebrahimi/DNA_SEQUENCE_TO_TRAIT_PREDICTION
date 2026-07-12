# HIV blind-spot position-novelty flag — is out-of-catalog resistance at known positions? (2026-07-11)

**Verdict: FLAG_RECOVERS_BLINDSPOT** — median flag sensitivity on the blind spot = 0.604, median lift = 3.98 (bars: sens>=0.3, lift>=1.2).

Do catalog-MISSED resistant HIV isolates carry an uncatalogued mutation at a KNOWN DRM position (a deterministic position-novelty flag can catch), or is the resistance elsewhere?

Flag: position_novel = a non-WT AA at a catalogued DRM position whose specific substitution is NOT catalogued (a 'catalog call may be incomplete' self-awareness flag, not a resistance predictor). Label = Stanford HIVDB PhenoSense fold-change (free, independent wet-lab).

| drug | catalog-neg | blind-spot (true R) | **flag sens** | flag FP (on S) | base R | flag R | lift |
|---|---|---|---|---|---|---|---|
| efavirenz | 1111 | 53 | **0.604** | 0.105 | 0.048 | 0.224 | 4.691 |
| nevirapine | 1089 | 99 | **0.515** | 0.077 | 0.091 | 0.402 | 4.417 |
| etravirine | 440 | 32 | **0.688** | 0.132 | 0.073 | 0.289 | 3.98 |
| rilpivirine | 126 | 21 | **0.571** | 0.267 | 0.167 | 0.3 | 1.8 |
| doravirine | 37 | 21 | **0.714** | 0.375 | 0.568 | 0.714 | 1.259 |

## Honest caveats
- RISK-FLAGGED family: the LEARNED/likelihood blind-spot rescue is a CLOSED negative (ESM below chance; blind spot is pocket-mediated per the ddG probe). This tests only the DETERMINISTIC position-novelty angle.
- A BLINDSPOT_NOT_POSITION_LOCAL verdict CONFIRMS the closed negative (out-of-catalog resistance is not at known positions) — a valid finding, not a failure to build.
- The flag is a self-awareness DIAGNOSTIC ('catalog may be incomplete here'), NOT a resistance call.
- fold>=3 is the illustrative NNRTI cutoff (not a per-drug clinical breakpoint).

Citation: Rhee 2003 Nucleic Acids Res 31:298-303.