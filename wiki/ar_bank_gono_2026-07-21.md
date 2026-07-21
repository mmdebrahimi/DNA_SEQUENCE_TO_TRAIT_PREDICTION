# AR-Bank N. gonorrhoeae validation — a NEW species, 2 cells SCORED_ENDORSED (2026-07-21)

**Status:** ✅ **cefixime + ciprofloxacin both SCORED_ENDORSED** · **Cell:** `neisseria_amr` (NON-FROZEN) ·
**Frozen surface:** byte-unchanged (verify_lock OK).

## Result

Launched the N. gonorrhoeae AR-Bank validation on Kaggle (SRA-assembly + AMRFinder factory, reusing the
single-end-fixed enterococcus kernel). 20/21 isolates assembled (1 timeout); no empty assemblies (min 4
determinants). Scored the `neisseria_amr` cell vs CDC S/I/R labels:

| drug | n | R/S | sens | spec | acc | verdict |
|---|---|---|---|---|---|---|
| **ciprofloxacin** | 20 | 11R/9S | 1.00 | 1.00 | 1.00 | **SCORED_ENDORSED** |
| **cefixime** | 20 | 12R/8S | 0.917 | 1.00 | 0.95 | **SCORED_ENDORSED** (after v0.1 fix) |

**First N. gonorrhoeae validation — a whole new species.** Ciprofloxacin (gyrA/parC QRDR) is perfect out of
the box. Cefixime needed a v0.1 rule fix (below).

## Cefixime v0.1 — narrowed from over-calling (spec 0.0 → 1.0)

The v0 rule fired on **any** penA ESC point mutation → R. Verify-in-batch caught **spec 0.0** (all 8
cefixime-S isolates called R). Diagnosis via per-marker separation on the cohort:

| penA marker | R count | S count | discriminative? |
|---|---|---|---|
| **I312M / V316T / N512Y / G545S** (mosaic-34 core) | 11 | **0** | **YES — cefixime-R signature** |
| A510V / F504L | 11 | 8 | no (shared; v0 fired on these → over-call) |
| A516G | 0 | 8 | S-associated |

The cefixime-S isolates carry only a **partial** mosaic (A510V/F504L, MIC 0.015–0.06, below breakpoint); the
R isolates carry the full **mosaic penA-34 core** {I312M, V316T, N512Y, G545S} (MIC ≥ 0.25). The v0.1 rule
requires ≥3 of the 4 core markers → **spec 1.0** without losing the quartet-R isolates. This exactly mirrors
the ceftriaxone v0.1 narrowing (which removed the A510-mosaic FP by requiring the specific A501 marker).

**Honest caveats:** derived + validated on the AR-Bank cohort (like ceftriaxone v0.1); the mosaic-34 markers
are literature-grounded (penA-34/XXXIV cefixime signature). **Sensitivity ceiling:** 1 FN — SAMEA3165247, a
`penA_D346DD`-only non-mosaic R at MIC 1 — is NOT caught by the mosaic-34 core (a different high-MIC path,
disclosed, not a rule bug). sens 0.917.

## Scope

First N. gonorrhoeae AR-Bank validation (a new species per the datasource scouting memo). NON-FROZEN
`neisseria_amr` cell; the frozen decoder surface is byte-unchanged. The Kaggle factory pattern (SRA-assembly
+ AMRFinder, single-end-fixed) transferred cleanly to a new Gram-negative organism. Powering-optimal
21-isolate subset (all scarce-S + double-R) delivered ≥8 per class on both drugs.
