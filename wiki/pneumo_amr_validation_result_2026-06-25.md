# S. pneumoniae gene-presence AMR cell — validation result (2026-06-25)

Built + validated the **gene-presence** half of the pneumococcus AMR cell (the clean, tractable sub-cell the
go/no-go identified). The β-lactam PBP half stays deferred (separate engine; see
`wiki/pneumo_amr_cell_go_no_go_2026-06-25.md`).

## The cell
- Rule module: `dna_decode/organism_rules/pneumo_amr.py` (NON-FROZEN, `KNOWLEDGE_BASELINE` / `organism_routed`).
  - Macrolide (erythromycin): R iff `erm(B)` OR `mef(A)`/`msr(D)` present.
  - Tetracycline: R iff `tet(M)` OR `tet(O)` present.
- Validator: `scripts/pneumo_amr_validate.py` → `wiki/pneumo_amr_genepresence_validation.json`.
- Tests: `tests/test_pneumo_amr.py` (9, pure-logic). FROZEN E. coli AMR surface byte-unchanged.

## Result — vs WET-LAB measured AST (GPS Poland cohort, CLSI disc-diffusion zone breakpoints)
| Drug | n | accuracy | sensitivity | specificity | confusion |
|---|---|---|---|---|---|
| **Macrolide (erm/mef)** | 127 | **0.961** | 0.968 | 0.954 | TP60 FP3 TN62 FN2 |
| **Tetracycline (tetM)** | 74 | **0.932** | 0.964 | 0.913 | TP27 FP4 TN42 FN1 |

## Honesty tier — FULLY-INDEPENDENT measured-label (confirmed 2026-06-25)
Two arms, both vs the same wet-lab measured AST; the second closes the gap the first left open.
- **Label = WET-LAB measured AST** (disc/agar zone diameters), independent of any caller → clears the
  circularity gate (G1/G3). A real measured phenotype, not a gene-call.
- **Arm 1 (rule-validation):** the headline 0.961/0.932 above used **GPS pipeline determinant calls** as the
  genotype (Supplementary Data 2) — external, not our caller. On its own that arm only **supports rule plausibility**.
- **Arm 2 (fully-independent — below):** OUR AMRFinder end-to-end on the assemblies AGREES with GPS determinants
  1.0/1.0 across n=127, so the same 0.961/0.932 holds with NO external-genotype dependency. **The headline is
  therefore a fully-independent measured-label number.**

### Fully-independent arm — COMPLETE (2026-06-25, full cohort n=127)
The "AMRFinder ≈ GPS BLAST" premise is **measured and confirmed.** `scripts/pneumo_amrfinder_swap.py` ran OUR
AMRFinder (`ncbi/amr` 4.2.7) end-to-end on all 127 GPS assemblies → our own erm/mef/tet determinant calls → re-score.
- **our-AMRFinder-vs-GPS determinant agreement = 1.0 across the full n=127** (both drugs) — perfect; our caller
  and GPS's pipeline call identical erm/mef/tet determinants.
- **FULLY-INDEPENDENT re-score (our caller end-to-end vs measured AST): macrolide acc 0.961 (TP60/FP3/TN62/FN2,
  n=127), tetracycline acc 0.932 (TP27/FP4/TN42/FN1, n=74)** — **identical** to the GPS-determinant rule-validation
  numbers, because the determinant agreement is perfect.
- **→ TIER UPGRADED: this is a FULLY-INDEPENDENT measured-label number** (our caller end-to-end, no GPS genotype
  dependency), not just rule-validation. The C1 honesty gap is fully closed. Artifact:
  `wiki/pneumo_amr_amrfinder_swap_pilot.json` (full n=127).

## What this establishes
- The pneumococcus AMR cell's **gene-presence component is a clean GO** (0.93–0.96 vs measured), confirming
  the go/no-go's split: gene-presence drugs tractable now; β-lactams need the deferred PBP engine.
- A NEW near-independent measured-label AMR validation, on a NEW organism (S. pneumoniae) — distinct from the
  10 frozen E. coli/Klebsiella/etc. cells and from the negative-results map's bacterial-AMR saturation (that
  was NCBI-PD AST on the frozen-engine organisms; this is GPS disc/agar AST on a new organism via a new
  non-frozen rule).

## Deferred (named, not done)
1. **Independent AMRFinder run** (Docker, post-ktype) → fully-independent number.
2. **β-lactam PBP-typing engine** + pneumo MIC breakpoints (meningitis/non-meningitis/oral) — the larger build.
3. **Multi-cohort** beyond Poland n=263 (G4/G8 powering + clonality).
4. **Co-trimoxazole / FQ** — folA/folP + gyrA/parC are point-mutations (not gene-presence); FQ had 0 R in
   this cohort (no signal). Out of the gene-presence v0.
