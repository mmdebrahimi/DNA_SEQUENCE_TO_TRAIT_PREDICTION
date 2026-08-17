# Nitrogen conditional essentiality — RESULTS (2026-08-17)

**Pre-registration:** `wiki/fba_nitrogen_prereg_2026-08-17.md` (frozen before any solve; Amendment 1
recorded before any verdict was read from the amended gate).
**Artifact:** `wiki/fba_conditional_nitrogen_2026-08-17.json` · **Runner:**
`scripts/fba_conditional_nitrogen.py` · **Probe:** `scripts/fba_nitrogen_determinism_probe.py`

Panel: 13 nitrogen sources × 155 two-sided genes (of 1,339 with complete rows), iML1515, Fitness Browser
RB-TnSeq `orgId=Keio`, `fit < -2.0`, `FRAC = 0.01`, `processes=1`.

## Determinism first (the gate that failed, and why)

The first run FAILED its determinism gate: 147 of 2,015 cells. **That failure was real but misdiagnosed
by its own gate.** Probing before changing anything measured the largest disagreement anywhere at
**3.24e-11**, **zero** cells crossing the call line, and the headline metric identical to six decimals
across three passes including one on a freshly loaded model. The differing-cell count was itself unstable
(147 → 58 → 66) — the signature of a 1e-12 cutoff sitting inside float noise. See Amendment 1 for the
replacement gate and why it is net stricter.

**The amended run PASSES**, and the margin design earned itself immediately:

| | |
|---|---|
| call flips across `FRAC=0.01` | **0** |
| largest numerical drift | **1.03e-8** |
| nearest cell to the threshold | **1.00e-2** |
| safety factor (margin ÷ drift) | **9.7e5** (bar 1,000) |
| headline metric, both passes | 0.6799 / 0.6799 |

Note the drift on this run (1.03e-8) is ~300× larger than the probe's (3.24e-11). A **fixed** `max_delta
< 1e-9` bar would have FAILED this run. The derived margin passed it correctly, because what changed was
noise, not any claim.

## The three pre-registered predictions

| # | prediction | carbon | nitrogen | verdict |
|---|---|---|---|---|
| **P1** | ratio distribution bimodal; `FRAC` in a near-empty band | 0 of 16,676 in [0.001, 0.05] | **0.0 %** of 2,015 (bar < 1 %) | **REPLICATES** |
| **P2** | missed essential cells are mostly FLAT (deletion changed nothing) | 76.9 % | **91.07 %** of 448 (bar ≥ 50 %) | **REPLICATES** |
| **P3** | per-cell agreement beats the best-constant null | lift positive | **0.6799 vs 0.5355, lift +0.1444** | **REPLICATES** |

Per-condition MCC spans 0.115 (L-Arginine) to 0.489 (L-Glutamine).

## What this does and does not establish

**P1 is the load-bearing one and it holds on an independent substrate axis.** The emptiness around the
essentiality threshold is therefore a property of **the readout** — a threshold on an FBA growth ratio —
not of the carbon panel. That generalises the mechanism behind the E-Flux null (`a880001`) beyond the
substrate it was measured on, and correspondingly generalises the closure of the graded escape hatch
(`wiki/fba_eflux_graded_RETRACTION_2026-08-17.md`).

**P1 and P2 are NOT independent, and "three independent replications" would overstate this.** Given P1,
every missed-essential cell must have ratio ≥ 0.05, so P2 was always going to be high; its actual content
is narrower — *of the cells the model declines to call essential, 91.07 % are perfectly flat (≥ 0.999)
and the remaining 8.93 % lie in [0.05, 0.999)*. P2 adds "flat, not merely under-impaired"; it does not
add a second independent confirmation of P1.

**The actionable reading of P2:** for 91 % of the cells the model gets wrong, deleting the gene changes
the growth rate by nothing at all. There is no route in the model by which those genes matter under
those conditions. That is a **model-structure** deficit, not a threshold deficit and not a constraint
deficit — consistent with the E-Flux null, and it says that tightening or reweighting flux bounds cannot
address the residual. Missing coupling/regulation is where the remaining error lives.

## Caveats (carried from the pre-registration, unchanged)

1. Several N sources (alanine, serine, aspartate, glutamine) also supply **carbon**. True of the real
   assay too (glucose minimal medium + test compound as sole N source), so these are not nitrogen-only
   perturbations and must not be described as such.
2. 13 of 16 assay sources; two dipeptides and casamino acids have no iML1515 exchange.
3. 2 replicates per source, averaged — thinner than carbon's ~2.5.
4. `best_constant_null` is a weak baseline (majority class). Beating it is necessary, not impressive.

## The control that was missing

The first artifact wrote `"deterministic": false` and then reported `per_cell_agreement`,
`per_condition` and all three P1/P2/P3 verdicts **beside it** — quotable numbers the pre-registration
said may not be reported. A flag standing next to the numbers it invalidates is not a control. The
numbers are now **removed** on a gate failure (`redact_unverified`), pinned by
`test_redact_unverified_removes_every_solve_derived_claim`, which reproduces that exact payload.
