# S. pneumoniae β-lactam AMR — rec-1 assessment + wrap-vs-build decision (2026-06-25)

Completed the β-lactam half of the pneumococcus AMR cell (the deferred half from the go/no-go), as far as it
goes WITHOUT a PBP→MIC engine + WITHOUT Docker (ktype occupies it). Built the breakpoint foundation + the
rigorous breakpoint-stratified ceiling, and made the wrap-vs-build decision explicit.

## Built (reusable, no Docker)
- `dna_decode/data/pneumo_breakpoints.py` (NON-FROZEN) — CLSI S. pneumoniae β-lactam breakpoints, keyed by
  **breakpoint CONTEXT** (meningitis / non-meningitis / oral). `classify(drug, context, mic)` → R/I/S. This is
  the breakpoint-discipline foundation any pneumo β-lactam cell needs (the frozen `mic_tiers` is E. coli).
- `tests/test_pneumo_breakpoints.py` (4) — pins the ambiguity (same MIC → R meningitis / S non-meningitis).

## Breakpoint-stratified ceiling (GPS PBP→MIC predicted R/S vs measured MIC)
| Drug @ context | n | acc | sens | spec | R | verdict |
|---|---|---|---|---|---|---|
| **Penicillin @ meningitis** | 255 | **0.973** | 0.979 | 0.968 | 97 | **clean + well-powered** |
| Penicillin @ non-meningitis | 244 | 0.881 | 0.5 | 0.887 | 4 | underpowered R (4 R total — can't assess) |
| Meropenem @ either | 229 | 0.782 | 1.0 | 0.779 | 3 | GPS over-calls (50 FP, 3 R) |

**Finding:** the β-lactam cell is really ONE clean drug-context — **penicillin @ the meningitis breakpoint
(0.973, 97 R)**. Non-meningitis is underpowered (4 R); meropenem is over-called. The acc swing on the SAME
genomes/predictions (0.973 → 0.881) IS the breakpoint-ambiguity, now made explicit by the context-keyed module.

## The wrap-vs-build decision (the cell's crux)
The 0.973 ceiling uses GPS's PBP→MIC predictions. To turn that into OUR cell, two paths:
- **(A) WRAP the CDC/GPS PBP→MIC method** → a `KNOWLEDGE_BASELINE` cell (the tier of TB CRyPTIC + SARS CoV-RDB:
  honest, in-distribution, not an independent baseline). BUT a deployable wrap needs the CDC PBP-type→MIC
  MODEL in-repo — which we do NOT have; using GPS's per-isolate predictions directly adds literally nothing
  (we'd ship GPS's answer). So a pure wrap is not a real decoder.
- **(B) REIMPLEMENT the CDC PBP-type→MIC model** (Li et al. 2017 random-forest / the PBP-type lookup) → a real
  engine that takes pbp1a/2b/2x types → MIC → R/S. THIS is the value-add, but it is the **multi-session build**
  the go/no-go named: a PBP-typer + the regression model + a training/validation split.

## Decision: DEFER the β-lactam cell to a scoped PBP-engine session
- The valuable β-lactam cell = path (B), a real PBP→MIC engine. That is a genuine multi-session build, not a
  tonight increment, and is **Docker-gated** for the PBP-typing step.
- Tonight delivered the **foundation** (breakpoint module + the breakpoint-stratified ceiling proving
  penicillin@meningitis is the one clean well-powered target at 0.973) — so the future build starts from a
  grounded target, not a guess.
- The **gene-presence cell** (macrolide 0.961 / tet 0.932, shipped this session) remains the real near-independent
  pneumo-AMR win; the β-lactam cell is the deferred PBP-engine extension.

## The wall (named)
- **Code-closable (large):** the CDC PBP-type→MIC engine (path B) — needs the published model + a PBP-typer +
  Docker. A scoped future session.
- **Not a tonight task.** Continuing to analyze β-lactams without the engine hits diminishing returns — the
  conclusion has converged: penicillin@meningitis is the target, the PBP engine is the wall. Banked.
