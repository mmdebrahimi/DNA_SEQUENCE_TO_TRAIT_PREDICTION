# Independent verification — J2 Phase-1 ESM2-650M result (2026-07-05)

Second-reviewer (mosfaer session) READ-ONLY check of DNA-11's committed north-star result
`wiki/j2_phase1_esm2_result_2026-07-05.md`. Recomputed the headline metric from the committed
per-assay JSON (`wiki/j2_phase1_esm2_result_2026-07-05.json`, 201 assays) — no ESM re-run needed.

## Verdict: PASS — exact match, no drift
| metric | claimed (md) | recomputed (JSON, independent) | match |
|---|---|---|---|
| median \|Spearman\| (ESM2-650M vs wet-lab DMS) | 0.4914 | **0.4914** | ✅ |
| shuffled negative control (median) | 0.0136 | **0.0136** | ✅ |
| n assays | 201 | **201** | ✅ |
| PASS ≥ 0.45 | PASS | **True** | ✅ |
| field-match ≥ 0.48 | YES | **True** | ✅ |

## Significance (the through-line, now complete on the real model)
- 2026-07-04 (proxy): AlphaMissense vs DMS median |ρ| **0.417** — substrate proven.
- 2026-07-05 (real): **ESM2-650M itself** median |ρ| **0.491** — the learned/JEPA protein-representation
  thesis is empirically confirmed on our own model, matching the published field number. **J2 Phase 1 DONE.**
- Scope honest + unchanged: molecular-variant-effect layer only; does NOT rescue the complex-organismal
  0-for-5 de-confounded negative.

## Method / honesty
- Read-only; DNA-11's `main` tree untouched (verified). Recompute is a pure median over the committed JSON —
  confirms the md faithfully renders its own data (drift guard), NOT a re-run of ESM (no GPU here).
- What this does NOT re-verify: the ESM inference itself (that ran on Kaggle P100; the JSON is trusted as
  DNA-11's produced artifact) — only that the reported headline == the data it was computed from.
