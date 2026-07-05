# Independent verification audit — DNA-decode load-bearing claims (2026-07-04)

> A **second-reviewer** pass by a parallel session (Soraya), run **read-only** against the DNA-11 session's
> committed work — zero collision (no writes to the shared `main` working tree; every number recomputed from
> committed code/data, bytecode + pytest-cache suppressed). Purpose: independently confirm the load-bearing
> claims reproduce and the **self-graded** certification capstone is faithful to its code, before those claims
> harden. NOT a re-validation of the underlying biology — each domain card keeps its own honest tier.

## Verdict: PASS — no discrepancies found

The decoder's load-bearing claims reproduce independently, and the self-graded capstone is faithful to the
code that generates it (no hand-edit drift) and honest by construction.

## Checks run + results

| # | Check | Method (read-only) | Result |
|---|---|---|---|
| 1 | Load-bearing test suite | `pytest` on 5 headline files (certification_capstone, dms_learned_model_falsifier, abo_blood, cipro_bounded_falsifier, amr_rules) | **59/59 pass** (9.2s) |
| 2 | DMS learned-signal number | Re-ran `scripts/dms_learned_model_falsifier.py --max-assays 40` from committed data | **median \|Spearman\| = 0.417** — exact match to committed `wiki/jepa_dms_learned_signal_result_2026-07-04.md`; field context reproduced (ESM2-650M 0.484, GEMME 0.484, EVE 0.466) |
| 3 | Certification-capstone census (drift guard) | Called `dna_decode.data.cell_registry.cells()` directly; recomputed track + tier counts | **EXACT match** to committed `wiki/certification_capstone.md`: total **67**; track amr25/finder4/pgx3/typing6/viral29; tiers faithful_to_tool11 / independent_measured25 / knowledge_baseline4 / near_independent15 / not_censused1 / no_free_source11 |
| 4 | Capstone honesty posture | Read `certification_capstone.{md,json}` | Honest **by construction**: `no_aggregate_verdict:true`, per-cell tiers preserved, caveats disclosed. The self-validation grading risk does **not** materialize here. |

## Method / honesty notes
- **Zero collision:** no branch switch, no write to the shared `main` tree; `PYTHONDONTWRITEBYTECODE=1` + `-p no:cacheprovider` so no cache/bytecode landed in DNA-11's tree. This audit file lives on the `mosfaer` branch only.
- **Deliberately NOT run:** `scripts/build_certification_capstone.py` (it would OVERWRITE DNA-11's in-flight capstone files) — instead the census was recomputed by calling the registry directly and diffing against the committed artifact.
- **Scope bound (no silent truncation):** verified a **targeted load-bearing subset**, not all **225** test files (full suite exceeds a 10-min run — that is DNA-11's CI concern, out of scope for a claim-verification pass).
- **What this does NOT assert:** the underlying biological validity of each claim (that stays with each domain card's own honest tier); only that the claims **reproduce** from code + data and the capstone **faithfully renders** what the code computes.

## Residual / recommended
- None blocking. If desired, a scheduled full-suite run (all 225 files, chunked to beat the 10-min wall) would upgrade check #1 from "load-bearing subset" to "whole suite".
