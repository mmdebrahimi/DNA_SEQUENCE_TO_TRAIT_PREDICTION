# Brainstorm review — calibrate_organism + wiring + promotion gate (2026-06-08)

> Design review of `dna_decode/eval/calibrate_organism.py`, its opt-in wiring into `call_resistance`, the
> committed registry, and the independent-cohort promotion gate. All findings verified grounded against the
> code. Fixes applied 2026-06-09 (ledger rows 41, 43; commits 96bffc4, efcf4dd). Persisted for reference.

## Critical issues (fixed)

**C1 — Pseudomonas registry entry was a degenerate one-class verdict (n_R=0, n_S=7).** Two root causes:
`build_calibrated_registry._load_cohort` read only `data/raw/<slug>/amrfinder_runs` (7 main.tsv there)
while the validator resolves runs via `reuse_glob = data/raw/<genus>_*/amrfinder_runs` (the other 23, cached
in a sibling cohort dir) → under-load to 7; and those 7 were all-S, with no min-class guard in `calibrate()`
→ a bogus `EXPRESSION_FLOOR` on one-class data. **Fix:** builder loads via `reuse_glob`; `calibrate()`
returns a new `INSUFFICIENT_EVIDENCE` verdict when either scored class < `MIN_CLASS_COUNT` (5), BEFORE the
abstain-floor check; registry rebuilt (all 5 entries now on valid 15R/15S; Pseudomonas → valid
EXPRESSION_FLOOR at LOO 0.667).

**C2 — `loo_balanced_accuracy` was plain accuracy.** `loo_correct/len`, misnamed. **Fix:** compute true
balanced accuracy from held-out LOO per-class predictions; `None` when a class is absent.

## Medium issues (fixed, except M2 deferred)

**M1 — LOO estimated a different procedure than the deployed rule** (per-fold intrinsic + modal config vs
full-cohort intrinsic; `Counter.most_common` tie-break ≠ the documented `_select_best_config` order).
**Fix:** deployed config = deterministic full-cohort `_select_best_config`; LOO separately estimates the
selection procedure (no modal/tie-break ambiguity).

**M3 — promotion gate lacked a specificity floor + config-match was too brittle.** **Fix:** gate =
OOS acc + sens + **spec** ≥ 0.80 AND ≥ 10 scored/class; config-match is a **flag, not a gate**
(non-inferior OOS performance suffices). Empirically confirmed by the overnight run: Salmonella's
re-calibration picked a different equally-good config yet the deployed rule scored 1.0 OOS.

**M2 — AMRFinder `Method` column dropped (point status via `"_"` proxy). DEFERRED.** The `"_"` proxy holds
on every observed AMRFinder symbol (acquired genes carry no `_`); documented in `is_point_mutation`. Future
hardening: carry `Method`=POINT through `features_from_main_tsv` (a feature-shape change).

## Verdict taxonomy (adopted)

`CALIBRATED` (powered + LOO bal-acc ≥ 0.70) / `EXPRESSION_FLOOR` (powered but sub-floor — a mechanistic
abstain) / `INSUFFICIENT_EVIDENCE` (a class < MIN_CLASS_COUNT). Keeps `EXPRESSION_FLOOR` a real biological
claim rather than a dumping ground for under-powered cohorts.

## What was solid (survived review)

- The core thesis (counter + threshold + intrinsic-exclusion are organism-specific; abstain where presence
  can't decode) — unchanged.
- The opt-in architecture — every issue was confined to the opt-in registry path; the default `DRUG_RULE`
  was never touched (tests stayed green throughout), which is exactly why opt-in-until-validated was right.
- Family-level intrinsic grouping + excluding point mutations from intrinsic detection — correct.

## Open tradeoffs (resolved)

- The 0.70 abstain floor is post-hoc on 5 cohorts; the 3-verdict split resolves it (floor applies only once
  a cohort is balanced + powered; weak cohorts route to INSUFFICIENT_EVIDENCE).
- Promotion-to-default of the 3 OOS-validated cipro configs: **NOT done** — once the CLI forwards
  `--organism`, "promotion to default" has no clean target (forcing `DRUG_RULE` itself to be organism-aware
  would corrupt the validator's baseline-vs-calibrated comparison). The CLI forwarding delivers the value;
  the registry stays opt-in pending a different-lab cohort.
