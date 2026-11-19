# Fix The Pool Rule Producer Detection And Diagnose The Two Zeros As One Defect

## Lens status
Repo-index: absent (`repo-index-status.sh` not installed) — source-only reading, no degradation.
sentrux: available; plan is NON-structural (additive pure function in one existing module) → architecture gate skipped.
codex: absent (verified exit 127 this session) → no cross-engine review available in-repo.
gh: available.
Wiki replay: 2 HIGH-salience priors in this exact module loaded (see Codebase Context).
project-rules.md: absent → Project Rules Check section omitted.

## Problem Statement

The 2026-11-19 pool-restriction result (`wiki/fba_pool_restriction_generality_2026-11-19.json`,
commit `a79f9f6f`) found that restricting the candidate pool turns an unsearchable 2266-reaction
space into an exhaustive search that finds real growth-coupled designs — 3 of 5 products. Two
products returned zero (acetate, formate) and the artifact attributed formate's zero to the rule
while recording acetate's as **UNDIAGNOSED**.

Diagnosis during planning shows the two zeros are **ONE defect with two faces**, and the defect is
worse than the artifact says:

1. The rule excludes every reaction that "produces" the target, which for a product whose production
   runs through one dominant fermentative route removes the intervention instead of protecting it.
2. **Producer detection itself is invalid.** It was done by stoichiometric sign
   (`r.metabolites[met] > 0`). For `ac_c` that returns 14 reactions and `ACKr` / `PTAr` / `PTA2` are
   **NOT among them**, because they are reversible and written in the consuming direction. So the
   rule both over-excludes (formate: cut `PFL`, the only fermentative route) and mis-identifies
   (acetate: the 4 reactions it excluded are not the 14 it would have detected).

Goal: replace sign-based producer detection with a **flux-role** test — a reaction is load-bearing
for the target iff deleting it collapses the achievable product — then re-validate across all five
products and pin the result.

Non-goals: changing any existing function signature or call path; shipping a per-product pool
default into the CLI; touching the frozen five or `uv.lock`; re-running the MILP arm.

## Codebase Context


See: ./codebase-context.md

## Pre-Change Baseline

- Full suite: **12 failed / 8983 passed / 191 skipped** (686.06s). The 12 failing files are the
  standing absent-data set: cyp4f2 7, essentiality-report-card 3, readme-quickstart 1, pneumo 1.
- `tests/test_fba_design.py`: **21 passed** (78.28s).
- Pool-rule generality, committed: lactate 63 designs / ethanol 30 / succinate 2 / **acetate 0** /
  **formate 0**. Best lactate design `PTAr,ACALD,PFL` at min_product_flux **13.42495**, growth
  **0.115583**; succinate best `LDH_D,ALCD2x,PFL` at **9.263835**.
- Producer detection today: sign-based, `ac_c` → 14 reactions, `for_c` → 11 reactions.
- Active prospective lock: `verify_lock → ok=True, drifted=[], missing=[], incomplete_pin=[]`.

## Verification Signal

- A pure `load_bearing_reactions` / `build_candidate_pool` pair exists in `design.py` and is
  unit-testable without a model where the logic is pure, and against the real model where it is not.
- Running the fixed constructor on all five products re-measures the generality table, and the result
  lands in ONE of three pre-registered branches (Step 3) rather than being interpreted after the fact.
- `tests/test_fba_design.py` grows tests that pin (a) the reversible-reaction trap directly — a
  sign-based test must be shown to MISS `ACKr` for `ac_c` while the flux-role test catches it — and
  (b) that no existing call path changed.
- Full suite failure SET unchanged from the 12-file baseline; `verify_lock` still `ok=True`; the five
  frozen files byte-identical by sha256.

## Implementation Steps

### Step 1: Measure the two zeros as one defect and record the invalid detection
Files: scripts/fba_pool_rule_probe.py, wiki/fba_pool_rule_diagnosis_2026-11-19.json
Depends on: none

**What changes:**
- New read-only probe. For each of the five targets, report: the sign-based producer set, the set my
  committed rule actually excluded, and the symmetric difference between them.
- For each sign-detected producer AND each excluded reaction, run the flux-role test via
  `product_range_at_growth` inside `apply_knockout`: does single deletion collapse `max_product` to
  <= `FLUX_TOL`? Report which reactions are load-bearing by that test.
- Record the measured acetate finding explicitly: `ACKr` / `PTAr` / `PTA2` are absent from the
  sign-detected `ac_c` producer set because they are reversible and written consuming-direction.
- Artifact re-parsed in-script with asserts (the invalid-JSON prior).

**Test strategy:**
- Self-check mode exercising the pure set-difference reporting on synthetic inputs.
- Non-vacuity: assert the symmetric difference is NON-EMPTY for acetate, since an empty one would
  mean the detection and the exclusion agreed and there would be nothing to fix.

### Step 2: Add the flux-role pool constructor to design.py
Files: dna_decode/fba/design.py
Depends on: Step 1

**What changes:**
- Add `load_bearing_reactions(model, target_rxn, candidates, growth_frac, wt_growth)` — pure w.r.t.
  the model (`with model:`); returns the subset whose single deletion drops `max_product` to
  <= `FLUX_TOL`. Built on `product_range_at_growth` + `apply_knockout`, not a new FBA call.
- Add `build_candidate_pool(model, target_rxn, routes, *, growth_frac, wt_growth)` — returns
  `routes` minus `load_bearing_reactions`, plus a report dict naming what was excluded and why.
- **ADDITIVE ONLY**: no existing signature changes, no change to `find_coupled_designs` behaviour.
  The constructor's output is passed BY THE CALLER as the existing `gene_ids` argument.
- Refuse rather than guess: raise if `routes` is empty, and if EVERY route is load-bearing (that pool
  would be empty, and an empty pool silently returns zero designs — the exact shape of the bug).

**Test strategy:**
- Covered by Step 4; Step 2 ships with the `with model:` non-mutation property exercised inline
  (model state identical before/after, checked on reaction bounds).

### Step 3: Re-validate all five products under the fixed constructor
Files: scripts/fba_pool_rule_probe.py, wiki/fba_pool_rule_fix_2026-11-19.json
Depends on: Step 2

**What changes:**
- Extend the Step-1 probe with a `--validate` mode that builds each pool via `build_candidate_pool`
  and runs `find_coupled_designs` with params read VERBATIM from the committed generality artifact
  (anaerobic, level reaction, max_knockouts 3, growth_frac 0.9) — not from this plan's prose.
- **Pre-registered branches, fixed before the run** (the FBA-remediation prior):
  - `fix-recovers-designs` — formate and/or acetate now yield >= 1 coupled design. The rule was the
    wall for those products.
  - `fix-changes-pools-but-not-outcomes` — pools provably differ (reported set diff) yet both still
    yield zero. Then the zeros are a GENUINE negative for those products, because the fermentative
    route was preserved this time, and that is a stronger result than the original zero.
  - `no-pool-change` — the constructor returns the same pools. Then the diagnosis was wrong and must
    be retracted; do not reinterpret.
- Assert the three already-working products (lactate / ethanol / succinate) do not REGRESS.

**Test strategy:**
- Assert the chosen branch is one of the three literals; a fourth value raises.
- Assert the lactate design count stays >= 10 and succinate still returns `LDH_D,ALCD2x,PFL`.

### Step 4: Pin the reversible-reaction trap and the no-regression property as tests
Files: tests/test_fba_design.py
Depends on: Step 3

**What changes:**
- A test that reproduces the DEFECT and shows the fix catches it: for `ac_c`, the sign-based producer
  set must NOT contain `ACKr`, while `load_bearing_reactions` on the fermentative routes MUST return
  a non-empty load-bearing set. Without this the fix has no regression guard against reverting to
  sign detection.
- A test that `build_candidate_pool` raises when every route is load-bearing, rather than returning
  an empty pool that would read as "no designs exist".
- A test pinning the measured Step-3 branch and the three non-regressing products, bars written loose
  against measured values the way the lactate tests are.
- `@pytest.mark.slow` + `pytest.importorskip("cobra")` per file convention.

**Test strategy:**
- `uv run pytest tests/test_fba_design.py -q` green, count strictly above the 21 baseline.
- Non-vacuity on the trap test: confirm it FAILS if the flux-role call is swapped back for the
  sign test.

### Step 5: Reconcile the committed generality artifact with the new finding
Files: wiki/fba_pool_restriction_generality_2026-11-19.json
Depends on: Step 3

**What changes:**
- The committed artifact records acetate as `UNEXPLAINED` and attributes only formate to the rule.
  Add a `superseded_by` block naming `wiki/fba_pool_rule_fix_2026-11-19.json` and state that both
  zeros trace to ONE defect with the sign-detection half unrecorded at the time.
- Do NOT rewrite the original measured fields — the numbers stand; only the interpretation is
  augmented, in a namespace-separate key (the shared-key / augment-only discipline this repo uses
  for the report card).

**Test strategy:**
- Re-parse the artifact and assert the original `results` block is byte-equal to its committed value
  (`git show HEAD:<path>`), so the augmentation provably changed no measured number.

## Execution Preview

- Wave 0: Step 1
- Wave 1: Step 2
- Wave 2: Step 3
- Wave 3: Step 4, Step 5 (different files, both depend only on Step 3)
- Total waves: 4 · max parallelism: 2 · critical path: Step 1 → 2 → 3 → 4
- Recommended mode: **sequential**. Steps 1 and 3 edit the same script, Step 5 edits a committed
  evidence artifact, and `data/` is gitignored so a worktree has no substrate — matching the
  HIGH-salience Oxford precedent that forced sequential mode for contract-bearing files.

## Risk Flags

- **The fix may not recover any design.** Pre-registered branch 2 covers it and treats it as a
  stronger negative, not a failure. This is explicitly modelled on the in-module prior where a
  pre-committed verdict mislabelled a zero.
- **`design.py` is a shipped module.** Mitigated by additive-only: no signature changes, no behaviour
  change to `find_coupled_designs`, pool construction stays caller-driven. It is NOT one of the
  frozen five, verified against the active lock manifest's `surface_sha256`.
- **The flux-role test costs one LP per candidate route.** ~14 routes x 5 targets = ~70 solves, which
  is seconds at the scale measured this session. No MILP, no overnight budget.
- **No unconfirmed external surfaces.** Every surface this plan relies on was read from the code or
  measured during planning, never inherited from prose: the `find_coupled_designs` signature
  (`design.py:391`), `load_model(None)` yielding iML1515 at 2712/1516, the `ac_c` 14-reaction and
  `for_c` 11-reaction sign-detected producer sets, and every baseline number in Pre-Change Baseline.
- **CLAUDE.md is at its word ceiling** (35,459 of 35,460), so this plan does NOT document itself
  there; the artifacts and tests are the record.

## Open Questions

- Whether a documented per-product pool should eventually ship in the CLI remains a USER authority
  call and is explicitly out of scope. This plan makes that decision better-informed by fixing the
  rule first; it does not take it.
- Whether `L_LACD2` / `L_LACD3` (L-lactate, not D-) belong in a D-lactate pool at all. They were in
  the committed run and are retained for comparability; flagged rather than silently changed.

## Verification

1. `uv run pytest tests/test_fba_design.py -q` — green, count > 21.
2. `uv run pytest tests/ -q` — failure SET identical to the 12-file baseline (reconcile the SET and
   its cause, never the count).
3. `python -c "from dna_decode.eval import prospective_lock as pl; _,m=pl.resolve_active_lock(); print(pl.verify_lock(m))"` — `ok=True`.
4. sha256 of the five frozen files unchanged; `uv.lock` unchanged (no dependency added).
5. Every new JSON artifact re-parsed in the script that writes it, with non-vacuity asserts.

<!-- toolkit: check=clean waves=clean gate=fired:open-questions -->
