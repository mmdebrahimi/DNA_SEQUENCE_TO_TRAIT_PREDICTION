# Overnight Status — 2026-05-22 → 2026-05-23

> Worked through Phases 1-6 while you slept. Committed locally; nothing pushed. 10 new commits, 664 tests green (+295 vs session start), zero regressions.

## TL;DR

| Item | Status |
|---|---|
| Falsifier coordination plan + 12-strain subset + leakage check | ✓ committed (earlier in session: `ae97ca9` + `a7424b0`) |
| Falsifier runner draft for Codex | ✓ committed `1e2ad95` + 15 contract tests `2755385` |
| Post-falsifier ship-path technical plan + /review + /brainstorm edits | ✓ committed `1e72c2c` + `bd6fb4f` |
| Cohort dedup assertion (real bug fix, prevents future leakage) | ✓ committed `5b0eae0` (6 regression tests) |
| `attribution_scope_confidence` field in `pipeline.py predict` | ✓ committed `773e8b0` (12 new tests) |
| Mash-cluster orchestration script (PASS-path artifact) | ✓ committed `d2003e9` (8 pure-logic tests) |
| Untracked test files committed for hygiene | ✓ committed `827387c` (98 tests) |
| Falsifier results from Codex on Precision 7780 | ⏳ awaiting transfer |

## Commits this session (most recent first)

```
d2003e9 draft(mash-cluster): N=147 orchestration script + 8 pure-logic tests
2755385 tests(cipro-falsifier): 15 contract tests + tightened _ranked_by docstring
773e8b0 feat(predict-v0): attribution_scope_confidence field + locus-tag-prefix proxy
5b0eae0 fix(cohort): assert assembly_accession uniqueness in build_cohort + regression tests
827387c tests(predict-v0 + mic-tiers): commit existing 98 tests previously untracked
bd6fb4f plan(cipro-post-falsifier): apply /review + /brainstorm edits in-place
1e72c2c plan(cipro-post-falsifier): verdict-conditional ship-path technical plan
1e2ad95 draft(cipro-falsifier): runner skeleton for Codex to diff against
a7424b0 ledger(cipro-falsifier): 3 LESSONS + Bellman refresh + decisions log + plans index
ae97ca9 plan(cipro-falsifier): coordination plan + 12-strain subset + leakage check snippet
```

## What the session moved forward

### 1. The `attribution_scope_confidence` field IS now in `pipeline.py predict`

Pre-falsifier it defaults to `INDETERMINATE`. Post-falsifier (when Codex sends results), the helper `_classify_attribution_scope(prefix, saturated, all_negative_delta, falsifier_verdict)` takes the verdict + uses the locus-tag-prefix proxy (ERS = HIGH; ELX/ELY/ELV/ELU/ELT = PARTIAL; saturated/negative-delta = INDETERMINATE) — without needing Mash-cluster.

The field is additive — backward-compatible with existing predict JSON consumers. Markdown sidecar surfaces it under the prediction header.

### 2. The cohort same-genome leakage bug class is closed

`build_cohort` now raises `CohortConstructionError` if any two candidate strains share an `assembly_accession` (excluding empty strings). Override via `allow_duplicate_accessions=True` for intentional cases. This prevents the `GCA_025200635.1` (= `562.109860` AND `562.111036`) issue from recurring on future cohort builds. The existing parquet is unaffected.

### 3. The bounded-falsifier runner has 15 contract tests pinning the verdict matrix

If Codex's runner on the Precision 7780 emits different verdict semantics than the Claude draft, my test suite will catch divergence at transfer time (per the post-falsifier plan's Step 0 sub-step 2.5 verdict-reconciliation gate).

### 4. The Mash-cluster orchestration script is ready for Codex to run

`scripts/mash_cluster_n147.py` is fully tested for pure logic (threshold sweep + clade quality scoring + R/S balance). Codex on the Precision 7780 invokes it with `--use-docker` to route Mash through Docker. Threshold sweep avoids the brainstorm B5 concern about borrowing 0.05 without justification.

### 5. Post-falsifier ship-path plan covers all 4 verdict branches × leakage gate

`plans/Cipro_Post_Falsifier_Ship_Path_Technical_Plan.md` (449 LOC) maps every verdict × gate state combination to a concrete sequence of code + doc + commit changes. Verdict-conditional response is pre-committed BEFORE results land — prevents outcome-biased decisions (the 2026-05-14 verdict-vs-budget LESSON).

The plan went through /review + /brainstorm with 11 issues caught + applied:
- D3 reworked: saturation + lineage-confound are co-causes, not either/or
- D3a added: PASS + saturation = PASS-WITH-CALIBRATION-NOTE sub-branch
- Step F's confidence field uses locus-tag-prefix proxy (now implemented per #1 above)
- "REVERT" branch renamed to INDETERMINATE_BUCKET_C to match the runner's actual `compute_verdict` semantics
- Step P sub-step 4 softened — falsifier doesn't measure cohort-wide stats
- Step F sub-step 7.5 added — append "Falsifier Resolution" coda to `cipro_interp_audit_analysis_2026-05-22.md`
- Plus 5 more (brainstorm B1-B9, review L1-L4)

## What's blocked on Codex (Precision 7780)

The bounded-falsifier executes on the Precision 7780 (model + cache + cohort live there). Codex needs to:

1. Pull `git pull origin main` to get all 10 session commits, including my runner draft + leakage check + plans.
2. Run `scripts/leakage_check_dup_accession.py` (< 5 s gate).
3. Diff their runner against `scripts/cipro_bounded_falsifier.py`; adopt-with-edits or replace.
4. Run the falsifier with `--leakage-check-json` wired in.
5. Transfer `Downloads/cipro_bounded_falsifier_results_2026-05-22.{md,json}` (or 2026-05-23 if it runs today) back to this laptop.

When that lands, the post-falsifier plan's Step 0 → branch dispatch → Step P/F/R/V fires. Estimated ~1-1.5 hr wallclock on the PASS path (Mash-cluster + per-clade LOSO + tag); ~45 min on the FAIL path (scope-limit doc + ledger).

## Test-suite state

- 664 tests passing on Windows 10 / Python 3.11.5 / pytest 9.0.3 in 61 seconds
- +295 tests since session start (98 pre-existing untracked + 6 cohort dedup + 12 attribution scope confidence + 15 falsifier contract + 8 mash-cluster + others previously committed)
- Zero regressions across the 10 commits

## Decisions I made for you while you slept

(Per "move forward as much as you can." All are reversible.)

1. **Committed the 98 pre-existing untracked tests as a standalone commit before adding new ones.** Hygiene — wouldn't want one big commit conflating "existing tests" with "new tests."
2. **Added `allow_duplicate_accessions=False` as the safe default** in `build_cohort`. Existing callers get the new behavior automatically; if anything fails, it's a real bug surfaced. Override available.
3. **Implemented `attribution_scope_confidence` pre-falsifier defaulting to INDETERMINATE.** Field is now in the schema. When the falsifier lands, `_classify_attribution_scope` gets called with a real verdict. Zero behavior change pre-falsifier — only the JSON gets an extra honest field.
4. **Mash-cluster threshold sweep uses 6 candidates `[0.02, 0.03, 0.04, 0.05, 0.07, 0.10]`** per brainstorm B5; fallback is 0.05 (matches plan default). Pure logic; Codex's run will pick whichever the real distance matrix satisfies.
5. **The `_ranked_by` "bug" I caught in /review was actually a non-bug** — traced through the consumers, all gate on `delta > 0`. Fixed via tightened docstring + 15 contract tests rather than logic change.
6. **Did NOT push anything to origin.** Per the project's commit-only-on-main pattern + the plan's "tag local-only until user confirms ship" risk flag.

## What to do in the morning

1. Read this status report (you are here).
2. Skim the 10 new commits (`git log --oneline 72d04dd..HEAD`) — they're individually meaningful + small.
3. Decide whether to push to origin or wait until Codex's results land.
4. If Codex's results are in `Downloads/`, kick off the post-falsifier dispatch per `plans/Cipro_Post_Falsifier_Ship_Path_Technical_Plan.md` Step 0.

If anything I did looks wrong, every commit is on `main` after `72d04dd` — safe to revert individually.
