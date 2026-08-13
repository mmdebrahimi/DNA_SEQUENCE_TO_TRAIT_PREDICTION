# FBA Conditional-Essentiality Remediation — Solver-Status Audit, Condition-Parametric Metrics, and the Stratified Commit Claim

## Lens status

- **Correctness:** primary lens. One published claim (pFBA regulatory lift) has an identified mechanism by which it could be a solver artifact; one published sentence is factually false.
- **Risks/Regressions:** secondary. All changes are additive to non-frozen research code; the frozen AMR decoder surface is untouched.
- **External-tool surface:** verified — cobrapy 0.31.1 `single_gene_deletion` frame columns confirmed from upstream source, not inherited.
- **Complexity:** low. One new pure module, three signature extensions, four script wirings, two re-runs.
- **Test gaps:** addressed — every pure change is unit-testable offline without `feba.db` or a solver.

## Problem Statement

An adversarial review of the corrected FBA conditional-essentiality claim set found that two published claims rest on unaudited foundations and one published sentence is false:

1. **The pFBA regulatory lift may be a solver artifact.** `scripts/fba_regulatory_conditional_test.py:106` codes a NaN growth as ratio `0.0`, and `0.0 <= ESSENTIAL_FRAC (0.01)` means "essential". cobrapy returns NaN exactly when a deletion solve is non-optimal or infeasible. The restricted arm forces off ~1,865 of 2,712 gene-associated reactions (69%) before deleting genes — the condition in this codebase most likely to drive LPs infeasible — and never reads the `status` column. The observed result has precisely the shape that artifact would produce: cells called essential 16 → 96, TP 10 → 56, FP 6 → 40, precision *down* 0.625 → 0.583, threshold-free AUROC *down* 0.612 → 0.576. The existing rate-matched null cannot detect this, because it scores labels only and never invokes the solver.

2. **A published honest-limits sentence is false.** `wiki/fba_conditional_carbon_2026-08-12.md:89-91` states "the t-statistic that the loader reads is not used." `dna_decode/fba/fitness_browser.py:87` and `:122` select `g.sysName, f.expName, f.fit` — the loader never reads it. `PRAGMA table_info(GeneFitness)` on the live 7.4 GB `feba.db` returns `['orgId','locusId','expName','fit','t']`, so the column exists and is one word away.

3. **Two more metric helpers carry the hardcoded-condition assumption that caused the retraction, and they fail silently.** `continuous_readout` (`conditional_essentiality.py:164`) and `deployable_threshold` (`:236`) both do `keys = sorted(CONDITIONS)` and accept no `conditions` parameter. Both guard with `ratios.get(c, {})`, so on a 25-source ratios dict every 4-media key misses, zero cells accumulate, and the function returns `{"auroc": None, "note": "degenerate: one class only"}` — which reads as "your data was degenerate", not "you passed the wrong keys". `rate_matched_null` (`fba_regulatory_conditional_test.py:59`) is hardcoded twice over: explicitly at `:59`, and again by falling through `switch_accuracy`'s 4-media default at `conditional_essentiality.py:371`.

4. **The commit-rate headline is over-compressed, and the non-optimal cells are uncrossed.** 39/5,425 non-optimal cells is a footnote for aggregate Claim A, but one non-optimal cell can create or destroy a 1-of-25 exact-set match, and 25 of the 33 commitments are exactly that shape.

Scope boundary: this plan remediates and re-derives. The full label-threshold sensitivity sweep is enabled here (Step 3 makes `GeneFitness.t` reachable) but the sweep itself is a follow-on research run, not part of this plan.

## Codebase Context

**Module:** `dna_decode/fba/` — 13 modules, non-frozen research code. The frozen AMR decoder surface (`amr_rules.py`, `calibrated_amr_rules.json`, `mic_tiers.py`, `shipped_decoder_surface.py`, `cohort_manifest.py`) is not touched by any step.

**Metric helpers in `conditional_essentiality.py`** — 4 of 6 are condition-parametric, 2 are not:

| helper | line | condition-parametric? |
|---|---|---|
| `pattern_distribution` | 291 / 310 | yes (`conditions` param, default 4-media) |
| `constant_baselines` | 340 / 351 | yes |
| `switch_accuracy` | 360 / 371 | yes |
| `confusion_from_calls` / `mcc` | 122 / 132 | n/a — condition-free |
| `continuous_readout` | 139 / 164 | **no — no param at all** |
| `deployable_threshold` | 220 / 236 | **no — no param at all** |

**Deletion scripts and their status-audit state** (`"status"` reference counts):

| script | conditions | audits status |
|---|---|---|
| `scripts/fba_conditional_carbon_validate.py` | 25 carbon | yes (counts only, `:105-110`) |
| `scripts/fba_conditional_essentiality_validate.py` | 4 media | yes (counts only) |
| `scripts/fba_regulatory_conditional_test.py` | 4 media | **no** |
| `scripts/fba_gapfill_carbon_recheck.py` | 25 carbon | **no** |
| `scripts/fba_gapfill_conditional_test.py` | 4 media | **no** |

**NaN-to-essential coding sites** (both grounded by direct read):
- `fba_regulatory_conditional_test.py:106` — `d[gid] = 0.0 if g != g else g / wt`, then `calls[c] = {g: v <= ESSENTIAL_FRAC}` at `:108`.
- `fba_gapfill_carbon_recheck.py:71` — `d[gid] = (g != g) or (g < FRAC * wt)`.

**External tool surface — verified, not inherited.** cobrapy `0.31.1`; `cobra.flux_analysis.deletion._multi_deletion` constructs its result with `columns=["ids", "growth", "status"]` and documents `status : str — The solution's status.` Confirmed by reading the installed package source via `inspect.getsource`, not from docs or a prior reference. The `ids` entry is a **frozenset**, not a bare id — existing consumers use `next(iter(row["ids"]))`.

**Test infrastructure:** pytest, 11 FBA test files. `tests/test_fba_conditional_essentiality.py` (31 tests) covers the metric helpers including the two regression tests pinning the hardcoded-4 fix. `tests/test_fba_regulatory_conditional.py` (4 tests) covers `rate_matched_null` as a pure function. Real-data tests skip when `feba.db` is absent — every change below is unit-testable offline.

### Reusable-Code Survey

- **`scripts/fba_conditional_carbon_validate.py:105-110`** — the existing status-audit block. This is the pattern to extract into a shared helper; it currently counts per-condition only and does not record which cells, which is why the M4 enrichment cross-check cannot be run against the committed artifact.
- **`conditional_essentiality.py:310`** — `keys = sorted(conditions) if conditions is not None else sorted(CONDITIONS)` is the established condition-threading idiom (backward-compatible default). Steps 2 and 5 replicate it rather than inventing a second convention.
- **`confusion_from_calls` (`:122`) + `mcc` (`:132`)** — condition-free, reused unchanged by the abstention re-derivation.
- **`constant_baselines` (`:340`)** — reused for the abstained-denominator null in Step 8.
- Searched: `dna_decode/fba/`, `scripts/fba_*.py`, `tests/test_fba_*.py`. No `graphify-out/GRAPH_REPORT.md` present; no `src/utils`-class shared directory in this repo (package is flat under `dna_decode/`).

## Pre-Change Baseline

Committed values that every change below must either preserve or explicitly supersede.

**Test suite:** 3,338 passed, 10 skipped, 0 failures.

**Claim A — `wiki/fba_conditional_carbon_2026-08-12.json`** (25 carbon sources, 217 conditionally-essential genes, 5,425 cells):

| metric | value |
|---|---|
| exact-set match | 23/217 = 10.6% |
| per-cell agreement | 0.7368 |
| best constant null | 0.6623 |
| constant-pattern fraction | 184/217 = 84.8% |
| commits (varying predictions) | 33 |
| distinct true patterns | 141 |
| non-optimal solves | 39 across 15 of 25 conditions |

**Claim C — `wiki/fba_regulatory_conditional_test_2026-08-12.json`** (4 media, 67 genes, 268 cells):

| metric | baseline | pFBA-restricted |
|---|---|---|
| exact-set | 3/67 | 5/67 |
| per-cell agreement | 0.5709 | 0.6157 |
| mean per-condition MCC | 0.0611 | 0.217 |
| threshold-free AUROC | 0.6115 | 0.5755 |
| cells called essential | 16 | 96 |
| TP / FP | 10 / 6 | 56 / 40 |
| reactions forced off | 0 | 1,863–1,868 per condition |
| rate-matched null (200 draws) | — | mean 0.5172, sd 0.0284, max 0.5933, 0/200 ≥ observed |

**Claim B — `wiki/fba_gapfill_carbon_recheck_2026-08-12.json`:** 154 binary call flips of 5,425; exact-set −1; per-cell +0.0003.

**Frozen surface:** `amr_rules.py` + `calibrated_amr_rules.json` sha256 must be byte-unchanged (guarded by `tests/test_tb_leak_guard.py`).

## Verification Signal

The plan succeeds when all of these hold:

1. **`pytest tests/ -q` reports ≥ 3,338 passed, 0 failures.** New tests raise the count; no existing test may be deleted or weakened to pass.
2. **Every deletion script audits solver status.** `grep -c 'audit_deletion_frame' scripts/fba_*.py` returns ≥ 1 for all five deletion scripts.
3. **No metric helper silently assumes 4 conditions.** `continuous_readout`, `deployable_threshold`, and `rate_matched_null` each accept a `conditions` argument, and a test passes a 25-key ratios dict to each and asserts a non-degenerate result (the current behavior returns `{"auroc": None, "note": "degenerate: one class only"}`, which the test pins as the *old* bug).
4. **Claim C carries a pre-committed verdict.** `wiki/fba_regulatory_conditional_recheck_<date>.json` exists with a `verdict` field drawn from the frozen three-value set in Step 8, and the memo quotes it.
5. **The false sentence is gone.** `grep -n "t-statistic that the loader reads" wiki/` returns nothing.
6. **The commit claim is stratified.** `wiki/fba_conditional_carbon_2026-08-12.json` carries a `commit_strata` object with `predicted_constant` / `predicted_1_of_25` / `predicted_2plus` counts, each with exact-set, per-cell, and non-optimal overlap.
7. **Frozen surface byte-unchanged** — `tests/test_tb_leak_guard.py` green.

## Implementation Steps

### Step 1: Shared solver-status audit helper
Files: dna_decode/fba/solver_audit.py, tests/test_fba_solver_audit.py
Depends on: none

**What changes:**
- New module `dna_decode/fba/solver_audit.py`, pure and solver-free.
- `DeletionAudit` dataclass: `n_rows`, `n_nonoptimal`, `n_nan_growth`, `statuses: dict[str, int]`, `nonoptimal_cells: set[str]` (gene ids), `nan_cells: set[str]`.
- `audit_deletion_frame(res, condition: str) -> DeletionAudit` — iterates a cobrapy deletion DataFrame, extracts `next(iter(row["ids"]))` as the gene id, reads `row["status"]` defensively (`row.get("status") if hasattr(row, "get") else None`, matching the existing idiom at `fba_conditional_carbon_validate.py:107`), and records **which cells**, not just counts. Recording the cell ids is the whole point — the committed 39-cell count cannot be crossed against the commit set because only counts were kept.
- `merge_audits(dict[str, DeletionAudit]) -> dict` — JSON-serializable rollup for artifact sidecars, with per-condition and total counts plus the `(gene, condition)` cell list.
- Docstring records the verified cobrapy 0.31.1 column contract and the `ids`-is-a-frozenset gotcha.

**Test strategy:**
- Unit tests build a synthetic `pandas.DataFrame` with `columns=["ids","growth","status"]` — no solver, no model, no `feba.db`.
- Cases: all-optimal frame → `n_nonoptimal == 0`; a mixed frame → correct gene ids in `nonoptimal_cells`; a NaN-growth-but-optimal row → counted in `nan_cells` and not in `nonoptimal_cells` (they are distinct failure modes); a frame with no `status` column → returns `status_available=False` rather than raising, so an older cobrapy cannot crash a run.
- One test asserts `merge_audits` output is `json.dumps`-able.

### Step 2: Make continuous_readout and deployable_threshold condition-parametric
Files: dna_decode/fba/conditional_essentiality.py, tests/test_fba_conditional_essentiality.py
Depends on: none

**What changes:**
- `continuous_readout(records, ratios, conditions=None)` — replace `keys = sorted(CONDITIONS)` at `:164` with the established idiom `keys = sorted(conditions) if conditions is not None else sorted(CONDITIONS)`.
- `deployable_threshold(records, ratios, n_folds=5, conditions=None)` — same replacement at `:236`.
- Both keep the 4-media default, so every existing caller is unchanged by construction.
- Add a module-level comment at each site matching the one already at `:307-309`, naming this as the third and fourth instances of the retracted bug class.

**Test strategy:**
- A test builds a 25-key ratios dict and calls each helper **without** `conditions`, asserting the degenerate return (`auroc is None` / `n_cells == 0`) — this pins the *old* silent-failure behavior as a known trap rather than leaving it undocumented.
- A paired test passes the same dict **with** `conditions=<25 keys>` and asserts a non-degenerate result (`n_cells > 0`, `auroc is not None`).
- A regression test asserts the 4-media default path returns values identical to the committed baseline, so the threading cannot perturb any published number.

### Step 3: Make the Fitness Browser t-statistic reachable
Files: dna_decode/fba/fitness_browser.py, tests/test_fba_fitness_browser_t.py
Depends on: none

**What changes:**
- `load_records` and `mean_fitness_matrix` SQL becomes `SELECT g.sysName, f.expName, f.fit, f.t` (column verified present on the live db: `['orgId','locusId','expName','fit','t']`).
- `load_records` gains `min_abs_t: float | None = None`. When set, a per-(gene, condition) measurement is admitted only if `abs(mean_t) >= min_abs_t`; a gene falling below in **any** condition is dropped, preserving the existing complete-row rule at `:104-105`.
- New `mean_t_matrix(conn, conditions, genes)` mirroring `mean_fitness_matrix`, so a sweep can consume t without re-querying.
- Default `min_abs_t=None` — **no existing number moves**; this step only makes the axis reachable.
- Docstring records that the inherited `fit < -2` cutoff is unchanged and that t is now available but unused by default.

**Test strategy:**
- Tests build an in-memory SQLite with the real `GeneFitness` / `Gene` / `Experiment` schema (5-column `GeneFitness`), so no 7.4 GB db is needed.
- Cases: `min_abs_t=None` returns exactly the same records as before (equivalence pin); a gene with a low-|t| measurement in one condition is dropped when `min_abs_t=3.0`; `mean_t_matrix` averages replicates the same way `mean_fitness_matrix` does.
- A test asserts the SQL selects four columns, guarding against a future edit silently dropping `f.t` again.

### Step 4: Correct the false honest-limits sentence
Files: wiki/fba_conditional_carbon_2026-08-12.md
Depends on: none

**What changes:**
- Replace "the t-statistic that the loader reads is not used" with the true statement: the per-measurement t-statistic (`GeneFitness.t`) is present in the source table and was never read by the loader; as of this plan it is selectable via `min_abs_t` but unused by default.
- Add a one-line correction note in the same bullet recording that this was a factual error inside the honest-limits section, found by adversarial review — consistent with how the 84.8% correction is recorded in the same file.

**Test strategy:**
- Manual: `grep -n "t-statistic that the loader reads" wiki/` returns nothing.
- No unit test — this is a prose correction. The behavioral guarantee is Step 3's four-column SQL test.

### Step 5: Audit and abstain in the pFBA regulatory script
Files: scripts/fba_regulatory_conditional_test.py, tests/test_fba_regulatory_conditional.py
Depends on: Step 1, Step 2

**What changes:**
- `score_model` calls `audit_deletion_frame(res, c)` per condition and returns the merged audit alongside its metrics.
- Record `pfba_status` per condition from the `pfba(model)` solution object at `:89`, so an inspected-but-clean pFBA solve is provable rather than assumed.
- **Abstention arm:** add `abstain_nonoptimal: bool` to `score_model`. When True, a `(gene, condition)` cell whose deletion row was non-optimal or NaN-growth is **excluded from scoring entirely** rather than coded as ratio `0.0` (= essential). The current `:106` coding is retained as the default so the baseline number stays reproducible; abstention is the new arm.
- Recompute `constant_baselines` and `rate_matched_null` **on the abstained cell set** — a null computed on the full denominator against a metric computed on a reduced one is not a control. This is the single most error-prone part of the step.
- Thread `conditions` through `rate_matched_null(records, n_called_essential, conditions=None, ...)`: replace the hardcoded `keys = sorted(CONDITIONS)` at `:59` and pass `conditions=keys` into `switch_accuracy` at `:68`, closing both the explicit and the fall-through hardcode.
- `continuous_readout` call at `:112` passes `conditions=keys`.

**Test strategy:**
- Existing 4 `rate_matched_null` tests must pass unchanged (the `conditions` default preserves behavior).
- New: `rate_matched_null` with an explicit 25-key condition set produces cells over 25 keys, not 4 (pins the closed hardcode).
- New: an abstention unit test with a synthetic records/calls/audit triple asserts that excluded cells appear in neither the numerator nor the denominator of `per_condition_agreement`, and that the constant null is recomputed on the reduced set.
- The live solver run is Step 8; this step is offline-testable.

### Step 6: Audit solver status in the two gap-fill scripts
Files: scripts/fba_gapfill_carbon_recheck.py, scripts/fba_gapfill_conditional_test.py, tests/test_fba_solver_audit.py
Depends on: Step 1

**What changes:**
- Both `score` functions call `audit_deletion_frame` and persist the merged audit into their result JSON under `solver_audit`.
- `fba_gapfill_carbon_recheck.py:71`'s `(g != g) or (g < FRAC * wt)` keeps its NaN-to-essential coding (changing it would silently move the published 154-flip number); the audit makes the NaN count **visible** so the flip count can be read against it.
- Add a caveat line to each script's result `caveats` list stating how many flips involve a non-optimal cell — the flip count is only meaningful net of solver noise.

**Test strategy:**
- A **static regression guard** in `tests/test_fba_solver_audit.py`, mirroring the repo's existing "no unmigrated consumer" guard pattern: parse all five `scripts/fba_*` deletion scripts with `ast` and assert each one references `audit_deletion_frame`. This is the test that stops a sixth deletion script from shipping unaudited — a coverage rule the wirings themselves cannot enforce.
- A second guard asserts every script that reads `row["growth"]` also reads a status, so the NaN-to-essential coding can never again appear without an audit beside it.
- No behavioral unit tests for the wirings themselves (they are thin calls into a helper tested in Step 1); the behavioral check is the re-run in Step 9.

### Step 7: Per-cell non-optimal recording and commit stratification in the carbon validator
Files: scripts/fba_conditional_carbon_validate.py
Depends on: Step 1, Step 2

**What changes:**
- Replace the count-only audit block at `:105-110` with `audit_deletion_frame`, so the 39 non-optimal cells are recorded **by (gene, condition)** and can be crossed against the commit set.
- New `commit_strata` section in the result JSON, partitioning the 217 genes into `predicted_constant` / `predicted_1_of_25` / `predicted_2plus`, each carrying `n_genes`, `n_exact_set_match`, `per_cell_agreement`, and `n_genes_touching_a_nonoptimal_cell`.
- New `nonoptimal_enrichment` field: how many of the 33 commit genes and how many of the 23 exact matches involve at least one non-optimal cell. This is the M4 cross-check — one non-optimal cell can create or destroy a 1-of-25 exact-set match, so a global 0.7% rate does not clear the concentrated subset.
- `continuous_readout` / `deployable_threshold` calls (if present) pass `conditions=keys`.

**Test strategy:**
- A pure unit test for the stratification function with synthetic records + calls, asserting the three strata partition the gene set exactly (counts sum to `n_conditionally_essential`, no gene in two strata).
- A test asserting exact-set matches can only occur in the two non-constant strata (a constant prediction cannot match a two-sided gene) — this is the arithmetic invariant that the retracted bug violated.

### Step 8: Re-derive Claim C under abstention and issue the pre-committed verdict
Files: wiki/fba_regulatory_conditional_recheck_2026-08-13.json, wiki/fba_regulatory_conditional_recheck_2026-08-13.md
Depends on: Step 5

**What changes:**
- Run `uv run python scripts/fba_regulatory_conditional_test.py` with both arms (default coding + abstention) and persist a new artifact. The original `wiki/fba_regulatory_conditional_test_2026-08-12.json` is **not overwritten** — a superseding artifact is written alongside it, matching the project's shared-key-overwrite discipline.
- **Verdict is pre-committed before the run** (the project's own verdict-vs-budget lesson: author the branches before the result lands). Let `A` = abstained per-cell agreement, `N_max` = max of the rate-matched null recomputed on the abstained denominator, `f` = fraction of restricted-arm cells that are non-optimal or NaN:
  - `REGULATORY_LIFT_CONFIRMED` — `f < 0.05` **and** `A > N_max` **and** `A` within 0.02 of the published 0.6157.
  - `REGULATORY_LIFT_IS_A_SOLVER_ARTIFACT` — `A <= N_max`, i.e. the lift does not survive removing non-optimal cells.
  - `REGULATORY_LIFT_PARTIALLY_SURVIVES` — anything else; the memo then reports the abstained number as the headline and the original as superseded.
- Sidecar records `f`, the per-condition non-optimal counts, the pFBA solution statuses, and the recomputed null.

**Test strategy:**
- Requires a solver; `feba.db` is **not** needed (this is the 4-media Orth substrate, labels come from the committed TSV).
- Post-run assertion: the default-coding arm reproduces the committed baseline (per-cell 0.6157, TP 56, FP 40) exactly. If it does not, the wiring changed behavior and the run is invalid — stop and fix before interpreting the abstention arm.

### Step 9: Re-run the carbon panel with per-cell audit and stratification
Files: wiki/fba_conditional_carbon_2026-08-12.json, wiki/fba_gapfill_carbon_recheck_2026-08-12.json
Depends on: Step 7, Step 6

**What changes:**
- Re-run `scripts/fba_conditional_carbon_validate.py` (5,425 knockouts) and `scripts/fba_gapfill_carbon_recheck.py`, refreshing both sidecars in place with the new `solver_audit`, `commit_strata`, and `nonoptimal_enrichment` fields.
- In-place refresh is correct here (unlike Step 8) because no headline number is expected to change — only new audit fields are added. **If any committed metric moves, that is a finding, not a refresh:** stop, and treat the delta as a defect in Step 7's wiring.

**Test strategy:**
- Requires `feba.db` on `D:` plus a solver.
- Post-run assertion: exact-set 23/217, per-cell 0.7368, constant 184/217 are byte-identical to the committed values. Any drift invalidates the run.

### Step 10: Restratify the memos and record both verdicts
Files: wiki/fba_conditional_carbon_2026-08-12.md, wiki/fba_gapfill_conditional_answer_2026-08-12.md, wiki/fba_regulatory_conditional_test_2026-08-12.json
Depends on: Step 4, Step 8, Step 9

**What changes:**
- Carbon memo: replace the "accurate when it commits" H1 and summary with the stratified table from Step 7 (predicted-constant / 1-of-25 / 2+, each with count, exact-set, per-cell, non-optimal overlap). The commit-rate decomposition stays, subordinated to the table. **"Commit" is NOT redefined to mean multi-condition only** — that would rest the headline on 8 genes and swap one over-claim for an under-powered statistic.
- Carbon memo honest-limits: replace the "39 non-optimal, same gap in three other scripts" bullet with the resolved state plus the enrichment cross-check result.
- `fba_regulatory_conditional_test_2026-08-12.json`: add a `superseded_by` field pointing at the Step 8 artifact and the verdict. The original numbers stay byte-unchanged — superseded, not rewritten.
- Gap-fill memo: add the solver-audit numbers to the 154-flip claim.

**Test strategy:**
- Manual read-through against the Verification Signal list.
- `grep` guards from Verification Signal items 2, 5.

## Execution Preview

| wave | steps | parallelism |
|---|---|---|
| 0 | Step 1, Step 2, Step 3, Step 4 | 4 |
| 1 | Step 5, Step 6, Step 7 | 3 |
| 2 | Step 8, Step 9 | 2 |
| 3 | Step 10 | 1 |

- **Total waves:** 4
- **Max parallelism:** 4 (wave 0)
- **Critical path:** Step 1 → Step 5 → Step 8 → Step 10 (4 waves)
- **No intra-wave file overlap.** Wave 0's four steps touch four disjoint files. Wave 1's three steps touch three disjoint scripts (Step 5 → regulatory, Step 6 → two gap-fill, Step 7 → carbon validator). Wave 2's two runs write disjoint artifacts. Step 10 depends on Step 4 solely because both edit `wiki/fba_conditional_carbon_2026-08-12.md`.

## Risk Flags

- **Step 8 is the load-bearing step and can invalidate a published claim.** That is the intent. The pre-committed verdict set exists specifically so the interpretation is not authored after seeing the number.
- **Abstention changes the denominator — the most likely place to introduce a new error.** If `constant_baselines` or `rate_matched_null` are computed on the full cell set while the metric is computed on the abstained set, the comparison is invalid and would likely *flatter* the result. Step 5's test explicitly pins the recomputation.
- **Step 9 is long-running** (5,425 knockouts, plus a second run over the 1,125-reaction augmented model) and requires the external `D:` drive. Prior sessions recorded USB disconnects corrupting long runs on this host. Run detached with flushed output; if `D:` drops, the run is void, not partial.
- **Steps 8 and 9 need a working solver and, for Step 9, `feba.db`.** Every other step (1–7, 10) is fully offline. If the solver environment is unavailable, waves 0–1 still complete and deliver the correctness fixes; only the re-derivations block.
- **`[unverified]` — the actual non-optimal rate in the restricted pFBA arm is unknown.** The whole C1 concern is a mechanism argument from the shape of the result (recall-only gain, precision and AUROC both degrading) plus the 69% forced-off rate. It is plausible the rate is near zero and Claim C is confirmed intact. The plan is written to make either outcome cheap to report.
- **Do not "fix" the NaN-to-essential coding as a cleanup.** It is load-bearing for reproducing every committed number. It is quarantined behind the new abstention flag deliberately.
- No `project-rules.md`, `DESIGN.md`, or `AUDIT_REPORT.md` in this repo — those sections are omitted rather than fabricated. No UI surface is touched, so `/visual-design-consultation` does not apply.

## Open Questions

1. **Should the abstention arm become the default once Step 8 lands?** If `f` turns out to be small and the verdict is CONFIRMED, keeping the default coding is harmless and preserves reproducibility. If `f` is large, the default coding is actively misleading and arguably should flip — but flipping it changes every previously published number from these scripts. Recommendation: leave the default alone regardless of verdict, and let the memo carry the abstained number as the headline. Decide at Step 10.
2. **Does the label-threshold sweep belong in this plan?** Step 3 makes `GeneFitness.t` reachable but runs no sweep. The sweep can move the 217-gene denominator and therefore every number in Claim A, which makes it a larger change than everything here combined. Scoped out deliberately; it should be its own plan once the metrics it would consume are trustworthy.
3. **Is `min_abs_t` the right filter shape?** Dropping a gene that fails |t| in any single condition preserves the complete-row rule but could shrink the 217-gene set sharply. A per-cell mask (admit the gene, abstain the cell) is the alternative and composes better with Step 5's abstention machinery — but it changes the switch-pattern semantics. Not decidable without seeing the t distribution; Step 3 ships the simple version and `mean_t_matrix` makes the distribution inspectable.

## Verification

```bash
cd /c/Users/Farshad/PythonProjects/dna_decode

# 1. Full suite -- must be >= 3,338 passed, 0 failures
uv run pytest tests/ -q

# 2. Focused FBA suite
uv run pytest tests/test_fba_solver_audit.py tests/test_fba_conditional_essentiality.py \
              tests/test_fba_regulatory_conditional.py tests/test_fba_fitness_browser_t.py -v

# 3. Every deletion script audits status (expect 5 files, each >= 1)
grep -c 'audit_deletion_frame' scripts/fba_conditional_carbon_validate.py \
    scripts/fba_conditional_essentiality_validate.py scripts/fba_regulatory_conditional_test.py \
    scripts/fba_gapfill_carbon_recheck.py scripts/fba_gapfill_conditional_test.py

# 4. No metric helper silently assumes 4 conditions (expect only the CONDITIONS
#    definition + apply_condition lookup + backward-compatible defaults)
grep -n 'sorted(CONDITIONS)\|tuple(CONDITIONS)' dna_decode/fba/conditional_essentiality.py

# 5. The false sentence is gone
grep -rn "t-statistic that the loader reads" wiki/ || echo "OK: corrected"

# 6. Frozen AMR surface byte-unchanged
uv run pytest tests/test_tb_leak_guard.py -q

# 7. Re-derivations (need solver; step 9 also needs D:/feba.db)
uv run python scripts/fba_regulatory_conditional_test.py
uv run python scripts/fba_conditional_carbon_validate.py
uv run python scripts/fba_gapfill_carbon_recheck.py

# 8. Verdict landed
python -c "import json; print(json.load(open('wiki/fba_regulatory_conditional_recheck_2026-08-13.json'))['verdict'])"
```

<!-- toolkit: check=clean waves=clean gate=fired:open-questions,unverified -->
