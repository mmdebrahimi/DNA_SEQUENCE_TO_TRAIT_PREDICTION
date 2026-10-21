# Close the Pre-Grid Defects in the Sampling-Sweep Instrumentation — Technical Plan

> Four defects in the guards that gate the temperature grid, all found by review before the grid ran; two would have silently voided or gutted it.

## Lens status
Inputs: conversation (review findings against the in-flight `Discharge_The_Queued_Sampling_Sweep_G_D_Unblock_Trigger_Plan/`, Steps 3–6)
Degradations: repo-index absent (`scripts/repo-index-status.sh` not installed) — direct file reading only; sentrux unavailable; external second-opinion surfaces unavailable (review was single-pass); `python3` is a non-functional WindowsApps stub, toolkit run under `python`

## Problem Statement

The live sampling-sweep plan is mid-execution: Steps 1–2 shipped, the Step 3 reproduction anchor is in flight detached, and Steps 4–6 (the two downward grid points, the result artifact, the ledger discharge) have not run. A review of the instrumentation built for Steps 1–3 found four defects, two of them capable of destroying the grid's value without raising an error.

- **The geometry guard false-refuses a legitimate run.** The artifact records `positional_bins` as `a.bins` unconditionally (`scripts/glm_generator_falsifier.py:258`); the `--compare-to` projection records `a.bins if "positional" in a.modes or "both" in a.modes else None` (`:129`). Measured: a `--modes kmer` run compared against its own kmer-only baseline yields `{'positional_bins': (10, None)}` and is refused. It cannot bite the current `kmer,both` runs — both branches take `a.bins` — which is exactly why it survived the author's own check.
- **The grid's exclusion rule, as written in Step 4, would exclude the T=1.0 anchor** — the comparator the other two points are measured against. The plan's own Verification Signal gate 2 is correctly scoped (`knob_verified_end_to_end == true` at every **non-default** grid point, `:110`), but Step 4's bullet (`:194`) and test-strategy line (`:199`) drop the scoping. Verification is gated on `requested_non_default` (`:170`), which is `False` at T=1.0, and the real in-flight artifact confirms the consequence: `knob_contract = {requested_non_default: false, knob_wired: false, knob_verified_end_to_end: false}`. A default point is not *unverified*; it is *nothing-to-verify*.
- **The knob verification exercises a path the sweep does not use.** `_at` calls `gen.generate(...)` (`:178`) — unconditional — while the sweep runs `generate_conditional`. Both pass `logits_processor=self._sampling_processors()` (`generate.py:232`, `:277`), so the risk is low; but in an episode whose entire subject is a knob that looked wired and was inert, proving liveness on a different path is a proxy.
- **The detached runner has no lock against a double launch.** Two concurrent launches means two torch processes, this host's documented access-violation hazard, and a job's liveness was already misjudged once during this execution.

Scope is the instrumentation only. This plan does **not** run the grid (that is the live plan's Step 4), does not touch the pre-registered comparator, and does not widen the grid.

## Codebase Context

- `scripts/glm_generator_falsifier.py` — the runner. Two independent dicts describe the same run's geometry: the `--compare-to` projection (`:126-135`) and the artifact (`:255-277`). They were hand-written in parallel and have already diverged on one field. `knob_record` is built at `:168-170`, updated at `:181`, and stamped into the artifact at `:268`. Prompt construction for the conditional path happens at `:195+`, *after* the knob-verification block at `:170-188`.
- `dna_decode/glm/knob_contract.py` — pre-run refusals: `validate_sampling` / `assert_knob_wired` / `verify_knob_live` / `validate_comparison_mode` / `TOP_K_KMER_FLOOR`, plus the geometry contract (`COMPARISON_GEOMETRY_FIELDS`, `geometry_of`, `geometry_diff`, `assert_geometry_comparable`). Every function is pure. The module's stated job is *preconditions on a run*, which is why post-run grid gates do not belong in it.
- `dna_decode/glm/generate.py` — `torch.manual_seed(self.seed)` heads both `generate_conditional` (`:222`) and `generate` (`:266`), so verification generations cannot perturb the main run's RNG stream. `generate_conditional(prompts, length, *, ...)` takes a prompt list.
- `run_glm_anchor_detached.bat` — repo-root detached launcher, modelled on the existing `run_mechanism_audit_detached.bat`. Deletes its sentinel at start; writes `EXIT=<code>` as the last line.
- Test infrastructure: 24 `tests/test_glm_*.py` files. The relevant four are `test_glm_knob_contract.py` (26 cases), `test_glm_sampling_knob.py` (6), `test_glm_mode_contract.py` (15), and `test_glm_geometry_contract.py` (22 cases **written but never executed** — a detached torch job holds a model and the full suite against that is the documented access-violation hazard). Established patterns in these files: a stub-reference + tripwire pair for ordering proofs, and an explicit non-vacuity partner for every guard.

### Reusable-Code Survey
- `dna_decode/glm/knob_contract.py` — already the home of the pure pre-run contract functions and the `geometry_*` helpers Step 1 must single-source against; `COMPARISON_GEOMETRY_FIELDS` is the field list Step 1 reuses rather than redeclaring.
- `dna_decode/glm/generate.py` — `generate_conditional` is the existing primitive Step 3 points the verification at; no new generation code is needed.
- `tests/test_glm_mode_contract.py` — the stub-reference + `_CorpusReached` tripwire pattern (and its non-vacuity partner) is copied rather than reinvented for any new ordering proof.
- `run_mechanism_audit_detached.bat` — the pre-existing detached-runner precedent the new launcher follows.
- None — searched: `graphify-out/GRAPH_REPORT.md` (absent), `src/`, `lib/`, `utils/`, `common/` (absent in this repo; package root is `dna_decode/`) for a generic locking/mutex helper; nothing exists, so Step 5 writes its own.

## Pre-Change Baseline

Captured 2026-10-20/21 on a clean tree at `53b8c1f`:

- **Suite:** 15 failed / 6726 passed / 185 skipped. Every failure is accounted for and none is a code defect (15 absent-data + 4 uninstalled heavy extras, per `CLAUDE.md`).
- **The false refusal, measured:** `geometry_diff` on a `--modes kmer` run against its own kmer-only baseline returns `{'positional_bins': (10, None)}`.
- **The guard's correct behaviour, measured:** `{'positional_bins': (6, 10)}` on the real `wiki/glm_sweep_anchor_bins10_UNCOMPARABLE_2026-10-20.json` vs `wiki/glm_generator_falsifier_REAL_2026-10-07.json`; empty on self-comparison; empty on the planned `--bins 6` re-run.
- **The anchor's knob record, from the real artifact:** all three of `requested_non_default` / `knob_wired` / `knob_verified_end_to_end` are `false` at T=1.0.
- **Unrun coverage:** `tests/test_glm_geometry_contract.py` — 22 cases, 0 executed.
- **Recorded comparator (unchanged by this plan):** `kmer` GLM 0.7237 vs null 0.5990; `both` GLM 0.7706 vs null 0.6679, at `positional_bins=6`.

## Verification Signal

1. `geometry_diff` on a `--modes kmer` run against its own kmer-only baseline returns `{}` — the measured `{'positional_bins': (10, None)}` is gone.
2. A test asserts the `--compare-to` projection and the artifact's geometry fields come from **one** expression, and fails if a second expression is reintroduced.
3. `verdict_eligibility` classifies the real anchor record (all-false at T=1.0) as eligible-by-default-sampling, and a requested-but-unverified record as excluded. Both directions pinned.
4. The live plan's Step 4 bullet and test-strategy line cite the same scoping as its own gate 2, and a test pins that the plan text and the function agree.
5. The knob verification runs `generate_conditional`, demonstrated by a test that fails if it is pointed back at the unconditional path.
6. A second concurrent launch of the detached runner exits non-zero without starting a second torch process.
7. Suite after: 15 failed / ≥6748 passed (6726 + 22 newly-executed geometry cases, plus the cases added here), zero committed-artifact churn.

## Implementation Steps

### Step 1: Single-source the run geometry so the two descriptions cannot diverge
Files: scripts/glm_generator_falsifier.py, tests/test_glm_geometry_contract.py
Depends on: none

**What changes:**
- Add `geometry_from_args(a) -> dict` to the runner, returning exactly the `COMPARISON_GEOMETRY_FIELDS` keys, with `positional_bins` taken as `a.bins` — matching what the artifact already records.
- Replace the hand-written `--compare-to` projection (`:126-135`) with a call to it.
- Spread the same call into the artifact dict (`:255-277`), removing the duplicate literals for those nine fields and leaving the artifact's other keys (`natural_sets`, `fit_eval_split`, `device`, `sampling`, `knob_contract`, `results`, `honest_limits`) untouched.
- **The fix is structural, not a value edit.** Changing `None` to `a.bins` in one place would fix today's divergence and leave the mechanism that produced it — two parallel hand-written views of the same run — fully intact. One expression is what makes a recurrence impossible.
- Keep the `context_bp` conditional (`a.context if a.conditional else None`) exactly as the artifact writes it; it is consistent across both sites today and must stay that way.

**Test strategy:**
- Unit: the measured regression — a `--modes kmer` run against its own kmer-only baseline diffs empty.
- Unit: `set(geometry_from_args(ns)) == set(COMPARISON_GEOMETRY_FIELDS)`, so a field added to the tuple without a value fails loudly rather than silently reporting `<absent>`.
- Unit (non-vacuity, the load-bearing one): the artifact's geometry fields and the `--compare-to` projection are asserted **equal for the same args**, across several arg shapes including `kmer`-only and `positional`-only. This is the test that would have caught the original divergence, so it must fail if the two sites are un-merged.
- Regression: the existing `{'positional_bins': (6, 10)}` real-artifact case still refuses.

### Step 2: Code-own the grid's validity gates in a post-run module
Files: dna_decode/glm/sweep_gates.py, tests/test_glm_sweep_gates.py
Depends on: none

**What changes:**
- New module `dna_decode/glm/sweep_gates.py` for **post-run** grid gates, kept separate from `knob_contract.py` because that module's job is preconditions on a run and mixing the two makes neither name true.
- `verdict_eligibility(knob_record) -> str`, four outcomes, first match wins: `eligible:default-sampling` (no non-default knob was requested — nothing to verify) · `eligible:verified` · `excluded:requested-but-unverified` · `excluded:knob-dead` (`requested_non_default` true, `knob_wired` false). Returning a reason string rather than a bool is deliberate: the bool is what collapses *nothing-to-verify* into *unverified*, which is the defect.
- `null_reproduced(point_null, anchor_null, recorded_null, *, tol=0.0) -> dict` — the cross-point invariant. Temperature touches only the GLM candidate, so every grid point's Markov null must equal the anchor's. **Honest correction to the review that prompted this:** the live plan already specifies a null cross-check, against the *recorded* `0.5990 / 0.6679` (`:102-103`); what is missing is that it is prose rather than a coded gate, and that it compares against the recorded baseline only, never anchor-vs-grid-point. This function covers both, and `tol=0.0` is the default because the null is empirically deterministic — it reproduced to four decimals across two independent runs.
- `assemble_grid(points) -> dict` — applies `verdict_eligibility` and `null_reproduced` per point and returns the eligible set plus the excluded set **with reasons**, never dropping a point silently.

**Test strategy:**
- Unit: `verdict_eligibility` on the **real** anchor record (all three fields false) returns `eligible:default-sampling`; a test comment records that a bool-returning version would have excluded the comparator.
- Unit: all four outcomes reachable, and the ladder's order pinned (a dead knob that also requested non-default classifies `excluded:knob-dead`, not `excluded:requested-but-unverified`).
- Unit: `null_reproduced` fails on a moved null and passes on an identical one; a moved null is reported as a **validity failure of that point**, not averaged in.
- Unit: `assemble_grid` never returns an excluded point without a reason string.

### Step 3: Verify the sampling knob on the conditional path the sweep actually uses
Files: scripts/glm_generator_falsifier.py
Depends on: Step 1

**What changes:**
- Hoist the conditional prompt construction (currently `:195+`) above the knob-verification block so one real prompt is available to it.
- Point `_at` at `gen.generate_conditional(prompts[:1], a.length, ...)`, falling back to `gen.generate` only if `a.conditional` is false — unreachable today, since `validate_comparison_mode` leaves `promoter + conditional` as the single valid pairing, but kept so the helper is not silently wrong if that contract ever widens.
- Wrap the `gen.temperature` mutation in `try/finally` so an exception inside verification cannot leave the generator at the probe's temperature (`:177-181`).
- Depends on Step 1 only because both edit this file; there is no logical coupling.

**Test strategy:**
- Unit: a stubbed generator records which method the verification called; the test fails if it is `generate` rather than `generate_conditional`.
- Unit (non-vacuity): the same harness proves the assertion can fail — pointing `_at` back at the unconditional path makes it red.
- Unit: an exception raised inside the probe leaves `gen.temperature` at its original value.
- Note for the executor: this reorders real work in `main()`, so confirm the genome is loaded once, not twice, after the hoist.

### Step 4: Align the live plan's Step-4 prose with its own correctly-scoped gate
Files: plans/Discharge_The_Queued_Sampling_Sweep_G_D_Unblock_Trigger_Plan/technical-plan.md, tests/test_glm_sweep_gates.py
Depends on: Step 2

**What changes:**
- Amend the Step 4 bullet (`:194`) and its test-strategy line (`:199`) to scope exclusion to points that **requested** a non-default knob, matching gate 2 (`:110`) which was already right. Record it as a save-time amendment rather than a silent edit, since the plan is executing.
- Cite `sweep_gates.verdict_eligibility` as the single definition, so the rule lives in one place instead of three.
- **Deliberately not changed:** the comparator, the outcome branches, the downward-only grid. Those are pre-registered, and editing them mid-execution after seeing the anchor is the thing pre-registration exists to prevent.

**Test strategy:**
- Unit: the plan's Step-4 text contains the non-default scoping and does **not** contain the unscoped wording — the same shape as the existing prose-pinning test that keeps the report card from naming a retired lock.
- Unit (non-vacuity): restoring the unscoped sentence makes it fail.

### Step 5: Make the detached runner refuse a second concurrent launch
Files: run_glm_anchor_detached.bat, scripts/glm_run_lock.py, tests/test_glm_run_lock.py
Depends on: none

**What changes:**
- `scripts/glm_run_lock.py` — acquire/release a PID-bearing lockfile under the scratch/run directory; `acquire` exits non-zero when a live holder exists, and **reclaims a stale lock whose PID is gone**, so a killed job cannot block the next launch forever (the failure mode that already cost this execution two restarts).
- The `.bat` calls `acquire` before `uv run` and releases in all paths; a failed acquire exits without starting python.
- Liveness is checked by PID existence, not by the sentinel, because the sentinel is written only at exit and so says nothing while a job is running.

**Test strategy:**
- Unit: second `acquire` against a live PID fails; against a dead PID it succeeds and reclaims.
- Unit: release is idempotent and a release by a non-owner does not remove the owner's lock.
- Manual: two launches in quick succession start exactly one `glm_generator_falsifier` process — verified by process listing, which is the only check that covers the `.bat` wiring rather than the Python helper.

### Step 6: Run the deferred suites and reconcile against the baseline
Files: tests/test_glm_geometry_contract.py, tests/test_glm_sweep_gates.py, tests/test_glm_run_lock.py
Depends on: Step 1, Step 2, Step 3, Step 5

**What changes:**
- Run the four knob/geometry/sweep/lock files first, then the full suite, and reconcile the delta against 15 failed / 6726 passed / 185 skipped — naming every moved count. An unexplained delta is a finding, not noise.
- Confirm zero churn in committed artifacts after a full run (the autouse restore fixtures were added for exactly this).
- Confirm the frozen five are byte-unchanged and the prospective lock still verifies at 2026-08-31.
- **Blocked while a torch job holds a model.** This step additionally requires no detached generation in flight — running the suite against one is this host's documented access-violation hazard, which is why the 22 geometry cases are unrun today. Check for a live `glm_generator_falsifier` process before starting; this is an environmental precondition the dependency graph cannot express.

**Test strategy:**
- This step *is* the test execution; its own signal is the reconciled count plus a clean `git status`.
- Any newly-red test is diagnosed before the step closes — a red suite that is merely explained is not a passing suite.

## Execution Preview

Wave 0 (3 parallel):  Step 1 — single-source the run geometry, Step 2 — post-run grid gates module, Step 5 — detached-runner lock
Wave 1 (2 parallel):  Step 3 — verify on the conditional path, Step 4 — align the live plan's prose
Wave 2 (1):           Step 6 — run the deferred suites and reconcile

Critical path: Step 1 → Step 3 → Step 6 (3 waves)
Max parallelism: 3 agents

Note: Parallel execution requires a git repository with a configured remote. If unavailable, /execute-plan falls back to sequential mode.

## Risk Flags

- Severity: medium — File overlap: Steps 1 and 3 both edit `scripts/glm_generator_falsifier.py`. Resolved by an explicit `Depends on: Step 1` on Step 3, which serialises them across waves; there is no logical coupling, only a shared file.
- Severity: medium — Step 3 reorders real work inside `main()` (hoisting prompt construction above the verification block). Transitive risk not captured by the `Files:` list: it moves a `load_fasta` call, so a careless hoist loads the 4.6 MB genome twice. Named in the step.
- Severity: medium — Step 6 has an environmental precondition the wave graph cannot express: no torch job in flight. If the executor ignores it, the symptom is an access violation that reads like a test failure.
- Severity: low — `tests/test_glm_sweep_gates.py` is written by Step 2 and extended by Step 4, so the two steps share a file across waves. Serialised by `Depends on: Step 2`.
- Severity: low — Step 4 edits a plan that is mid-execution. Mitigated by recording it as an amendment and by scoping the edit to the two lines that contradict the plan's own gate; the pre-registered comparator and grid are untouched.
- Severity: low — External tool integrations: `uv run` (verified by the two successful detached launches this session); `pytest` (verified in use); the `plan-toolkit.py` surface is verified by `--help`. `python3` is a **non-functional** WindowsApps stub on this host — `python3 --version` prints a Store-install message and exits 0, so any script invoking `python3` fails open; the toolkit was run under `python`. Windows `.bat` PID-liveness checking in Step 5 is `[unverified]` against an authoritative source and must be confirmed empirically in that step, not assumed.
- Severity: low — Restructuring applied: the review's findings were split so that the two file-sharing runner edits land in different waves, and the plan-prose fix was separated from the module that defines the rule it cites, so the citation cannot be written before the function exists.

## Open Questions

- **Spend ~22 min on a seed replicate at T=1.0 to get a generation-variance floor?** Every grid point uses seed 7, so the only noise figure available is the discriminator's natural-vs-natural self-control (max 0.5447 over five splits), which is not generation variance. Without a replicate, a small gap between temperatures is uninterpretable. Against: the pre-registered decision rule is a verdict plus the sign of a gap, not a tolerance, so a small gap was never going to be claimed, and the arm already costs 3 × 22 min.
- **Add T=1.2 to the grid?** The downward-only restriction rests on the uniform control (0.9408) being far worse than the GLM (0.7237), which rules out the *limit* of flattening, not mild flattening. If both sharpened points come back worse, "sharpening does not help" and "T=1.0 is a local optimum" will be indistinguishable. Against: widening a pre-registered grid after seeing the anchor is what pre-registration exists to prevent; cleaner to record it as an untested direction in the result memo.
- Should `null_reproduced` default to `tol=0.0`? It is justified — the null reproduced to four decimals across two independent runs — but it makes any floating-point drift a hard validity failure rather than a warning.

## Verification

1. `uv run pytest tests/test_glm_geometry_contract.py tests/test_glm_sweep_gates.py tests/test_glm_knob_contract.py tests/test_glm_run_lock.py -q` — all green, with the 22 previously-unrun geometry cases now executed.
2. Re-run the measured baseline by hand: a `--modes kmer` run against its own kmer-only baseline diffs empty, and the real `{'positional_bins': (6, 10)}` case still refuses.
3. `verdict_eligibility` on the in-flight anchor artifact's `knob_contract` block returns `eligible:default-sampling`.
4. Full suite reconciles against 15 failed / 6726 passed / 185 skipped with every moved count named, and `git status` is clean afterwards.
5. Frozen five byte-unchanged; prospective lock verifies at 2026-08-31.
6. Two rapid launches of the detached runner yield exactly one `glm_generator_falsifier` process.


## Save-time amendments

Captured at: 2026-10-21
Source: `/save-plan` arguments

These are **audit notes for human readers, not executable instructions**. `/execute-plan`
reads only `## Implementation Steps` for work. If an amendment ever changes a step's file
list, dependencies or structure, re-run `/technical-plan` before executing.

- corrected the review's over-claim that the null cross-check was unspecified (the live plan already specifies it against the recorded values; what was missing was a coded gate and the anchor-vs-point comparison)
- the positional_bins divergence was measured rather than inferred before planning
- the fix was made structural (one shared expression) rather than a one-line value edit
- grid widening and the seed replicate left as open questions rather than folded into scope

<!-- toolkit: check=clean waves=clean gate=fired:open-questions,unverified -->
