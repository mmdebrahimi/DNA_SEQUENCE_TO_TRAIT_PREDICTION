# GLM Ceiling-Arm Accounting Correction — Technical Plan

> Correct two self-contradictions in the pyramid-arm memo, close the test gap the new constant opened, record the substrate-infeasibility as a revivable scoped negative, and stage the one arms-ledger move that needs user authority.

## Lens status
Inputs: conversation (the 2026-10-12 adversarial review of the 2026-10-11 `--advance` conclusions)
Degradations: repo-index unconfigured (no `.repo-index/` — Step 1 was direct file reading, no Layer 2 used); sentrux unavailable (not on PATH — no architecture gate); `project-rules.md` absent (Project Rules Check section omitted); `DESIGN.md` absent but plan touches no UI, so no risk flag needed

## Problem Statement

The 2026-10-11 run measured a new arm against the GLM ceiling criterion and drew six conclusions. An adversarial review verified three defects, all in the reporting rather than the measurement:

1. `wiki/glm_pyramid_arm_2026-10-11.md:5` claims the run *"closes one of them with a number"*, while `:70` says *"the arm measured here is not either declared arm"* and the criterion stays `blocked:untried-arm`. Both cannot hold — a NEW arm was measured and NEITHER declared arm was closed.
2. `:31` claims *"more scale is not the missing lever"* — a lever-level negative from one encoder at one capacity on one substrate, which contradicts the same artifact's own `honest_limits` block.
3. `multiresolution-pyramid-4-decades-UNRUN` sits at `status: planned` in `soraya_runs/arms-ledger.json` while `dna_decode/glm/encoder.py:148 validate_widths` provably refuses it on this substrate — so the criterion is held in a state the repo already knows is impossible.

Step 1 research found a fourth, unreported issue: `PYRAMID_WIDTHS` (added 2026-10-11) has **zero** test references, so neither the constant's span nor the `validate_widths` refusal that makes the 4-decade arm infeasible is pinned by a test.

Defects 1, 2 and 4 are fully decided and reversible. Defect 3 is **authority-gated**: resolving it means either killing an expensive arm or splitting an acceptance bar, and this repo treats both as the user's call.

## Codebase Context

- `dna_decode/glm/encoder.py` — `DEFAULT_WIDTHS = (6, 20, 75)` (line 46), `PYRAMID_WIDTHS = (3, 12, 48, 150)` (line 62), `validate_widths(widths, length)` (line 148) computes `too_wide = [x for x in w if conv_output_length(length, x) < 1]` and raises `EncoderDataError`. Executed check: `validate_widths((10,100,1000,10000), 150)` raises `"kernel width(s) [1000, 10000] exceed the sequence length 150"`; `validate_widths((3,12,48,150), 150)` returns cleanly. This is the mechanism that makes the 4-decade arm infeasible, and nothing tests it.
- `scripts/glm_encoder_gate.py` — imports `PYRAMID_WIDTHS` (line 50) and appends two arms under `--pyramid` (lines 317-319). Additive: without the flag the three published arms are byte-identical and the reconcile gate passes (28 tests).
- `soraya_runs/arms-ledger.json` — three `measured` arms (gc 0.3112, conv-3width 0.3234, pyramid-4scale-1.70decades 0.3167) and two `planned`. `arms_ledger` exposes exactly three statuses (`_STATUSES = ('planned','measured','killed')`); `Arm` fields are `name, cost_tier, status, metric, evidence, declared_at, resolved_at` — **no revival field**. `kill_arm` refuses only redundancy reasons (`_REDUNDANCY_REASON` matches bar-met / already-cleared / redundant / sibling / not-needed / goal-is-met), so a substrate-infeasibility reason is mechanically accepted.
- `soraya_runs/negative-results.json` — exists, currently 1 record. `NegativeLog.record(candidate, *, recorded_at, project, token=None, scope='', revival_condition='')`; `NegativeRecord` carries `scope` and `revival_condition`. This is the only surface that can hold a revival trigger, since the arm record cannot.
- Tests: `tests/test_glm_encoder.py` and `tests/test_glm_encoder_gate.py` exist (28 passing in the gate file). Neither references `PYRAMID_WIDTHS`.
- Prior-work replay (`wiki/decisions-log.md`, 35 HIGH-salience entries): one lesson applies directly — *before scheduling a "replace X with Y in file F" step, verify F currently contains X.* Both target lines were read and confirmed verbatim before this plan was written. A second lesson warns that running another planning skill on a closeout-ready artifact is a decision-avoidance pattern; see Risk Flags.

### Reusable-Code Survey
- `~/.claude/skills/soraya/scripts/negative_accrual.py` — `NegativeLog.record(..., scope=, revival_condition=)` already models a scoped negative with a revival trigger; use it rather than hand-writing JSON into the store.
- `~/.claude/skills/soraya/scripts/arms_ledger.py` — `ArmsLedger.kill_arm(criterion, name, reason=, evidence=)` enforces the redundancy refusal; use it rather than editing `arms-ledger.json` by hand.
- `~/.claude/skills/soraya/scripts/ce_ledger.py` — `append_action(ledger_path, action_class, description, outcome, date)` is the enum-validated AC9 row writer already used for this family's two existing rows.
- `dna_decode/glm/encoder.py::validate_widths` — call it directly in the new test rather than reimplementing the width-vs-length arithmetic.
- None of the plan's own logic is new code; searched: `graphify-out/GRAPH_REPORT.md` (absent), `src/`, `lib/`, `utils/`, `common/`, `helpers/` (all absent), `dna_decode/glm/`, `scripts/`, `tests/`, `~/.claude/skills/soraya/scripts/`.

## Pre-Change Baseline

Captured 2026-10-12 before any edit:

- `wiki/glm_pyramid_arm_2026-10-11.md` line 5 contains `This closes one of them with a number`; line 31 contains `scale is not the missing lever`; line 70 contains `the arm measured here is not`. Verified by `sed -n '5p;31p;70p'`.
- `grep -rn "PYRAMID_WIDTHS" tests/` → **0 matches** (only `scripts/glm_encoder_gate.py:50,317,319` and `dna_decode/glm/encoder.py:62`).
- `tests/test_glm_encoder_gate.py` → 28 passed.
- `soraya_runs/negative-results.json` → 1 record.
- `arms_ledger.ArmsLedger('.').all_ceiling_verdicts()` → `verdict='blocked:untried-arm'`, `measured=['conv-3width-13k-params-150bp','gc-single-feature','pyramid-4scale-1.70decades-150bp']`, `planned=['multiresolution-pyramid-4-decades-UNRUN','substrate-change-constructed-eukaryotic-or-human-panel-UNRUN']`, `best='conv-3width-13k-params-150bp'`.
- Frozen five (`amr_rules.py`, `calibrated_amr_rules.json`, `mic_tiers.py`, `shipped_decoder_surface.py`, `cohort_manifest.py`) → `git diff --stat` empty; prospective lock `verify_lock` → `ok=True`.

## Verification Signal

- `grep -c "closes one of them" wiki/glm_pyramid_arm_2026-10-11.md` → `0`; the replacement sentence names "a third measured arm" and no longer claims an arm was closed.
- `grep -c "scale is not the missing lever" wiki/glm_pyramid_arm_2026-10-11.md` → `0`; the substituted claim is scoped to 150 bp tiles and the registered protocol.
- `grep -rn "PYRAMID_WIDTHS" tests/` → **≥1 match**, and the new test fails if `PYRAMID_WIDTHS` is widened past 150 or if `validate_widths` stops refusing an over-wide kernel (both proven by temporary mutation, then reverted).
- `tests/test_glm_encoder.py` and `tests/test_glm_encoder_gate.py` → all pass, count ≥ the 28 baseline plus the new cases.
- `soraya_runs/negative-results.json` → 2 records; the new one has non-empty `scope` and `revival_condition`, and `negative_accrual.scope_warning(rec)` returns no warning.
- Frozen five still byte-unchanged and `verify_lock` still `ok=True` after every step.
- Step 5 only: `all_ceiling_verdicts()` no longer reports `blocked:untried-arm` **for the reason the user authorized**, and `best` is unchanged at `conv-3width-13k-params-150bp`.

## Implementation Steps

### Step 1: Correct the two memo self-contradictions
Files: wiki/glm_pyramid_arm_2026-10-11.md
Depends on: none

**What changes:**
- Line 5: replace `This closes one of them with a number` with wording that says the run ADDS a third measured arm and RECLASSIFIES one of the two unrun arms — because neither declared arm was closed, which `:70` already states correctly.
- Line 31: replace `more scale is not the missing lever` with the supported claim — within 150 bp tiles, widening 1.10 → 1.70 decades under the registered small-CNN protocol did not help.
- Add one sentence recording that the earlier headline over-reached its own `honest_limits` block, so the correction is auditable rather than silent.
- Do NOT touch the results table, the paired statistics, or the substrate census — the measurement was not in question.

**Test strategy:**
- `grep -c` both retired phrases → expect `0` each (the Verification Signal greps).
- Re-read `:70` and the `honest_limits` list to confirm the corrected headline no longer contradicts either.
- No automated test; this is prose in a dated artifact. Correctness is the two greps plus the contradiction re-read.

### Step 2: Pin PYRAMID_WIDTHS and the over-wide-kernel refusal with tests
Files: tests/test_glm_encoder.py
Depends on: none

**What changes:**
- Add a test asserting `PYRAMID_WIDTHS == (3, 12, 48, 150)`, that every width is ≤ 150, and that its span is ~1.70 decades — so widening it past the substrate silently becomes a failing test rather than a runtime surprise.
- Add a test asserting `validate_widths((10, 100, 1000, 10000), 150)` raises `EncoderDataError` naming the offending widths, and that `validate_widths(PYRAMID_WIDTHS, 150)` returns cleanly. This is the mechanism that makes the 4-decade arm infeasible; it is currently untested.
- Add a test asserting `encoder_feature_width(PYRAMID_WIDTHS)` is 4 branches' worth (256 at the default 32 channels), so a branch silently dropped from the tuple is caught.
- Docstring each test with WHY it exists, naming the 4-decade arm it bounds.

**Test strategy:**
- Run the encoder test module plus the gate test module; expect all green and a count above the 28 baseline.
- Prove non-vacuity by temporary mutation: widen `PYRAMID_WIDTHS` to include 1000 and confirm the new test FAILS; relax the `too_wide` check in `validate_widths` and confirm the refusal test FAILS. Revert both and confirm green. A guard never shown to fail is not known to work.

### Step 3: Record the substrate infeasibility as a revivable scoped negative
Files: soraya_runs/negative-results.json
Depends on: none

**What changes:**
- Via `NegativeLog.record(...)` — not by hand-editing the store — add one record for the 4-decade pyramid arm with:
  - `scope` = the narrowest wording the executed evidence supports: this arm is unrunnable on GSE144621 150 bp tiles because `validate_widths` refuses kernels wider than the sequence; it is NOT evidence against multi-decade pyramids on a long-context substrate.
  - `revival_condition` = a substrate with ≥10,000 bp of context (or any substrate where `validate_widths` admits four decades) becomes available for this criterion.
  - `evidence` = the executed `validate_widths` refusal plus the substrate length census (tile max 150, neg_control 150, random 150, pos_control 48-75).
- Confirm `scope_warning(rec)` returns empty — a missing scope or revival trigger is exactly what that helper exists to catch.

**What this buys:** it makes Step 5 safe under EITHER user decision. If the arm is later killed, the idea survives with a dated revival trigger; if it is held, the reason it is held is on the record instead of in prose. This step is purely additive and needs no authority.

**Test strategy:**
- Re-read the store and assert it has exactly one more record than the baseline 1, that the new record's `scope` and `revival_condition` are non-empty, and that `scope_warning` is silent.
- Confirm the pre-existing record is byte-unchanged (append-only, not a rewrite).

### Step 4: Commit the corrections and record the batch
Files: project_state/glm-learned-representation-2026-10-07.md
Depends on: Step 1, Step 2, Step 3

**What changes:**
- Append one AC9 Action Log row (`ce_ledger.append_action`, class `edit-local-code`) describing the three corrections, naming both retired phrases verbatim so the record says what was wrong rather than only what is now right.
- Commit on the existing side branch `soraya/f3-calibration-and-cabbage-cell` — do NOT push, and do NOT commit into any branch another session holds.
- Re-verify the frozen five are byte-unchanged and the prospective lock still verifies, and state both in the commit body.

**Test strategy:**
- `git diff --stat HEAD` against the frozen five → empty.
- `verify_lock` → `ok=True`, no drift.
- Both encoder test modules green; `git status --short` shows no unintended files staged.

### Step 5: Resolve the arms-ledger state — AUTHORITY-GATED, do not execute without a decision
Files: soraya_runs/arms-ledger.json, project_state/glm-learned-representation-2026-10-07.md
Depends on: Step 3, Step 4

**What changes:** exactly one of the following, per the user's answer to the single narrowed question. **`/execute-plan` must STOP here if no decision has been recorded** — this step is not auto-runnable.

- **Option A — kill within this criterion.** `kill_arm(criterion, 'multiresolution-pyramid-4-decades-UNRUN', reason=<substrate-infeasibility wording>, evidence=<the validate_widths refusal + census>)`. Mechanically permitted (the reason matches no redundancy token). Also move or kill `substrate-change-constructed-eukaryotic-or-human-panel-UNRUN`, which the review found MISFILED: the criterion names the GSE144621-derived `0.941` ceiling, and a eukaryotic panel would have its own different ceiling. Result: criterion closes on three measured arms.
- **Option B — split the criterion.** Keep the current criterion scoped to 150 bp GSE144621 tiles (closes now, three measured arms, no planned arm) and declare a separate long-context criterion carrying the 4-decade arm. Nothing is killed. This redefines terminal closure, which is why it needs ratification.
- **Option C — hold.** Leave both planned and accept that the criterion stays blocked indefinitely. Step 3 has already made the reason legible, so this is now a recorded choice rather than drift.

Whichever is taken, append one AC9 Action Log row via `ce_ledger.append_action` naming the option, the authorizing decision, and the resulting `all_ceiling_verdicts()` output.

**Test strategy:**
- Before: snapshot `all_ceiling_verdicts()`. After: re-read and assert the verdict changed only in the way the chosen option predicts, and that `best` is still `conv-3width-13k-params-150bp` — no option may alter which arm is best.
- Assert no `measured` arm's `metric` or `evidence` changed.
- Back up `arms-ledger.json` before the write and diff after, so the change is reviewable and revertible.

## Execution Preview

Wave 0 (3 parallel):  Step 1 — Correct the two memo self-contradictions, Step 2 — Pin PYRAMID_WIDTHS and the over-wide-kernel refusal with tests, Step 3 — Record the substrate infeasibility as a revivable scoped negative
Wave 1 (1 parallel):  Step 4 — Commit the corrections and record the batch
Wave 2 (1 parallel):  Step 5 — Resolve the arms-ledger state (AUTHORITY-GATED)

Critical path: Step 1 → Step 4 → Step 5 (3 waves)
Max parallelism: 3 agents

Note: Parallel execution requires a git repository with a configured remote. If unavailable, /execute-plan falls back to sequential mode.

## Risk Flags

- Severity: high — **Step 5 is authority-gated and must not auto-execute.** It changes a governance record that decides when a ceiling criterion is terminal, and two of its three options are hard to reverse in spirit (killing an expensive arm, or redefining an acceptance bar). `/execute-plan` reaching Wave 2 without a recorded user decision must stop, not pick a default. Mitigation: Step 3 runs first and is purely additive, so the expensive idea is preserved with a revival trigger under every option.
- Severity: medium — **This plan may itself be the decision-avoidance pattern.** A HIGH-salience lesson in `wiki/decisions-log.md` warns that invoking another planning skill on a closeout-ready artifact is a known failure mode, and that the honest move is "apply N obvious in-place edits + commit". Steps 1-3 are exactly that: three small, reversible edits that did not need a plan. The part that genuinely needed planning is the Step 5 gate. Mitigation: Steps 1-3 are executable immediately and independently; nothing is blocked waiting on the plan.
- Severity: medium — **File overlap between steps:** none within any wave. Wave 0's three steps touch three disjoint files (`wiki/glm_pyramid_arm_2026-10-11.md`, `tests/test_glm_encoder.py`, `soraya_runs/negative-results.json`). Steps 4 and 5 both touch `project_state/glm-learned-representation-2026-10-07.md`. The first draft had them in ONE wave sharing that file — caught by the wave checker, not by review. Resolved by ORDERING the commit (Step 4) BEFORE the gated step (Step 5) rather than the reverse: making the commit depend on the gate would have blocked three reversible corrections behind an authority decision.
- Severity: medium — **External tool surfaces, all confirmed by direct introspection, with no surface left unconfirmed:** `arms_ledger.ArmsLedger.kill_arm(self, criterion, name, reason, evidence)` and `declare_arm(self, criterion, name, cost_tier='unknown')` and `measure_arm(self, criterion, name, metric, evidence)` — verified via `inspect.signature`; `_STATUSES = ('planned','measured','killed')` and `_REDUNDANCY_REASON` regex — verified by reading module constants. `NegativeLog.record(self, candidate, *, recorded_at, project, token=None, scope='', revival_condition='')` and `NegativeRecord` fields — verified via `inspect.signature` and `__dataclass_fields__`. `ce_ledger.append_action(ledger_path, action_class, description, outcome, date)` — verified via `inspect.signature` and already used twice on this family. `validate_widths` behavior — verified by EXECUTING both the passing and raising cases.
- Severity: low — **`python3` is not on this host; only `python` resolves.** The toolkit invocation in this skill's own instructions specifies `python3`, which exits 127 here. Every command in this plan uses `python`. Mitigation: noted in Lens status degradations and in the toolkit footer handling.
- Severity: low — **Steps depending on external services or env vars:** none. Every step is local file I/O plus local pytest. No network, no Docker, no GPU, no paid surface.
- Severity: low — **Files modified by build tools or code generators:** none. `soraya_runs/*.json` are written by the helper APIs this plan calls, not by a generator; both are backed up or append-only.
- Severity: low — **Restructuring applied:** (a) the commit was split out of Steps 1-3 so one reviewable commit carries the three corrections; (b) Step 3 is ordered BEFORE the gated step so the revival trigger exists before any destructive option can run; (c) the commit and the gated step were SWAPPED after the wave checker found them sharing a ledger file in one wave — the commit now runs first so it is not hostage to the gate.

## Open Questions

- **The one narrowed question from review:** may the current GSE144621-150 bp criterion close by retiring or moving substrate-incompatible arms, while a separate long-context criterion carries the 4-decade idea forward? This is Step 5's gate; Options A/B/C above are the concrete branches.
- Was `substrate-change-constructed-eukaryotic-or-human-panel-UNRUN` ever intended to measure *bacterial* promoter activity, or was it attached to this criterion by mistake? The criterion names a GSE144621-derived `0.941` ceiling, which a eukaryotic panel cannot share — evidence favours "misfiled", but the arm's author's intent is not recorded.
- Should `arms_ledger` gain a non-destructive `deferred` / `blocked_by_substrate` status? That would make Step 5 unnecessary in future, but the module is shared tooling another session owns, so it is out of scope here and noted only.

## Verification

1. Run `tests/test_glm_encoder.py` and `tests/test_glm_encoder_gate.py` — all green, count above the 28 baseline.
2. Run the two Verification-Signal greps on the memo — both retired phrases return `0`.
3. Prove the new tests non-vacuous by temporary mutation (widen `PYRAMID_WIDTHS`; relax `validate_widths`), confirm each fails, revert, confirm green.
4. Re-read `soraya_runs/negative-results.json` — 2 records, new one has non-empty `scope` + `revival_condition`, `scope_warning` silent, prior record unchanged.
5. Re-run `python scripts/glm_encoder_gate.py --seeds 2 --no-download` WITHOUT `--pyramid` and confirm the three published arms and `reconcile_ok` are unchanged — the additive flag must stay additive.
6. `git diff --stat HEAD` against the frozen five → empty; `verify_lock` → `ok=True`.
7. Step 5 only, after a decision: `all_ceiling_verdicts()` reflects the authorized option and `best` is still `conv-3width-13k-params-150bp`.

<!-- toolkit: check=clean waves=clean gate=fired:open-questions,severity-high mode=conversation -->
