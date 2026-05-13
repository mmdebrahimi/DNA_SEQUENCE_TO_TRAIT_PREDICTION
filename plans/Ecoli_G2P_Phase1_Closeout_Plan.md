# Ecoli G2P Phase 1 Closeout Plan

> Wrap up the stalled `/execute-plan` epilogue for `Ecoli_G2P_Phase1_Ship_Path_Plan.md` — toolchain restore, doc reconciliation, first authoritative test pass, archive, state cleanup, push, final report.

---

## Problem Statement

`/execute-plan` on `plans/Ecoli_G2P_Phase1_Ship_Path_Plan.md` ran to code-complete (all 7 ship-path steps committed: `3ca785a`, `0316cd7`, `2aba7d0`, `4c72157`, `401223c`, `95ce043`, `45f9447`, plus `39ca362` for Step 11.5 tests), then stalled mid-epilogue when the prior session ran out of context. Closeout is partially done and inconsistent:

- Two `.claude/execute-plan-state/*.json` files are stale (mark Steps 3.5/13/14/15/16/17/11.5 as `pending` or `in_progress` despite commits landing).
- Untracked: `CLAUDE.md`, `FUTURE_FEATURES.md`, `LESSONS_LEARNED.md`, `TODOS.md`, `wiki/`, `plans/Ecoli_G2P_Phase1_Ship_Path_Plan.md`, `.claude/execute-plan-state/Ecoli_G2P_Phase1_Ship_Path_Plan.json`.
- `TODOS.md` has open `[ ]` checkboxes for shipped work; `CLAUDE.md` has "not yet implemented" annotations for Steps 13/14 and "(pending)" on Wave 3.5 hardening — committing as-is would encode wrong state.
- `test_baseline: null` in both state files; only test-epilogue commit is `39ca362` covering Step 11.5 only.
- `uv` is not on PATH; `.venv` is gone — `uv run pytest tests/ -v` cannot run as-is.
- `/documentation` skill never ran.
- `/retrospective` skill ran in prior session, output lost to compaction.
- `main` is ahead of `origin/main` by 1 commit.
- No `executed_plans/` precedent exists in the repo.

This plan finishes closeout; it does not re-run any code wave.

## Design Decisions

### D1: Selective expansion over hold scope

**Decision:** Add four concrete corrections to the original 6-step closeout (toolchain restore, reorder docs-before-commit, explicit doc deliverables, record test outcome non-gating).

**Rationale:** Engineering review surfaced a hard blocker (`uv` missing) and a step-ordering bug (commit before `/documentation` would commit stale docs). Product review surfaced a credibility framing concern (unit-tests-only ship gate vs. G2P-platform claims).

### D2: Real-data validation = Phase 2 entry criterion

**Decision:** Document as Phase 2 entry criterion in the final report; do not add a real-data step inside closeout.

**Rationale:** Adding it to closeout strengthens "Phase 1 shipped" credibility for any external consumer (SME, future-self) but expands scope. Leaving it to Phase 2 ships closeout cleanly but inherits silent debt unless explicitly named in the report.

### D3: Toolchain restore — uv with Python 3.11 pin

**Decision:** Install `uv` (Windows: `irm https://astral.sh/uv/install.ps1 | iex`) then `uv python pin 3.11` + `uv sync`.

**Rationale:** uv matches project convention (`pyproject.toml` is uv-shaped). Explicit Python 3.11 pin prevents uv from picking 3.12+ which may have compatibility issues with project deps.

### D4: Archival convention — status header + git tag

**Decision:** In-file `Status: archived 2026-05-12` header + git tag `phase-1-shipped`.

**Rationale:** Solo-project preference is status-header + git tag (cheapest, no convention precedent to set). Tag creation is conditional on green test gate.

### D5: Step ordering — `/documentation` before commit

**Decision:** Run `/documentation` skill BEFORE committing the untracked files.

**Rationale:** TODOS.md and CLAUDE.md contain factually wrong state (open checkboxes for shipped work, "not yet implemented" annotations on shipped steps). Committing then doc-skill-fixing would create two commits where one suffices; doc-skill-fixing then committing produces a single clean commit.

### D6: Retrospective re-derivation, not skip

**Decision:** Spend 15 min re-deriving lessons from git log + plan diff, append to `LESSONS_LEARNED.md`.

**Rationale:** Retrospective output from the prior session is lost to compaction. Cost is low (15 min); value is high (carries highest-density learnings of any closeout artifact). Trade-off was "accept loss" — rejected as false economy.

### D7: Test outcome recorded, not gated

**Decision:** Step 3 records `{passed, failed, errors, collected, pytest_exit_code, env}` regardless of color; non-green produces follow-up TODOS rows before archive — does NOT block archive.

**Rationale:** 195 tests have never been run end-to-end under real deps; first authoritative pass will likely surface integration bugs. Treating green as table-stakes risks indefinite closeout loop. Recording outcome + filing follow-ups preserves the audit trail.

### D8: State-file cleanup conditional on test gate

**Decision:** Delete state files only if test gate passes (green). If red, preserve state files OR create `closeout-status.json` with `verification_deferred: true`.

**Rationale:** Preserves machine-readable "not verified" signal for future sessions when tests fail.

### D9: `.claude/execute-plan-state/` added to `.gitignore`

**Decision:** Same commit as Step 4 (the untracked-files commit) adds `.claude/execute-plan-state/` to `.gitignore`.

**Rationale:** Internal scaffolding; future runs will recreate it. Avoids future untracked-file noise.

### D10: Git tag after report, conditional on green

**Decision:** Create `phase-1-shipped` tag in Step 7 (after report commit), not Step 5, and only if Step 3 test gate passed.

**Rationale:** Tag should point at the commit containing the ship report, not an earlier intermediate state. Red-path ships should not create a "shipped" tag.

---

## Implementation Plan

### Step 0: Toolchain restore
Files: pyproject.toml, .python-version
Depends on: none

Install `uv` with explicit Python 3.11 lock:
```powershell
# Install uv if not on PATH
irm https://astral.sh/uv/install.ps1 | iex
# Pin Python 3.11 and sync
uv python pin 3.11
uv sync
```
**Verify:** `uv run pytest --collect-only -q` reports ≥195 collected, 0 collection errors.

### Step 1: Run /documentation skill (main-orchestration-only)
Files: TODOS.md, CLAUDE.md, README.md, LESSONS_LEARNED.md, FUTURE_FEATURES.md
Depends on: Step 0

> **Note:** `/documentation` skill runs in main orchestration context, not worktree agent.

- `TODOS.md`: archive Phase 1 done sections (Wave 3.5 C7/C8/M4/M5/cleanup, Steps 15/14/13/17/11.5/16); keep post-Phase-1 polish + known limitations only.
- `CLAUDE.md`: strike "not yet implemented" annotations on Steps 13/14; flip Wave 3.5 hardening from "(pending)" → "shipped".
- `README.md`: spot-check the `# 4. Populate the embedding cache (deferred...)` note is still accurate.
- `LESSONS_LEARNED.md` + `FUTURE_FEATURES.md`: append-only; no rewrite needed.

### Step 2: Retrospective re-derivation
Files: LESSONS_LEARNED.md
Depends on: Step 1

15-min pass over `git log --oneline` + the two plans.
Append 3-5 bullets to `LESSONS_LEARNED.md` capturing Phase-1-shipping lessons.

### Step 3: Full test pass + outcome recording
Files: test_baseline.json, TODOS.md
Depends on: Step 0

Run `uv run pytest tests/ -v` from `dna_decode/`.

**Record structured outcome to `test_baseline.json`:**
```json
{
  "date": "2026-05-12",
  "passed": <int>,
  "failed": <int>,
  "errors": <int>,
  "collected": <int>,
  "pytest_exit_code": <int>,
  "env": {"python": "<version>", "uv": "<version>"},
  "command": "uv run pytest tests/ -v"
}
```

**Gate condition for green:** `pytest_exit_code == 0 && collected >= 195 && failed == 0 && errors == 0`

**If red (gate fails):** file specific failures as new `TODOS.md` rows; continue to Step 4 anyway (per D7).

### Step 4: Triage untracked files & single commit
Files: .gitignore, CLAUDE.md, FUTURE_FEATURES.md, LESSONS_LEARNED.md, TODOS.md, wiki/, plans/Ecoli_G2P_Phase1_Ship_Path_Plan.md, test_baseline.json
Depends on: Step 1, Step 2, Step 3

- Add `.claude/execute-plan-state/` to `.gitignore`.
- Commit all untracked/modified files.
- Single commit, message: `docs: Phase 1 closeout — docs reconciled, ship-path plan + closeout-state captured`.

### Step 5: Archive the ship-path plan (commit only, no tag)
Files: plans/Ecoli_G2P_Phase1_Ship_Path_Plan.md
Depends on: Step 4

Add in-file status header to `plans/Ecoli_G2P_Phase1_Ship_Path_Plan.md`:
```markdown
> **Status:** archived 2026-05-12 — Phase 1 shipped
```
Commit: `docs: archive ship-path plan with status header`

**Note:** Git tag `phase-1-shipped` deferred to Step 7 (after report, conditional on green gate per D10).

### Step 6: State-file cleanup (conditional)
Files: .claude/execute-plan-state/
Depends on: Step 3, Step 5

**If Step 3 gate passed (green):**
- Delete `.claude/execute-plan-state/Ecoli_G2P_Phase1_Ship_Path_Plan.json` and `.claude/execute-plan-state/Ecoli_G2P_Platform_Technical_Plan.json`.
- The `.gitignore` entry from Step 4 means future runs of `/execute-plan` recreate these without untracked-file noise.

**If Step 3 gate failed (red):**
- Preserve state files OR create `closeout-status.json` with:
```json
{
  "verification_deferred": true,
  "reason": "pytest gate failed",
  "test_baseline": "test_baseline.json"
}
```
- This preserves machine-readable "not verified" signal for future sessions.

### Step 7: Push + final report + conditional tag
Files: wiki/phase1_ship_report.md
Depends on: Step 5, Step 6

1. **Push main:** `git push origin main`
2. **Write report to `wiki/phase1_ship_report.md`:**
   - Header: `## Phase 1 Ship Report (GREEN)` or `## Phase 1 Ship Report (RED — tests deferred)`
   - **Waves shipped:** Waves 0-3 + 1.5/2.5/3.5 hardening (technical plan) + Steps 3.5/13/14/15/17/11.5/16 (ship-path plan).
   - **Baseline test count + env:** from Step 3 `test_baseline.json`.
   - **Doc updates:** what `/documentation` skill changed.
   - **Archive status:** in-file status header + git tag (if green).
   - **What Phase 1 ship gate validated:** code compiles, unit tests pass on real deps (if green), synthetic smoke (`scripts/smoke_pipeline.py`) runs end-to-end.
   - **What Phase 1 ship gate did NOT validate:** prediction accuracy on real E. coli genomes; region-attribution biological correctness; multi-strain leaderboard performance.
   - **Phase 2 entry criterion:** real-data smoke (one E. coli genome end-to-end with prediction + attribution captured; no ground-truth comparison required at entry).
3. **Commit report:** `docs: Phase 1 ship report`
4. **Conditional tag (if green gate passed in Step 3):**
   - `git tag -a phase-1-shipped -m "Phase 1 shipped — all tests green, docs reconciled"`
   - `git push origin phase-1-shipped`
5. **Push main again (includes report commit):** `git push origin main`

---

## Verification

- `git status` shows clean working tree (no untracked, no modified).
- `git log --oneline` shows the closeout commit (Step 4), archive commit (Step 5), and report commit (Step 7) on top of `39ca362`.
- `git ls-files plans/` includes `Ecoli_G2P_Phase1_Ship_Path_Plan.md` with status header.
- `.claude/execute-plan-state/` is empty (if green) or contains `closeout-status.json` (if red).
- `wiki/plans-index.md` has an entry for this closeout plan.
- `origin/main` matches local `main` after Step 7 push.
- If green: `git tag -l phase-1-shipped` shows the tag pointing at the report commit.
- `wiki/phase1_ship_report.md` exists with all required sections.

## Wave Grouping (for /execute-plan)

| Wave | Steps | Rationale |
|------|-------|-----------|
| 0 | Step 0 | Toolchain must exist before anything else |
| 1 | Step 1, Step 3 | Docs and tests can run in parallel (both depend only on Step 0) |
| 2 | Step 2 | Depends on Step 1 (needs doc changes visible) |
| 3 | Step 4 | Depends on Steps 1, 2, 3 (commits all outputs) |
| 4 | Step 5 | Depends on Step 4 |
| 5 | Step 6 | Depends on Steps 3, 5 (needs test outcome + archive done) |
| 6 | Step 7 | Depends on Steps 5, 6 (final push + report) |

**Note:** Step 1 (`/documentation`) is main-orchestration-only — cannot run in worktree agent.
