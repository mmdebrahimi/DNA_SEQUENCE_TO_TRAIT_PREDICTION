# Return Decision Tree Patch + Bat Hardening — Technical Plan

> Apply the /review synthesis's 3 surgical correctness fixes + structural restructure to `wiki/return_decision_tree_2026-05-16.md`, plus the one-line `-u` flag fix to `run_stage1b_detached.bat`. Implementation-ready under 4 steps; max parallelism 2 (doc + bat are independent).

---

## Problem Statement

The just-shipped return decision tree has 3 verified correctness defects (malformed `/project-state` calls that would REFUSE; wrong `taskkill /F` flush claim; missing MSYS path-mangling caveat) plus 2-3 structural UX issues (buried lede, overbuilt Step 3 + Step 4, heartbeat dependence on a flaky signal). The current Stage 1b launch script is unbuffered-mode-unaware — next detached relaunch will repeat the silent-buffering failure mode.

Fix scope: surgical eng patches FIRST (commands must work as written), then optional CEO restructure for fatigue-resilience, plus one-line `.bat` hardening. All changes are reversible single-file edits with no module-coupling risk.

## Codebase Context

Verified via direct read during /review synthesis:

- **`wiki/return_decision_tree_2026-05-16.md`** — 95 lines; file paths + globs in §1 commands match the producing scripts (`scripts/stage1_n40_cipro.py:533` writes `wiki/stage1_n40_cipro_mean-plus-max_<date>.md`; `cipro_attribution_preflight.py:85,245` writes `wiki/cipro_attribution_preflight_<date>.md` + `.per_strain.json`). Glob patterns in §1 verified correct.
- **`run_stage1b_detached.bat`** — 4-line .bat; line 4 is `uv run python scripts\stage1_n40_cipro.py > stage1b_detached.log 2>&1`; no `-u` flag → block-buffered stdout → `stage1b_detached.log` stays 0 bytes until process exits.
- **`/project-state` SKILL.md** at `C:\Users\Farshad\.claude\skills\project-state\SKILL.md` — Step 1 mandates `<project-id>` as first positional arg + every mutation flag has required companion args. The slug for this project is `dna-decode-2026-05-11` (per CLAUDE.md + ledger filename).
- **CLAUDE.md gotcha** (added today's session in commit `4e3474c`): "Git Bash silently breaks Docker `-v <host>:/<container>` volume mounts via MSYS path conversion. For direct `docker run` from Git Bash: prefix every invocation with `MSYS_NO_PATHCONV=1`." Same applies to `tasklist /FI "PID eq <N>"` (verified).
- **Prior decisions-log entries surveyed**: HIGH-salience phase→Evidence-Packets reframing (today, 2026-05-15) — not directly applicable. HIGH-salience verdict-vs-action-separation lesson (2026-05-14) — already reflected in the doc structure. No prior decision conflicts with this patch plan.

No `.repo-index/repo-summary.md` present; skipped repo-index step.

## Design Decisions

### D1: Eng correctness patches come FIRST, restructure SECOND

**Decision:** Step 1 fixes the 3 verified correctness defects (malformed commands + wrong flush claim + missing MSYS caveat). Step 2 restructures for fatigue-resilience. Sequenced as Step 1 → Step 2 because both modify the same file.

**Rationale:** Eng fixes alone leave a wall-of-text but the commands work. Restructure alone (without eng fixes) would leave broken commands in the appendix. Together they cover both "commands work" + "user can read under fatigue." Restructure depending on patches ensures the patches don't get accidentally reverted by the restructure rewrite.

**Trade-off:** Could have merged into one bigger step. Splitting keeps each step's blast radius small + makes review easier mid-execution.

### D2: `-u` flag (not env var) for the .bat fix

**Decision:** Modify `run_stage1b_detached.bat` line 4 to add `-u` flag to `uv run python`. Don't use `set PYTHONUNBUFFERED=1`.

**Rationale:** `-u` is more explicit + targeted to this one invocation. Env var would leak to any nested python.exe. The .bat is a one-off launcher; explicit > implicit.

**Trade-off:** Either works. `-u` was chosen for tighter blast radius.

### D3: Merge enumerated sub-steps into 4 actual steps

**Decision:** Original enumeration had 8 sub-steps (one per fix); merged into 4 because adjacent doc-fix sub-steps all modify the same file + serialize anyway.

**Rationale:** Merging same-file consecutive steps reduces wave count without reducing safety. Spec says "merge tiny steps — if two adjacent steps touch the same files and one depends on the other, merging reduces wave count."

**Trade-off:** Bigger steps mean longer review-during-execution windows. Acceptable for this small-blast-radius docs work.

### D4: Don't touch the running Stage 1b process

**Decision:** All 4 steps avoid modifying the in-memory PID 4736 (currently running detached Stage 1b). Step 3's `.bat` fix only affects FUTURE relaunches.

**Rationale:** Stage 1b at +6 min mid-bootstrap-CI is expensive to lose. The `.bat` change is for the next session's relaunch, not this one.

**Trade-off:** Future relaunches benefit; current run still has buffered stdout in memory (acceptable since verdict packet writes via `write_text` regardless of stdout state).

## Implementation Plan

### Step 1: Apply 3 correctness patches to return_decision_tree
**Files:** `wiki/return_decision_tree_2026-05-16.md`
**Depends on:** none

Three Edit operations against unique strings (byte-preserving):

1. **Step 3 ledger commands (lines 60-65):** replace `/project-state --append-action --class run-tests with the Stage 1b outcome row` (descriptive prose, non-callable) with the full callable form: `/project-state dna-decode-2026-05-11 --append-action --class run-tests --description "Stage 1b N=38 cipro outcome" --outcome "<paste verdict text here>"`. Same fix for `--update-hypothesis H17` (add slug + retain existing companion args) and `--refresh-frame` (add slug + `--current-state "<text>"` companion).

2. **§2B/§2D `taskkill /F` flush claim (lines 49-50, 67-70):** replace `inspect stage1b_detached.log for partial output (may need to flush by killing process explicitly with taskkill /F /PID <pid>)` with: `taskkill /F` does NOT flush Python stdio buffers. Partial stdout in `stage1b_detached.log` is unrecoverable if python.exe was launched without `-u` / `PYTHONUNBUFFERED=1`. The verdict packet (file write via `write_text`, not stdout) is the authoritative signal — its presence/absence is reliable independent of buffering.

3. **§1 tasklist caveat (line 17):** keep `tasklist | grep python.exe` but add a parenthetical: `(grep form works in Git Bash; if you switch to PID-specific filter use MSYS_NO_PATHCONV=1 tasklist /FI "PID eq 4736" — MSYS otherwise mangles the quoted arg per CLAUDE.md gotcha)`.

**Key details:**
- All 3 edits are byte-preserving Edit operations against unique strings; no risk of overlap.
- /project-state slug for THIS project: `dna-decode-2026-05-11` (verified against ledger filename + CLAUDE.md).
- Edit ordering within the step: doesn't matter (3 disjoint regions).

**Test strategy:**
- Visual: `grep -c '/project-state --append-action' wiki/return_decision_tree_2026-05-16.md` should return 0 occurrences without slug.
- Visual: `grep 'flush by killing process explicitly'` should return zero matches.
- Visual: `grep 'MSYS_NO_PATHCONV'` should return exactly one new occurrence in §1.
- Mock-run: paste one of the patched `/project-state` commands into shell with a test slug; confirm syntax matches SKILL.md Step 1 arg validation.

### Step 2: Restructure for fatigue-resilience (CEO scope reduction)
**Files:** `wiki/return_decision_tree_2026-05-16.md`
**Depends on:** Step 1

Structural rewrite (same file, larger blast radius than Step 1's targeted edits):

1. **Cut entire Step 4** (EP2 cef + tet smoke prereq) — EP2 has its own design plan at `plans/EP2_Cef_Tet_Smoke_Design_Plan.md`; doesn't belong in triage.
2. **Cut Step 3's 6-item ledger checklist** down to one line: "Log outcome via `/project-state dna-decode-2026-05-11 --append-action --class run-tests --description '<X>' --outcome '<verdict>'`, commit, push."
3. **Add 5-line bash triage block at top** (above existing §1) that runs all 4 presence checks + prints a 4-letter state code (`S1+PF+` / `S1+PF-` / `S1-PF+` / `S1-PF-`) so first command resolves quadrant.
4. **Add 2×2 verdict table** mapping (Stage 1b stage2_action × preflight verdict) → next-action one-liner. Replaces the verbose "Combined call" bulleted prose in §2A.
5. **Move §2B + §2C + §2D into a "Recovery procedures" appendix below the fold.** Keep them complete but de-prioritize.
6. **Add timestamp gate** at top of §1: "If you're returning >7 hr after launch (verdict expected ~19:43 EDT), Stage 1b should be DONE — packet present OR process dead, no third state. If process is still alive past 7 hr, suspect hang."
7. **Add Databricks-burst sleep-on-it gate** in the new 2×2 table: any cell that triggers `BURST_STAGE_2` carries a "if returning <2 hr energy, defer firing burst to next session" note.
8. **Remove heartbeat dependence from decision branches.** Keep heartbeat log mention as informational appendix line only ("heartbeat log shows process aliveness signal at ~2 lines/sec until packet lands; never load-bearing for decisions").

**Key details:**
- Above-the-fold target: ~30 lines total (was ~80).
- Appendix preserves the original Quadrant B/C/D recovery procedures (post-Step-1 patched).
- No new external file references beyond what's already in the doc.

**Test strategy:**
- Visual: confirm above-the-fold (top 30 lines) contains: bash triage block, 4-letter state-code emission, 2×2 verdict table, timestamp gate.
- Visual: confirm appendix exists with §B/§C/§D recovery details.
- Cognitive: read top 30 lines as if returning fatigued; can you make a decision without scrolling? If no, iterate.
- Counter-check: confirm Step 1's 3 patches are PRESERVED through the restructure (the restructure must not silently revert the patch fixes).

### Step 3: Add `-u` unbuffered flag to detach helper
**Files:** `run_stage1b_detached.bat`
**Depends on:** none

Single-line .bat edit:

- `run_stage1b_detached.bat` line 4 — replace `uv run python scripts\stage1_n40_cipro.py > stage1b_detached.log 2>&1` with `uv run python -u scripts\stage1_n40_cipro.py > stage1b_detached.log 2>&1`. The `-u` flag forces stdin/stdout/stderr to be unbuffered.

**Key details:**
- Effect: future detached relaunches will write `stage1b_detached.log` in near-real-time instead of staying empty until process exit. Current Stage 1b run (PID 4736, launched before this patch) is unaffected — it already has buffered stdout in memory.
- Risk: `-u` has negligible perf overhead (<1%) on long-running CPU-bound jobs.
- Alternative considered: `set PYTHONUNBUFFERED=1` env var. Equivalent effect; `-u` is more explicit + targeted to this one invocation. Chose `-u` for tighter blast radius (D2).

**Test strategy:**
- Quick: `uv run python -u -c "import sys; print(sys.flags.unbuffered)"` → expect `1`.
- End-to-end: would require relaunching Stage 1b to verify the log fills in real-time. Defer this validation; don't kill the current healthy run.

### Step 4: Commit + push patches
**Files:** `wiki/return_decision_tree_2026-05-16.md`, `run_stage1b_detached.bat`
**Depends on:** Step 2, Step 3

- Stage both modified files via `git add wiki/return_decision_tree_2026-05-16.md run_stage1b_detached.bat`.
- Commit with message: `docs+chore: harden return decision tree (review patches) + unbuffered detach launcher`. Body: enumerate the 3 correctness fixes (Step 3 slug + taskkill flush + MSYS caveat) + the structural restructure (Step 4 cut, 2×2 table, appendix demote, heartbeat de-emphasis) + the `-u` flag rationale.
- Push to `origin/main`.

**Key details:**
- No tests to run; this is docs + a .bat config tweak.
- Pre-push: `git status --short` should show only these 2 files staged + nothing else surprising.

**Test strategy:**
- Confirm `git log -1` shows the new commit.
- Confirm `git status --short --branch` reports `## main...origin/main` (no `[ahead N]`) after push.
- Confirm the file's `wiki/` path resolves from anywhere via `git show HEAD:wiki/return_decision_tree_2026-05-16.md | head -10`.

## Execution Preview

```
Wave 0 (2 parallel):  Step 1 — Apply 3 correctness patches,   Step 3 — Add -u flag to detach helper
Wave 1 (1):           Step 2 — Restructure for fatigue-resilience
Wave 2 (1):           Step 4 — Commit + push
```

- **Critical path:** Step 1 → Step 2 → Step 4 (3 waves)
- **Max parallelism:** 2 agents (Wave 0)
- **Total wall-time estimate:** ~20 min sequential / ~15 min parallel (Step 2 is the biggest single chunk; restructure work)

Note: Parallel execution requires a git repository with a configured remote — both exist (origin = `mmdebrahimi/DNA_SEQUENCE_TO_TRAIT_PREDICTION`). If unavailable, /execute-plan falls back to sequential mode.

## Risk Flags

- **File overlap:** Steps 1 + 2 both modify `wiki/return_decision_tree_2026-05-16.md`. Dependency enforced (Step 2 depends on Step 1). Step 2's restructure MUST preserve Step 1's correctness patches — flagged in Step 2's test strategy as a counter-check.
- **No transitive imports:** This plan is docs + a .bat tweak; no code modules + no test dependencies.
- **No build-tool / code-gen interactions:** Plain markdown + plain CMD batch; no generators touch either file.
- **External services:** Step 4 push requires GitHub origin reachability. If push fails, treat as recoverable (retry later); the local commit is durable.
- **Restructuring applied during planning:** Originally enumerated 8 sub-steps (one per fix). Merged into 4 steps because the doc-fix sub-steps all modify the same file + serialize anyway; merging reduces wave count without reducing safety (D3). The `.bat` fix was originally going to bundle into the commit step but split off as Step 3 because it touches a different file + can fire in parallel with Step 1.
- **Stage 1b is running.** None of these steps touch the running process. Step 3 only affects FUTURE relaunches, not the in-memory PID 4736. Step 4's git push lands separately from runtime artifacts (verdict packets land via `write_text` regardless of git state).
- **Pre-departure context:** the user is outside. This plan exists to be executed when they return OR by Claude on explicit "go." Either way, no time pressure.

## Verification

After all 4 steps complete:

1. **Doc passes correctness audit:**
   - `grep -c '/project-state --' wiki/return_decision_tree_2026-05-16.md` returns 0 (all uses have slug).
   - `grep -c 'flush by killing' wiki/return_decision_tree_2026-05-16.md` returns 0.
   - `grep -c 'MSYS_NO_PATHCONV' wiki/return_decision_tree_2026-05-16.md` returns ≥1.
2. **Doc passes UX audit:**
   - Above-the-fold (top ~30 lines) contains the bash triage block + 4-letter state-code emission + 2×2 verdict table.
   - Step 4 (EP2) section no longer exists.
   - Heartbeat is not load-bearing for any decision branch.
3. **`.bat` fix lands:**
   - `grep 'python -u scripts' run_stage1b_detached.bat` returns 1 match.
4. **Git state clean:**
   - `git status --short --branch` reports `## main...origin/main` (no ahead/behind).
   - `git log --oneline -1` shows the new commit hash + descriptive subject line.
5. **Stage 1b unaffected:**
   - `tasklist | grep python.exe` still shows PID 4736 (or whatever Stage 1b's PID is) alive.
   - `wiki/stage1_n40_cipro_*.md` file list unchanged (no new verdict packet means Stage 1b is still computing; verdict packet means it completed during the patch work — handle per the (now-improved) return decision tree).
