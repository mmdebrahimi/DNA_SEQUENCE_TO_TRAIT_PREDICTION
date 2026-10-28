# Amend The Saved Correction Plan So Its Two Defects Cannot Ship

## Lens status

Dispatch: conversation-driven (no args) — scope is the pre-exec review findings C1, C2 and M1–M5 against the saved plan.
Repo index: absent (`repo-index-status.sh` not installed) — direct file reading only, no degradation.
Tools: `sentrux` present, `gh` present, `codex` ABSENT — the review that produced these findings was single-voice, which is itself a Risk Flag below.
Project rules: `project-rules.md` absent — Project Rules Check omitted.
Wiki replay: `wiki/` present; `wiki/plans-index.md` is one of the two files under change.
Scope discipline: this is an amendment to **one markdown file plus one index line**. Three steps is the honest size; a seven-step plan to edit a plan would be motion.

## Problem Statement

The plan saved at `plans/Correct_Two_Wrong_Claims_I_Committed_And_Make_The_Next_Check_Runnable_Plan/technical-plan.md` is sound in intent and carries two defects that would bite during execution, plus five smaller ones.

1. **Its Step 3 hides a refactor of a test-pinned function.** It requires `unparsed_rows()` to share `published_incumbents()`'s consumption predicates rather than restate them — correct, but those predicates are **stateful** (`cur` and `methods` mutate inside the loop and the numeric branch reads both), so sharing them means extracting a classifier and rewriting that loop. `Files:` and the test strategy acknowledge neither, and the equivalence proof that would make such a refactor safe is absent. As written the step is either unsatisfiable or a silent rewrite of a function pinned at `tests/test_lentimpra.py:182`.
2. **Its Step 7 contradicts itself in a way that breaks the execute-plan contract.** `Files:` declares read-only and the plan's Verification item 8 asserts read-only "by contract", while the step body authorises *"Correct only prose surfaces found."* If it finds one it edits a file absent from `Files:` — the line that drives wave grouping and overlap detection.

The five medium findings: a `superseded_claim` field that turns a runtime return value into a changelog; an unsourced "~0.01–0.03" interval inserted into the very memo arguing that unsourced values must be refused; a shared-global-path constraint stated only in Execution Preview rather than at the step that writes it; an index entry missing the `**Status:**` field that 38 other entries carry; and a docstring correction that could read as a weakened refusal.

**Success criteria.** The amended plan contains no step whose `Files:` disagrees with its body; Step 3's refactor is explicit and carries an equivalence proof; no correction history lives in a runtime return value; no unsourced quantitative claim is introduced; and the index entry declares a status.

**Non-goals.** Deciding the narrow-vs-full scope question (whether `unparsed_rows()` ships now or waits for the volume) — that is a user decision and is carried forward as an Open Question, with this amendment making **both** paths safe rather than picking one. Executing the amended plan. Touching `dna_decode/` or the arms ledger at all: every file in scope here is a planning artifact.

## Codebase Context

- `plans/Correct_Two_Wrong_Claims_I_Committed_And_Make_The_Next_Check_Runnable_Plan/technical-plan.md` — **20,206 bytes**, written this session. 7 steps; waves `[[1,4,5,7],[2,6],[3]]`; max parallelism 4; footer `check=clean waves=clean gate=fired:open-questions,unverified`.
- Its **Step 3** names `Files: dna_decode/glm/lentimpra.py, tests/test_lentimpra.py` and says "Derive the consumption rule from the same predicates `published_incumbents` uses — do not restate them as a second copy".
- Its **Step 7** names `Files: (read-only — soraya_runs/arms-ledger.json, dna_decode/eval/regime.py, wiki/)` and its body says "Correct only **prose** surfaces found".
- **The statefulness that grounds C1**, read from `dna_decode/glm/lentimpra.py` around line 320: `blocks, cur, methods = {}, None, None` then, inside the row loop, the `CELL_TYPES` branch assigns `cur` and resets `methods`, the `test fold` branch assigns `methods`, and the numeric branch is guarded by `if cur and methods and isinstance(cells[0], (int, float))`. The row classification is therefore not a pure function of the row.
- **The toolkit's step-heading contract** requires `### Step N:` with a numeric N — so splitting Step 7 cannot be done as `7a`/`7b`; it renumbers the plan to **8 steps**.
- `wiki/plans-index.md` — **109,801 bytes**, pure LF, canonical `/save-plan` header. **38 entries carry `**Status:**`**; observed values include `executed — …` and `active (…)`. The entry for this plan family carries none.
- `tests/test_lentimpra.py:182` — `assert inc["metric"] == "UNSPECIFIED_IN_SOURCE"`, the pin that makes Step 3's refactor risky.

### Reusable-Code Survey

- **`plan-toolkit.py check` / `waves` / `validate`** — the deterministic validator; Step 3 re-runs it rather than hand-computing the amended wave table, because a hand-written wave grouping is a second expression of the dependency graph.
- **The prior plan's own Step 3 pattern** (`Make A Verdict Reachable…`, this session) — it paired a refactor with a byte-identical before/after proof plus a source-shape assertion. That is the pattern C1 says is missing here; reused rather than reinvented.
- **The 38 existing `**Status:**` entries in `wiki/plans-index.md`** — the format to match for M4, rather than inventing a vocabulary.
- Searched: `graphify-out/GRAPH_REPORT.md` (absent), `plans/`, `wiki/`, `dna_decode/glm/`, `tests/`, `~/.claude/skills/technical-plan/scripts/`.

## Pre-Change Baseline

Measured in-session, 2026-10-28:

- Target plan = **20,206 bytes**, 7 steps, waves `[[1,4,5,7],[2,6],[3]]`, max parallelism **4**, critical path `1→2→3`, footer `gate=fired:open-questions,unverified`.
- Step 7's `Files:` line declares read-only; its body authorises prose edits — **the contradiction is present in the saved file**, not hypothetical.
- Step 3 contains **no** equivalence test and **no** mention of extracting a classifier; `grep` for `equivalence`/`byte-identical` in the saved plan: no hits.
- `wiki/plans-index.md` = **109,801 bytes**; `grep -c "^\*\*Status:\*\*"` = **38**; the entry for this plan has none.
- `dna_decode/glm/lentimpra.py` = 17,365 bytes, unchanged since `b28fa8dc`; `uv run pytest tests/test_lentimpra.py -q` → 4 passed / 18 skipped.
- Frozen five + `uv.lock`: byte-unchanged.

## Verification Signal

- **No step's `Files:` disagrees with its body** — checked by reading each amended step's `Files:` line against what its bullets authorise. Step 7 either stops at reporting, or the file it edits is declared.
- **Step 3 is executable as written**: it names the classifier extraction, lists the function it rewrites, and carries an equivalence check — `published_incumbents()` output identical before and after against a recorded fixture, plus a source-shape assertion that the predicates exist in exactly one place.
- **`grep` over the amended plan returns no hit for `superseded_claim`** and no unsourced numeric interval (`0.01`, `0.03`) in the Step 4 bullets.
- **The toolkit re-validates after amendment**: `check` reports no `missing`/`errors`, `waves` reports no `errors`, and the recomputed wave table is written from the toolkit's output rather than by hand.
- **The index entry declares `**Status:** candidate`**, matching the format of the 38 existing entries.
- The amended plan's footer is the toolkit's own output for the amended body — not the old footer carried over.

## Implementation Steps

### Step 1: Close the two execution-blocking defects
Files: plans/Correct_Two_Wrong_Claims_I_Committed_And_Make_The_Next_Check_Runnable_Plan/technical-plan.md
Depends on: none

**What changes:**
- **Step 3 (C1):** state that sharing the consumption rule requires extracting the row classification out of `published_incumbents()`'s loop into a helper both callers use, and that this REWRITES that loop. Add the missing proof to its test strategy: `published_incumbents()` output byte-identical before and after against a recorded fixture, plus a source-shape assertion that the predicates appear in exactly one place (value equality cannot see a reintroduced duplicate). Note that the state the classifier needs (`cur`, `methods`) must be passed in, not re-derived.
- **Step 7 (C2):** split into two numbered steps — a read-only reporting step, and a correcting step whose `Files:` declares the specific prose files it may edit. Because the toolkit requires a numeric `### Step N:`, this renumbers the plan to **8 steps**; the ledger-asserts-it case stays a surface-and-stop, so the arms ledger never appears in any `Files:` line.
- Leave every other step's text alone — this step fixes exactly the two findings that block execution.

**Test strategy:**
- Read each touched step and confirm its `Files:` line covers every file its bullets authorise writing.
- `grep` the amended plan for `read-only` and confirm no step claims it while authorising an edit.
- Defer wave/footer recomputation to Step 3 — hand-editing the wave table here would be a second expression of the dependency graph.

### Step 2: Apply the five medium corrections
Files: plans/Correct_Two_Wrong_Claims_I_Committed_And_Make_The_Next_Check_Runnable_Plan/technical-plan.md
Depends on: Step 1

**What changes:**
- **M1:** drop `superseded_claim` from the planned return-value additions. The correction history goes to the memo and `CHANGELOG.md`; the return value keeps **one** pointer to the memo. Record the reason in the step so the precedent is explicit: a runtime payload that accumulates a field per correction grows without bound.
- **M2:** remove the "~0.01–0.03" interval from Step 4's bullets and state the point qualitatively (the two coefficients are often close on standardized, near-linear data, so a two-decimal match discriminates poorly). An unsourced number inside a provenance-discipline memo undercuts the memo.
- **M3:** move the shared-global-path constraint into Step 5's own text — `~/.claude/soraya-resume/` is outside the repo, gets no overlap protection, and a worktree gives zero isolation there, so that step must not run concurrently with another execution.
- **M5:** require the docstring edit to keep the readout-type-mismatch rationale adjacent and primary, so the scoping sentence cannot be read as conceding the refusal was premature.
- M4 is deliberately **not** here — it edits a different file and lands in Step 3.

**Test strategy:**
- `grep` for `superseded_claim`, `0.01`, `0.03` in the amended plan — expect no hits.
- Confirm Step 5's body now carries the no-concurrency constraint and that Execution Preview still agrees with it.

### Step 3: Re-validate the amended plan and declare its index status
Files: plans/Correct_Two_Wrong_Claims_I_Committed_And_Make_The_Next_Check_Runnable_Plan/technical-plan.md, wiki/plans-index.md
Depends on: Step 2

**What changes:**
- Re-run `plan-toolkit.py check` and `waves` on the amended body; fix any `missing`/`errors` before proceeding. Write the recomputed wave table, total waves, max parallelism and critical path from the toolkit's output into Execution Preview — the step count changed from 7 to 8, so the old table is stale by construction.
- Re-run `validate` and replace the trailing footer with its output, preserving the footer-is-final-line invariant.
- **M4:** add `**Status:** candidate` to the index entry, matching the 38 existing entries' format, and note in the entry that the plan was amended after a pre-exec review.
- Record the save-order inversion in one line: `/save-plan` ran before the pre-exec review, so this amendment is an in-place edit rather than a second save — which is what the save-once rule prescribes.

**Test strategy:**
- `tail -c 500` the amended plan and confirm exactly one `<!-- toolkit:` line, as the final line.
- Assert the plan's stated wave table equals the toolkit's `waves` output verbatim.
- `grep -c "^\*\*Status:\*\*"` on the index returns **39**, and `wiki/plans-index.md` stays pure LF.

## Execution Preview

Waves (topological, from the `Depends on` graph):

- **Wave 0:** Step 1
- **Wave 1:** Step 2
- **Wave 2:** Step 3

Total waves: 3. Max parallelism: 1. Critical path: Step 1 → 2 → 3.

Strictly sequential by construction: all three steps edit the same plan file, and Step 3's revalidation is only meaningful once Steps 1 and 2 have landed. No parallelism is available and none should be manufactured — a worktree here would isolate nothing and would cost a `uv.lock` reconciliation.

## Risk Flags

- **Amending a saved plan edits an artifact `/execute-plan` may already be pointed at.** If execution starts against the pre-amendment file, it executes the two defects. The amendment should land before any `/execute-plan` on that path.
- **The findings driving this plan came from a single-voice review** (`codex` absent on this host, verified). C1 and C2 rest on reading the parser's loop and the saved plan's own `Files:` lines — both quoted in Codebase Context so a reader can check them rather than take them on trust.
- **Step 1 renumbers the plan to 8 steps.** Any external reference to "Step 7" of that plan — including the index entry's bullets and this conversation — becomes ambiguous. Step 3's index note is the mitigation; nothing else in the repo references those step numbers (`grep` for the plan's folder name outside `plans/` and `wiki/plans-index.md`: no hits).
- **Recursion risk, named.** This is a plan amending a plan. It is justified only because the two defects are execution-blocking and the edit is two files; a further round of planning on top of this one would be motion, not signal. If the next review of this amendment finds only wording issues, apply them directly rather than planning again.
- **M1 removes a field the earlier review asked for.** The prior plan added `superseded_claim` to make the correction legible; dropping it trades call-site visibility for a bounded payload. The memo must then actually carry the superseded wording, or the correction becomes invisible — that obligation moves to Step 2's text.
- Architecture gate: `sentrux` present, but this plan touches only planning artifacts — no source, no module boundary. A structural compare would be vacuous.

## Open Questions

- **Narrow or full scope for the amended plan?** The narrow version is Steps 1, 2, 4, 5, 7 of the original — pure wording corrections, zero refactor, all executable today. The full version keeps the `unparsed_rows()` work, which C1 shows means refactoring a test-pinned function for a check that cannot run until the volume is attached and may then find nothing. This amendment deliberately makes **both** safe rather than choosing; the choice is yours and it changes which steps exist.
- **Should the equivalence fixture be committed or generated at execution time?** A committed fixture is a stable reference across sessions; a generated one cannot drift out of sync with the code. The prior plan's analogous step left this implicit.
- **Does `**Status:** candidate` or `active` apply once `/execute-plan` starts?** The convention lists both; nothing in the index says who flips it or when, and 38 entries do not settle it by example.
- Whether to record the save-order inversion as a one-off or as a process note worth carrying into the pipeline rules — this is the second time this session that a plan was saved before its mandatory pre-exec review.

## Verification

1. Each amended step's `Files:` line is read against its own bullets; no step claims read-only while authorising an edit.
2. Step 3 of the amended plan names the classifier extraction, the function it rewrites, and an equivalence proof; `grep` for `equivalence` in the amended plan returns a hit where previously it returned none.
3. `grep` for `superseded_claim`, `0.01`, `0.03` in the amended plan returns no hits.
4. `plan-toolkit.py check` reports no `missing`/`errors`; `waves` reports no `errors`; the written wave table equals the toolkit's output verbatim.
5. `tail -c 500` shows exactly one `<!-- toolkit:` line and it is the final line.
6. `grep -c "^\*\*Status:\*\*" wiki/plans-index.md` returns 39; the file stays pure LF (no CRLF introduced).
7. `dna_decode/`, `tests/`, `soraya_runs/arms-ledger.json` and the frozen five are byte-unchanged — this plan touches planning artifacts only, and `git diff --stat` on those paths must be empty.

<!-- toolkit: check=clean waves=clean gate=fired:open-questions -->
