# Close The Verdict-Vocabulary Defect The Review Found And Scope The Positive Finding

## Lens status

```
Repo-index:        absent (no .repo-index/, no repo-index-status.sh) — source-only reading
sentrux:           n/a — not installed
codex:             n/a — not installed (the review that produced this plan ran Claude-only)
gh:                present (/c/Program Files/GitHub CLI/gh) — unused, nothing to PR
graphify:          absent — no god-node reuse signal
DESIGN.md:         absent — plan touches no UI
project-rules.md:  absent — Project Rules Check section omitted
Baseline tests:    77 pass / 0 fail across the two affected files, measured at plan time
Prediction tables: present, 10 seeds + manifest — every step re-derives, none recomputes
```

## Problem Statement

A post-ship review of the stacking measurement (uncommitted, 14 modified + 3 untracked paths) found
one defect that matters and two claims that are scoped wrongly. The review ran without an independent
model, so every finding below was verified by executing against the code rather than by reasoning
about it.

**The defect.** `stacking_verdict_for` returns `WEAK_STACK` for any statistic with `median > 0`. No
seed-count condition, no CI condition. Verified by direct call:

```
POSITIVE NOISE 4/10, CI spans zero     -> WEAK_STACK
REAL 10/10, CI excludes zero           -> WEAK_STACK
```

So the machine-readable verdict cannot separate a gain that reproduced on every seed from one that
landed positive by chance. That is the same collapse the file deliberately prevents one function away
— `INDETERMINATE_STACKING_NOT_MEASURED` exists precisely so a missing measurement cannot read as a
negative one — and it is load-bearing for the shipped result, whose entire argument is that +0.0145 at
10/10 with a CI excluding zero is a real effect that merely fell under a margin. The prose detail
string says so; the verdict field cannot.

**And the defect is in the PRE-REGISTRATION, not only in the code.**
`PREREGISTERED["stacking"]["verdicts"]["WEAK_STACK"]` reads *"gain > 0 but below the margin or the
seed count"* — the code implements the frozen definition faithfully. This is what shapes the whole
plan: the fix edits a pre-registered bar after the result is known, which is the exact retro-fit the
original plan existed to prevent. The only thing that makes it legitimate is a **neutrality gate**:
the change must be provably outcome-neutral on the data already in hand, and it must introduce **no
new threshold**. Both hold, and Step 5 proves the first rather than asserting it.

**Two scoping errors.**
- The aggregate `stacking_verdict_detail` names only the combos that cleared. Read from the shipped
  artifact: `stacking_per_combo` records `hybrid_conv_maxtx_gc -> NO_STACK`, and the detail string
  names only the two `WEAK_STACK` combos. The markdown per-combo table does include the negative, so
  the human surface is honest and the JSON field a programmatic reader consumes is not. The artifact's
  own `honest_limits` item 5 instructs the reader to "read the per-combo block, not only the
  aggregate" — and the detail string is where that instruction is violated.
- The published finding "the orthogonality argument got the direction backwards" is true but
  over-general. `scripts/glm_tile_predictions.py:357-358` passes the SAME `dec` dict to both
  `tile_features()` and `max_tx_scalar()`, so `max_tx` is a lossy scalar reduction of exactly what the
  `decomposition` channel sees in full. The +0.7812 Spearman is therefore mechanically expected, and
  the counterexample is a **reduction/superset pair, not two modalities**. The memo states this
  relationship — but only inside the quoted `negative_scope`, which scopes the gc NEGATIVE. The
  positive headline carries no caveat, and the headline is what gets quoted.

**One finding is deliberately recorded and NOT fixed.** `promcalc_decomposition_arm.paired` builds a
normal-approximation interval on `fmean(d)` with `z=1.96` while reporting `median(d)` as the point
estimate — an estimand mismatch, and `t(9)=2.262` is the correct critical value at n=10, so published
intervals are ~15% narrower than their nominal level. I re-derived all three under t: decomposition
`[0.0103, 0.0177]` and gc+decomposition `[0.0030, 0.0165]` still exclude zero, gc still includes it;
sign tests corroborate independently (p=0.00098 / 0.0107 / 0.828). Nothing published is wrong. But the
function is **shared with the published parent arm**, so correcting it moves parent numbers and
belongs to its own pre-registered change.

### Success criteria

1. A reproducible gain and positive noise receive DIFFERENT verdict strings.
2. No verdict in the published artifact changes — proven by diff, not asserted.
3. No new numeric threshold is introduced; the narrowed condition reuses the frozen
   `min_seeds_winning`.
4. Every combo appears in the aggregate detail string, with its own verdict.
5. The reduction/superset caveat sits beside the positive headline, not only beside the negative.
6. The CI estimand finding is recorded with its re-derivation and left unfixed.
7. Every published 2-channel value stays byte-identical to `git show HEAD:`; top-level verdict stays
   `HYBRID_ADDS`.

### Non-goals

- No change to `promcalc_decomposition_arm.paired` (see above).
- No change to `margin` (0.02) or `min_seeds_winning` (8). The margin question stays an open
  pre-registration for a different question, exactly as the memo states.
- No change to the rung, the split, the combine rule, the 0.941 ceiling, or any prediction table.
- No history rewrite. The pushed commit `e96b5591` stays as it is; corrections NAME it.
- No new `planned` arm. Open Question 1 carries that decision; this plan does not take it.

## Codebase Context


See: ./codebase-context.md

## Pre-Change Baseline

Measured at plan time, not quoted from the prior run:

| quantity | value | how obtained |
|---|---|---|
| tests, two affected files | **77 passed / 0 failed** (63 + 14) | `uv run pytest tests/test_glm_tile_hybrid_arm.py tests/test_tile_chain_substrate_wall_amendment5.py -q` |
| `stacking_verdict` | `WEAK_STACK` | read from `wiki/glm_tile_hybrid_arm_2026-10-26.json` |
| `stacking_per_combo` | `hybrid_conv_maxtx_gc` → `NO_STACK`; `hybrid_conv_maxtx_decomp` → `WEAK_STACK`; `hybrid_conv_maxtx_gc_decomp` → `WEAK_STACK` | same artifact |
| `stacking_verdict_detail` | names 2 of 3 combos (the two `WEAK_STACK`) | same artifact |
| `honest_limits` | 6 items | same artifact |
| `verdicts` map | 4 keys | `PREREGISTERED["stacking"]["verdicts"]` |
| verdict on noise `_p(0.0001, 4)` | `WEAK_STACK` — the defect | direct call to `stacking_verdict_for` |
| top-level `verdict` | `HYBRID_ADDS` | same artifact |
| published 2-channel paired median | `+0.0339`, CI95 `[0.0238, 0.0409]` | same artifact; byte-identical to `git show HEAD:` |
| git state | `main` at `e96b5591`, 14 modified + 3 untracked, nothing committed | `git status --short --branch` |

## Verification Signal

What would show this change worked, and what would show it failed:

1. **Vocabulary separates the two cases.** `stacking_verdict_for(_p(0.0001, 4))` returns the new
   `INDETERMINATE_*` key while `stacking_verdict_for(_p(0.0145, 10))` returns `WEAK_STACK`. Failure
   signal: both still return the same string.
2. **NEUTRALITY GATE — the decisive one.** All three entries in `stacking_per_combo` carry the same
   verdict after the re-run as before (`NO_STACK` / `WEAK_STACK` / `WEAK_STACK`), and
   `stacking_verdict` is still `WEAK_STACK`. Failure signal: any verdict moves — which would mean a
   frozen bar was relaxed after the fact, and the change must be reverted rather than argued for.
3. **No new threshold.** `git diff` on `scripts/glm_tile_hybrid_arm.py` introduces no numeric literal
   in the stacking bar. `margin == 0.02` and `min_seeds_winning == 8` unchanged, still pinned by
   `:249-250` and `:253-254`.
4. **Regression gate.** Every published 2-channel value in the regenerated artifact is byte-identical
   to `git show HEAD:wiki/glm_tile_hybrid_arm_2026-10-26.json`, and top-level `verdict` is
   `HYBRID_ADDS`. Failure signal: any drift — the stacking numbers stop being trustworthy and nothing
   downstream should be believed.
5. **Every combo named.** `stacking_verdict_detail` contains all three combo names and all three
   verdict strings. Failure signal: a reader of the JSON alone still cannot see the negative.
6. **Caveat placed.** The memo states the reduction/superset relationship above line 48 (before the
   `WHY THE VERDICT IS STILL WEAK_STACK` section), not only in the line-114 quote.
7. **Test delta is accounted for.** 77 → expected 83-85 in the two files, with exactly one assertion
   CHANGED (`:399`) and the rest added. Failure signal: any assertion silently deleted.
8. **Frozen five + `uv.lock` byte-unchanged** — `git status --porcelain` on exactly those six paths
   returns empty.

## Implementation Steps

### Step 1: Add a fourth stacking verdict and narrow WEAK_STACK to reproducible gains
Files: scripts/glm_tile_hybrid_arm.py
Depends on: none

**What changes:**
- Add `INDETERMINATE_GAIN_NOT_REPRODUCIBLE` to `PREREGISTERED["stacking"]["verdicts"]`, described as
  *a positive paired gain that did not reproduce on enough seeds — neither a measured stack nor a
  measured negative*. Joins the existing `INDETERMINATE_*` family.
- Narrow the `WEAK_STACK` entry's text to *gain > 0 AND >= min_seeds_winning seeds win, but below the
  margin*, and record IN THAT ENTRY that the condition was tightened after the run, that it reuses the
  already-frozen `min_seeds_winning` rather than introducing a threshold, and that Step 5 proves no
  measured verdict moved. A silent edit to a pre-registered bar is the failure being avoided; a
  disclosed, outcome-neutral one is not.
- In `stacking_verdict_for`, insert the reproducibility condition: `median > 0 and wins >=
  s["min_seeds_winning"]` → `WEAK_STACK`; `median > 0` (seeds short) → the new key, with a detail
  string that says the gain is positive but did not reproduce and must not be read as either a stack
  or a negative. `NO_STACK` remains the `median <= 0` fallback with its text unchanged.
- Extend the aggregate precedence tuple in `stacking_verdict` from `("STACKS", "WEAK_STACK")` to
  `("STACKS", "WEAK_STACK", "INDETERMINATE_GAIN_NOT_REPRODUCIBLE")`, ordered. Without this, an arm in
  which every combo is unreproducible-positive would fall through to the `NO_STACK` return — reporting
  a measured negative for a set of positive gains, which is the original defect re-created at the
  aggregate level.
- Do not touch `margin`, `min_seeds_winning`, `comparator`, `mirrors`, `predicted_outcome`,
  `negative_scope`, or any top-level `PREREGISTERED` key.

**Test strategy:**
- Unit: `stacking_verdict_for(_p(0.0001, 4))` → the new key; `_p(0.0145, 10)` → `WEAK_STACK`;
  `_p(margin, min_seeds_winning - 1)` → the new key (the moved boundary); `_p(margin/2,
  min_seeds_winning)` → `WEAK_STACK`; `_p(-0.01, 2)` → `NO_STACK`.
- Unit: the four pre-existing verdict keys are all still present and their STRINGS for `STACKS` /
  `NO_STACK` / `INDETERMINATE_STACKING_NOT_MEASURED` are unchanged, so only `WEAK_STACK` moved.
- Unit: `margin == 0.02` and `min_seeds_winning == 8` still assert, and the new condition references
  `min_seeds_winning` rather than a literal (grep the source for a second numeric literal in the
  stacking branch and assert there is none).
- Unit: an all-unreproducible-positive arm aggregates to the new key, NOT to `NO_STACK` — the
  precedence bug, pinned.

### Step 2: Make the aggregate detail string name every combo, not only the winners
Files: scripts/glm_tile_hybrid_arm.py
Depends on: Step 1

**What changes:**
- In `stacking_verdict`, replace the hit-path join over `hits` with a join over `sorted(per)`, each
  entry annotated with its own verdict, so the returned detail reads e.g.
  `hybrid_conv_maxtx_decomp [WEAK_STACK]: ...; hybrid_conv_maxtx_gc [NO_STACK]: ...`. The aggregate
  verdict itself is unchanged — still best-of — but the string stops presenting best-of as if it were
  all-of.
- Reuse `sorted(per)`, already used by the `NO_STACK` fallback four lines below, so both paths iterate
  identically.
- Leave `render_md` alone: its per-combo table already lists losers, verified.

**Test strategy:**
- Unit: an arm with one `WEAK_STACK` and two `NO_STACK` combos returns `WEAK_STACK` and a detail string
  containing ALL THREE combo names and all distinct verdict strings present.
- Unit, non-vacuity: assert the detail names a combo whose verdict differs from the aggregate — this
  is the assertion the pre-change code fails. Confirm it fails against the pre-change function before
  accepting the fix, so the test is proven to bite.
- Unit: the aggregate verdict returned is byte-identical before and after for the three published
  statistics (the string changed; the verdict did not).

### Step 3: Scope the positive finding and record the CI estimand finding in honest_limits
Files: scripts/glm_tile_hybrid_arm.py
Depends on: Step 2

**What changes:**
- Append a 7th `honest_limits` item scoping the direction claim: `max_tx` and `decomposition` are two
  readings of the SAME Promoter Calculator output — `scripts/glm_tile_predictions.py:357-358` passes
  one `dec` dict to both `tile_features()` and `max_tx_scalar()` — so the pair whose redundancy is
  highest and which nonetheless stacks is a **reduction/superset pair**. Pairwise prediction-vector
  Spearman is refuted as a screen ACROSS THIS PAIR; it is not shown to be invalid between genuinely
  different modality families, which remain untested.
- Append an 8th item recording the CI finding: the interval is a normal approximation on the MEAN with
  `z=1.96` while the reported point estimate is the MEDIAN, and `t(9)=2.262` is the correct critical
  value at n=10, so published intervals are narrower than nominal. State that the conclusion was
  re-derived under `t` and survives on all three combos (`[0.0103, 0.0177]`, `[0.0030, 0.0165]`, gc
  still including zero) with sign tests corroborating, and that the function
  (`promcalc_decomposition_arm.paired`) is SHARED with the published parent arm, so correcting it
  moves parent numbers and is deliberately out of scope here.
- Do not modify items 1-6.

**Test strategy:**
- Unit: `honest_limits` has 8 items; items 1-6 are byte-identical to their current text.
- Unit: item 7 names both `max_tx` and `decomposition` and the word `reduction`; item 8 names both
  `mean` and `median` and the `t` critical value — so neither can be satisfied by vague prose.
- Unit against the SHIPPED artifact (after Step 4): both new items are present in the published
  `honest_limits` and both reach the rendered markdown via `render_md`'s limits section.

### Step 4: Re-run the arm and prove the change moved no measured verdict
Files: wiki/glm_tile_hybrid_arm_2026-10-26.json, wiki/glm_tile_hybrid_arm_2026-10-26.md
Depends on: Step 3

**What changes:**
- Capture the pre-change artifact to the scratchpad first, so the neutrality comparison is against a
  file that exists rather than against a remembered value.
- Re-run `uv run python scripts/glm_tile_hybrid_arm.py --seeds 10` against the existing 10 prediction
  tables. No recompute: the tables and `manifest.json` are present and untouched.
- Preserve the artifact's CRLF line endings and `indent=1` — `indent=2` would produce a whole-file
  reformat diff that buries the real change, and these two files are CRLF while the source is LF.
- **NEUTRALITY GATE, blocking:** diff the regenerated artifact against the captured pre-change copy and
  assert the changed field set is exactly `{stacking_verdict_detail, honest_limits, preregistered}` —
  plus nothing else. If `stacking_verdict`, any `stacking_per_combo` verdict, or any `arms` value
  moved, STOP and revert Steps 1-3 rather than argue the result forward.
- **REGRESSION GATE, blocking:** compare every published 2-channel value against
  `git show HEAD:wiki/glm_tile_hybrid_arm_2026-10-26.json` — `hybrid_conv_maxtx`,
  `hybrid_conv_maxtx_gc`, `max_tx`, `conv`, `control_shuffled_conv`, the four legacy `_paired_*` keys,
  the whole `orthogonality` block, and top-level `verdict` — and require byte equality.
- Confirm `_short_arm` and `_skipped_combos` are still both absent, as on the original run.

**Test strategy:**
- Integration: the two blocking gates above, run as explicit comparisons whose output is recorded, not
  as assertions inside a test that could pass vacuously if a key were renamed.
- Unit: re-assert `uv run pytest tests/test_glm_tile_hybrid_arm.py -q` green against the regenerated
  artifact (the `:807-829` tests read it directly).
- Manual: `git diff --numstat` on both artifacts should show a small additive change, not a reformat.

### Step 5: Update and extend the test suite for the moved boundary
Files: tests/test_glm_tile_hybrid_arm.py
Depends on: Step 4

**What changes:**
- Rewrite the one assertion the change breaks, `:399`
  (`_p(margin, min_seeds_winning - 1)` → `WEAK_STACK`), to expect
  `INDETERMINATE_GAIN_NOT_REPRODUCIBLE`, and comment WHY in place: a gain that clears the effect-size
  bar on 7 of 10 seeds is positive-but-unreproducible, and the old expectation is precisely the
  conflation this change removes. Change it; do not delete it.
- Add the Step-1, Step-2 and Step-3 unit tests named above, including the two non-vacuity checks (the
  detail-string test must be shown to fail against the pre-change function; the no-new-literal check
  must be shown to fail if a literal is substituted for `min_seeds_winning`).
- Add one test asserting the three verdict strings OTHER than `WEAK_STACK` are byte-unchanged, so a
  future edit cannot quietly widen a second branch.
- Do not touch `tests/test_tile_chain_substrate_wall_amendment5.py` (CRLF, unrelated to this defect).

**Test strategy:**
- `uv run pytest tests/test_glm_tile_hybrid_arm.py tests/test_tile_chain_substrate_wall_amendment5.py -q`
  — expect 83-85 passed / 0 failed against a 77 baseline, with exactly one assertion changed.
- Each added test proven non-vacuous by mutation: revert the specific line it guards, confirm the test
  fails, restore.
- Full suite deferred to Step 8's verification; never run while a torch job holds a model.

### Step 6: Move the reduction/superset caveat up beside the positive headline
Files: wiki/glm_tile_stack_arm_2026-10-26.md
Depends on: Step 4

**What changes:**
- Add a short scoping paragraph inside `## THE PREDICTION WAS FALSIFIED, AND THAT IS THE RESULT`
  (currently line 9), immediately after the direction claim, stating that `max_tx` and `decomposition`
  are two readings of one Promoter Calculator output with `max_tx` a scalar reduction of it — citing
  `scripts/glm_tile_predictions.py:357-358` — so what is refuted is prediction-vector Spearman as a
  screen across a reduction/superset pair, NOT across genuinely different modality families.
- Amend the document headline so it does not read as a general result about orthogonality screening.
- Add the CI estimand note to the memo's `## Honest limits`, with the re-derived `t` intervals and the
  sign-test p-values, and state plainly that the conclusion survives and the shared function is out of
  scope.
- Add a note that the verdict vocabulary was tightened after the run, that the tightening was
  outcome-neutral on these three statistics, and point at Step 4's recorded neutrality gate as the
  evidence rather than asserting it.
- Preserve the file's existing line endings; write with `open(..., newline="")` and verify byte counts.
- Do NOT restate any measured number from memory — read each from the regenerated artifact, as the
  memo's original build script did.

**Test strategy:**
- Manual: the caveat appears before line 48 (`## WHY THE VERDICT IS STILL WEAK_STACK`).
- Manual: every number in the new text matches the regenerated artifact, checked by reading the
  artifact rather than the memo.
- Automated: a grep assertion that the memo's headline no longer makes the unscoped direction claim,
  and that the new paragraph names `max_tx`, `decomposition` and `reduction`.

### Step 7: Record the residuals in TODOS and the umbrella ledger
Files: TODOS.md, project_state/glm-umbrella-2026-10-07.md
Depends on: Step 4

**What changes:**
- Amend the existing `WEAK_STACK is an AGGREGATE over three superset combos` TODOS residual — it
  already instructs the reader not to trust the aggregate alone, and the detail-string fix partly
  discharges it; say which part remains (one combo per bar, or a stated correction before the run).
- Add a TODOS residual for the CI estimand finding: the mean/median mismatch, `z` vs `t(9)`, the
  re-derivation that shows the conclusion survives, and the fact that
  `promcalc_decomposition_arm.paired` is shared with the published parent arm so a fix moves parent
  numbers. Name it as its own pre-registered change, not a cleanup.
- Add a TODOS residual recording that the `WEAK_STACK` condition was tightened AFTER the result was
  known, that it was outcome-neutral on the three measured statistics, and that this is the one case
  where editing a frozen definition was accepted — with the neutrality gate as the reason, so the
  precedent is bounded rather than open.
- Append one `## Action Log` row (class `propose`) to the umbrella naming the defect the review found,
  the neutrality gate, and the scoping correction. Verify the row count increments by exactly 1, the
  row appears once, and the `## Action Log` / `## Last Evaluation` markers each appear once. LF.

**Test strategy:**
- Automated: umbrella row count delta is exactly 1; new row present exactly once; both section markers
  still unique; no CRLF introduced.
- Manual: the three TODOS entries each name a file or function a reader can open.

### Step 8: Verify the whole change and commit
Files: CHANGELOG.md, LESSONS_LEARNED.md
Depends on: Step 5, Step 6, Step 7

**What changes:**
- Run the full suite and compare against the recorded 135-pass baseline from the prior run plus this
  plan's additions. Report the real numbers; do not cite the earlier run as verification.
- Confirm the frozen five and `uv.lock` are byte-unchanged via `git status --porcelain` on exactly
  those six paths.
- Add a CHANGELOG entry under the existing `## [Unreleased]` and one LESSONS_LEARNED bullet: a verdict
  vocabulary that cannot separate a reproducible gain from positive noise is the same defect class as
  one that cannot separate a missing measurement from a negative one, and the second was guarded while
  the first was not — in the same function.
- Commit the whole working tree as one commit naming the review that found the defect. **This step
  needs explicit user authorization** (see Open Question 3); if withheld, every other step still
  stands and the tree stays uncommitted.

**Test strategy:**
- `uv run pytest tests/ -q` — full suite, real output quoted.
- `git status --short --branch` before and after, so what entered the commit is auditable.
- `git show --stat HEAD` after committing, to confirm the file set matches the plan.

## Execution Preview

| wave | steps | parallelism |
|---|---|---|
| 0 | Step 1 | 1 |
| 1 | Step 2 | 1 |
| 2 | Step 3 | 1 |
| 3 | Step 4 | 1 |
| 4 | Step 5, Step 6, Step 7 | **3** |
| 5 | Step 8 | 1 |

Total waves: **6**. Max parallelism: **3**. Critical path:
Step 1 → 2 → 3 → 4 → (5 | 6 | 7) → 8 — six of eight steps.

**Run this SEQUENTIALLY, not in parallel worktrees.** Steps 1-3 all edit one file and are serial by
construction. Wave 4's three steps touch disjoint files and could parallelise, but
`data/promcalc/tile_predictions/` is gitignored, so a worktree cannot see the prediction tables and its
baseline for the Step-4 artifact would be unreconcilable. The whole chain also costs minutes, not
hours — the re-run reads 10 cached tables.

## Risk Flags

1. **Editing a pre-registered bar after the result is known is the thing this project exists to
   prevent, and this plan does it.** The mitigation is not an argument, it is Step 4's blocking
   neutrality gate plus the reuse of `min_seeds_winning` instead of a new number. If the gate shows any
   measured verdict moved, the correct response is to revert Steps 1-3 — not to justify the move. This
   flag should survive into the commit message.
2. **The narrowed condition keys on seed count, not on the CI, and that is a compromise.** "CI excludes
   zero" is the more direct statement of reproducibility, but the test helper `_p` synthesises `ci95`
   as `[median-0.01, median+0.01]`, so a CI-keyed condition would be untestable through the existing
   fixtures and would introduce a second judgement call. Seed count reuses a frozen number. Recorded as
   a known approximation, not as the ideal.
3. **`_short_arm` remains `[inferred]` in production** and this plan does not change that — both guards
   stayed silent on the real run. Pre-existing, already in TODOS, non-blocking.
4. **The Step-6 memo edit is the step most likely to introduce a wrong number**, because it restates
   measured values in prose. Mitigated by reading each from the regenerated artifact. This is the same
   class as the defect the prior run hit, where a write-up step ran before the measurement it
   characterised — which is why Step 6 depends on Step 4 rather than carrying `Depends on: none`.
5. **EOL hazard is live on this host.** `Path.write_text` converts `\n` to `\r\n`. The source and tests
   are LF; the two published artifacts are CRLF. Every write must use `open(..., newline="")` with a
   byte-count check. A silent EOL flip is the class that broke the TB barcode checksum in this repo.
6. **Nothing has been committed across this entire run** — 14 modified + 3 untracked paths, including
   two untracked wiki memos that exist nowhere else. A crash before Step 8 loses them. This argues for
   running the chain soon rather than for widening its scope.
7. **No independent review of THIS plan is available** (codex absent). The plan it implements came from
   a Claude-only review, so a single model found the defect and is now proposing the fix. The
   neutrality gate is what substitutes for a second opinion; it is a weaker substitute.
8. **No unverified external surface.** Every external surface named (`uv run python`, `uv run pytest`,
   `git`, and the arm's own `--in-dir` / `--seeds` / `--out`) was read from the repo's own
   `add_argument` calls or is already in routine use here. No flag or output format is inherited from
   documentation. One host note worth carrying: `python3` on this machine resolves to the Windows Store
   alias and cannot run the plan toolkit — `uv run python` is the working interpreter, and this plan was
   validated through it.

## Open Questions

1. **Should the untested new-modality-family direction become a `planned` arm?** It is currently a
   field on a resolved arm, so it never appears in `blocking()` output and is discoverable only by
   someone already reading that arm. Declaring an arm re-blocks the ceiling criterion on a scope
   question deliberately left open. A third option neither side has taken: a `planned` arm explicitly
   tagged scope-undecided, so it surfaces without claiming to be runnable. **Deliberately not folded
   into scope** — it changes what `blocking()` reports.
2. **Was a bar inherited from a different question ever the right bar here?** The inheritance rationale
   is "so the two results stay comparable", which presumes the two questions are the same question —
   and the measurement partly contradicts that (one asks whether a learned encoder beats a scalar, the
   other whether finer-grained access to an existing channel pays for itself). The memo's stated path
   is a new pre-registered bar for the new question. **Not folded in** — it is a new pre-registration,
   not an edit.
3. **Is committing authorized?** Step 8 commits the whole tree. The review recommended it and the tree
   has been uncommitted across the run, but no commit has been asked for in this session. If the answer
   is no, Steps 1-7 stand unchanged and the tree stays dirty.
4. **Should the new verdict key also be applied to `verdict()`'s parent branch?** The parent
   `WEAK_GAIN` branch has the identical shape — positive median, no seed condition — so the same
   conflation exists for the 2-channel question. **Not folded in**: that branch produced the PUBLISHED
   `HYBRID_ADDS` result, so touching it risks the regression gate for no measured benefit, and it
   deserves its own neutrality analysis.

## Verification

- [ ] `stacking_verdict_for(_p(0.0001, 4))` and `_p(0.0145, 10)` return DIFFERENT verdict strings.
- [ ] **Neutrality gate:** all three `stacking_per_combo` verdicts and `stacking_verdict` unchanged
      after the re-run; changed field set is exactly `{stacking_verdict_detail, honest_limits,
      preregistered}`.
- [ ] No numeric literal added to the stacking bar; `margin == 0.02` and `min_seeds_winning == 8`
      still pinned.
- [ ] An all-unreproducible-positive arm aggregates to the new key, not to `NO_STACK`.
- [ ] `stacking_verdict_detail` names all three combos with their own verdicts; the test for this is
      shown to fail against the pre-change function.
- [ ] **Regression gate:** every published 2-channel value byte-identical to `git show HEAD:`; top-level
      verdict `HYBRID_ADDS`; `_short_arm` and `_skipped_combos` both still absent.
- [ ] `honest_limits` has 8 items, items 1-6 byte-identical, both new items reaching the markdown.
- [ ] Memo carries the reduction/superset caveat above line 48 and a headline that does not overclaim;
      every number in it read from the regenerated artifact.
- [ ] Exactly one existing assertion changed (`tests/test_glm_tile_hybrid_arm.py:399`), none deleted;
      every added test proven non-vacuous by mutation.
- [ ] Two affected files: 77 → 83-85 passed / 0 failed. Full suite run and its real numbers reported.
- [ ] Frozen five + `uv.lock`: `git status --porcelain` on exactly those six paths returns empty.
- [ ] EOL preserved per file (source/tests LF, artifacts CRLF, umbrella LF); byte counts verified.
- [ ] Umbrella Action Log row count +1 exactly; row present once; both section markers unique.
- [ ] Three TODOS residuals added/amended, each naming an openable file or function.

<!-- toolkit: check=clean waves=clean gate=fired:open-questions -->
