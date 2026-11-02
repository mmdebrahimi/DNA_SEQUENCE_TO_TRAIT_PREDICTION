# Correct the Amendment Plan's Two Defects and Settle Its One Open Choice — Technical Plan

> Deliberately small: all three corrections are edits to one plan document, the chain is now two planning cycles deep with no code changed, and the only live defect in it stays red until that plan is executed.

## Lens status
Inputs: conversation (the review's recommended next step), grounded in `plans/Clear_The_Three_Blockers_The_Frame_Conditioning_Review_Found_Plan/technical-plan.md` and the preceding review of it.
Degradations: the repo-index helper script is absent on this host, so understanding comes from reading files directly; `project-rules.md` is absent, so the Project Rules Check section is omitted; `DESIGN.md` is absent and nothing here touches a user interface; the graphify dependency-graph report is absent, so the reuse survey was built by direct search; `codex` is absent; `sentrux` and `gh` resolve but no architecture gate applies, because this plan changes no source module.

## Problem Statement

A review of the amendment plan found two defects in it and left one choice open. All three are corrections to a saved plan document rather than to code.

**The first defect is a test that would fail on prose.** That plan's Step 1 specifies a test asserting the control-loader module "no longer contains the substring `gzip`". The module's docstring at line 12 contains the phrase `embeds genomewide.py as base64+gzip`, and that sentence survives the fix it describes, so the test fails while the code is correct. This is the same confusion between executable instructions and human commentary that the ownership guard's own text pattern suffers from — reproduced one step after diagnosing it there, which is why the correction is to change the *kind* of check rather than to reword the docstring around it.

**The second defect leaves a committed surface restating a claim the plan removes.** That plan's Step 4 rewrites the frame plan's Step 1 so it stops claiming exact reproduction of the upstream extraction function. But `wiki/plans-index.md` line 869 carries a key-decision bullet asserting that Step 1 "replicates the authors' `extract_promoter_regions` exactly, INCLUDING two apparent defects", and Step 4's file list does not include the index. After the amendment that bullet is false, on the plan whose stated purpose is correcting false claims.

**The open choice concerns where a measurement's output goes.** That plan's Step 3 extends an existing probe script and states that the existing artifact is left untouched. The script writes its artifact unconditionally and the record includes an elapsed-time field, so re-running it always produces a diff against a committed file that two other documents cite. The plan's intent and its mechanism disagree, and the choice is between giving the new comparison its own script or gating the original emission behind a flag.

**Scope, stated plainly because it bears on whether this plan should exist.** Every change here is an edit to one markdown file plus a record of what was accepted. Nothing is built, nothing is measured, and no test changes state. The chain that produced it is now two planning cycles deep: the exposure work shipped code and ran, the frame plan produced no code, the amendment plan produced no code, and this corrects the amendment plan. Meanwhile the one genuine code defect in the chain — a red ownership guard caused by a second reader of a data substrate — remains red, and the fix for it is a small module change plus the widening of seven test fixtures, currently sitting unexecuted in the plan this one corrects. The proportionate alternative was to apply these three edits through the save-time amendment mechanism rather than author another plan; that alternative is recorded in Open Questions rather than taken, because a plan was asked for.

**Non-goal.** This plan does not execute the amendment plan, does not touch the control loader, and does not clear the red guard.

## Codebase Context

**The pattern the corrected test should use already exists here, in eight or more test files.** `tests/test_amr_uncounted.py` parses a module with `ast.parse` at line 145 and then walks the tree collecting called names and imported names at lines 146 to 158, asserting on the resulting sets. The same approach appears in `tests/test_codon_consumers_behavioural.py`, `tests/test_codon_consumers_canonical.py`, `tests/test_codon_table_single_source.py`, `tests/test_data_inventory.py`, `tests/test_eval_transfer.py`, `tests/test_fba_solver_audit.py` and `tests/test_glm_finetune_kernel_builder.py`. So the correction adopts an established in-repo convention rather than introducing one, and the convention is specifically for asserting on the *shape* of source rather than on its text.

**Where the stale claim lives.** `wiki/plans-index.md` line 869 is a key-decision bullet under the frame plan's entry. The index is the file a reader scans to see what a plan decided without opening it, which is what makes a false bullet there worse than a false sentence buried in the body.

**Why the probe's artifact cannot be left untouched as written.** `scripts/frame_plan_feasibility_probe.py` defines its output path as a module constant at line 38, includes a `wall_seconds` field in the record at line 262, and writes the file unconditionally at line 264. There is no flag and no conditional. Any run of a modified version of that script rewrites the committed artifact, and the timing field guarantees the rewrite differs even when every input is identical.

**What cites that artifact.** Both `plans/Frame_Conditioning_Does_Corpus_Window_Position_Explain_PromoterAtlas_0_5076_Plan/technical-plan.md` and its sibling `brainstorm.md` reference it as the source of their pinned figures, which is why silently regenerating it is a correctness concern for those documents rather than only an untidy diff.

**Current state of the documents being edited.** The amendment plan is 215 lines and carries the toolkit footer `check=clean waves=clean gate=fired:open-questions,severity-high` as its final line. Its wave grouping is three waves at a maximum parallelism of two, with steps one and three in the first wave. Adding a path to one step's file list is an input to that computation, so the grouping has to be recomputed rather than assumed unchanged.

**Version-control state.** The working branch is `promoteratlas-exposure-2026-11-01`, three commits ahead of `main`, with the amendment plan folder and the modified index currently uncommitted on it. Nothing has been pushed.

### Reusable-Code Survey
- `tests/test_amr_uncounted.py` lines 145 to 158 provide the exact syntax-tree assertion pattern the corrected test should follow: parse the module, walk for import nodes, assert on the name set
- Seven further test files use the same `ast.parse` and `ast.walk` convention, so the pattern is established rather than novel
- `scripts/plan-toolkit.py` provides the `check`, `waves` and `validate` subcommands that re-validate an amended plan and emit its replacement footer
- `scripts/frame_plan_feasibility_probe.py` is the script whose output location is in question, and its existing structure is what makes the flag-versus-new-script choice concrete
- None additionally — searched: `graphify-out/GRAPH_REPORT.md` (absent), `src/`, `lib/`, `utils/`, `helpers/` (none exist at the repository root), `tests/`, `scripts/`

## Pre-Change Baseline

- The amendment plan is **215 lines**, its footer reads `check=clean waves=clean gate=fired:open-questions,severity-high`, and its wave grouping is **three waves, maximum parallelism two**, with the first wave holding steps one and three.
- Its Step 1 test strategy currently specifies a substring assertion on `gzip`. The target module contains that substring **three times**: once in prose at line 12, once as the import at line 27, once as the call at line 59. After the fix only the prose occurrence survives, so the specified test fails.
- `wiki/plans-index.md` line **869** asserts exact replication including two apparent defects. The amendment plan's Step 4 file list contains **one** path and does not include the index.
- `scripts/frame_plan_feasibility_probe.py` writes its artifact unconditionally at line **264**, with the output path fixed at line **38** and a timing field at line **262**.
- The ownership guard `tests/test_substrate_loader_ownership.py` reports **1 failed, 7 passed**, unchanged by anything in this plan.
- The full suite reports **14 failed, 7227 passed, 180 skipped**, also unchanged by this plan.
- The syntax-tree assertion convention appears in **8** existing test files.

## Verification Signal

Success here is narrow and worth stating precisely, because this plan builds nothing: the amendment plan becomes executable without a reviewer finding the same two defects again. Concretely, its Step 1 specifies a check that passes when the code is correct and fails when it is not; its Step 4 names every file carrying a claim that step changes; and its Step 3's stated intent matches its mechanism, so running it does not rewrite a committed artifact that two documents cite.

Failure is equally clear. If recomputing the wave grouping after the file-list change produces a different grouping than the three waves currently recorded, the amendment plan's Execution Preview is stale and must be corrected in the same pass rather than left to contradict its own machine-readable output.

What this plan explicitly does **not** move: the red ownership guard, the suite's failure count, and any measured number.

## Implementation Steps

### Step 1: Apply the three corrections to the amendment plan and re-validate it
Files: plans/Clear_The_Three_Blockers_The_Frame_Conditioning_Review_Found_Plan/technical-plan.md
Depends on: none

**What changes:**
- Step 1's test strategy in that plan replaces the substring assertion with a syntax-tree assertion following the convention in `tests/test_amr_uncounted.py`: parse the module, walk for import nodes, and assert that none names `gzip` and that no attribute access targets it. The reason is recorded in the step text, so a future reader sees that a text search was rejected rather than overlooked.
- Step 4's file list in that plan gains `wiki/plans-index.md`, and its body gains the requirement to rewrite the line-869 key-decision bullet so the index stops asserting exact replication. The existing requirement to re-label the three pinned figures with their ordering is unchanged.
- Step 3's body in that plan resolves the output-location question by giving the ordering comparison its own script rather than extending the committed probe, which removes the conflict entirely instead of adding a flag that a future caller can forget. Its file list changes accordingly, and its stated promise that the existing artifact is untouched becomes true by construction.
- Because a file list changed, the wave grouping is recomputed with the toolkit rather than reasoned about, and the Execution Preview section is corrected to whatever the toolkit reports. The footer is regenerated and replaces the old one as the final line.

**Test strategy:**
- The toolkit's `check` must report no missing sections and no errors, and `waves` must report no errors and no unintended intra-wave file overlap.
- The recomputed grouping is compared against the three waves currently recorded. If it differs, the Execution Preview is corrected; if it matches, that is asserted rather than assumed, because a changed file list is exactly the input that can silently move it.
- The footer must be the final line of the file after the edit, since consumers read it from the file tail.
- A scope check by reading: the amendment plan still has five steps, still clears all three original blockers, and the only things weakened are the two defective specifications. No step is removed and no acceptance requirement is relaxed.

### Step 2: Record which review findings were accepted, and how
Files: plans/Clear_The_Three_Blockers_The_Frame_Conditioning_Review_Found_Plan/brainstorm.md
Depends on: Step 1

**What changes:**
- A dated entry records the two critical findings as accepted with the correction applied for each, the medium finding about the artifact as accepted and resolved by the new-script route, and the two remaining medium findings as acknowledged but not applied, with the reason stated.
- The entry also records that one finding the review raised was checked and refuted: the concern that translating the owner's error type would miss an absent-file case does not hold, because the owner raises its own typed error with the message about downloading being disabled.
- Nothing already written is deleted. The original findings stay as they were so the chain from finding to resolution is readable in one place.
- The entry notes the single-model limitation that applied to the review, so a later reader knows the findings came from a review without an independent second opinion.

**Test strategy:**
- Manual read-through that each of the five findings in the review has exactly one disposition recorded, and that no disposition claims an edit the first step did not make.
- Assert the file still parses as the Markdown it was: headings intact, no truncated section.
- No test asserts the dispositions are *correct*, which is the honest limit of a record: it captures what was decided, not whether the decision was sound.

## Execution Preview

Wave 0 (1):  Step 1 — apply the three corrections and re-validate
Wave 1 (1):  Step 2 — record the dispositions

Critical path: Step 1 → Step 2 (2 waves)
Max parallelism: 1 agent

There is no parallelism to exploit, and that is a property of the work rather than of the decomposition: both steps edit documents in the same plan folder and the second describes what the first did. Sequential execution is the right mode, and a fresh worktree would offer nothing here because no step reads data excluded from version control.

## Risk Flags

- Severity: high — **The planning chain is recursing while the only live defect stays red.** This is the third review-then-plan cycle, no code has changed since the exposure work, and the fix for the red ownership guard is sitting unexecuted one plan away. The mitigation is not technical: it is to execute the amendment plan after this correction rather than review it again. A fourth cycle on the same subject should be treated as a signal to stop planning.
- Severity: medium — **Changing a file list changes the wave-grouping input, so the Execution Preview can silently go stale.** The correction recomputes the grouping with the toolkit and compares it against the recorded three waves rather than assuming it held. This is the recorded lesson that a split or an edit resets a step's ordering and the computation must be re-run rather than hand-reasoned.
- Severity: medium — **Editing a plan to weaken two of its specifications could quietly weaken more.** Step 1's test strategy includes a read-through confirming all five steps survive and no acceptance requirement is relaxed, because "correct a defect" and "soften a bar" look similar in a diff.
- Severity: medium — **Resolving the artifact question by adding a script rather than a flag leaves two scripts that share most of their logic.** That is a real maintenance cost and the alternative was considered: a flag keeps one script but relies on every future caller passing it, and a forgotten flag silently rewrites a committed artifact that two documents cite. The duplication is the cheaper failure.
- Severity: low — **Two documents this plan edits are uncommitted on a working branch.** The amendment plan folder and the modified index sit uncommitted on `promoteratlas-exposure-2026-11-01`, which is three commits ahead of `main` and unpushed. Editing uncommitted files is fine, but there is no committed state to diff the corrections against, so the read-through in Step 1 is the only scope check available.
- Severity: low — **External tool surfaces.** The toolkit's `check`, `waves` and `validate` subcommands were confirmed by direct invocation repeatedly during this session, and on this host the name `python3` resolves to a Microsoft Store placeholder so `python` must be used. `sentrux` resolves but no architecture gate applies, because this plan changes no source module and adds no dependency edge. No step contacts a network or a service, so no external schema is left unconfirmed.
- Severity: low — **No restructuring was applied during wave grouping.** Both steps edit files in the same folder and the second depends on the first, so there is nothing to parallelise. The one choice worth naming is that the three corrections were kept as a single step rather than split into three, because all three edit the same file and splitting them would have produced three steps that must run in strict sequence with no independent verification between them.

## Open Questions

- **Should these three edits have gone through the save-time amendment mechanism instead of a new plan?** That mechanism exists precisely for appending amendments to a saved plan and re-validating its footer, which is most of what Step 1 does. The case for a plan is that two of the three corrections change step contracts rather than annotating them, and the amendment mechanism's own contract says amendments are audit notes and not executable instructions, so a contract change arguably belongs in the body. The case against is proportionality. This plan takes the plan route because one was asked for, and the question is recorded rather than settled.
- **Should the ordering comparison live in its own script, or behind a flag on the existing one?** This plan picks the new script, on the grounds that a flag relies on every future caller remembering it while duplication fails loudly. The counter-argument is that two scripts sharing a corpus-window replica will drift, and a drifted replica produces two different window sets under the same name. If the shared logic is extracted into the package rather than copied, that objection mostly dissolves — which is a third option this plan does not take.
- **Is the single-model limitation on the preceding review material enough to warrant redoing it?** Three consecutive reviews of this chain ran without an independent second opinion, and each reviewed work the same author produced. Two findings per round were established by measurement rather than argument, which is the partial mitigation, but the accumulated anchoring risk is real and unquantified.

## Architecture Impact

**Expected structural changes:**
- None. No module is added, removed or modified. No import edge changes. The only source-adjacent effect is that a future execution of the corrected plan will create one new script, and that creation is not part of this plan.

**Risk areas:**
- None in the dependency graph, because no Python file under `dna_decode/` is touched.
- The documents edited are inputs to the planning toolchain rather than to the runtime, so the regression surface is "a plan that no longer validates" rather than "code that no longer works". Step 1's test strategy covers exactly that.

**Sentrux gates:**
- Baseline: not needed. One already exists at `.sentrux/baseline.json` and this plan changes no module.
- Compare: no. A compare would measure an unchanged graph and report noise.
- Acceptable delta: none expected, since no source file changes.
- Regression-response stance: investigate. Any reported structural change on a documentation-only plan would mean something outside this plan's scope moved, which is worth knowing.

**Anti-drift checks:**
- Every Python file under `dna_decode/` and every file under `scripts/` shows an empty diff at the end of this plan.
- The five frozen decoder files and the dependency lockfile are byte-unchanged, and the active prospective lock re-verifies.
- The ownership guard's result is unchanged at 1 failed and 7 passed, confirming this plan did not accidentally touch the code path it describes.

## Verification

1. The corrected amendment plan passes the toolkit's `check` with no missing sections and no errors, and `waves` with no errors and no unintended overlap.
2. The recomputed wave grouping is compared against the recorded three waves, and the Execution Preview matches whatever the toolkit reports.
3. The toolkit footer is the final line of the corrected plan.
4. The corrected plan's Step 1 specifies a syntax-tree assertion rather than a substring search, and its Step 4 file list includes `wiki/plans-index.md`.
5. The corrected plan's Step 3 gives the ordering comparison its own script, so its claim that the existing artifact is untouched is true by construction.
6. A read-through confirms all five steps survive, all three original blockers are still addressed, and no acceptance requirement was relaxed.
7. The review sidecar records exactly one disposition for each of the five findings, including the one that verification refuted.
8. Every file under `dna_decode/` and `scripts/` shows an empty diff, the frozen files and lockfile are byte-unchanged, and the prospective lock re-verifies.
9. The ownership guard still reports 1 failed and 7 passed, and the full suite still reports 14 failed, confirming this plan moved no test.

<!-- toolkit: check=clean waves=clean gate=fired:open-questions,severity-high -->
