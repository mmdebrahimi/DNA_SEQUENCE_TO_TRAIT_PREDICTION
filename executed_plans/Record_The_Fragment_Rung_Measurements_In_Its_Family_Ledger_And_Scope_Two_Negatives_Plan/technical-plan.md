## Lens status

Repo-index: absent (`repo-index-status.sh` produced no output) — source-only reading, no degradation.
Reuse signals: graphify absent; no generic util dirs. sentrux present / gh present / **codex absent**.
Prior review: `/brainstorm` ran Claude-only (codex 127, LM Studio unreachable, auditor script absent) and
its accepted findings are this plan's input. Optional context present: `wiki/` only — no `project-rules.md`
(so `## Project Rules Check` is omitted), no `DESIGN.md`, no `AUDIT_REPORT.md`.

## Problem Statement

The fragment-pyramid rung produced three measurements over 2026-10-28/29 — a 10-seed Kaggle T4 run, an
independent server1 replication, and a 40-seed null extension that pooled to 50 seeds — and **none of them
is recorded in any family ledger.** They exist only as `wiki/` artifacts and commit messages. The owning
family is G-A (`glm-learned-representation-2026-10-07`), whose refined goal is a learned multi-resolution
encoder plus a trivial head; the fragment rung is a second rung of that family alongside the tile rung.

That matters for a specific, measured reason: `advance_ranker` and every fresh session read
`project_state/`, not commit messages. A measurement that corrected a twice-published number
(the null's apparent +0.05 bias) is therefore invisible to the mechanism that chooses what to work on next.

Two smaller scope notes from the same review are folded in because they are the same class of defect —
a true result whose boundary is not written down where a later reader will meet it:

1. G-C's `NO_GAIN` verdict is richly documented mechanically but never states that the eukaryotic
   mechanism which would make condition-conditioning pay is *structurally absent* from a prokaryotic
   substrate. The ledger's originating goal does say "bacterial analogue", so the scope cue is partly
   there; what is missing is the consequence.
2. `scripts/promcalc_decomposition_arm.py:14` already names sigma70 physics, but nowhere records that this
   makes the decomposition non-transferable to eukaryotic promoters.

Non-goals: no eukaryotic work; no change to any frozen acceptance bar; no edit to the frozen five; no new
family ledger; no re-run of any arm.

## Codebase Context


See: ./codebase-context.md

## Pre-Change Baseline

Measured 2026-10-30, all read-only:

| quantity | value |
|---|---|
| files in `project_state/` + `features/` naming `fragment_pyramid` | 0 |
| files in `project_state/` containing `0.0242` (the pooled null median) | 0 |
| hits for `kaggle` / `INVALID_NULL` / `learned_contribution` in G-A's ledger | 0 |
| G-A ledger Action-Log rows | 11 |
| G-C ledger Action-Log rows | 6 |
| occurrences of `eukaryot` in G-C's ledger | 0 |
| occurrences of `eukaryot` in `scripts/promcalc_*.py` | 0 |
| occurrences of `sigma70` in `scripts/promcalc_*.py` | 1 (`promcalc_decomposition_arm.py:14`) |

Behavioural baseline to preserve: `advance_ranker.rank` on the GLM umbrella currently reports 3 eligible /
2 blocked:predecessor, and the frozen five are byte-unchanged. No step may alter either.

## Verification Signal

- `grep -c 0.0242 project_state/glm-learned-representation-2026-10-07.md` goes 0 -> at least 1.
- G-A Action-Log rows 11 -> 12; G-C 6 -> 7; each file grows by exactly one line.
- `eukaryot` goes 0 -> at least 1 in both G-C's ledger and `scripts/promcalc_decomposition_arm.py`.
- Every target file's CRLF count rises by exactly the number of appended lines and its bare-LF count
  stays 0.
- `uv run pytest tests/test_promcalc_decomposition_arm.py -q` still passes after Step 3.
- `advance_ranker.rank` output on the umbrella is UNCHANGED (same eligible set, same blocked set) — these
  are Action-Log and comment edits, so any change to eligibility is a defect.
- `git diff --stat` touches exactly the three intended files; the frozen five appear nowhere.

## Implementation Steps

### Step 1: Record the fragment-pyramid rung's measured history in G-A's ledger
Files: project_state/glm-learned-representation-2026-10-07.md
Depends on: none

**What changes:**
- Append ONE Action-Log row (row 12, class `run-tests`) at the
  `<!-- project-state:end:action-log -->` marker, per the AC9 sanctioned cheap path.
- Every figure is read from the committed artifacts, not retyped from prose: the 10-seed Kaggle run
  (`kmer4_median 0.112`, `pyramid_median 0.1001`, `delta_median -0.0127`, 2/10 wins, null median 0.0517)
  and the 40-seed extension (`seeds 10..49`, null median 0.0191, sd 0.0649, delta -0.0061, 13/40 wins,
  ci95 `[-0.0133, -0.0037]`, `learned_contribution` +0.0965 with 37/40 positive).
- Record the three load-bearing readings: (a) the pooled 50-seed null median is **+0.0242**, which FAILS
  the `|median| < 0.02` bar, so quote the pooled figure and not the fresh-seed 0.0191; (b) both earlier
  runs drew `seed_range {start: 0, n: 10}`, so their agreement at ~+0.05 was ONE sample reported twice and
  the apparent bias was an artifact of those ten initialisations; (c) `learned_contribution` is
  `pyramid_spearman - null_spearman` on the same init and split, so a lower null mechanically raises it —
  its rise and the null's fall are the SAME fact, not two independent wins.
- Record that the gate may be unsatisfiable as specified — sd is initialisation-driven and stable across
  all three samples (0.0661 / 0.0649 / 0.0654), so SE(median) shrinks only as 1/sqrt(n) while the point
  estimate sits ABOVE the bar — and mark re-specifying the bar as an ACCEPTANCE-BAR decision that is
  surfaced, NOT taken.
- Note that `kmer4` holds ρ 0.1009–0.1241 across all ten seeds while the pyramid holds 0.0839–0.1249, so
  the task, splits, labels and scoring are all validated by the comparator and the pyramid is learning
  substantially rather than failing — what is unvalidated is the torch path specifically.
- Touch nothing else: no criterion, no hypothesis row, no candidate row, no frame row, no progress proxy.

**Test strategy:**
- Byte-level: assert the file grows by exactly one line, CRLF 257 -> 258, bare LF stays 0.
- Re-read every quoted figure out of `wiki/glm_fragment_pyramid_arm_kaggle_full_2026-10-28.json`,
  `wiki/glm_fragment_pyramid_arm_kaggle_ext_2026-10-29.json` and
  `wiki/fragment_null_ext_summary_2026-10-29.json` and fail on any mismatch — no figure typed from memory.
- Re-run `advance_ranker.rank` and diff against the baseline eligible/blocked sets; any change is a defect.
- Confirm the appended row's class is one of the 7-value enum.

### Step 2: Scope the G-C NO_GAIN verdict to its substrate class
Files: project_state/glm-condition-conditioning-2026-10-07.md
Depends on: none

**What changes:**
- Append ONE Action-Log row (row 7, class `research`) recording the boundary the ledger does not yet
  state: the one-sequence-to-many-outputs expansion that makes cell-state conditioning pay in eukaryotes
  is implemented by machinery ABSENT from this substrate (spliceosome composition, cell-type-specific
  snRNA expression, chromatin state), so `NO_GAIN` is evidence from a system where the mechanism cannot
  operate — not evidence against condition-conditioning generally.
- State plainly what does NOT change: the verdict, the +0.0001 figure, the 0.02 bar, `null_clean=true`,
  and the family's terminal status recorded in row 3 all stand. This row adds a boundary, not a revision.
- Acknowledge in the row that the ledger's originating goal already says "bacterial analogue", so the
  scope cue is partly present; the row supplies the consequence, which is the part a later reader needs.

**Test strategy:**
- Byte-level: exactly one line added, CRLF 199 -> 200, bare LF stays 0.
- Assert `eukaryot` goes 0 -> at least 1 in this file.
- Assert the strings `NO_GAIN`, `+0.0001` and `null_clean=true` still appear with unchanged counts —
  the scope note must not have edited the verdict.
- Confirm no candidate row, criterion or hypothesis row changed (`git diff` shows a single insertion).

### Step 3: Name the eukaryotic non-transfer consequence beside the existing sigma70 claim
Files: scripts/promcalc_decomposition_arm.py
Depends on: none

**What changes:**
- Extend the module docstring at the existing sigma70 statement (line 14) with the consequence: the
  Promoter Calculator models sigma70 bacterial promoter architecture — the field set is
  `dG_10, dG_35, dG_ext10, dG_spacer, dG_disc, dG_ITR, dG_UP`, i.e. -10 and -35 boxes and the spacer
  between them — so the decomposition does NOT transfer to eukaryotic promoters, whose initiation runs
  through a TATA box, roughly six general transcription factors, Mediator, distal enhancers and
  chromatin accessibility. A different object, not a scaled one.
- Label the organism attribution honestly as INFERRED from the tool's identity plus its field names
  rather than read from this repo, since no promcalc script states it beyond the existing line 14.
- Add the one-line consequence a future session actually needs: reusing promcalc for the north star's
  "climb to animals" step would be a category error, and the mirror of the already-recorded finding that
  AlphaGenome's architecture does not transfer DOWN to a 100-300 bp bacterial promoter.
- Comment-only change. No code, no signature, no behaviour.

**Test strategy:**
- `uv run pytest tests/test_promcalc_decomposition_arm.py -q` passes (comment-only, so it must).
- `python -c "import ast, pathlib; ast.parse(...)"` confirms the file still parses.
- Byte-level: bare LF stays 0; CRLF rises by exactly the number of appended comment lines.
- Assert `eukaryot` goes 0 -> at least 1 in `scripts/promcalc_decomposition_arm.py` (the file this step
  edits), and that the pre-existing sigma70 sentence is still present verbatim (extended, not replaced).

## Execution Preview

- **Wave 0** (parallel, 3): Step 1, Step 2, Step 3

1 wave, max parallelism 3. Critical path: any single step. No intra-wave file overlap — the three steps
touch three distinct files.

**Run it SEQUENTIALLY anyway.** Parallel mode isolates via `git worktree`, and the recorded lesson here is
that worktrees are the wrong shape when the targets are small edits in a shared tree with gitignored
siblings; three minute-scale edits gain nothing from isolation and the reconciliation cost is real. Cost:
Step 1 ~10 min (the figure re-read is the work), Steps 2-3 ~5 min each. Total under 25 minutes.

## Risk Flags

- **Step 1 is the one that can go wrong, and the failure mode is a wrong number.** The row quotes about a
  dozen figures. The mitigation is in its test strategy — every figure is re-read from the committed JSON
  and mismatches fail — but a figure that is correct yet mislabelled (a null median presented as a delta,
  say) would pass that check. This repo has a measured instance of exactly that: a sensitivity and a
  specificity of equal value swapped without detection. Read the row back against the artifacts once.
- **Quoting 0.0191 instead of the pooled 0.0242 would be cherry-picking the flattering subsample.** The
  fresh 40 seeds clear the bar as a point estimate; the pooled 50 do not. Step 1 must carry the pooled
  figure as the headline.
- **The sigma70 organism claim is `[inferred]`, not `[grounded]`.** It rests on the Promoter Calculator
  being a sigma70 model plus the field names being sigma70 architecture — sound domain reasoning, but not
  stated in this repo beyond line 14. Step 3 labels it as inferred rather than asserting it.
- **Step 2's marginal value is smaller than the review implied** — the ledger's originating goal already
  reads "bacterial analogue", and rows 2-3 give a thorough mechanistic account. If a reviewer judges the
  consequence sentence redundant, dropping Step 2 costs little.
- **No step clears any advisory.** G-A's five "listed as pending" candidate warnings are unclearable by
  any safe action (recorded in that ledger's row 10) and remain after this plan. Do not expect the ranker
  output to improve.
- **No external surface is in scope**, so there is no flag, config key, exit code or output format to
  ground. No verification-gap marker is claimed either, because none applies: such a marker asserts a
  surface existed and resisted checking, and here nothing external is touched at all. `sentrux` is
  present but these are doc and ledger edits, so the architecture gate is not structural for this plan.

## Open Questions

- **Should the fragment rung get its own family ledger rather than rows in G-A?** G-A's refined goal is
  pinned to TILES (ceiling 0.9410) while this rung is FRAGMENTS (321,123 at 48-499 bp, a different noise
  ceiling the feature record explicitly distinguishes). Two rungs under one family is defensible and
  cheaper; a separate family would make the rung visible to the ranker in its own right but spends a
  self-init slot and adds an umbrella flow-down edge. Not resolved here; Step 1 takes the cheap option.
- **Is the `|median| < 0.02` null bar satisfiable at all?** The plan RECORDS the arithmetic and refuses to
  re-specify the bar. That is an acceptance-bar decision and therefore the user's.
- **Should a test pin the Step-3 scope note?** This repo pins prose elsewhere (help text, report-card
  prose), but those are trust surfaces a user reads. A module docstring is not, so no test is proposed —
  deliberately, as the alternative is ceremony with no consumer.
- Whether `features/glm-learned-representation/feature.md` should carry the same fragment-rung record as
  the ledger. Left out to avoid two copies of one claim drifting apart; the ledger is what the ranker reads.

## Verification

1. `git diff --stat` shows exactly three files changed, with +1 / +1 / +N insertions and zero deletions.
2. Per-file byte check: bare LF stays 0 in all three; CRLF rises by exactly the appended line count.
3. `uv run pytest tests/test_promcalc_decomposition_arm.py tests/test_masked_pyramid.py -q` passes.
4. Figure audit: re-read every number in the Step-1 row out of the three committed `wiki/` artifacts and
   confirm each matches, including that the pooled median reads 0.0242 and not 0.0191.
5. `advance_ranker.rank` on `project_state/glm-umbrella-2026-10-07.md` returns the same eligible and
   blocked sets as the baseline.
6. Frozen five unchanged: `git diff --stat` names none of `dna_decode/eval/amr_rules.py`,
   `dna_decode/data/calibrated_amr_rules.json`, `dna_decode/data/mic_tiers.py`,
   `dna_decode/data/shipped_decoder_surface.py`, `dna_decode/eval/cohort_manifest.py`.
7. Post-change greps match the Verification Signal table: `0.0242` at least 1 in G-A's ledger,
   `eukaryot` at least 1 in G-C's ledger, and at least 1 across `scripts/promcalc_*.py` as a repo-wide
   post-check (the per-step assertion in Step 3 is scoped to the single file it edits).

## Save-time amendments

Captured at: 2026-10-30
Source: `/save-plan` arguments

These are AUDIT NOTES for human readers, not executable instructions: `/execute-plan` reads only
`## Implementation Steps`. If an amendment changes a Step contract (file list, dependency,
structure), the toolkit fires on re-validation and the plan should go back through
`/technical-plan` before execution.

- scope corrected during planning: the fragment rung IS owned by G-A (found only after a vocabulary fix — `fragment_pyramid` returns 0 hits while the family says 'pyramid arm'), so the plan records the rung in G-A rather than claiming no family owns it
- the sigma70 claim already exists at promcalc_decomposition_arm.py:14 so Step 3 extends it in place rather than starting a second copy
- Step 2's marginal value downgraded because the G-C ledger's originating goal already reads 'bacterial analogue'
- two toolkit gate findings fixed before emit (a test-strategy glob outside Step 3's Files list, and a misused verification-gap marker)
- a separate fragment-rung family ledger and re-specifying the null bar both left as open questions rather than folded in

<!-- toolkit: check=clean waves=clean gate=fired:open-questions,reuse-missing -->
