# Execution Log — Clear_The_Three_Blockers_The_Frame_Conditioning_Review_Found_Plan
Date: 2026-11-01
Waves: 3 (max parallelism 2) — executed SEQUENTIALLY; see rationale below
Files changed: dna_decode/glm/tile_controls.py, scripts/frame_window_ordering_probe.py, tests/test_tile_controls.py, tests/test_frame_window_ordering.py, tests/test_substrate_loader_ownership.py, plans/Frame_Conditioning_Does_Corpus_Window_Position_Explain_PromoterAtlas_0_5076_Plan/technical-plan.md, plans/Frame_Conditioning_Does_Corpus_Window_Position_Explain_PromoterAtlas_0_5076_Plan/brainstorm.md, wiki/plans-index.md, wiki/frame_plan_window_ordering_2026-11-01.json, wiki/frame_plan_amendments_2026-11-01.md, project_state/glm-learned-representation-2026-10-07.md
Sentrux verdict: n/a — the plan's own Architecture Impact specifies no compare for a single added import edge; the importer-set assertion was used instead and passed
Commit: n/a — branch only (uncommitted on promoteratlas-exposure-2026-11-01)

## Mode

Sequential, although the toolkit permitted parallel (wave 0 holds 2 steps and a remote is
configured). All three required data files are untracked, so a fresh worktree would have been unable
to run Step 3 at all and would have silently skipped Step 1's real-data assertion. Parallel mode also
opens pull requests, which were not requested.

## Outcome

THE RED OWNERSHIP GUARD IS GREEN. tests/test_substrate_loader_ownership.py moved from 1 failed /
7 passed to 8 passed. dna_decode/glm/tile_controls.py no longer parses the GSE144621 substrate; it
calls genomewide.load_peak_tiles and filters the result, making it a consumer of the registered owner
rather than a second reader.

Routing changed no data: 936 random + 470 neg_control at 150 bp, total 1,406, byte-identical to what
the deleted parser produced. The cost was elsewhere — the owner enforces MIN_COLUMNS=8 and requires
peak_start/peak_end/strand, so seven test fixtures had to widen from five columns. That is why an
earlier estimate of "roughly a ten-line change" was wrong.

THE FIX NEEDED BOTH ACTS, IN ORDER, AND THAT WAS MEASURED. With the parser deleted the guard STILL
FAILED, on the one surviving docstring mention, because its criterion is the bare text pattern
\bGSE\d+ over whole files. The ALLOWED_MENTIONS row was added only after observing that failure,
which is what makes it evidence rather than assumption. An exemption alone would have covered a live
reader; a deletion alone was insufficient.

THE ORDERING DIVERGENCE IS REAL. File-order threading gives 3,490 windows over 624,812 bases against
sorted ordering's 3,496 over 625,255 — a symmetric difference of 26 windows, with 146 of 43,061
uniquely-located tiles changing in-window status. Small (0.74% and 0.34%) but non-zero, so the frame
plan's claim of exact reproduction could not stand. A zero difference had been pre-committed as an
equally reportable finding, so there was no incentive to tune the comparison until it found something.

## Tests

- Baseline 14 failed / 7227 passed / 180 skipped.
- After the five plan steps: 13 failed / 7238 passed / 180 skipped. The +11 reconciles exactly as
  2 new control-loader tests + 8 new ordering tests + 1 guard test flipping to pass.
- After the test epilogue's additions: 13 failed / 7255 passed / 180 skipped. The +17 is exactly the
  new tests, and the failure SET was diffed rather than only counted — identical.
- The one failure that disappeared is the ownership guard. The remaining 13 are pre-existing.

## Corrections this execution produced

CLAUDE.md's standing claim that "all 14 remaining are absent DATA, NONE a code defect" is wrong
twice: the count is now 13, and one of the 13 is a registry count drift rather than missing data
(tests/test_tile_chain_constraint_registered.py asserts 14 == 12 against a registry holding 14, with
its eleven named constraints all passing). Deliberately NOT re-pinned: a one-character bump is
exactly the tidying that rail exists to interrupt, so it is a guard-structure decision for the user.

RETRACTED 2026-11-01, same day. This log originally claimed "CLAUDE.md is NOT one word from its
ceiling; wc -w gives 34,661 against a recorded ceiling of 35,460, so the 35,459 figure in NEXT.md is
stale." THAT CLAIM WAS WRONG AND THE RECORDED FIGURE IS CORRECT.

The ceiling is DEFINED by tests/test_claude_md_growth_ceiling.py, whose measure is
_words(text) = len(text.split()). That measure reports 35,459 against WORD_CEILING = 35_460, i.e.
exactly one word of headroom, and the test passes. A sibling test,
test_the_ceiling_is_not_vacuously_slack, exists precisely to catch a ceiling drifting away from the
live file (it asserts slack <= 4,000) and it passes at slack = 1. The guard is working.

wc -w reports 34,661 on the same file -- a 798-word disagreement with the defining measure whose
cause is NOT established (the file contains no non-ASCII whitespace, and Python's whitespace split
and an ASCII-only split agree exactly, so the divergence is on wc's side and is left unexplained
rather than guessed). Do not use wc -w against this ceiling.

THE ERROR WORTH KEEPING: a number is only comparable against the measure that DEFINES it. The same
session had already applied that rule correctly hours earlier -- pinning ridge-protocol parity at the
reference's own rounding precision rather than widening a tolerance -- and then broke it here, using
a convenient tool instead of the authoritative one and telling the user their own correct record was
stale prose.

## Known hole, recorded not closed

load_control_tiles(categories=()) returns [] — the empty value the module exists never to return,
because the refusal is constructed from the REQUESTED categories. No call site does this, so it is
latent rather than live. Closing it is a one-line guard but is a behaviour change, so it was left to
the user and pinned as a known hole.

## Not done

The frame-conditioning experiment itself. This made its plan executable and cleared a regression; the
three cells, the acceptance bar and every measured question remain ahead.
