# Pre-execution review of this plan, and what was done about it

<!-- doc-status: active -->
<!-- measured-on: 2026-11-01 -->

<!-- Review of `technical-plan.md` in this folder, run after it was saved and before any step
     executed. SINGLE-MODEL REVIEW: the second-opinion tool was unavailable (`codex` exited 127, not
     on PATH; LM Studio UNREACHABLE and `lmstudio-auditor.sh` absent), and the reviewer wrote the
     plan, so anchoring-bias risk is real. Two of the three accepted findings were established by
     measurement rather than argument, which is the partial mitigation.

     FILE PROVENANCE, stated because the plan got this wrong: the plan's Step 2 specified appending
     a dated entry to this file and promised that "nothing already written is deleted". This file did
     not exist. The review had been emitted to chat and never saved, so there was nothing to append
     to and nothing to preserve. The step was executed as CREATE rather than APPEND, which satisfies
     its intent and additionally rescues a review that would otherwise have been lost with the
     conversation. Recorded as a deviation rather than silently absorbed. -->

## DECISIONS NEEDED

Two findings were acknowledged and deliberately **not** applied, and whether to apply them is yours.

1. **M2 — the guard exemption's comment records a count but cannot enforce it.** The plan's Step 2
   adds a row to the ownership guard's allow-list with a comment recording how many prose mentions
   survive. A comment cannot fail, so an increase is visible only if someone reads it. A test
   asserting the exact occurrence count in `tile_controls.py` would convert the observation into a
   guard. Not applied because it adds a test the plan did not ask for.
2. **M3 — the matched-strength baseline control is not specified concretely enough to implement.**
   The plan requires the frame plan's bar to name "a series of baselines of increasing standalone
   strength" without saying what they are. The feature registry in `expression.py` offers GC content,
   four-mer and six-mer frequencies, and positional one-hot encoding, so a concrete ladder is
   constructible. Leaving it unnamed invites the series to be chosen after seeing which composition
   gives the preferred answer, which is the failure pre-registration exists to prevent. Not applied
   because naming the ladder is a measurement-design decision rather than a defect repair.

## CORRECTIONS

**A correction introduced during execution, caught only by re-validating.** Applying the first fix
below added a new advisory: citing `tests/test_amr_uncounted.py` inside a Test strategy block leaked
a path that was not in Step 1's file list, so the toolkit footer went from
`gate=fired:open-questions,severity-high` to `…,test-strategy-leak`. Resolved by adding a
`Test-strategy-ignore:` line to Step 1 at line start, since that path is inspection-only — the plan
reads it as a pattern to copy and never modifies it. **The lesson is the ordering:** the footer was
regenerated rather than assumed still valid, which is the only reason the new firing was seen at all.

## Dispositions

| # | Finding | Disposition |
|---|---|---|
| C1 | The structural test would fail on prose | **Accepted, applied** |
| C2 | The index restates the claim Step 4 removes | **Accepted, applied** |
| M1 | Extending the probe rewrites a committed artifact | **Accepted, applied** |
| M2 | The exemption comment enforces nothing | Acknowledged, not applied |
| M3 | The matched-strength control is under-specified | Acknowledged, not applied |

### C1 — accepted and applied

The plan's Step 1 specified a test asserting the control-loader module "no longer contains the
substring `gzip`". The module's docstring at line 12 contains `embeds genomewide.py as base64+gzip`,
and that sentence survives the fix it describes, so the test would fail while the code is correct.

Replaced with a syntax-tree assertion following the convention at `tests/test_amr_uncounted.py`
lines 145 to 158, which is used in seven further test files: parse the module, walk it, and assert
that no import node names `gzip` and no attribute access targets it. The step text now records that
a text search was **rejected** rather than overlooked, so it cannot be reinstated as a
simplification. The underlying reason is the same one that makes the ownership guard's own
`\bGSE\d+` pattern fire on docstrings: a text search over source cannot distinguish an instruction
from a comment.

### C2 — accepted and applied

The plan's Step 4 rewrites the frame plan's Step 1 so it stops claiming exact reproduction of the
upstream extraction function, but `wiki/plans-index.md` line 869 asserts that Step 1 *replicates the
authors' `extract_promoter_regions` exactly, INCLUDING two apparent defects*. Step 4's file list did
not include the index, so the amendment would have left a committed surface restating the claim it
removes — on the plan whose purpose is correcting false claims.

`wiki/plans-index.md` was added to Step 4's file list and its body gained the requirement to rewrite
that bullet. The wave grouping was recomputed afterwards rather than assumed: it is unchanged at
`[[1,3],[2,4],[5]]` with maximum parallelism two and zero file overlaps, because Step 4's wave
partner, Step 2, touches only the guard test.

### M1 — accepted and applied

The plan's Step 3 extended `scripts/frame_plan_feasibility_probe.py` while promising that the
existing artifact was left untouched. That promise was false as written: the script fixes its output
path as a module constant at line 38, includes a `wall_seconds` field at line 262, and writes the
file unconditionally at line 264, so any run of a modified version rewrites
`wiki/frame_plan_feasibility_2026-11-01.json` — and the timing field guarantees the rewrite differs
even when every input is identical. That artifact is committed and cited by both the frame plan and
its review.

Resolved by giving the ordering comparison its own script, `scripts/frame_window_ordering_probe.py`,
which makes the promise true by construction. A flag on the existing script was considered and
rejected, because it would rely on every future caller remembering to pass it. The accepted cost is
duplication: both scripts carry a corpus-window replica and two copies can drift into producing
different window sets under the same name. If that drift becomes real, the recorded fix is to
extract the replica into the package rather than to re-merge the scripts.

### A concern that verification refuted

The review raised that translating the owner's error type into this module's own might miss the
absent-file case. It does not. `ensure_cached` raises its own `GenomeWideDataError` with the message
about downloading being disabled when the file is missing and downloading is off, so the translation
covers the case the existing test exercises. Recorded because a refuted concern is part of the
record, not something to drop.

## Scope note

This review and its dispositions changed two documents and no code. The chain that produced them is
three planning cycles deep on one subject, and the live defect it all concerns — a red ownership
guard caused by a second reader of a data substrate — is still red. The plan's own high-severity risk
flag says a fourth cycle should be read as a signal to stop planning and execute, and that remains
the standing recommendation.
