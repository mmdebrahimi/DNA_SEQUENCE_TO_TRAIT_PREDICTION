# Correct Two Wrong Claims I Committed, And Make The Next Check Runnable

## Lens status

Dispatch: conversation-driven (no args) — scope is the review findings C1, C2, M1, M3 against commit `b28fa8dc`.
Repo index: absent (`repo-index-status.sh` not installed) — direct file reading only, no degradation.
Tools: `sentrux` present, `gh` present, `codex` ABSENT — the review that produced these findings was single-voice.
Project rules: `project-rules.md` absent — Project Rules Check omitted.
Wiki replay: `wiki/` present; the subject memo is itself one of the files under correction.
Substrate: `D:` and `E:` detached — every step here is deliberately drive-independent.

## Problem Statement

Commit `b28fa8dc` recorded an honest refusal — the lentiMPRA published metric stays `UNSPECIFIED_IN_SOURCE` — but wrapped it in two claims that are wrong, plus one over-claim inherited from before. All three are in the surfaces a future session reads first.

1. **The "decisive next check" cannot work as written.** It says to read the workbook's *header/units row*. The parser shows the header is `test fold | <method names>` and carries no metric. What a future session should read is the rows the parser **discards** — every row whose first cell is not a cell type, not `test fold`, and not numeric falls through silently, and a caption or units note is exactly the kind of row that lives there. The real check survives; my stated version of it does not.
2. **"Blocked on provenance, not on compute" conflates cannot-COMPARE with cannot-RUN.** Nothing gates running the eukaryotic arm on the published metric. The acquisition memo's own wording is *"no number here may be compared against it"*. The arm can be trained on the 159,795 verified pairs and report its own coefficient; the incumbents are context. Calling the metric that arm's blocker invents a wall — the same shape as the `kaggle.json` precondition named earlier in this session for a route nobody was taking.
3. **The docstring over-claims.** *"the sheet gives numbers and method names and never says whether they are Pearson, Spearman or something else"* is a statement about the **rows the parser reads**. The parser never looks at the rest, so the sheet has not been shown to be silent.

Two softer findings are in scope because they are cheap and bound the committed claim: the table-identity mismatch may be a **version artifact** (the text read was the author-accepted version, whose supplementary numbering can differ from the published one), and the publisher's own supplementary descriptions — the direct route to table identity — were never fetched.

**Success criteria.** No surface claims a check that cannot work; no surface says the eukaryotic arm is blocked on the metric; the docstring's claim is scoped to what the parser reads; the real next check is **runnable code** rather than a prose instruction; and the version hypothesis is either settled or recorded as attempted-and-open.

**Non-goals.** Writing a metric value (the refusal stands, and `tests/test_lentimpra.py:182` pins it). Running any arm. Editing the arms ledger — if a ledger field turns out to assert the metric blocker, that is surfaced, not rewritten, because arm records change only on their own executed evidence.

## Codebase Context

- `dna_decode/glm/lentimpra.py` (17,365 bytes, LF). `published_incumbents()` at line 304; docstring over-claim at lines 307–310; the `metric_resolution_attempt` block added in `b28fa8dc` with `decisive_next_check` at lines 361–364. `CELL_TYPES = ("HepG2", "K562", "WTC11")` at line 59; `S6_FILE = "Table S6 - folds and performance.xlsx"` at line 51.
- **The parser's consumption rule, read from the code** (this is the grounding for Step 1 and the design of Step 3): a row is consumed only if its first cell is in `CELL_TYPES`, or equals `test fold` (case-insensitive), or — with a block open and methods known — is numeric. Every other row is skipped with no record that it existed.
- `wiki/lentimpra_metric_resolution_2026-10-28.md` (3,965 bytes). The "blocked on provenance, not on compute" sentence is at line 63; the header-row instruction at line 61.
- `~/.claude/soraya-resume/dna_decode-77502f867d.json` (3,571 bytes, written 12:23 today) — the live resume checkpoint. **Outside the repo, so not version-controlled**, and the first thing a fresh session reads. Written only via `resume_state.write_checkpoint`.
- `tests/test_lentimpra.py` — 22 tests; `test_the_published_metric_is_declared_UNSPECIFIED_rather_than_guessed` asserts `inc["metric"] == "UNSPECIFIED_IN_SOURCE"` at line 182. Data-dependent tests sit behind a `needs_data` skipif.
- `wiki/lentimpra_acquisition_2026-10-27.md` — line 157 is the comparison-scoped wording that C2 contradicts; it also refers to "Table S6" in a way that invites the paper-vs-deposit conflation.

### Reusable-Code Survey

- **`dna_decode/data/ecoff_catalog.py`** — the provenance gate this refusal already follows (`source_url` + `verbatim_quote`, accessor raises rather than defaulting). Step 6 reuses its discipline for anything the publisher route returns.
- **`resume_state.write_checkpoint`** (`~/.claude/skills/soraya/scripts/resume_state.py`) — the only sanctioned writer for the checkpoint; Step 5 calls it rather than editing the JSON by hand.
- **`dna_decode/glm/lentimpra.py::_CACHE`** — the existing tuple-backed cache pattern; Step 3's new reader must not bypass or duplicate it.
- Searched: `graphify-out/GRAPH_REPORT.md` (absent), `dna_decode/glm/`, `dna_decode/data/`, `tests/`, `wiki/`, `~/.claude/skills/soraya/scripts/`.

## Pre-Change Baseline

Measured in-session, 2026-10-28, at `HEAD b28fa8dc`:

- `dna_decode/glm/lentimpra.py` = **17,365 bytes**, pure LF. Docstring asserts the sheet "never says" the metric; `decisive_next_check` names the header/units row.
- `wiki/lentimpra_metric_resolution_2026-10-28.md` = **3,965 bytes**; line 63 asserts the arm "stays blocked on provenance, not on compute".
- `~/.claude/soraya-resume/dna_decode-77502f867d.json` = **3,571 bytes**; its `next[]` item 4 repeats the header-row instruction and its `done[]` repeats the provenance framing.
- `uv run pytest tests/test_lentimpra.py -q` → **4 passed, 18 skipped** in 1.89 s (substrate detached, by design).
- `published_incumbents()` returns `metric="UNSPECIFIED_IN_SOURCE"`; **no** function exists that reports the sheet's unparsed rows (`grep` for `unparsed`/`skipped` in the module: no hits).
- Frozen five + `uv.lock`: byte-unchanged (`git diff --stat` empty for all six).

## Verification Signal

- **Nothing left asserting the wrong check**: a grep for `header/units row` and `header row` across `dna_decode/`, `wiki/lentimpra_metric_resolution_2026-10-28.md` and the checkpoint returns no hit that presents it as the resolving step.
- **Nothing left saying the eukaryotic arm is metric-blocked**: a grep for `blocked on provenance` returns nothing, and the arm's own ledger entry is read and reported (not edited).
- **The docstring's claim is scoped**: it names the parser's consumption rule rather than asserting what the sheet contains.
- **The next check is runnable, not prose**: `unparsed_rows()` exists, is covered by a synthetic-workbook test that needs no `D:`, and its test proves it RETURNS a planted caption row that `published_incumbents()` drops — the non-vacuity that distinguishes a real reader from a stub.
- **The refusal is untouched**: `test_the_published_metric_is_declared_UNSPECIFIED_rather_than_guessed` still passes and no `"metric": "Pearson"`/`"Spearman"` literal appears anywhere.
- Step 6 passes in either direction: table identity settled **with** a locator and quote, or recorded as attempted-and-open naming the routes tried.

## Implementation Steps

### Step 1: Scope the docstring's claim to what the parser actually reads
Files: dna_decode/glm/lentimpra.py
Depends on: none

**What changes:**
- Rewrite the `published_incumbents` docstring claim. Replace "the sheet ... never says whether they are Pearson, Spearman or something else" with the parser-scoped truth: the rows this function consumes are a cell-type marker, a `test fold` header carrying method names, and numeric fold rows — and **no row it reads states a coefficient**.
- Add one sentence stating the limit explicitly: rows outside that shape are skipped, so the sheet has NOT been shown to be silent, and a caption is exactly where a statement would sit.
- Keep the readout-type-mismatch rationale **adjacent and primary** — it is the reason the field refuses, and the scoping sentence must not displace it. A reader who meets the scoping first could conclude the refusal was premature and go guess, which is the failure this whole plan exists to prevent.

**Test strategy:**
- No behaviour change; assert the module still imports and `published_incumbents`' source no longer contains the string `never says`.
- Run the existing 22-test file to confirm 4 passed / 18 skipped is unchanged.

### Step 2: Replace the decisive-next-check text with the check that can work
Files: dna_decode/glm/lentimpra.py
Depends on: Step 1

**What changes:**
- Rewrite `decisive_next_check` in the `metric_resolution_attempt` block: the resolving read is the **rows the parser discards**, not the header row. Name why (the header is `test fold | <methods>`), and point at the `unparsed_rows()` helper Step 3 adds.
- Remove the "blocked today because the substrate volume (D:) is detached" framing as the *reason the check fails* — the volume is a prerequisite for running it, not the reason the stated check was wrong.
- **Do NOT add a `superseded_claim` field.** The attempt block already returns five fields to every caller, and a runtime payload that grows one field per correction is unbounded. Keep **one** pointer to the memo and let the memo carry the history — which makes it Step 4's obligation to record the superseded wording verbatim, or the correction becomes invisible instead of merely relocated.

**Test strategy:**
- Assert `header/units row` no longer appears in the module and `unparsed_rows` is referenced by the attempt block.
- Assert `metric` is still `UNSPECIFIED_IN_SOURCE` and no coefficient literal was introduced.

### Step 3: Make the real check runnable instead of a prose instruction
Files: dna_decode/glm/lentimpra.py, tests/test_lentimpra.py
Depends on: Step 2

**What changes:**
- Add `unparsed_rows(data_dir=DEFAULT_DIR, sheet="Benchmarking results") -> list[dict]` returning, for every row `published_incumbents()` does **not** consume, its 1-based row index and its non-empty cell values.
- **This REWRITES `published_incumbents()`'s row loop, and that is the step's main risk.** The consumption rule must be shared, not restated — but the predicates are **stateful**: `blocks, cur, methods = {}, None, None` is initialised before the loop, the `CELL_TYPES` branch assigns `cur` and resets `methods`, the `test fold` branch assigns `methods`, and the numeric branch is guarded by `if cur and methods and isinstance(cells[0], (int, float))`. Row classification is therefore **not a pure function of the row**, so sharing it means extracting a classifier that takes `(first, cells, cur, methods)` and returns the row kind, then rewiring the existing loop to call it. The carried state is passed IN, never re-derived inside the classifier.
- Raise the module's existing `LentiMpraDataError` when the workbook is absent, matching the sibling accessor rather than returning an empty list (an empty list would read as "no unparsed rows" when the real state is "could not look").

**Test strategy:**
- **Equivalence proof FIRST, because this rewrites a function pinned at `tests/test_lentimpra.py:182`:** capture `published_incumbents()`'s output on a recorded fixture before the refactor and assert it is identical after. Without this, the extraction is an unguarded rewrite of the function that decides what the incumbents are.
- **Source-shape assertion:** the classification predicates appear in exactly ONE place. Value equality cannot see a reintroduced duplicate, so assert on the source — `published_incumbents` contains no `in CELL_TYPES` / `== "test fold"` test of its own once the classifier exists.
- Synthetic workbook via `openpyxl`, no `D:` needed: a sheet containing a planted caption row, a `test fold` header, a cell-type marker and numeric fold rows.
- **Non-vacuity, the load-bearing test:** `published_incumbents()` drops the planted caption while `unparsed_rows()` RETURNS it. A reader that returns everything, or nothing, fails this pair.
- Absent-workbook test asserts `LentiMpraDataError`, not an empty list.

### Step 4: Correct the memo's two wrong claims in place
Files: wiki/lentimpra_metric_resolution_2026-10-28.md
Depends on: none

**What changes:**
- Replace the "stays blocked on provenance, not on compute" sentence: the arm is blocked on compute and a detached volume like the others; the metric blocks **comparison against the published incumbents**, which is a different and smaller claim. Quote the acquisition memo's own comparison-scoped wording as the evidence.
- Replace the header-row instruction with the unparsed-rows check and name `unparsed_rows()`.
- Add a short CORRECTIONS section at the top (per the front-load-by-reader-action convention) naming both corrections and the commit they correct, so a reader who stops after one screen gets the retraction rather than the superseded version.
- Soften "suggestive of Pearson": the two coefficients are often close on standardized, near-linear data, so a two-decimal match discriminates poorly between them. State this **qualitatively and without a numeric interval** — an unsourced number inside a memo whose argument is that unsourced values must be refused would undercut the memo. This strengthens the refusal while weakening the stated reason.
- Record the superseded `decisive_next_check` and "blocked on provenance" wording verbatim here (the obligation M1 moved out of the return value), so the correction is legible rather than a silent overwrite.

**Test strategy:**
- Manual read-through; assert the file no longer contains `blocked on provenance` or the header-row instruction, and that the CORRECTIONS section precedes the findings.

### Step 5: Correct the live resume checkpoint, which the next session reads first
Files: (none in repo — writes ~/.claude/soraya-resume/dna_decode-77502f867d.json)
Depends on: none

**What changes:**
- Re-call `resume_state.write_checkpoint` with the same goal/done/facts, amending the two wrong items: the provenance framing in `done[]` and the header-row instruction in `next[]` item 4.
- Add a `facts[]` entry recording the parser's consumption rule, so the next session does not have to re-derive why the header read is the wrong one.
- Keep every other item verbatim — this is a correction, not a rewrite, and the Kaggle items are still live.
- **This step must not run concurrently with another execution of any plan.** `~/.claude/soraya-resume/` is OUTSIDE the repo, so it gets no overlap protection from `Files:` and a git worktree gives it zero isolation — two concurrent runs would clobber one checkpoint, and the checkpoint is the artifact a fresh session reads first.

**Test strategy:**
- Re-read via `resume_state.render(cwd)` and confirm the corrected strings are present and the Kaggle `next[]` items survived.
- Confirm the file is still a single JSON per cwd (not duplicated) and that `write_checkpoint` overwrote rather than appended.

### Step 6: Test the version-artifact hypothesis on the publisher route
Files: wiki/lentimpra_metric_resolution_2026-10-28.md
Depends on: Step 4

**What changes:**
- Fetch the publisher's supplementary-file descriptions for the article and check whether a file named for folds/performance is listed, and under what number — the direct test of whether the table-identity mismatch is a **version artifact** (author-accepted vs published numbering) rather than two different tables.
- **If settled:** record the finding with `source_url` plus a short verbatim locator, following the ECOFF gate's discipline, and downgrade or retire the mismatch finding accordingly.
- **If not settled:** record the route as attempted-and-open, naming exactly what was fetched and what it did not say. A failed read is evidence about a route, not a property of the source — and this is now the third route, after the search index and the author-version text.
- Either way, demote "table identity" from the load-bearing finding to one of three, since the paper-reports-both finding is independent of it and survives regardless.

**Test strategy:**
- Manual: the memo states which routes were tried and what each returned; no claim about table identity is left unlabelled as to its evidence.
- If a value or locator is recorded, assert it carries both a URL and a quote — a bare assertion is refused.

### Step 7: Report which surfaces say the arm is metric-blocked
Files: soraya_runs/arms-ledger.json, dna_decode/eval/regime.py
Depends on: none

**What changes:**
- **Read-only by construction — this step writes nothing.** Grep the arms ledger, `regime.py` and `wiki/` for any field or sentence asserting that `substrate-change-constructed-eukaryotic-or-human-panel-UNRUN` is blocked on the published metric, and emit the hit list with its searched paths.
- Report the arm's actual recorded blocker **verbatim**, so the VOI ranking rests on what the ledger says rather than on a session's recollection.
- If a LEDGER field asserts the metric blocker, **surface it and stop there**: arm records change on their own executed evidence, and rewriting one to match a reasoning correction is the move the redundancy refusal exists to prevent. That case ends in a report, never an edit.
- The two files are named in `Files:` because this step READS them; it is the correcting step below that declares what it may write.

**Test strategy:**
- The grep output is the artifact. A zero-hit result is reported as zero-hit, not as confirmation, and the searched paths are named so the scan's extent is auditable.
- Assert the arms ledger is byte-unchanged after this step — the read-only claim is checked, not asserted.

### Step 8: Correct the prose surfaces the report named
Files: wiki/lentimpra_acquisition_2026-10-27.md, wiki/lentimpra_metric_resolution_2026-10-28.md
Depends on: Step 6, Step 7

**What changes:**
- For each PROSE surface Step 7 named, correct the claim that the metric blocks the arm: it blocks **comparison against the published incumbents**, not running the arm.
- `Files:` declares the two memos that can plausibly carry such wording. **If Step 7 names a prose file outside this list, stop and surface it** rather than editing an undeclared file — the `Files:` line is what the executor's wave and overlap logic reads, so silently widening it defeats the contract.
- If Step 7 found no prose surface, this step is a no-op and says so. A no-op is the expected outcome and must not be padded into an edit.

**Test strategy:**
- Re-run Step 7's grep after the edit; the asserting sentences are gone and no new file was touched.
- `git diff --stat` names only files from this step's `Files:` line.

## Execution Preview

Waves (topological, from the `Depends on` graph):

- **Wave 0:** Step 1, Step 4, Step 5, Step 7 — four-way, no file overlap (module / memo / out-of-repo checkpoint / read-only report).
- **Wave 1:** Step 2, Step 6
- **Wave 2:** Step 3, Step 8

Total waves: 3. Max parallelism: 4. Critical path: Step 1 → 2 → 3.

Recomputed by `plan-toolkit.py waves` after the 7-to-8-step split, not by hand — a hand-written wave table is a second expression of the dependency graph. **Step 8 depends on Step 6 as well as Step 7** because the split first produced a real intra-wave collision: both write `wiki/lentimpra_metric_resolution_2026-10-28.md`, which the toolkit's `overlaps` field caught before it could ship.

Sequential execution is recommended despite the computed 4-way parallelism: Steps 1, 2 and 3 all edit `dna_decode/glm/lentimpra.py`, Steps 4, 6 and 8 share the memo, **Step 5 writes a shared path outside the repo that a worktree cannot isolate**, and the whole change is under an hour — a worktree's isolation buys nothing here and costs a `uv.lock` reconciliation.

## Risk Flags

- **This plan corrects a commit that is already pushed** (`b28fa8dc`). The corrections are additive text changes, so no rewrite of history; the superseded wording is recorded in **Step 4's memo** rather than silently replaced — relocated there from a planned `superseded_claim` return-value field, which was dropped because a runtime payload that grows one field per correction is unbounded.
- **`[unverified]` — the publisher route has not been fetched.** No claim about whether the table mismatch is a version artifact appears anywhere in this plan; Step 6 is that test, and it is allowed to fail.
- **Step 5 writes outside the repo** and is therefore not covered by git. If it is skipped, the next session reads the uncorrected checkpoint first and the other six steps are largely wasted — this is the highest-leverage step and the easiest to drop.
- **Step 3 risks a second copy of the consumption rule.** If the predicates are restated instead of shared, the two expressions drift and `unparsed_rows()` silently stops matching what `published_incumbents()` skips. The step requires sharing, not restating.
- **Step 7 may find nothing**, which is a clean result and must be reported as such. The named paths bound the scan; a zero-hit grep over the wrong paths is the failure shape to avoid.
- **The review behind these findings was single-voice** (`codex` absent), so C1/C2 rest on my own reading of the parser and the acquisition memo — both quoted in Codebase Context so a reader can check them rather than take them.
- Architecture gate: `sentrux` present; one module plus docs, no boundary crossed — a structural compare is cheap but unlikely to be informative.

## Open Questions

- **Should `unparsed_rows()` ship as a module function or a script?** A module function is testable and importable; a script is more discoverable for a future session that only reads `scripts/`. Defaulting to the module function with the attempt block pointing at it, but discoverability is a real argument the other way.
- **Is the memo's `verdict` string still right?** `UNRESOLVED_TABLE_IDENTITY_UNVERIFIED` names only the table question, while the paper-reports-both finding is independent and would survive the table question being settled. `UNRESOLVED_NO_SHEET_LEVEL_STATEMENT` is truer. Changing it is cosmetic today and would matter the moment anything branches on the string.
- **Does the eukaryotic arm's re-ranking follow from C2 in this plan, or separately?** C2 says the metric is not that arm's blocker, which raises its rank. This plan deliberately only corrects the claim; actually re-ranking and running it is separate work and needs the volume attached.
- Whether to pre-commit now to the outcome where the caption check finds nothing — i.e. agree in advance that the metric is then unanswerable from our artifacts, rather than discovering that twice.

## Verification

1. `uv run pytest tests/test_lentimpra.py -q` — counts reported, not inferred; the 4 passed / 18 skipped baseline must not regress, and Step 3 adds passing synthetic-workbook tests that need no `D:`.
2. The non-vacuity pair holds: `published_incumbents()` drops the planted caption, `unparsed_rows()` returns it.
3. `grep -rn "never says\|header/units row\|blocked on provenance"` over `dna_decode/`, `wiki/lentimpra_metric_resolution_2026-10-28.md` and the rendered checkpoint returns no hit presenting the superseded claims as current.
4. `test_the_published_metric_is_declared_UNSPECIFIED_rather_than_guessed` still passes; no `"metric": "Pearson"` or `"metric": "Spearman"` literal exists.
5. `git diff --stat` empty for the frozen five and `uv.lock` after every step.
6. EOL discipline: `dna_decode/glm/lentimpra.py` and `tests/test_lentimpra.py` stay pure LF; each patch asserts its own file's ending unchanged rather than absolute.
7. `resume_state.render(cwd)` shows the corrected items and the surviving Kaggle items.
8. The arms ledger is byte-unchanged, asserted inside Step 7 rather than claimed by it; Step 8 touches only files named in its own `Files:` line, and a prose surface outside that list stops the step instead of widening it.

<!-- toolkit: check=clean waves=clean gate=fired:open-questions,unverified -->
