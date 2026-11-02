# Clear the Three Blockers the Frame-Conditioning Review Found — Technical Plan

> One of the three is live code and the other two are amendments to an unexecuted plan, and the first one is larger than the review estimated because routing a loader through its owner means satisfying the owner's stricter column contract in seven test fixtures.

## Lens status
Inputs: conversation (the review's recommended-next-step directive), grounded in `plans/Frame_Conditioning_Does_Corpus_Window_Position_Explain_PromoterAtlas_0_5076_Plan/brainstorm.md` and `wiki/frame_plan_feasibility_2026-11-01.json`, both committed.
Degradations: the repo-index helper script is absent on this host, so all understanding comes from reading files directly rather than from generated summaries; `project-rules.md` is absent, so the Project Rules Check section is omitted; `DESIGN.md` is absent and this plan touches no user interface; the graphify dependency-graph report is absent, so the reuse survey was built by direct search; `codex` is absent; `sentrux` and `gh` resolve.

## Problem Statement

A pre-execution review of the frame-conditioning plan found three blocking issues. One is a live defect in committed code; the other two are defects in a plan document that has not been executed, so fixing them costs nothing but attention. This plan clears all three, and corrects two claims the review itself made.

**The first blocker is live and is a regression.** `tests/test_substrate_loader_ownership.py` is red because `dna_decode/glm/tile_controls.py` opens the GSE144621 peak-tile data file itself rather than calling `dna_decode/glm/genomewide.py`, which that test registers as the file's sole owner. The guard exists because the sibling fragment files in the same data series are separated by spaces while the peak-tile file is separated by tabs, so a second loader that does not inherit the owner's delimiter-sniffing can parse one of them into zero rows — a result that reads as "no signal" rather than "no data". The control loader hardcodes a tab split, which is correct for this one file and wrong by exactly the mechanism the guard names.

**The review's own sizing of that fix was wrong, and this plan corrects it.** The review called it roughly a ten-line change. The module change is indeed small, but the owner enforces a minimum of eight columns and requires `variant`, `expn_med`, `peak_start`, `peak_end`, `strand`, `RNA_exp_sum_1` and `RNA_exp_sum_2`, whereas the control loader's synthetic test fixtures supply five columns and lack three of those. Seven of the nine tests in `tests/test_tile_controls.py` use those fixtures, so all seven break until the fixtures are rebuilt to the owner's contract. That rebuild is mechanical but it is the bulk of the work.

**The review's claim about the right remedy was also half wrong.** It stated that adding a row to the guard's `ALLOWED_MENTIONS` list would be the wrong fix. That is true only of adding such a row *instead of* deleting the parser, which would exempt a genuine second reader. It is false as a blanket claim, because the guard matches the bare text pattern `\bGSE\d+` across the whole file including docstrings, and the control loader names the accession in a module docstring sentence that will survive the fix. Clearing the guard therefore requires both acts in order: delete the parser so the module genuinely stops reading the file, then add a row for the surviving prose mention — which is precisely what that list is for, since its own comment describes each row as "a module that TALKS ABOUT a substrate it does not read" and calls the row "the intended friction".

**The second blocker is a contradiction in the frame plan's Step 1.** That step says to reproduce the upstream corpus-extraction function exactly. The upstream function threads a running `last_end` variable through gene features in file order, and the annotation file's coding-sequence rows are measurably not in ascending coordinate order on either strand. The feasibility probe sorts before threading, so the three figures the plan pins as acceptance values come from a replica that diverges from the thing it claims to reproduce. The divergence is currently asserted rather than quantified, which this plan fixes by measuring both orderings.

**The third blocker is that the frame plan's two additive measurements lack the control its subtractive predecessor had.** Adding an informative feature to a baseline strengthens that baseline, which shrinks the encoder's margin over it whether or not the encoder's advantage is positional, so margin shrinkage alone is not evidence of absorption. Separately, the two strata the plan would compare within are not exchangeable: the in-window interquartile range is 0.4735 against the out-of-window 0.2710, while the standard deviation inverts at 3.380 against 5.535, so a rank correlation computed inside the narrower stratum is attenuated for both channels and a smaller margin there is expected from range restriction alone.

**Non-goal.** This plan does not run the frame-conditioning experiment and does not freeze its acceptance bar. It makes the plan safe to execute and clears the red guard. Nothing here touches the measured headline metric or the benchmark's ceiling verdict.

## Codebase Context

**The guard and what it will and will not accept.** `tests/test_substrate_loader_ownership.py` compiles the pattern `\bGSE\d+` at line 46 and applies it to the full text of every Python file under `dna_decode/`, so a docstring mention counts as a match exactly as an import does. `OWNERS` registers two accessions, GSE108535 owned by `glm/expression.py` and GSE144621 owned by `glm/genomewide.py`, each with the specific trap a second reader would re-derive. `ALLOWED_MENTIONS` currently exempts `glm/genomewide.py` for GSE108535, and `glm/encoder.py` together with `glm/lentimpra.py` for GSE144621 — so a prose-only exemption row has clear precedent and is already used twice for this accession. The file's own standard for adding a row is that the hits were verified by reading them rather than assumed.

**What the owner requires of a file it parses, and why that reaches the fixtures.** `genomewide.py` sets `MIN_COLUMNS = 8` and `load_peak_tiles` requires the seven named columns listed in the Problem Statement. It resolves its own path through `ensure_cached`, which joins a cache *directory* with the file name, so it takes a directory and not a file path. The current `load_control_tiles` accepts either, which means the public signature narrows slightly under this change.

**What the owner already returns, measured rather than assumed.** Calling `load_peak_tiles` on the cached directory returns all 46,713 rows across four categories, of which 936 `random` and 470 `neg_control` survive at 150 base pairs. Those are byte-identical to the counts the second loader produced, which is why this fix dissolves the violation instead of trading accuracy for compliance. The owner also applies its own filters for non-ACGT characters and unparseable numbers, so the wrapper inherits them rather than repeating them.

**What the control loader keeps and what it loses.** The module is 138 lines. The parts that survive are the `ControlSeq` dataclass, the `CONTROL_CATEGORIES` and `EXPECTED_AT_150` constants, the `ControlCohortError` exception, the zero-rows refusal, and `control_counts`. The parts that go are the `gzip` import, the `PEAK_TILE_FILE` constant that names the accession in code, and the private `_rows` parser. `ControlSeq` is deliberately retained rather than replaced by the owner's `Tile` type, because `Tile` carries a `peak` field that is the split unit for tiles and meaningless for controls, and because the existing tests assert on `ControlSeq`'s fields.

**The ordering question, measured.** `parse_gff3` returns 9,523 rows from the cached annotation, of which 4,340 are coding sequences on the single contig `NC_000913.3`, split 2,118 on the plus strand and 2,222 on the minus, and every one of them carries a `locus_tag` — which matters because the upstream extractor filters on that qualifier's presence. Checked per strand in the order the parser returns them, the start coordinates are **not** ascending on either strand. The upstream code also reads GenBank files rather than this GFF3, so the file order it relies on is a different file's order; faithful replication is not merely unverified here but not well-defined.

**Two modules must not be edited.** `dna_decode/glm/genomewide.py` is embedded as a compressed payload inside `scripts/kaggle/build_promoteratlas_tile_arm_kernel.py`, and `scripts/promoteratlas_arm_dryrun.py` asserts the embedded copy is byte-identical to the repository copy, so an edit silently invalidates that check. `dna_decode/glm/expression.py` holds the ridge-regression protocol every comparison on this benchmark is anchored to. Both are read-only throughout, and both were verified byte-unchanged at the end of the preceding session.

**Test infrastructure.** pytest, configured through `pyproject.toml`, invoked in this project as `uv run pytest`. Tests that need data which is excluded from version control skip rather than fail, which is the pattern `tests/test_tile_controls.py` already uses for its single real-data assertion. `tests/test_contract_number_audit.py` uses the `wiki/` directory as a pool of decoy numbers, so any new artifact with a large quantity of decimal values can break it; a 102,755-value dump did exactly that on 2026-11-01.

### Reusable-Code Survey
- `dna_decode/glm/genomewide.py` provides `load_peak_tiles`, which is the owner this plan routes through, plus `gc_features` and `spearman`; read-only reuse
- `dna_decode/glm/genomewide.py` also provides the `Tile` dataclass whose fields `ControlSeq` is mapped from, so no new parsing logic is written anywhere
- `scripts/frame_plan_feasibility_probe.py` already contains the corpus-window replica and the locate sweep, so the ordering comparison extends an existing script rather than starting a new one
- `scripts/promoteratlas_subspace_removal.py` provides the tie-averaging `auroc` implementation and the seeded random-control pattern that the two missing controls follow
- `tests/test_substrate_loader_ownership.py` provides the `ALLOWED_MENTIONS` mechanism and its documented standard for adding a row
- None additionally — searched: `graphify-out/GRAPH_REPORT.md` (absent), `src/`, `lib/`, `utils/`, `helpers/` (none exist at the repository root), `dna_decode/glm/`, `dna_decode/data/`, `scripts/`

## Pre-Change Baseline

Every figure below was measured during this planning session or read from a committed artifact.

- The guard file `tests/test_substrate_loader_ownership.py` reports **1 failed, 7 passed**. The failure names `glm/tile_controls.py` alongside the owner `glm/genomewide.py` for accession GSE144621.
- The full suite reports **14 failed, 7227 passed, 180 skipped** in 797.9 seconds. Thirteen failures are the recorded absent-data cases; the fourteenth is the guard above.
- `dna_decode/glm/tile_controls.py` is **138 lines**. `tests/test_tile_controls.py` holds **9 tests**, of which **7** build the five-column synthetic fixture.
- `ALLOWED_MENTIONS["GSE144621"]` is currently `{glm/encoder.py, glm/lentimpra.py}`; `OWNERS` covers GSE108535 and GSE144621.
- The owner requires a minimum of **8** columns and the seven named peak-tile columns, and returns **46,713** rows including **936 random and 470 neg_control at 150 base pairs**.
- Coding-sequence rows in the cached annotation are **not** ascending by start coordinate on either strand; there are 4,340 such rows, all carrying a `locus_tag`.
- The frame plan pins, all derived under sorted ordering: **3,496** corpus windows covering **625,255** bases, which is **13.47%** of the genome; **43,061** of 43,417 distinct tile sequences uniquely located, a rate of **0.9918**, with 325 ambiguous and 31 not found; **28,622** tiles overlapping a window and **8,239** entirely inside one.
- Per-stratum expression spread: in-window interquartile range **0.4735** with standard deviation **3.380**; out-of-window **0.2710** and **5.535**.
- The inherited decision margin is **0.02**, and the committed metrics it is compared against are **0.5076** for the encoder and **0.3112** for the GC baseline.

## Verification Signal

This plan succeeds when three things hold simultaneously. The ownership guard passes without the parser being exempted, which is demonstrated by the guard failing after the parser is deleted and passing only after the prose row is added. The frame plan no longer claims exact replication of the upstream extractor, and the magnitude of the divergence it does carry is a measured number in an artifact rather than an assertion. And the frame plan's acceptance bar specification names both missing controls, so that freezing it later cannot produce a verdict graded without them.

Failure is equally well defined. If the owner-routed loader returns control counts other than 936 and 470, the routing changed the data and the change is rejected rather than re-pinned. If the two orderings turn out to produce identical window sets, the second blocker was not real and the plan should say so instead of adding a parameter nobody needs.

## Implementation Steps

### Step 1: Route the control loader through the substrate's owner
Files: dna_decode/glm/tile_controls.py, tests/test_tile_controls.py
Depends on: none

**What changes:**
- `load_control_tiles` calls `genomewide.load_peak_tiles(cache_dir, allow_download=False)` and filters its result by category and sequence length, instead of opening the file. The `gzip` import, the `PEAK_TILE_FILE` constant and the private `_rows` parser are deleted, which removes every code-level naming of the accession.
- The signature narrows to a cache **directory**, because the owner resolves its own path. The docstring states this and says why, so the narrowing is not discovered by a caller at runtime.
- `GenomeWideDataError` from the owner is translated into `ControlCohortError`, keeping this module's public failure type stable for its existing callers.
- `ControlSeq` is retained and mapped from the owner's `Tile`, rather than returning `Tile` directly, because `Tile` carries a `peak` field that is the split unit for tiles and has no meaning for a control sequence.
- The zero-rows refusal and the geometry-incompatibility message are preserved verbatim, since they encode the measured fact that no `pos_control` row survives at 150 base pairs and that padding was measured to contaminate this encoder.
- The seven synthetic fixtures are rebuilt to satisfy the owner's contract: at least eight columns, including `peak_start`, `peak_end` and `strand`. This is the bulk of the change and the reason the review's ten-line estimate was wrong.

**Test strategy:**
- All nine existing tests must pass unchanged in intent, with only their fixtures widened. The behavioural assertions are the contract and must not be relaxed to accommodate the new path.
- The real-data test continues to assert exactly 936 `random` and 470 `neg_control`. If routing through the owner changes either number, the routing changed the data and must be rejected rather than re-pinned.
- A new test parses the module with `ast.parse` and walks it, asserting that no import node names `gzip` and that no attribute access targets it. This follows the convention already used in `tests/test_amr_uncounted.py` lines 145 to 158 and seven other test files. **A substring search was specifically rejected and must not be reinstated:** the module docstring at line 12 contains the phrase `embeds genomewide.py as base64+gzip`, and that sentence survives the fix it describes, so a text search would fail while the code is correct. A text search over source cannot distinguish an instruction from a comment, which is the same weakness the ownership guard's own `\bGSE\d+` pattern has.
- A new test asserts that passing a file path rather than a directory now raises, pinning the deliberate signature narrowing instead of leaving it as undocumented drift.
Test-strategy-ignore: tests/test_amr_uncounted.py

### Step 2: Admit the surviving prose mention, and prove the admission was necessary
Files: tests/test_substrate_loader_ownership.py
Depends on: Step 1

**What changes:**
- `glm/tile_controls.py` is added to `ALLOWED_MENTIONS["GSE144621"]`, joining `glm/encoder.py` and `glm/lentimpra.py`, with a comment recording what was read: the module docstring explains that the owner parses the peak-tile file and why this module is separate from it, and after Step 1 the module opens no file and names no accession outside that prose.
- The comment follows the file's own stated standard, which is that the hits were verified by reading them rather than assumed, and it records the count of surviving mentions so a future increase is visible in a diff.
- Nothing else in the file changes. `OWNERS` is untouched, and no other accession gains a row.

**Test strategy:**
- The ordering is the test. Before this step the guard must still FAIL on the prose mention alone, which proves Step 1's deletion did not clear it and that this row is doing real work rather than papering over a live reader. After this step the guard must PASS. Both states are recorded.
- The full guard file must report 8 passed, up from 7 passed and 1 failed.
- Deliberately not done: no mechanism is added to distinguish prose from code automatically. The guard's visible one-line exemption is its design, and replacing it with inference would remove the friction that makes the admission reviewable.

### Step 3: Quantify the window-ordering divergence instead of asserting it
Files: scripts/frame_window_ordering_probe.py, wiki/frame_plan_window_ordering_2026-11-01.json
Depends on: none

**What changes:**
- A NEW script builds the corpus-window set twice, once threading `last_end` through coding-sequence rows in the order the parser returns them and once in sorted ascending order, and reports both window counts, both covered-base totals, the size of the symmetric difference between the two window sets, and how many uniquely-located tiles change their in-window status.
- **A new script rather than an extension of `scripts/frame_plan_feasibility_probe.py`, and the reason is mechanical rather than stylistic.** That script fixes its output path as a module constant at line 38, includes a `wall_seconds` field at line 262, and writes the file unconditionally at line 264. So any run of a modified version rewrites `wiki/frame_plan_feasibility_2026-11-01.json`, and the timing field guarantees the rewrite differs even when every input is identical. That artifact is committed and cited by both the frame plan and its review, so the original plan's promise to leave it untouched was false as written. A new script makes the promise true by construction, where a flag on the existing script would rely on every future caller remembering to pass it.
- The accepted cost is duplication: both scripts carry a corpus-window replica, and two copies can drift into producing different window sets under the same name. If that drift becomes real, the fix is to extract the replica into the package rather than to re-merge the scripts.

**Test strategy:**
- A unit test on a small synthetic annotation deliberately shuffled out of coordinate order, asserting that the two orderings produce different window sets. If they do not differ on a case constructed to make them differ, the comparison is not measuring what it claims.
- The real-data run is the deliverable. A symmetric difference of zero is a legitimate and informative outcome, and would mean the second blocker is not real; it must be reported as such rather than treated as a failed run.
Test-strategy-ignore: wiki/frame_plan_feasibility_2026-11-01.json

### Step 4: Amend the frame plan for the ordering contradiction and the two missing controls
Files: plans/Frame_Conditioning_Does_Corpus_Window_Position_Explain_PromoterAtlas_0_5076_Plan/technical-plan.md, wiki/plans-index.md
Depends on: Step 3

**What changes:**
- Step 1 of that plan stops claiming exact reproduction of the upstream extractor. It instead makes the row ordering an explicit named parameter defaulting to sorted ascending, states that this diverges from the upstream reliance on file order, cites the Step 3 artifact for the size of that divergence, and scopes the resulting claim to the corpus rule as specified at 200 base pairs with a 5 base-pair intergenic gate rather than as the authors ran it.
- The three pinned figures are re-labelled with the ordering they were derived under, so a future reader cannot mistake them for order-independent facts.
- **`wiki/plans-index.md` line 869 is rewritten in the same step, because it restates the claim this step removes.** That key-decision bullet currently asserts that Step 1 *replicates the authors' `extract_promoter_regions` exactly, INCLUDING two apparent defects*, which becomes false the moment the body stops claiming exact reproduction. The index is what a reader scans to see what a plan decided without opening it, so a false bullet there outlives a false sentence in the body. This is the standing rule that a correction is finished only when every surface restating the claim carries it.
- Step 5 of that plan, which freezes the acceptance bar, gains both controls by name so that the bar cannot be frozen without them. The first is a matched-strength comparison: the encoder's margin is reported over a series of baselines of increasing standalone strength, including the frame, so that absorption is judged against how much any equally strong feature would shrink the margin rather than against the unstrengthened baseline. The second is a size-matched random stratification, splitting tiles into strata of the same sizes at random and measuring the within-stratum margins, which gives the null for how much stratification alone moves the margin.
- Step 6 of that plan gains the requirement that each stratum's expression spread is reported beside its margin, because the measured interquartile ranges differ by a factor of about 1.7 and the standard deviations invert, so a reader cannot otherwise tell range restriction from absorption.
- The review's two corrected figures replace the superseded ones: stated-coordinate agreement becomes 1193 of 43061, which is 2.77%, measured across every uniquely-located tile rather than the first two thousand rows, and the trivial-overlap share becomes 2,359 or 8.2%.

**Test strategy:**
- The plan toolkit's `check` and `waves` subcommands must both remain clean on the amended document, and the wave grouping must be recomputed rather than hand-reasoned, because an amendment that adds a dependency changes the graph.
- The amended document must still carry its toolkit footer as its final line, since consumers read that from the file tail.
- Manual review that the amendment adds scope and does not remove any: the three cells and seven named outcomes all survive, and only the ordering claim is weakened.

### Step 5: Record what changed, and re-verify the surfaces that restate it
Files: wiki/frame_plan_amendments_2026-11-01.md, project_state/glm-learned-representation-2026-10-07.md, plans/Frame_Conditioning_Does_Corpus_Window_Position_Explain_PromoterAtlas_0_5076_Plan/brainstorm.md
Depends on: Step 1, Step 2, Step 4

**What changes:**
- A short memo records the three blockers as cleared, the two corrections this plan makes to the review itself, and the measured ordering divergence. It is ordered with decisions first and corrections second.
- The project ledger gains one Action-Log row through the sanctioned append path, using a class from the seven-value enumeration, recording that the loader regression is fixed and that the review's sizing and remedy claims were both corrected.
- The review document gains a dated amendment noting that its "roughly a ten-line change" estimate and its "an ALLOWED_MENTIONS entry would be the wrong fix" claim were both corrected during planning, without deleting either original statement.
- The three pre-existing untracked files from earlier sessions are left untouched, as they were in the preceding commits.

**Test strategy:**
- Run the ownership guard, the control-loader tests, and the contract-number audit, and assert the audit still reports nothing to adjudicate and that the `wiki/` decoy pool has not gained a large quantity of decimal values.
- Assert the five frozen decoder files, the dependency lockfile, `genomewide.py` and `expression.py` all show an empty diff.
- Run the affected-file subset rather than the full suite, and state the before and after counts for the guard specifically, since that is the number this plan exists to move.
Test-strategy-ignore: scripts/contract_number_audit.py, tests/test_substrate_loader_ownership.py, tests/test_tile_controls.py

## Execution Preview

Wave 0 (2 parallel):  Step 1 — route the loader through its owner, Step 3 — quantify the ordering divergence
Wave 1 (2 parallel):  Step 2 — admit the prose mention, Step 4 — amend the frame plan
Wave 2 (1):           Step 5 — record and re-verify

Critical path: Step 1 → Step 2 → Step 5 (3 waves)
Max parallelism: 2 agents

Note: Parallel execution requires a git repository with a configured remote. If unavailable, /execute-plan falls back to sequential mode. For this plan sequential is likely the right mode regardless, because two of the five steps need data files that are excluded from version control and therefore absent from a fresh worktree.

## Risk Flags

- Severity: high — **Deleting the parser does not by itself clear the guard, and assuming it would leaves the work looking done while the suite stays red.** The guard matches the bare text pattern `\bGSE\d+` across whole files, and a docstring sentence in the control loader names the accession. Steps 1 and 2 are therefore ordered and the ordering is itself asserted: the guard must fail between them.
- Severity: high — **Routing through the owner breaks seven of nine existing tests via their fixtures, not via their assertions.** The temptation is to relax an assertion to make a fixture pass. The fixtures are what must change, and the test strategy says so explicitly.
- Severity: medium — **The signature narrows from accepting either a file or a directory to accepting only a directory.** Any existing caller passing a file path breaks. The only current callers are this plan's own scripts and tests, which is checkable by search, but the narrowing is pinned by a new test rather than left implicit.
- Severity: medium — **`genomewide.py` is embedded as a compressed payload in a Kaggle kernel that asserts byte-identity against it.** This plan routes *through* that module and must not edit it. Verified by an empty diff at the end.
- Severity: medium — **Amending a plan document after a review could quietly narrow its scope.** The amendment is specified as additive: the three cells and seven outcomes survive and only the ordering claim is weakened. A manual scope check is in Step 4's test strategy.
- Severity: medium — **A zero symmetric difference in Step 3 would mean the second blocker was not real.** That outcome is pre-committed as reportable rather than as a failed run, so there is no incentive to tune the comparison until it finds something.
- Severity: low — **New artifacts enter the decoy pool that the contract-number audit draws on.** The ordering artifact carries a small number of values, far below the scale that broke two guards on 2026-11-01, and Step 5 re-runs the audit.
- Severity: low — **External tool surfaces.** The plan toolkit's `check`, `waves` and `validate` subcommands were confirmed by direct invocation during this session, and on this host the `python3` name resolves to a Microsoft Store placeholder so `python` must be used. `sentrux` resolves but no architecture gate is proposed here, because this plan adds no module and changes no dependency edge. No step contacts a network or service, so no external schema is left unconfirmed.
- Severity: low — **No restructuring was applied during wave grouping.** Steps 4 and 5 both touch files that Step 4 writes, which is why Step 5 depends on it; Steps 1 and 3 touch disjoint paths and so do Steps 2 and 4. The one ordering choice worth naming is that Step 2 was kept separate from Step 1 rather than merged, specifically so that the exemption row is a distinct reviewable act and the guard's intermediate red state can be asserted.

## Open Questions

- **Should the control loader's docstring stop naming the accession instead of taking an exemption row?** Rewording it to say "the peak-tile substrate" would clear the guard with no exemption at all, at the cost of precision about which substrate is meant. This plan takes the exemption row because the guard's own comment describes that row as the intended friction and a reviewable one-line admission, but the alternative is cheap and defensible and the choice is the user's.
- **Does the frame plan's Step 1 still need a corrected-ordering variant at all?** If Step 3 shows the two orderings differ only trivially, the honest move may be to drop the parameter and state the divergence as negligible rather than carry configuration nobody will use. That cannot be decided before Step 3 reports.
- **Is a matched-strength baseline series the right control, or should it be a single feature tuned to match the frame's standalone strength?** The series is more informative and avoids the problem of constructing an exactly matched feature, but it is also more work and produces a curve rather than a single comparison. This plan specifies the series.

## Architecture Impact

**Expected structural changes:**
- No module is added and none is deleted. One module, `dna_decode/glm/tile_controls.py`, shrinks and gains a dependency on `dna_decode/glm/genomewide.py`.
- That new dependency edge is the point of the change rather than a side effect: the whole purpose is to make the control loader a consumer of the substrate's owner instead of a second reader of the substrate.
- The edge direction is consistent with the existing package layout, since `genomewide.py` is already imported by `atlas_embed.py` and `encoder.py`. No module under `dna_decode/` gains a dependency on anything under `scripts/`.

**Risk areas:**
- `genomewide.py` is the shared abstraction this plan routes through, and it stays read-only for the Kaggle byte-identity reason recorded above.
- `expression.py` is untouched and gains no new caller.
- No cycle is possible: `tile_controls.py` becomes a leaf that imports `genomewide.py`, which imports nothing from it.

**Sentrux gates:**
- Baseline before the first step: not needed. A baseline already exists at `.sentrux/baseline.json` from the preceding session and this plan adds no module.
- Compare after each wave: no. The change is one import edge inside an existing package.
- Acceptable delta: the total import-edge count rising by one. A new cycle or a new edge into `genomewide.py` from outside `dna_decode/glm/` is not acceptable.
- Regression-response stance: investigate rather than accept or revert, because on a change this small a flagged regression most likely indicates an unintended import.
- Honest limitation, carried from the preceding session: `.sentrux/baseline.json` holds nine aggregate scalar values with no per-edge record, so a compare can report that the edge count moved but cannot identify which edge moved. The real protection for the read-only invariant is the anti-drift check below, not the gate.

**Anti-drift checks:**
- `dna_decode/glm/genomewide.py` and `dna_decode/glm/expression.py` show an empty diff at the end.
- `scripts/promoteratlas_arm_dryrun.py` still passes its embedded-payload byte-identity assertion.
- The set of files under `dna_decode/` carrying an actual import of `genomewide` is exactly `atlas_embed.py`, `encoder.py` and the newly-routed `tile_controls.py` — three, up from two, and the third is the intended addition.
- The five frozen decoder files and the dependency lockfile are byte-unchanged, and the active prospective lock re-verifies.

## Verification

1. The ownership guard moves from 1 failed and 7 passed to **8 passed**, and it is observed failing between Steps 1 and 2 so the exemption row is proven necessary rather than assumed.
2. The control loader's nine tests pass with widened fixtures and unrelaxed assertions, and the real-data test still asserts exactly 936 `random` and 470 `neg_control`.
3. The control loader contains no `gzip` usage and no accession name outside its module docstring.
4. Passing a file path where a directory is now required raises, pinning the signature narrowing.
5. The ordering comparison reports both window sets, their symmetric difference, and the number of tiles whose in-window status changes, into its own artifact; a zero difference is reported as a finding.
6. The amended frame plan passes the toolkit's `check` and `waves` cleanly, carries its footer as the final line, and names both controls in its bar-freezing step.
7. The frozen decoder files, the lockfile, `genomewide.py` and `expression.py` all show an empty diff, and the prospective lock re-verifies.
8. The contract-number audit still reports nothing to adjudicate, and the `wiki/` decoy pool has not gained a large quantity of decimal values.
9. The full suite's failure count drops from 14 to **13**, and the one that disappears is the ownership guard. No other failure count changes.

<!-- toolkit: check=clean waves=clean gate=fired:open-questions,severity-high -->
