# Contract-Number Provenance: Sound Decoy Control + Typed Provenance Kinds — Technical Plan

> Make the standing contract-number audit's clean verdict mean something: replace an underpowered random decoy sample with a deterministic full-pool control, then type the 34 unverifiable numbers by provenance kind instead of citing them blind.

## Lens status
Inputs: conversation (the banked follow-on from the 2026-09-27 `/soraya --advance` run: 34 unverifiable contract numbers in three named classes)
Degradations: sentrux unavailable (Architecture Impact gates n/a — section omitted, change is additive to one script plus contract prose, not structural); repo-index unconfigured (`update-index` CLI not installed) so Step 1 research was direct file reading; no `project-rules.md` (Project Rules Check omitted); no `DESIGN.md` (no UI surface in this plan)

## Problem Statement

`scripts/contract_number_audit.py` is the standing check that every number in a `cell_registry` contract traces to the artifact the contract cites. As of commit `42ffe1f` it reports `NOTHING_TO_ADJUDICATE 0/157 across 16 cells`, plus a coverage block naming 34 numbers in 10 cells as unverifiable and 2 cells as LOW discrimination.

Two problems, found while grounding this plan:

1. **The decoy control I shipped in `42ffe1f` is underpowered, so its grades are not trustworthy.** It samples 10 random unrelated artifacts with a fixed seed and asks how many contain all of a cell's numbers. Measured against the full 603-artifact pool the estimates are noisy in both directions, and one case is flatly wrong: `essentiality:any:essentiality`'s seven numbers match **7/7 against at least six completely unrelated artifacts** (a yeast benchmark, two TB baselines, two staleness snapshots, a pigment run) while the k=10 sample reported **0/10**. A `0/10` grade therefore does not mean "no decoy matches". Since the grade is the thing that says whether a citation verifies anything, an underpowered grade is worse than no grade.

2. **A second confound in the same function:** a cell citing N artifacts is scored against a haystack of N concatenated files while each decoy is a single file. Bigger haystack, more numbers, easier full match — so multi-citation cells get inflated discrimination for free. My own research probe hit this (a 33-file `klebsiella` bundle scored 6/6 own against single-file decoys).

Downstream of those, the 34 unverifiable numbers are **not one class and must not be treated as 34 defects.** Grounded in the contracts, they are at least four kinds:

- **Real measurements of our own cells, artifact exists** — `typing:Streptococcus_pneumoniae:pneumoserotype` (0.939 / 0.661 …), `hla:human:b5701` (0.979 / 0.61 …, artifact located at `wiki/hla_b5701_validation_2026-07-06.json`), `typing:klebsiella:kleb` (top-1 0.45 / top-5 0.60 / lift 0.49).
- **Enforced by test, not by artifact** — `typing:Escherichia_coli:pathotype` 0.833 / 1.0, pinned by exact-equality asserts in `tests/test_pathotype_expec_recall.py`. It has no wiki artifact by design, and inventing one would be worse than naming the test.
- **External reference values, not our measurements** — `pgx:human:nudt15` `9.5` is a `*3` EAS allele frequency from a population database. No artifact of ours can or should contain it.
- **Structural non-measurements the extractor mis-captures** — `typing:klebsiella:kleb` `10.5281` is a **Zenodo DOI prefix**; `typing:cat:catcolor` `5.1` and `typing:horse:horsecolor` `4.6` are **kb lengths** (`ARHGAP36 5.1-kb intron-1 deletion`, `STX17 4.6-kb dup`). These are not numbers a measurement artifact would ever carry.

`ADJUDICATED_BENIGN` cannot absorb the last two kinds: it holds **1 entry** today and `test_adjudicated_benign_is_small_and_every_entry_has_a_reason` pins `len(...) <= 5`. Growing it per-number would also repeat the hand-enumerated-exclusion-list trap this repo has hit five times. The fix is a small number of **typed patterns** plus a per-number provenance **kind**, so "unverifiable" stops meaning four different things at once.

## Codebase Context

- `scripts/contract_number_audit.py` — the audit. `audit()` iterates `cell_registry.cells()`, builds a prose blob from `PROSE_FIELDS = ("claim", "validation_slice", "label_provenance", "demotion_rule", "claim_status")`, extracts decimals with `DECIMAL_RE`, drops year-like tokens via `NON_MEASUREMENT = ^(19|20)\d{2}\.`, and matches each against the concatenated text of the artifacts the prose cites. Helpers to reuse unchanged: `_cited_artifacts` (with `_expand_braces` for the `{md,json}` form), `_load_artifact_text`, `_artifact_numbers`, `_number_variants`, `_matches_by_rounding`. `_decoy_full_match_count` is the function both defects live in. `main()` writes `wiki/contract_number_audit_<today>.json` — a NEW dated file per run day; these accrue and are committed.
- **Rejected approach already in the file, do not reintroduce:** tracing a number to the package source. The comment records that it accepted 17/17 tokens including 12 randomly generated ones. The same vacuity is what this plan is fixing one layer out, so the replacement must be measured, not assumed.
- `dna_decode/data/cell_registry.py` — 115 `CellContract` dataclasses; 18 fields; prose fields carry the numbers. Citations are plain `wiki/...` paths inside prose (added for four cells in `4d18399`). Not a frozen file, but it IS the live evidence surface, so edits must be provenance-only and proven so by diff.
- `tests/test_contract_number_audit.py` — 24 tests. Its own docstring states the load-bearing tests are the non-vacuity ones. Relevant pins: `test_adjudicated_benign_is_small_and_every_entry_has_a_reason` (`<= 5`), `test_live_audit_is_clean_and_non_trivial` (`n_numbers_checked >= 100`, `n_cells_audited >= 10`), `test_discrimination_actually_discriminates` (requires at least one HIGH cell and more than one grade), `test_discrimination_is_deterministic` (currently satisfied by the fixed seed; a full-pool scan satisfies it without one).
- `scripts/tier_evidence_audit.py` — sibling standing audit over the same registry (reports `under-claim 1 / over-claim 0` across 115 cells). Same read-only, artifact-comparing shape; its output must not move.
- Pool measured: **603 `wiki/*.json` (15 MB)** plus 619 `.md`. A full-pool scan over the json pool for four cells, with one shared load, ran in **16.5 s** — so the deterministic control is affordable for a manually-run audit.
- Frozen surface (must stay byte-unchanged): `dna_decode/eval/amr_rules.py`, `dna_decode/data/calibrated_amr_rules.json`.
- Test-run side effect to expect and revert: any broad pytest selection rewrites `wiki/certification_capstone.{json,md}` and `wiki/pgx_report_card.{json,md}` with date-only diffs.

### Reusable-Code Survey
- `scripts/contract_number_audit.py` — `_load_artifact_text` / `_artifact_numbers` / `_number_variants` / `_matches_by_rounding` reused verbatim by the full-pool control; only the sampling and haystack-construction change.
- `scripts/contract_number_audit.py` — `_cited_artifacts` + `_expand_braces` reused for the size-matched haystack (the citation count per cell is what sets decoy haystack size).
- `scripts/tier_evidence_audit.py` — sibling registry audit; reuse its read-only posture and per-cell row shape rather than inventing a second reporting convention.
- `tests/test_contract_number_audit.py` — existing non-vacuity idiom (plant drift, require it caught) reused for the new control's tests.
- None found for a generic artifact-corpus loader: `None — searched: graphify-out/GRAPH_REPORT.md (absent), src/utils, src/lib, src/common, src/helpers, utils, lib, common, helpers (all absent), scripts/`

## Pre-Change Baseline

Captured from the live audit at commit `bcfcbc2` before any edit:

- `verdict = NOTHING_TO_ADJUDICATE`, `n_candidate_drift = 0`, `n_numbers_checked = 157`, `n_cells_audited = 16`.
- `n_cells_with_numbers_but_no_citation = 10`, `n_numbers_unverifiable = 34`.
- `n_cells_low_discrimination = 2` (`pgx:human:ugt1a1` 8/10, `typing:Klebsiella:ktype` 5/10 under the k=10 seed-11 sample).
- k=10 grades to be superseded: `mlst` 0/10, `resfinder` 0/10, `ktype` 7/10, `pointfinder` 6/10 (the 7/10 and 6/10 figures are from the pre-commit probe; the committed run reports 5/10 and 3/10 for the same cells under a different decoy pool, which is itself evidence of the instability this plan removes).
- Full-pool truth measured during planning, for comparison after Step 1: `mlst` **6/602**, `resfinder` **1/602**, `ktype` **308/602**, `pointfinder` **355/602**; `essentiality:any:essentiality` matches 7/7 against at least 6 unrelated artifacts while its k=10 sample reported 0/10.
- `len(ADJUDICATED_BENIGN) == 1`.
- `tests/test_contract_number_audit.py` = **24 passed**.
- `scripts/tier_evidence_audit.py` = `under-claim 1 / over-claim 0 across 115 cells`.
- `dna_decode/eval/amr_rules.py` + `dna_decode/data/calibrated_amr_rules.json` sha — recorded via `git hash-object` so byte-equality is checkable at the end.

## Verification Signal

- The control is deterministic **without a seed**: two consecutive `audit()` calls return identical `decoy_full_match` for every cell, and `_DECOY_SEED` no longer exists.
- Every cell reports `decoy_full_match` out of the **full pool size** (`n_pool` in the hundreds, not 10), and the four measured cells land at the full-pool figures above (`mlst` 6, `resfinder` 1, `ktype` 308, `pointfinder` 355) — numbers I measured before writing this plan, so a mismatch means the implementation is wrong, not the expectation.
- Haystack size is matched: a cell citing N artifacts is compared against decoy haystacks of N artifacts, pinned by a test that a synthetic 1-number cell citing many artifacts does not grade HIGH.
- `essentiality:any:essentiality` grades **LOW**, not HIGH — the specific case that exposed the underpowering.
- Each number carries a provenance `kind`; `n_numbers_unverifiable` drops to only those genuinely lacking any provenance, and the sum of per-kind counts equals the old total (nothing is silently dropped).
- `len(ADJUDICATED_BENIGN) <= 5` still holds and the existing test is unmodified.
- `verdict`, `n_candidate_drift` and `n_numbers_checked` are unchanged by every kind/grade addition — proven by asserting the pre/post values, since a coverage or grading change must never read as drift.
- `scripts/tier_evidence_audit.py` output unchanged; frozen AMR surface byte-unchanged; `wiki/` churn limited to date-only diffs plus the new dated audit artifact.

## Implementation Steps

### Step 1: Replace the sampled decoy control with a deterministic full-pool, size-matched control
Files: scripts/contract_number_audit.py, tests/test_contract_number_audit.py
Depends on: none

**What changes:**
- Delete `_DECOY_SEED` and `_N_DECOYS`; `_decoy_full_match_count` no longer imports `random` or samples. It scans the entire `wiki/*.json` pool, excluding the cell's own cited artifacts, and returns `(full_match_count, n_pool_compared)`.
- Load the pool once per `audit()` call into the existing `_art_cache` and share it across all cells, so cost is one corpus read (~16.5 s measured for 4 cells; the load dominates, not the per-cell comparison).
- Size-match the haystack: when a cell cites N artifacts, each decoy trial concatenates N pool artifacts rather than 1, so haystack size is held constant between the cited and decoy conditions. Iterate consecutive non-overlapping N-sized groups over the sorted pool — deterministic, no sampling.
- Re-express the grade as a RATE against `n_pool_compared` rather than a count out of 10: HIGH when the rate is 0, LOW when >= 0.50, MODERATE otherwise. Keep the field names `decoy_full_match` / `n_decoys` so existing consumers and the `cells_low_discrimination` summary keep working; `n_decoys` now carries the pool size.
- `main()` printing gains the rate as a percentage so a reader sees `308/602 (51%)` rather than a bare count.

**Test strategy:**
- Determinism without a seed: two `audit()` calls agree on every cell's `decoy_full_match`, and `not hasattr(module, "_DECOY_SEED")`.
- Anchor on the four figures measured during planning (`mlst` 6, `resfinder` 1, `ktype` 308, `pointfinder` 355 out of 602) — these were measured before implementation, so they are a genuine expectation, not a snapshot of whatever the code does.
- Non-vacuity, the load-bearing test: `essentiality:any:essentiality` must grade LOW. Under the old k=10 sample it graded HIGH at 0/10 while 7/7 matched six unrelated artifacts, so this single assertion is what proves the underpowering is fixed.
- Size-matching: a synthetic cell with ONE round number citing many artifacts must not grade HIGH (its haystack advantage is removed).
- `test_discrimination_actually_discriminates` and `test_discrimination_is_deterministic` must pass unmodified; if either needs editing, say so explicitly rather than quietly relaxing a pin.
- Assert `verdict` / `n_candidate_drift` / `n_numbers_checked` equal the recorded baseline values.

### Step 2: Add a provenance `kind` per number, with pattern-typed structural non-measurements
Files: scripts/contract_number_audit.py, tests/test_contract_number_audit.py
Depends on: Step 1

**What changes:**
- Introduce a `kind` for every extracted number: `artifact` (matched in a cited artifact), `structural-non-measurement` (pattern-typed, below), `external-reference`, `enforced-by-test`, `threshold` (the existing `ADJUDICATED_BENIGN` role), `unverifiable` (the honest residual).
- Pattern-type the two structural families rather than listing numbers: a **DOI prefix** (`10.\d{4}` immediately followed by `/`, e.g. `Zenodo 10.5281/zenodo.14065540`) and a **kb/bp length** (a decimal immediately followed by `-kb`, `kb`, `-bp`, `bp`). Patterns are matched against the number's surrounding context in the prose, not the bare token, so `0.60` is never mistaken for a length.
- `external-reference` and `enforced-by-test` are per-number declarations carried in a small explicit table keyed by `(cell_id, number)` with a mandatory reason string — the same shape and discipline as `ADJUDICATED_BENIGN`, kept separate so the existing `<= 5` pin is untouched and each table stays individually small and auditable.
- Report per-kind counts in the summary, and make `n_numbers_unverifiable` mean only the `unverifiable` residual. Assert internally that per-kind counts sum to the total extracted, so no number can be silently dropped by a typing bug.

**Test strategy:**
- The DOI pattern types `10.5281` in the real `typing:klebsiella:kleb` prose and does NOT type an ordinary decimal; the kb pattern types `5.1` (catcolor) and `4.6` (horsecolor) and does not type `0.45`.
- Non-vacuity in the dangerous direction: a genuine measurement adjacent to the word `kb` elsewhere in the sentence is NOT typed as a length (guards over-broad context matching, which would silently exempt real numbers).
- Per-kind counts sum to the total extracted number count, for every cell.
- `len(ADJUDICATED_BENIGN) <= 5` unchanged; the new tables each carry a reason per entry.
- `verdict` / `n_candidate_drift` unchanged — a typing change must not read as drift.

### Step 3: Classify all 34 currently-unverifiable numbers by kind (research; no production code)
Files: wiki/contract_number_provenance_triage_2026-09-27.md
Depends on: Step 2

**What changes:**
- For each of the 34, record the kind, the evidence for that kind (the quoted prose context), and for the `artifact` kind the specific artifact path plus its **full-pool discrimination rate** from Step 1.
- Decide per cell whether a citation is worth adding, on one stated rule fixed before looking: cite only when the artifact resolves AND full-pool discrimination is better than LOW. A citation that cannot discriminate is recorded as such rather than added, because adding it would restate the `4d18399` over-claim this plan exists to correct.
- Known starting points from planning research, to be confirmed not assumed: `hla:human:b5701` -> `wiki/hla_b5701_validation_2026-07-06.json` (located, unverified); `pneumoserotype` -> one specific pneumo artifact, NOT the 7-file glob used in the probe; `typing:klebsiella:kleb` -> must be narrowed from the 33-file bundle to the actual source, since the bundle result was confounded by size.
- `essentiality:any:essentiality` is expected to land as not-worth-citing (its numbers match many unrelated artifacts); record that as a measured property of the numbers, not a defect of the cell.

**Test strategy:**
- Not a code step; verification is that every one of the 34 appears exactly once with a kind and evidence, and that the counts reconcile against the Step 1 audit output.
- Reconciliation is the check that matters: if the memo's per-kind totals disagree with the script's, one of them is wrong and both get re-read.

### Step 4: Add only the citations Step 3 admits, proven provenance-only
Files: dna_decode/data/cell_registry.py
Depends on: Step 3

**What changes:**
- Add the artifact citation to each cell Step 3 admitted, inline next to the script or cohort the prose already names (the `4d18399` shape).
- Dry-run each citation before adding it: compute would-be drift AND full-pool discrimination first. A citation that would surface drift is NOT added — it is reported as a finding, because drift in a shipped contract number is exactly what this audit exists to catch.
- Prove the edit is provenance-only rather than asserting it: the numeric-token multiset of the file must change only by date fragments from the new filenames, and `EvidenceTier` mention count must be unchanged.

**Test strategy:**
- Re-run the audit: `n_cells_audited` rises, `n_numbers_unverifiable` falls by exactly the count Step 3 admitted, `n_candidate_drift` stays 0.
- Registry diff check: only prose lines changed; numeric multiset delta is date fragments only; `EvidenceTier` count constant.
- Existing registry tests plus `tests/test_cell_registry.py` and the advertised-command guards stay green.
- `scripts/tier_evidence_audit.py` output unchanged (a prose edit must not move a tier).

### Step 5: Record the honest residual, and correct the two stale figures this run supersedes
Files: CLAUDE.md, project_state/evidence-surface-2026-08-31.md
Depends on: Step 4

**What changes:**
- CLAUDE.md: correct the standing audit figure and state the decoy control's real posture. The prior figure `0/128` was superseded to `0/157` during the run that produced this plan, and the k=10 grades it shipped with are now superseded again by full-pool rates — so the doc must carry the rates, the pool size, and the named limitation that a cell with few round numbers cannot be verified by this check at all.
- Record explicitly that `42ffe1f`'s k=10 control was underpowered and how it was caught (the `essentiality` 0/10-vs-7/7 case), so the next reader does not re-trust a sampled grade.
- AC9 row in the evidence-surface ledger via `/project-state --append-action`, superseding row 32's figures.

**Test strategy:**
- No code; verification is that the numbers in CLAUDE.md match the committed audit artifact exactly, re-derived from it rather than copied from this plan.
- Closing checks: frozen AMR surface byte-unchanged, prospective lock re-verifies, `tier_evidence_audit` unchanged, `wiki/` churn date-only plus the new dated artifact.

## Execution Preview

Wave 0 (1 parallel):  Step 1 — Replace the sampled decoy control with a deterministic full-pool, size-matched control
Wave 1 (1 parallel):  Step 2 — Add a provenance kind per number, with pattern-typed structural non-measurements
Wave 2 (1 parallel):  Step 3 — Classify all 34 currently-unverifiable numbers by kind
Wave 3 (1 parallel):  Step 4 — Add only the citations Step 3 admits, proven provenance-only
Wave 4 (1 parallel):  Step 5 — Record the honest residual and correct the stale figures

Critical path: Step 1 → Step 2 → Step 3 → Step 4 → Step 5 (5 waves)
Max parallelism: 1 agent

Note: Parallel execution requires a git repository with a configured remote. If unavailable, /execute-plan falls back to sequential mode.

This plan is **fully sequential by design, not by accident.** Steps 1 and 2 edit the same two files, and Step 4 edits the live evidence surface. The HIGH-salience prior decision on `[plan_file: Oxford_Cohort_External_Revalidation_Plan/technical-plan.md]` (2026-06-15) forced sequential mode for exactly this shape — a contract-coupled chain — specifically to keep autonomous worktree agents away from frozen and contract-bearing files. The same reasoning applies here; there is no parallelism to recover.

## Risk Flags

- Severity: high — **Step 1 will change numbers that are already committed and cited.** `ktype` and `pointfinder` move from MODERATE/LOW at k=10 to roughly 51% and 59% of the full pool, and `essentiality` flips HIGH to LOW. That is the point of the step, but it means the `42ffe1f` artifact and the CLAUDE.md figures become wrong the moment Step 1 lands, and Step 5 is the only thing that fixes them. If execution stops between Step 1 and Step 5, the repo carries a fresher audit artifact than its own documentation.
- Severity: high — **File overlap between Steps 1 and 2** (both `scripts/contract_number_audit.py` + `tests/test_contract_number_audit.py`), resolved by an explicit dependency and sequential execution. Do not let a wave-parallel executor pick these up together.
- Severity: medium — **The size-matched haystack is a design choice with no measured optimum.** Grouping the pool into consecutive N-sized blocks is deterministic and removes the size advantage, but a different grouping could give a different rate for multi-citation cells. Report the grouping rule in the artifact so the number is interpretable; do not present the rate as the only possible one.
- Severity: medium — **Pattern-typing can silently exempt a real measurement.** A context-matched `kb` or DOI pattern that is too greedy would type a genuine number as structural and remove it from the audit entirely — the exact failure mode (a check that quietly stops checking) this plan is correcting. Mitigated by the Step 2 test that a measurement near the word `kb` is not typed, but the residual risk is real and the pattern must stay narrow.
- Severity: medium — **`essentiality`'s 1/7 own-cover is unexplained.** Its candidate artifact contains only one of its seven numbers. That is either the wrong artifact or genuine drift in a shipped contract, and Step 3 must distinguish them. If it is drift, it is a finding that outranks the rest of this plan and should be surfaced rather than folded into a coverage statistic.
- Severity: low — **External tool integrations:** none. This plan shells out to nothing new. Surfaces verified by direct execution during planning: `scripts/contract_number_audit.py` (run, stdout + artifact shape read), `scripts/tier_evidence_audit.py` (run, `under-claim 1 / over-claim 0`), `uv run pytest` selections, `git hash-object` / `git show HEAD~N:<path>` for byte-equality. `sentrux` absent so the architecture gate records `n/a — sentrux not installed`.
- Severity: low — **Dated-artifact accrual.** `main()` writes `wiki/contract_number_audit_<today>.json`, so each run day adds a file. Not a shared-key overwrite (the prior HIGH-salience trap), but the directory grows and the test pinning `contract_number_audit_2026-09-24.json` validates an increasingly old file — worth noting, not fixed here.
- Severity: low — **Restructuring applied:** the originally-banked "three classes" framing was replaced by four kinds plus a control fix, because grounding showed the classes were not exhaustive (an external-reference kind exists) and that the control itself was the larger defect. Scope is not reduced — all 34 numbers are still addressed — but the ordering now puts the control first, since citations are only decidable once discrimination is measured correctly.

## Open Questions

- Should the `kb`/DOI pattern-typed numbers be dropped from extraction entirely, or retained with a `structural-non-measurement` kind and counted separately? The plan assumes retained-and-counted (nothing silently disappears), which is more auditable but keeps them visible in per-kind totals forever.
- `typing:Escherichia_coli:pathotype` carries `0.917`, which CLAUDE.md records as a **superseded over-rescue** value that must not be quoted, deliberately present in the prose as the rejected number. Should the audit gain a `superseded-value` kind so a deliberately-recorded wrong number is distinguishable from an unverifiable one? Treating it as `enforced-by-test` alongside 0.833 would be inaccurate.
- Is the LOW threshold of 0.50 of the pool the right bar? It is asserted, not derived. A cell at 49% would grade MODERATE while being nearly as unverifiable as one at 51%.

## Verification

1. `uv run python scripts/contract_number_audit.py` — confirm `verdict = NOTHING_TO_ADJUDICATE`, `n_candidate_drift = 0`, `n_numbers_checked >= 157`; confirm every cell reports a full-pool rate with `n_decoys` in the hundreds; confirm the four anchored figures (`mlst` 6, `resfinder` 1, `ktype` 308, `pointfinder` 355 of 602) and that `essentiality` grades LOW.
2. `uv run pytest tests/test_contract_number_audit.py -q` — all prior 24 tests plus the new ones green, with `test_discrimination_actually_discriminates` and `test_discrimination_is_deterministic` unmodified.
3. `uv run pytest tests/ -q -k "registry or contract or card or trust or tier"` — no regressions against the 717-passed baseline selection.
4. `uv run python scripts/tier_evidence_audit.py` — still `under-claim 1 / over-claim 0 across 115 cells`.
5. Registry provenance-only proof — numeric-token multiset of `dna_decode/data/cell_registry.py` differs only by date fragments; `EvidenceTier` mentions unchanged at 67.
6. Frozen surface — `git hash-object dna_decode/eval/amr_rules.py dna_decode/data/calibrated_amr_rules.json` equal to the pre-change hashes; prospective lock re-verifies.
7. `git status --porcelain` — only intended files plus the new dated audit artifact; `wiki/certification_capstone.*` and `wiki/pgx_report_card.*` reverted after any broad pytest run (verify non-date diff lines are zero first).
8. Read CLAUDE.md's audit figures back against the committed artifact and confirm they were re-derived from it, not transcribed from this plan.

## Save-time amendments

Captured at: 2026-09-27
Source: `/save-plan` arguments

> Audit notes only. `/execute-plan` reads executable work from `## Implementation Steps` alone;
> amendments here are provenance for human readers and are NOT executable instructions. If an
> amendment changes step contracts (file lists, dependencies, structure), re-run `/technical-plan`
> before `/execute-plan`.

- reordered to put the decoy-control fix first after grounding showed the k=10 control was underpowered; three classes became four kinds

<!-- toolkit: check=clean waves=clean gate=fired:open-questions,severity-high,test-strategy-leak -->
