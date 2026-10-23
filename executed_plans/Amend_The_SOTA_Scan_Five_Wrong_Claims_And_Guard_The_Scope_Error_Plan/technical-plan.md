## Lens status
Repo-index: absent — `.repo-index/` and `repo-index-status.sh` both missing on this host; source-only reading.
Adversarial review: Claude-only — `codex` absent (exit 1), LM Studio UNREACHABLE, both probe helpers absent.
Architecture gate: n/a — sentrux not installed.
Reuse signals: `graphify-out/GRAPH_REPORT.md` absent; no conventional util dirs.
Wiki replay: 48 HIGH-salience entries in `wiki/decisions-log.md`; two bear directly (see Codebase Context).
Pre-exec `/brainstorm`: ALREADY RUN — this plan exists because it returned three Critical findings.

## Problem Statement

Three wrong claims and one wrong accounting sit in two committed, pushed `wiki/` artifacts
(`wiki/ai_genomics_sota_scan_2026-10-22.{md,json}`, commits `cb1b2a05` / `614ed27f`) and in a pushed
commit-message body. An adversarial review of that scan found them; this plan amends them and guards
the one whose error class this project has already paid three months for.

What is wrong, in severity order:

1. **C-1 (scope).** The scan says two external 2026 benchmarks "reproduce" the project's supervised
   encoder-gate negative. Both score by **log-likelihood** — the Mechanistic Invariance Test's
   CSS/SCR are LL-derived and report an `AT-LL r` column; Nullsettes defines a detection as *"mutant
   log-likelihood significantly lower than its nonmutant original."* The encoder gate is not a
   likelihood probe: `wiki/glm_encoder_gate_result_2026-10-08.md` contains "likelihood" **zero
   times** (measured), and the comparison is a fitted conv encoder against one GC feature with
   per-seed deltas and a paired CI. CLAUDE.md names this exact conflation and its cost: *"the
   unscoped version cost three months… The body below is entirely about ZERO-SHOT likelihood."*
   The repo holds a measured instance where the two readouts give **opposite** verdicts (supervised
   0.8142 vs zero-shot 0.4750 on the same data, PASS vs FAIL).
2. **C-2 (wrong model, false alarm).** Decision D2 asks whether to re-download
   `GENERator-prokaryote-1.2b-base`. The code pins `GENERator-**v2**-prokaryote-1.2b-base` at five
   sites; the non-v2 id does not resolve publicly (HF API: `Invalid username or password`). And the
   local snapshot `8b2f768b0d293953518ff91d34600f9322ef1f94` is **byte-identical to upstream HEAD**
   for the pinned id (`lastModified 2026-07-05`). There is nothing to re-download and nothing to
   decide. D2 must be withdrawn, not answered.
3. **C-3 (incomplete read asserted as fact).** "The field left our question unasked" rests on a read
   covering characters 0–200,000 of a 209,274-character page; the final 9,274 were never read. The
   phrase is in the pushed commit message and stated as a finding.
4. **M-1 (accounting).** The scan says "roughly four of the nine `dG_*` components are computable
   there and three are structurally absent" — silently dropping two. `DG_FIELDS` is exactly nine, and
   the true split is **4 clean / 3 absent / 2 composites** (`dG_total`, `dG_bind`) that sum over the
   absent terms. `dG_total` is the headline feature.
5. **M-3 (category error).** "Freezing is not the bottleneck… retires 'we should have fine-tuned'"
   over-reads `arXiv 2608.05329`: its promoter tasks are classification (its own limitation says
   evaluation is classification-only) while the gate is expression regression, and the gate's encoder
   was trained from scratch rather than frozen, so frozen-vs-finetuned is not the gate's axis at all.

**Success criteria.** Both artifacts carry the corrected claims; the `.json` still parses; a guard
fails on the verbatim pre-fix prose and passes on the corrected prose; the frozen five, `uv.lock` and
`pyproject.toml` are byte-unchanged; the suite is no worse than baseline.

**Non-goals.** Do NOT reopen the encoder gate. Do NOT re-size the frozen pre-registered promcalc bar.
Do NOT predict the promcalc verdict. Do NOT adopt the NAR substrate (D1 stays an open, re-scoped
question). Do NOT touch the running stage-B job.

## Codebase Context


See: ./codebase-context.md

## Pre-Change Baseline

- **Artifacts:** both committed and pushed at `cb1b2a05` (+ ledger row 7 at `614ed27f`). The five
  defects are at the verified line/field anchors listed in Codebase Context. Pre-fix prose is
  recoverable via `git show cb1b2a05:wiki/ai_genomics_sota_scan_2026-10-22.md` (5 occurrences of
  "reproduce"), so Step 3's non-vacuity strings do not depend on capturing them before editing.
- **Test suite, measured this session:** 14 failed / 6943 passed / 186 skipped. All 14 failures are
  absent-data or host-dependent; none is a code defect (CLAUDE.md records the same 14).
- **Frozen surface:** `dna_decode/eval/amr_rules.py`, `dna_decode/data/calibrated_amr_rules.json`,
  `dna_decode/data/mic_tiers.py`, `dna_decode/data/shipped_decoder_surface.py`,
  `dna_decode/eval/cohort_manifest.py`, `uv.lock`, `pyproject.toml` — all byte-unchanged as of
  `614ed27f` (verified by `git diff --stat`).
- **Running job:** stage B mid-flight, ~20,000/43,611 tiles at ~55 tiles/min. This plan must not
  perturb it; no step touches `data/promcalc/`, `scripts/promcalc_decompose.py`, or the launcher.

## Verification Signal

The plan succeeds iff all of these hold after Step 5:

1. `uv run python -c "import json; json.load(open('wiki/ai_genomics_sota_scan_2026-10-22.json'))"`
   exits 0.
2. `uv run pytest tests/test_sota_scan_scope_guard.py -q` passes, **and** its non-vacuity cases
   demonstrably fail against the pre-fix strings pulled from `cb1b2a05` (the test asserts this
   internally, so a vacuous guard fails the run).
3. `grep -c "D2" wiki/ai_genomics_sota_scan_2026-10-22.md` returns a count consistent with D2 being
   present only as a **withdrawn** item naming the revision match, and
   `grep -c "GENERator-prokaryote-1.2b-base" ` returns 0 outside a quoted-as-wrong context.
4. The `.md` and `.json` agree: every amended claim appears in both, checked by Step 4's
   cross-consistency assertions rather than by eye.
5. Lines 58 and 92's correct "reproduces" usage (the 16,227-vs-16,226 fold check) is **still
   present** — a guard that removed it would have over-corrected.
6. `git diff --stat` over the frozen five + `uv.lock` + `pyproject.toml` is empty.
7. Full suite no worse than 14 failed / 6943 passed / 186 skipped, with the identical failure set and
   passes increased by exactly the number of new tests in Step 3.

## Implementation Steps

### Step 1: Amend the prose artifact
Files: wiki/ai_genomics_sota_scan_2026-10-22.md
Depends on: none

**What changes:**
- **C-1, the load-bearing edit.** Retitle F1 from "…results reproduce this project's encoder-gate
  negative" to a claim the evidence supports: the two benchmarks are **zero-shot likelihood probes**
  whose finding is **consistent in direction** with the supervised gate but is **not a replication**
  of it. Rewrite the "three teams, three benchmarks, one mechanism" synthesis to say what is actually
  shared — a finding about *what the pretrained representation carries*, not about what a supervised
  head can extract from it. Add an explicit scope rail in the section naming CLAUDE.md's standing
  rule and the measured opposite-verdict instance (supervised 0.8142 vs zero-shot 0.4750).
- Fix line 68's "independently reproduce them" the same way.
- **Leave lines 58 and 92 untouched** — their "reproduces" is the fold-range identity check and is
  correct.
- **C-2.** Replace the D2 block with a withdrawal: the code pins the v2 id at five sites, the non-v2
  id does not resolve, and the local snapshot equals upstream HEAD. Record the measured revision so
  the question is not re-asked. Keep it visible as a withdrawn decision rather than deleting it, so
  the record shows the claim was made and retracted.
- **C-3.** Downgrade the "unasked" finding to `[inferred]`, naming the unread tail (chars
  200,000–209,274) as the reason.
- **M-1.** Replace "roughly four… and three… absent" with the enumerated 4 / 3 / 2 split, naming
  `dG_total` and `dG_bind` as composites that sum over the absent terms, and noting that `dG_total`
  is the headline feature.
- **M-3.** Narrow the "freezing is not the bottleneck" claim to what the cited paper supports:
  classification-only evaluation, and an axis (frozen-vs-finetuned) that the from-scratch encoder
  does not sit on. Keep the surviving F2 point (that paper names GC/k-mer comparison as the thing it
  did not do).
- Add a short amendment header recording that these five changes came from an adversarial review on
  2026-10-22, so a future reader sees the correction rather than a silently-rewritten memo.

**Test strategy:**
- Mechanical: re-grep each anchor to confirm the intended string is gone and the replacement present.
- Confirm lines 58/92 survive (`grep -n "16,226"` still returns two hits).
- No automated test here; Step 3 supplies the standing guard.

### Step 2: Mirror the amendments into the JSON sidecar
Files: wiki/ai_genomics_sota_scan_2026-10-22.json
Depends on: Step 1

**What changes:**
- `findings[0].headline` and `.why_they_agree`: re-scoped per Step 1's decided wording, with an
  explicit `readout: "zero-shot log-likelihood"` on each external result and a
  `not_a_replication_of_supervised_gate: true` field so the distinction is machine-readable, not
  only prose.
- `findings[0].our_result`: add `readout: "fitted conv encoder vs one GC feature, paired CI over
  seeds"` and `likelihood_mentions_in_artifact: 0` (the measured fact that settles it).
- `decisions_needed`: D2 becomes a `withdrawn` entry carrying `pinned_id`,
  `local_snapshot_sha`, `upstream_sha`, `shas_identical: true`, `upstream_last_modified`, and
  `non_v2_id_resolves: false`. D1 is retained and re-scoped (see Step 1) — it stays open.
- `decisions_needed[0]`: add `components_composite_contaminated: ["dG_total", "dG_bind"]` and a
  `component_accounting: "4 clean / 3 absent / 2 composite"` field so the arithmetic closes to nine.
- `findings[1].consequence`: narrowed per M-3, with `cited_paper_evaluation: "classification only"`
  and `gate_encoder_was: "trained from scratch, not frozen"`.
- `measured_here.open_gap_the_field_left`: prefix with its tier and add
  `page_chars_read: "0-200000 of 209274"` plus `tier: "inferred"`.
- `verification_tiers`: unchanged. `read_only` / `frozen_surface_touched` / `running_job_touched`:
  unchanged.

**Test strategy:**
- `json.load` must succeed (the recorded gotcha: a hand-edited wiki JSON was invalid for a day
  because nothing read it).
- Assert the amended keys exist and the component arrays now sum to nine.

### Step 3: Add the standing scope guard
Files: tests/test_sota_scan_scope_guard.py
Depends on: Step 1

**What changes:**
- New test module modelled directly on `tests/test_essentiality_overclaim.py`, reusing all four of
  its parts: a surfaces tuple, `_segments_mentioning()` splitting on blank lines, `_violations()`,
  and a `PRE_FIX_STRINGS` non-vacuity tuple.
- The guarded **relationship**, stated generally rather than for this one memo: in any `wiki/*.md`
  segment that cites the encoder-gate delta (`+0.0122`) alongside either external benchmark
  (Mechanistic Invariance Test / Nullsettes / `2604.06549` / `2506.10271`), a zero-shot marker must
  appear in that same segment, and the words "reproduce"/"replicate" must not attach the external
  results to the gate. Scoped to `.md` because prose guards guard prose; the `.json` is covered by
  Step 2's structural assertions.
- `PRE_FIX_STRINGS` pulled from `git show cb1b2a05:wiki/ai_genomics_sota_scan_2026-10-22.md` — the
  F1 heading and the line-68 sentence verbatim — so the guard is proven to fail on the exact text
  that shipped.
- One premise-pinning test in the same spirit as the precedent's closing test: assert
  `wiki/glm_encoder_gate_result_2026-10-08.md` contains zero occurrences of "likelihood". If a future
  rewrite made the gate a likelihood probe, this guard's whole rationale would be void — so the
  rationale is checked, not trusted.

**Test strategy:**
- The module *is* the test. Non-vacuity is internal and mandatory: every pre-fix string must be
  flagged, so a guard matching nothing fails rather than passing forever.
- Run it alone first, then confirm it is collected by the full suite.

### Step 4: Cross-consistency and regression verification
Files: wiki/ai_genomics_sota_scan_2026-10-22.md, wiki/ai_genomics_sota_scan_2026-10-22.json, tests/test_sota_scan_scope_guard.py
Depends on: Step 2, Step 3

**What changes:**
- No content changes. Runs the Verification Signal checks and records their outputs.
- Cross-consistency `.md` vs `.json`: D2 withdrawn in both; the 4/3/2 split in both; the zero-shot
  scoping in both; the `[inferred]` tier on the NAR gap in both. Asserted, not eyeballed.
- Confirm lines 58/92's correct "reproduces" usage survived (over-correction check).
- `git diff --stat` over the frozen five + `uv.lock` + `pyproject.toml` must be empty.
- Full suite, reconciled against 14 / 6943 / 186: identical failure set, passes up by exactly the
  number of tests added in Step 3.
- **Precondition, checked before running the suite:** no torch job holds a model — CLAUDE.md records
  that running the full suite then dies with an access violation that reads like a test failure. The
  running stage-B job is not a torch job, but verify rather than assume.

**Test strategy:**
- This step is verification; its output is the evidence for Step 5's ledger row.
- If the suite count does not reconcile, stop and diagnose rather than recording a passing-looking
  result — the recorded discipline is that a flaw verification surfaces is diagnosed, not narrated.

### Step 5: Record the correction in the ledger and commit
Files: project_state/glm-umbrella-2026-10-07.md
Depends on: Step 4

**What changes:**
- One AC9 cheap-path Action Log row, class `write-plan`, at the `<!-- project-state:end:action-log -->`
  marker, with Step-6.5 post-mutation verification (re-read, row present, count incremented,
  marker still unique).
- The row records what was wrong rather than only what was done: the zero-shot/supervised conflation
  and that it is the error class CLAUDE.md prices at three months; that D2 named a model the code does
  not pin and dissolved on a two-call measurement; that a finding was asserted from a page whose tail
  was unread; and that the component split dropped two of nine.
- Commit the three touched files plus the ledger. Push.

**Test strategy:**
- Step 6.5 verification is the test, as in the two prior rows this session.
- Re-run the new guard post-commit to confirm the committed text passes it.

## Execution Preview

| Wave | Steps | Parallelism |
|---|---|---|
| 0 | Step 1 | 1 |
| 1 | Step 2, Step 3 | 2 |
| 2 | Step 4 | 1 |
| 3 | Step 5 | 1 |

Total waves: 4. Max parallelism: 2. Critical path: Step 1 → Step 2 → Step 4 → Step 5 (4 steps).

**Sequential execution is recommended over the computed 2-way parallelism.** Steps 2 and 3 touch
disjoint files, but both consume the *wording decided in Step 1*, and `data/` is gitignored
(`.gitignore:4`) so a fresh worktree cannot reconcile Step 4's suite baseline against the count
measured in this checkout — the same reason recorded for the previous plan in this repo. The
parallelism is real but worth ~one step of wall-clock against a real divergence risk.

## Risk Flags

- **Step 3 is the one item beyond the review's literal recommendation.** The review said amend the
  artifacts and leave D1 open; it did not ask for a guard. `decisions-log.md:132` is a standing
  caution against reaching for more scope after a review. The case for including it: this is a
  documented three-month error class, the pattern already exists in two places, and the marginal cost
  is one test file. The case against is in Open Questions, and dropping Step 3 leaves Steps 1/2/4/5 a
  complete plan.
- **Over-correction is a live risk.** Lines 58 and 92 use "reproduces" correctly about the fold-range
  identity check. A careless global replace would delete the one genuinely-reproduced result in the
  memo. Pinned explicitly in Verification Signal item 5.
- **`[unverified]` — the NAR page tail.** C-3's fix is to downgrade the claim, not to resolve it.
  Reading chars 200,000–209,274 would settle whether the biophysical-vs-learned comparison is truly
  absent; this plan deliberately does not, because the honest label is cheap and the resolution is
  not. If a later session wants the claim as fact, that read is the step.
- **Amending a pushed commit message is not possible without a rewrite**, which is forbidden here.
  The commit body for `cb1b2a05` will keep the "reproduce" framing and the "unasked" phrase. The
  amendment therefore lives in the artifacts plus a new commit that names the correction — the git
  history will show the claim and its retraction, which is the honest shape but worth stating.
- **Step 4's suite run is the only expensive step** (~40 min at the last measurement) and must not
  overlap a torch job.

## Open Questions

- **Should Step 3 ship at all?** Including it guards a documented-cost error class with an existing
  pattern; excluding it honors the review's literal scope and `decisions-log.md:132`. Recommended:
  include, because the guard targets a *relationship* that generalizes past this memo. Left open
  rather than folded in.
- **Should the guard cover all of `wiki/` or only this artifact?** Globbing `wiki/*.md` catches future
  artifacts that cite both the gate delta and the external benchmarks; scoping to one file is
  narrower and cannot false-positive. Recommended: glob, since the co-occurrence is specific enough
  that a false positive is unlikely and a missed future artifact is the failure that matters.
- **Should D1 be re-scoped in place or split into its own artifact?** It currently reads as "a second
  promcalc arm", which the review judged misleading. Renaming it inside this memo is cheap; moving it
  to a separate open-question artifact is cleaner but creates a second file to keep in sync. Not
  decided here.
- **Does the NAR substrate's GC/dG collinearity make D1 worth running?** Unmeasured. Settling it is
  one fit on the already-downloaded 112,955 rows, and it is the fact that would make D1 decidable —
  but it is new work, not an amendment, so it is out of this plan's scope.

## Verification

1. `uv run python -c "import json; json.load(open('wiki/ai_genomics_sota_scan_2026-10-22.json'))"`
   exits 0.
2. `uv run pytest tests/test_sota_scan_scope_guard.py -q` passes, including its internal non-vacuity
   cases against the `cb1b2a05` pre-fix strings.
3. `grep -n "16,226" wiki/ai_genomics_sota_scan_2026-10-22.md` still returns two hits
   (over-correction check).
4. Cross-consistency assertions between `.md` and `.json` pass for all five amendments.
5. `git diff --stat` over the frozen five + `uv.lock` + `pyproject.toml` is empty.
6. `uv run pytest tests/ -q` reconciles against 14 failed / 6943 passed / 186 skipped with an
   identical failure set and passes up by exactly the count of new tests.
7. Action Log row present, count incremented, end-marker unique (Step 6.5).
8. `data/promcalc/decompose_run.log` still advancing — the running job was not perturbed.

## Save-time amendments

Captured at: 2026-10-23
Source: `/save-plan` arguments

**Audit notes only.** `/execute-plan` reads executable work from `## Implementation Steps`
alone; nothing in this section is an instruction. If an amendment changes a step contract
(files, dependencies, structure) the toolkit footer will go non-clean on re-validation, which
is the signal to re-run `/technical-plan` before executing.

- withdrew D2 entirely rather than answering it — the code pins the v2 model id at five sites, the non-v2 id the scan named does not resolve publicly, and the local snapshot is byte-identical to upstream HEAD, so a two-call measurement dissolved the decision
- the zero-shot/supervised conflation identified as the load-bearing fix because CLAUDE.md prices that exact error class at three months and the repo holds a measured opposite-verdict instance
- the component split corrected from 'roughly four of nine' to an enumerated 4/3/2 that closes the arithmetic
- sequential execution recommended over the computed 2-way parallelism
- Step 3's guard flagged as the only item beyond the review's literal scope, with dropping it left as an open question rather than folded in
<!-- toolkit: check=clean waves=clean gate=fired:open-questions,unverified,reuse-missing -->
