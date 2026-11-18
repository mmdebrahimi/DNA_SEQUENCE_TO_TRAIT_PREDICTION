# Normalize the Active-Frame Result By Its Own Ceiling, And Make Link 4 Launch-Ready

## Lens status
Execution mode: to be decided by `/execute-plan` (2 waves available, max parallelism 3)
gh available: yes
sentrux available: yes — but the plan is artifact/doc-level over 7 files with no module graph change, so the structural-change heuristic excludes the Architecture Gate
Repo-index: absent (`$HOME/.claude/scripts/repo-index-status.sh` not found) — source-only reading, no degradation
Baseline tests: 12 failed / 8,904 passed / 205 skipped (`uv run pytest tests/ -q`)

## Problem Statement

Yesterday's scoring-frame work established that a pooled correlation must be read against its own frame's noise ceiling, because 0.941 is a pooled average over deciles whose own ceilings run 0.5013 to 0.9823. The two artifacts that work produced do not themselves apply that rule, and the omission makes a specific wrong inference available.

Measured during review: the same 13,217-parameter conv encoder scores 0.3234 on the pooled frame against a 0.941 ceiling (34.4% of it) and 0.1636 on the active frame against a 0.9690 ceiling (16.9% of it). The active-frame margin over GC grows from +0.0122 to +0.2190 — but it grows because GC falls 0.3576 while the encoder falls 0.1598. Both arms get worse; the comparator collapses. Without the normalization printed beside them, the growing margin reads as evidence that the learned representation is stronger than family G-A concluded, when relative to achievable signal it is weaker.

That inference was actually drawn during this session and had to be retracted, so the hazard is demonstrated rather than hypothetical. The fix is to record the normalization and the comparator-collapse reading in the artifacts that carry the margins, and to record the observation against the G-A arm as scope — explicitly not as a verdict change.

Second, independent problem: the finetune arm's final link is 28 epochs that need about 8.1 h of Kaggle GPU against 2.70 h remaining, refreshing 2026-11-21. Everything short of the push can be validated now, so that refresh day is a one-command day rather than a day of discovering build breakage.

## Codebase Context

`scripts/glm_genomewide_oracle.py` is the existing owner of ceiling normalization. It computes `a["fraction_of_ceiling"] = round(a["spearman"] / ceiling, 4)` at lines 134 and 169, and `_verdict(arms, ceiling)` at line 230 formats the phrase `rho {x} ({frac:.1%} of the {ceiling} ceiling)`. Its `PREREGISTERED["q1_bars"]` carries `fraction_of_ceiling_for_GENERALISES: 0.50`. The field name to reuse is `fraction_of_ceiling`, already present in `wiki/glm_tile_headroom_2026-10-08.json` and `wiki/glm_genomewide_oracle_2026-10-07.json`.

`dna_decode/glm/genomewide.py:270` owns `noise_ceiling(frags)`. No new ceiling derivation is needed: both ceilings are already recorded. The pooled 0.941 is `noise_ceiling` in `wiki/glm_learned_representation_2026-10-08.json`; the active 0.9690 is named in `wiki/glm_active_frame_gate_2026-11-16.json` under `preregistered.thresholds_are`, which states it is "DERIVED and contextualises them without setting them".

`wiki/glm_active_frame_gate_2026-11-16.json` already records the medians this plan needs — `phase2_active.gc_median` −0.0464, `encoder_multi_median` 0.1636, `margin_median` 0.219, `n_seeds_clearing_explains_bar` 9, `margin_ci` [0.1769, 0.2342]. It also already self-flags `conditions_on_the_outcome: true` and `not_a_ledger_candidate: true`. That artifact is committed at `eb0963eb` and is a recorded run result; this plan does not modify it.

`soraya_runs/arms-ledger.json` has established precedent for annotating an existing arm dict without adding a top-level key: the `fragment-rung-321k-rows-size-lever-LB-M9` arm carries a `metric_not_commensurable_with_siblings` field added on 2026-11-14 "to make a PROSE-ONLY caveat machine-visible". A top-level addition is structurally impossible there, because `len(criteria[crit])` is the arm count and a test pins `set(data) == {"schema_version", "criteria"}`.

`scripts/kaggle/build_promoteratlas_finetune_kernel.py` exposes `--out`, `--self-test`, `--smoke`, `--wall-hours`, `--ckpt-dataset`, `--horizon-epochs` (verified by parsing its `add_argument` calls). `scripts/kaggle_push_poll.py` exists. The Kaggle CLI surface used here was verified live this session: `kaggle quota`, `kaggle datasets list --mine`, `kaggle kernels list --mine` all return as documented, and `emanueleebrahimi/ft-ckpt-s3` exists at 186,587,892 bytes.

### Reusable-Code Survey

- `scripts/glm_genomewide_oracle.py` — owns `fraction_of_ceiling`; reuse the field name, the 4-decimal rounding and the `{frac:.1%} of the {ceiling} ceiling` phrasing rather than inventing a second spelling.
- `dna_decode/glm/genomewide.py:270::noise_ceiling` — owns ceiling derivation; not called here because both ceilings are already recorded in artifacts.
- `dna_decode/eval/preregistration.py` — this session's validator; its `assert_range_reported_with_aggregate` is the existing guard against quoting a stratified aggregate as a scalar, and is the natural place if a reusable check is wanted later. Not called in this plan.
- `soraya_runs/arms-ledger.json::fragment-rung-321k-rows-size-lever-LB-M9.metric_not_commensurable_with_siblings` — the precedent shape for Step 3's annotation.
- Searched: `scripts/`, `dna_decode/`, `wiki/*.json`, `wiki/*.md`; `graphify-out/GRAPH_REPORT.md` absent.

## Pre-Change Baseline

- Full suite: **12 failed / 8,904 passed / 205 skipped** under `uv run pytest tests/ -q`. Failure set is 7 in `tests/test_cyp4f2_getrm_concordance.py`, 3 in `tests/test_essentiality_report_card_ladder.py`, 1 in `tests/test_readme_quickstart_output.py`, 1 in `tests/test_pneumo_selection_probe.py` — all absent-data, all long-standing.
- `wiki/glm_active_frame_regrade_2026-11-17.json` verdict is `FRAME_EXPLAINS_THE_GAP` with `bar_diff_field_count: 1`.
- Neither active-frame artifact contains the string `fraction_of_ceiling`. Verified by grep: it appears in `wiki/glm_tile_headroom_2026-10-08.json`, `wiki/glm_genomewide_oracle_2026-10-07.json` and the eukaryotic-arm artifacts, and in neither active-frame artifact.
- `soraya_runs/arms-ledger.json` arm `conv-3width-13k-params-150bp` records `metric: 0.3234`, `status: measured`, `evidence: "wiki/glm_learned_representation_2026-10-08.json (WEAK; paired CI includes zero)"`.
- Link 4 readiness is unmeasured: no artifact records whether the kernel builds and self-tests for a `--ckpt-dataset ft-ckpt-s3` configuration.

## Verification Signal

1. `wiki/glm_active_frame_regrade_2026-11-17.json` contains `fraction_of_ceiling` for the encoder in both frames, with values rounding to 0.3437 pooled and 0.1688 active.
2. The same artifact states the comparator-collapse reading in words, including both deltas (encoder −0.1598, GC −0.3576).
3. The same artifact carries an explicit cross-frame scoping caution, and a test fails if it is removed.
4. `soraya_runs/arms-ledger.json` arm `conv-3width-13k-params-150bp` gains a scope field, and its `status`, `metric` and `evidence` values are **byte-identical** to baseline — asserted by test, because the whole point is that this is not a verdict change.
5. `wiki/glm_scoring_frame_resolution.md` carries the normalization beside its margin table.
6. Full suite failure **set** unchanged from the four files above; pass count rises by exactly the number of new tests.
7. `wiki/finetune_link4_readiness_2026-11-17.json` records a kernel `--self-test` result, the resolved checkpoint dataset, the chosen fresh slug, and an explicit `pushed: false`.

## Implementation Steps

### Step 1: Add cross-frame ceiling normalization to the re-grade
Files: scripts/active_frame_regrade.py, wiki/glm_active_frame_regrade_2026-11-17.json
Depends on: none

**What changes:**
- Add a `ceiling_normalization` block to the dict `regrade()` returns, carrying for each frame: the ceiling, the encoder median, the GC median, and `fraction_of_ceiling` computed as `encoder_median / ceiling` rounded to 4 places — reusing the field name and rounding from `scripts/glm_genomewide_oracle.py`.
- Read the pooled pair from `wiki/glm_learned_representation_2026-10-08.json` (`noise_ceiling`, `arms["encoder-multi"]["median"]`, `arms["gc"]["median"]`) and the active pair from the recorded `phase2_active` medians plus the active ceiling. Add a `RegradeRefused` branch for each missing field rather than defaulting, matching the existing `recorded_ci` refusal.
- Add the active ceiling as a module constant with its provenance string quoted from the gate's `preregistered.thresholds_are`, so the 0.9690 is sourced rather than typed.
- State the comparator-collapse reading explicitly: both arms fall, the encoder by 0.1598 and GC by 0.3576, so the margin grows by comparator collapse and the encoder reaches a *smaller* fraction of achievable signal on the active frame.
- Add a scoping caution noting that `scripts/glm_genomewide_oracle.py` documents "Do not compare the two fractions" for fractions taken across *substrates*; these two are across *frames of one substrate*, which is a weaker comparison but still requires naming which ceiling each fraction uses.
- Add a line stating this is **not** grounds to revisit G-A, citing the gate artifact's own `conditions_on_the_outcome: true` and `not_a_ledger_candidate: true`.
- Re-run `uv run python -m scripts.active_frame_regrade` to regenerate the artifact and inspect the printed output.

**Test strategy:**
- Covered by Step 2; this step's own check is the regenerated artifact plus reading the printed block.

### Step 2: Pin the normalization and its scoping caution
Files: tests/test_active_frame_regrade.py
Depends on: Step 1

**What changes:**
- Assert both fractions re-derive from the source artifacts rather than from literals: pooled from `glm_learned_representation_2026-10-08.json`, active from `glm_active_frame_gate_2026-11-16.json`.
- Assert the active fraction is strictly **less** than the pooled one — the finding that inverts the naive margin reading.
- Assert the comparator-collapse wording names both deltas, and that GC's fall exceeds the encoder's.
- Assert the scoping caution is present and mentions which ceiling each fraction uses.
- Assert the artifact says this is not grounds to revisit G-A.
- Add a non-vacuity test: with a synthetic artifact whose active encoder median is raised above the pooled fraction, the ordering assertion flips — so the ordering is measured, not hardcoded.

**Test strategy:**
- `uv run pytest tests/test_active_frame_regrade.py -q`; all new assertions pass and the non-vacuity test demonstrates the ordering check can fail.

### Step 3: Record the observation as SCOPE on the G-A arm, not as a verdict change
Files: soraya_runs/arms-ledger.json, tests/test_arms_ledger_ga_scope.py
Depends on: Step 1

**What changes:**
- Add one field to the existing `conv-3width-13k-params-150bp` arm dict, following the `metric_not_commensurable_with_siblings` precedent on the fragment arm: the pooled negative is scoped to the pooled frame; within-active ranking was not tested by G-A's bar; on the active frame the margin is +0.2190 while the encoder reaches a *smaller* fraction of its frame's ceiling (16.9% against 34.4%).
- State the revival condition per R7: what observation would re-open the question, namely a candidate that raises `fraction_of_ceiling` on either frame, not one that merely widens a margin over a collapsing comparator.
- Leave `status`, `metric`, `evidence` and `resolved_at` untouched.
- Do not write via `arms_ledger.declare_arm` or `kill_arm`; this is an annotation of an existing arm, and `register`-style helpers are idempotent-by-key and will not update an existing entry.

**Test strategy:**
- A new dedicated test file, owned solely by this step so no same-wave step shares it.
- Assert the arm's `status`, `metric` and `evidence` are byte-identical to their baseline values while the new scope field is present; assert the ledger's top-level key set is still exactly `{"schema_version", "criteria"}`; assert the arm count under the ceiling criterion is unchanged.

### Step 4: Carry the normalization into the human-facing memo
Files: wiki/glm_scoring_frame_resolution.md, tests/test_scoring_frame_artifact.py
Depends on: Step 1

**What changes:**
- Add the normalization beside the memo's existing margin table, so a reader seeing margins grow also sees that the encoder's fraction of achievable signal falls.
- State in one sentence that the margin grows by comparator collapse, and that this is not evidence the learned representation improved.
- Update the memo's pinning test for the added content.

**Test strategy:**
- `uv run pytest tests/test_scoring_frame_artifact.py -q`; assert the memo contains both fractions and the collapse sentence, and that the pre-existing assertions still hold.

### Step 5: Make link 4 launch-ready without pushing it
Files: wiki/finetune_link4_readiness_2026-11-17.json
Depends on: none

**What changes:**
- Build the link-4 kernel with `--ckpt-dataset` pointed at `emanueleebrahimi/ft-ckpt-s3`, `--wall-hours 9.5`, `--horizon-epochs 7`, and run `--self-test`. Capture the result verbatim.
- Choose a fresh slug and record it, asserting in the artifact that it is none of `pa-finetune-arm-v1`/`v2`/`v3`, `pa-finetune-ckpt-probe`/`2`/`3`, `pa-finetune-chain-s1`/`s2`/`s3`.
- Re-read `uv run kaggle quota` and record the reading with its timestamp, rather than carrying a figure forward.
- Record the remaining-work arithmetic from the measured rate: 28 epochs at the pooled 974.7 s/epoch, plus the per-session overhead the kernel incurs.
- Record `pushed: false` and the exact command the user would run, so refresh day is one command.
- Do **not** push. The push is outward-facing, spends the user's quota, and is blocked until 2026-11-21 regardless.

**Test strategy:**
- Manual: inspect the `--self-test` output and the artifact. A failing self-test is the finding — it is exactly what this step exists to surface early, and it would be reported rather than worked around.
- Add a light test asserting the readiness artifact parses, carries `pushed: false`, and that its slug is not on the forbidden list.

## Execution Preview

- **Wave 0** (parallelism 2): Step 1, Step 5
- **Wave 1** (parallelism 3): Step 2, Step 3, Step 4
- Total waves: 2. Max parallelism: 3. Critical path: Step 1 → Step 4 (length 2).
- No intra-wave file overlap, verified by the wave toolkit after Step 3 was moved onto its own test file. Wave 0 touches `scripts/active_frame_regrade.py` plus its artifact against a new readiness artifact; Wave 1 touches three disjoint sets — `tests/test_active_frame_regrade.py`, then `soraya_runs/arms-ledger.json` + `tests/test_arms_ledger_ga_scope.py`, then `wiki/glm_scoring_frame_resolution.md` + `tests/test_scoring_frame_artifact.py`.

## Risk Flags

- **RESOLVED, recorded because it was found by the wave toolkit rather than by reading.** A first draft of this plan had Step 3 and Step 4 both extending `tests/test_scoring_frame_artifact.py`, which the `waves` subcommand flagged as an intra-wave overlap at wave 1. The Risk Flag at that point only *named* a mitigation and left the collision in the step metadata. Step 3 now owns `tests/test_arms_ledger_ga_scope.py` and the overlap is gone from the re-run. The lesson is the reason it stays written down: flagging a defect is not fixing it.
- **The cross-frame fraction comparison sits next to a documented caution.** `wiki/glm_genomewide_oracle_2026-10-07.json` states "Do not compare the two fractions" — that caution is about fractions across *substrates* with different ceilings, not frames of one substrate, so it does not forbid this comparison. But the distinction is thin enough that Step 1 must name which ceiling each fraction uses. [grounded]
- **Do not rewrite `wiki/glm_active_frame_gate_2026-11-16.json`.** It is a recorded run result committed at `eb0963eb`. All additions land in the re-grade artifact, which this session produced. [grounded]
- **`fraction_of_ceiling_for_GENERALISES: 0.50` is a pre-registered bar belonging to a different family's question.** Both frames fall below it (34.4% and 16.9%). The plan records the fractions; it must not apply the oracle's bar to G-A's verdict, which would be importing a threshold across families. [grounded]
- **Step 5's self-test may fail.** That is a successful outcome for this plan, not a blocked step: it converts a quota-refresh-day surprise into a known defect four days early. Report it; do not fix the kernel inside this plan without re-planning.
- **The Architecture Gate is skipped.** sentrux is installed, but the change is artifact and test content across 7 files with no module-graph change, so the structural-change heuristic excludes it.

## Open Questions

- Should the scoped note also land in the G-A family ledger (`project_state/glm-learned-representation-2026-10-07.md`), or is the arms-ledger field sufficient? The review left this as a genuine tradeoff: recording it twice guards against a future session reading only one surface, but annotating a settled valid negative in the family ledger is closer to re-opening it.
- Should the re-grade compute the fractions, or should a small reusable helper live in `dna_decode/eval/preregistration.py` beside `assert_range_reported_with_aggregate`? The plan keeps it local because there is exactly one consumer; a second consumer would justify extracting it.
- The ~8.3 h of unaccounted Kaggle usage is unexplained. Step 5 re-reads the quota but does not diagnose that gap.

## Verification

1. `uv run python -m scripts.active_frame_regrade` regenerates the artifact; the printed block shows both fractions and the collapse reading.
2. `uv run pytest tests/test_active_frame_regrade.py tests/test_scoring_frame_artifact.py -q` passes, including the non-vacuity ordering test.
3. `uv run python -c "import json; d=json.load(open('soraya_runs/arms-ledger.json')); a=d['criteria'][[k for k in d['criteria']][0]]['conv-3width-13k-params-150bp']; print(a['status'], a['metric'])"` prints `measured 0.3234` — unchanged.
4. `uv run pytest tests/ -q` gives a failure set identical to the four baseline files, with the pass count up by the number of new tests.
5. `wiki/finetune_link4_readiness_2026-11-17.json` parses, carries `pushed: false`, a fresh slug not on the forbidden list, and a quota reading with its own timestamp.

<!-- toolkit: check=clean waves=clean gate=fired:open-questions -->
