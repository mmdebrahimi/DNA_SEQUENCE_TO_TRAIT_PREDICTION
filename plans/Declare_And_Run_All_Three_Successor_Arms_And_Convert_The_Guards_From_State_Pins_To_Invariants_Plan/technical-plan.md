# Declare And Run All Three Successor Arms, And Convert The Guards From State Pins To Invariants

## Lens status

Repo-index: ABSENT on this host (`~/.claude/scripts/repo-index-status.sh` does not exist) — direct
file reading only, no degradation.
Tools: `sentrux` present (structural gate eligible), `gh` present, **`codex` ABSENT** — so a
pre-execution `/brainstorm` will fall back to single-model and must be labelled as such.
`project-rules.md`: absent — the Project Rules Check section is omitted.
`DESIGN.md`: absent, and this plan touches no UI.
`AUDIT_REPORT.md`: absent.
Baseline: measured at `b0521d90`, not asserted — see Pre-Change Baseline.
Decision resolved: the user chose the elaborate branch — declare and pursue ALL THREE successor
arms rather than accept `met`. This plan does not re-litigate that.

## Problem Statement

`/execute-plan` on the positional arm finished with verdict `POSITIONAL_NO_GAIN` and then stopped on
one authority call: recording the arm removed the last `planned` arm, so `ceiling_verdict` flipped
`blocked:untried-arm` -> `met`. **That call is now answered: declare all three named successor arms
under R6 and run them, leaving no stone unturned.** This plan does that, and clears the four pieces
of residue the execution left behind.

**The answer changes the shape of the guard problem, and this is the plan's key structural point.**
Eight guard tests went red purely because the verdict flipped. Declaring any arm restores a
non-empty `planned` set, so `ceiling_verdict` returns to `blocked:untried-arm` and **four of the
eight go green on their own**. The remaining four are red for a legitimate, non-verdict reason: the
measured set gained the positional arm and the arm total moves 13 -> 16. With the verdict no longer
flipped, the documented precedent's condition is satisfied and re-pinning becomes an ordinary
re-pin — the sixth — rather than the escalation it would have been an hour ago.

**But a re-pin to exact counts would be red again the moment each new arm measures.** So the guards
are converted from EXACT-STATE pins to INVARIANT pins plus a historical-metric SUBSET pin. That
holds across declaration and across all three measurements, which is why it is done once here
instead of three more times later.

**Three arms, and the research settled that none is a drop-in.** The as-given variant is nearly free
(its features are already on disk). The expression heads are **geometry-locked to 86 bp** against a
150 bp substrate — measured, not assumed. The fine-tune is the long pole and needs a compute route
chosen by measurement.

**Non-goals.** No change to the frozen five, `uv.lock`, `genomewide.py` or `expression.py`. No edit
to `wiki/positional_arm_acceptance_bar.json` or `wiki/positional_arm_result_2026-11-01.json` — the
positional arm's verdict is settled and is never retro-edited. No re-baselining of sentrux (a
separate user call). No self-distillation of any kind: the project's ban stands and nothing here
needs it.

## Codebase Context


See: ./codebase-context.md

## Pre-Change Baseline

Every figure measured at `b0521d90`, not restated from memory:

- **Suite: 21 failed / 7514 passed / 180 skipped.** Diffed by NAME against the 13-failure baseline
  at `33b24322`: exactly 8 new, every one a `blocked:untried-arm` pin.
- `ceiling_verdict`: **`met`**, 12 measured / 0 planned / 1 killed, `best` =
  `promoteratlas-base-embeddings-frozen-150bp-padded-to-200` at **0.5076** (54% of the 0.941
  replicate-derived ceiling).
- Positional arm recorded at **0.4893** (`bins10_mean_oriented`); `autocorr` **0.3975**; the
  ineligible secondary `bins10_mean_asgiven` **0.5086** with paired CI **[-0.0079, +0.0017]**, which
  contains zero.
- Inherited margin **0.02**; shuffled-label null **+0.0053** over 10 seeds, median interval
  **[-0.0067, +0.0169]**.
- Constraint registry holds **14** rows against a pinned `len(rows) == 12`
  (`test_the_eleven_pre_existing_constraints_are_still_present`, red since before this work).
- `scripts/positional_arm_score.py::run` measures **cyclomatic ~58 over 304 lines** — the highest
  thing the last plan added. sentrux `complex_fn_count` 219 -> 223.
- **Heads geometry lock, measured:** all four `fc1.weight` are `(64, 11008)` and `11008 = 128 x 86`.
  A 150 bp tile flattens to `128 x 150 = 19200`. The heads cannot take a 150 bp latent.
- **Heads base weights are byte-identical to the frozen base** — 124/124 tensors, max |dev|
  **0.000e+00**, all four heads. So a heads arm measures FROZEN encoder + trained head.
- Compute: `torch.cuda.is_available()` **False**, 4 threads; local forward throughput **57.5 seq/s**
  measured. Kaggle CLI **2.2.5** present, and `scripts/kaggle/` already holds seven kernel builders.
- Stored positional block present and COMPLETE on `D:` at 44,106 x 5,120 float32 (0.90 GB).

## Verification Signal

The plan has succeeded when ALL of these hold, independent of which way any arm's verdict goes:

1. `ceiling_verdict` reads `blocked:untried-arm` immediately after Step 1 and again after Step 14,
   naming the arms that remain unrun — and the transition is reported, not inferred.
2. The suite's failure set, diffed BY NAME against the 13-failure baseline at `33b24322`, contains
   **zero** arms-guard entries and zero entries from the three files Step 2 touches.
3. `test_the_eleven_pre_existing_constraints_are_still_present` is green, and still fails when a
   NAMED constraint is deleted (proven non-vacuous, not merely passing).
4. Each of the three arms has a bar frozen BEFORE its score exists, and each bar's `verify()`
   refuses on a mutated anchor.
5. Each arm's reachable outcome set EQUALS its check order — not merely is contained in it.
6. The heads arm's windowing decision is frozen in its bar BEFORE any head output is scored, and
   the approximation it introduces is MEASURED against a true per-window recompute on a sample.
7. The fine-tune's compute route and seed/epoch budget are derived from a measured probe, not an
   estimate, and a capped run is stamped PARTIAL and refused by its scorer.
8. `scripts/positional_arm_score.py::run` drops below cyclomatic 25 with the gate ORDER still
   asserted by a test, and `wiki/positional_arm_result_2026-11-01.json` is byte-identical when
   re-derived from the refactored code.
9. The frozen five, `uv.lock`, `genomewide.py` and `expression.py` are byte-identical, verified by
   hash rather than by `git status`.
10. The arms-ledger diff is checked SEMANTICALLY — arms added, arms removed, which fields of which
    arm changed — because a `sort_keys` rewrite hides what moved.

## Implementation Steps

### Step 1: Declare the three successor arms under R6
Files: scripts/declare_successor_arms.py, soraya_runs/arms-ledger.json, tests/test_successor_arms_declared.py
Depends on: none

**What changes:**
- A committed, idempotent declaration script rather than a hand edit, so the three arms' rationales
  and cost tiers are auditable and re-derivable.
- Declares `asgiven-position-resolved-promoteratlas-150bp` (cost `cheap` — its features are already
  extracted), `shipped-expression-heads-zeroshot-150bp` (cost `medium`) and
  `finetuned-promoteratlas-encoder-150bp` (cost `large`).
- Writes each with `status=planned`, `metric=None`, `evidence=""` — the invariant
  `test_every_PLANNED_arm_still_has_an_empty_evidence_string` enforces, and the reason a rationale
  may NOT be parked in `evidence`.
- Each carries a `rationale` and a `declared_from` naming THIS plan. The as-given arm's rationale
  records that it is motivated by a MEASUREMENT (as-given 0.5086 vs oriented 0.4893) rather than a
  preference, and that it ties rather than beats the incumbent.
- Records on the as-given arm that it is NOT a re-read of the positional arm: it changes the
  DECLARED PRIMARY, which is what R6 makes a distinct arm.
- Uses a DIRECT dict write, not `measure_arm`/`declare_arm`'s return path, for the same measured
  reason as before: `Arm(**base)` has no field filter and raises on annotation keys.
- Reports `ceiling_verdict` before and after, and asserts it returns to `blocked:untried-arm`.

**Test strategy:**
- Assert all three arms exist with `planned`/`None`/`""` and carry a non-empty rationale.
- Assert the script is idempotent — a second run changes nothing (no duplicate arms, no status move).
- Assert `ceiling_verdict` is `blocked:untried-arm` and names exactly the three.
- Assert arm total 13 -> 16 and that NO pre-existing arm's `metric` or `status` changed, by semantic
  diff rather than line diff.

### Step 2: Convert the arms guards from exact-state pins to invariant pins
Files: tests/test_arms_ledger_stack_arm.py, tests/test_glm_tile_arms_declared.py, tests/test_tile_chain_substrate_wall_amendment5.py
Depends on: Step 1

**What changes:**
- This is a legitimate RE-PIN, the sixth, and the docstring says so with its reason: the verdict is
  back to `blocked:untried-arm` because Step 1 declared three arms, so the file's own standing rule
  ("had re-pinning flipped the verdict the correct action would have been to stop and escalate") is
  satisfied rather than bypassed. The escalation that blocked the previous plan is recorded as
  having been ANSWERED by the user, not as having been worked around.
- Replaces `assert set(verdict.planned) == EXPECTED_PLANNED` and `len(verdict.planned) == 1` with
  the INVARIANT the guard actually exists to protect: `planned` non-empty **iff** the verdict is
  `blocked:untried-arm`, asserted in both directions.
- Replaces `assert measured == MEASURED` with a SUBSET pin: every historical arm in `MEASURED` must
  still carry its exact metric, and no historical metric may move. New measured arms are permitted
  and counted. This keeps the real protection (a published number cannot drift) while surviving the
  three measurements this plan will add.
- Replaces `assert len(arms) == 13` with `planned + measured + killed == len(arms)` plus a floor of
  the historical count — arithmetic that reconciles the partition with the whole, which is what the
  original comment said it was for.
- Keeps `test_the_human_arms_metric_is_not_read_as_best` exactly as is: it is already an invariant
  and it already passes.
- Renames nothing to a count. The file's own lesson is that a guard whose name misstates what it
  guards is the same drift the project keeps finding in orientation strings.
- `test_a_planned_arm_blocks_the_criterion_in_the_json_itself` carries its own rewrite instruction
  ("replace this with an assertion about the settled verdict instead of deleting it"); that
  instruction is honoured in the invariant form rather than by deletion.

**Test strategy:**
- Each converted guard must FAIL on a synthetic ledger where the invariant is broken — a planned arm
  with a verdict of `met`, a historical metric moved by 0.0001, a partition that does not sum.
  Proven non-vacuous by construction, not asserted.
- Assert the invariant form holds under all three states the plan passes through: before Step 1
  (`met`), after Step 1 (`blocked:untried-arm`, three planned), and after Step 14 (fewer planned).
- Re-run the three files and assert zero failures.

### Step 3: Re-render the memo and amend the journal so no surface claims `met`
Files: wiki/positional_arm_result_2026-11-01.md, project_state/glm-learned-representation-2026-10-07.md
Depends on: Step 1

**What changes:**
- `wiki/positional_arm_result_2026-11-01.md` currently opens with "The ceiling criterion now reads
  `met`, and whether to accept that is yours." That is STALE the moment Step 1 lands. The memo is
  generated from the artifact plus the live verdict, so re-running its generator fixes it — the fix
  is a re-render, never a hand edit.
- Appends ONE Action-Log row recording the decision and its consequence: the authority call was
  answered in favour of the elaborate branch, three arms were declared, and the verdict returned to
  `blocked:untried-arm`, which converted the escalation into an ordinary re-pin.
- Amends the Open Questions entry rather than deleting it — it is ANSWERED, and the answer plus its
  date is the record. Deleting it would erase that a decision was made.
- The positional arm's own verdict, bar and result JSON are NOT touched.

**Test strategy:**
- Assert the rendered memo contains no claim that the criterion is `met` once a planned arm exists,
  and that the DECISIONS section instead names the unrun arms.
- Assert the memo on disk matches a fresh render (the existing guard, re-run).
- Assert the Action-Log row sits INSIDE the `project-state:end:action-log` anchor and that its class
  is from the 7-value enum.

### Step 4: Repair the constraint-count rail
Files: tests/test_tile_chain_constraint_registered.py
Depends on: none

**What changes:**
- `assert len(rows) == 12` is red against a registry of 14 and has been since before this work. It
  is the same defect class as Step 2's: an exact count standing where an invariant belongs.
- The 11 individually-named constraints stay exactly as they are — the docstring is right that "a
  count alone would pass if one were swapped", and the NAMES are the real guard.
- The trailing total becomes a floor tied to the named list (`>= len(NAMED)`) so it cannot go stale
  on the next append, and the two constraints appended since are NAMED so the file keeps its
  named-individually property rather than trading it for a looser bound.
- A drafted `>= 12` fix from an earlier session is superseded by the tied floor, which cannot drift.

**Test strategy:**
- Assert the test is green.
- Assert it still FAILS when any named constraint is removed from a synthetic store — the
  non-vacuity the floor must not cost.
- Assert it does NOT fail on a further append, by appending a synthetic row.

### Step 5: Decompose the scoring runner's gate chain
Files: scripts/positional_arm_score.py, tests/test_positional_arm_score.py
Depends on: none

**What changes:**
- `run` measures cyclomatic ~58 over 304 lines and is the highest thing the last plan added. The
  previous report recorded it as a fair follow-up and declined to touch it mid-execution.
- Extracts one pure-ish function per gate — `gate_split_equality`, `gate_parity`,
  `gate_shuffled_null`, `score_ladder` — each returning its gate dict, with `run` reduced to the
  ORDERED composition plus the stop-on-first-failure rule.
- The ordering stays in ONE place and stays asserted: the extraction must not scatter the
  "gates before the ladder" constraint across four call sites, which is the reason the function was
  left alone before. A module-level `CHECK_SEQUENCE` tuple drives the loop so the order is data.
- Behaviour-preserving by construction: no threshold, statistic or artifact field changes.

**Test strategy:**
- **Byte-identical re-derivation:** re-run the scorer against the committed block and assert the
  emitted artifact equals `wiki/positional_arm_result_2026-11-01.json` field-for-field except
  `wall_seconds` and the RAM readings. This is the test that makes the refactor safe.
- Assert `CHECK_SEQUENCE` order is what the frozen bar's `verdict_check_order` names, so a reordered
  refactor fails rather than silently scoring before a gate.
- Re-run the existing 43 scorer tests unchanged — they already drive every gate branch, so they are
  the regression suite for this refactor.
- Assert the measured cyclomatic complexity of `run` is below 25.

### Step 6: Freeze the as-given arm's bar
Files: scripts/freeze_asgiven_arm_bar.py, wiki/asgiven_arm_acceptance_bar.json, tests/test_freeze_asgiven_arm_bar.py
Depends on: Step 1

**What changes:**
- Freezes a bar for `asgiven-position-resolved-promoteratlas-150bp` BEFORE its score exists, even
  though the number is already visible as a secondary (0.5086). **That visibility is the hazard this
  bar exists to manage, and it is stated in the bar rather than hidden:** the arm is being declared
  AFTER its rung was observed, so this is not a blind pre-registration and must never be presented
  as one.
- The honest consequence, pre-committed: the arm's job is NOT to re-report 0.5086 as a win. Its
  interval already contains zero, so the bar's PREDICTED outcome is a TIE, and the bar names
  `ASGIVEN_TIES` as the expected terminal rather than `ASGIVEN_ADDS`.
- INHERITS the margin 0.02, the committed atlas comparator, the 10-seed `leave_peak_out(frac=0.25)`
  split and `ridge_spearman` unchanged — all re-read with a `verify()` that refuses on drift.
- Single primary (`bins10_mean_asgiven`), so NO multiplicity correction: one comparison at the plain
  95% level. Declaring a second rung here would re-open the forking path the positional bar closed.
- Closed outcome set with the reachable set required to EQUAL the check order:
  `ASGIVEN_ADDS` / `ASGIVEN_TIES` / `ASGIVEN_NO_GAIN` / `INDETERMINATE_SPLIT_MISMATCH` /
  `INDETERMINATE_PARITY_FAIL` / `INDETERMINATE_NULL_FAIL`.

**Test strategy:**
- Assert `verify()` refuses on a mutated comparator, a mutated arm status, and a missing arm.
- Assert the reachable-outcome set EQUALS the check order.
- Assert the bar records the observed-before-declared hazard explicitly — a bar that omitted it
  would read as a blind pre-registration.
- Assert the margin equals the ledger-derived inherited value.

### Step 7: Score the as-given arm from the existing block
Files: scripts/asgiven_arm_score.py, wiki/asgiven_arm_result_2026-11-02.json, tests/test_asgiven_arm_score.py
Depends on: Step 5, Step 6

**What changes:**
- Reuses the COMMITTED stored block on `D:` — no second forward pass — so this arm costs minutes.
- Runs the same three gates in the same order via the helpers Step 5 extracted, then scores the
  single primary with a plain 95% paired bootstrap.
- Emits the verdict from the closed set plus a re-derivation function, so the verdict is a pure
  function of the artifact.
- Writes per-tile predictions under `data/`, never `wiki/`, for the recorded decoy-pool reason.

**Test strategy:**
- Drive every `INDETERMINATE_*` branch with a stub scorer and assert each is reachable.
- Assert the parity gate trips on a perturbed block (non-vacuous).
- Assert the verdict is re-derivable from the emitted artifact alone, round-tripped through JSON.
- Assert the scored median reproduces the positional run's secondary value for the same rung and
  seeds — the two runs must agree on the same number, or one of them is wrong.

### Step 8: Expression-heads geometry and circularity probe
Files: scripts/expression_heads_probe.py, wiki/expression_heads_feasibility_2026-11-02.json, tests/test_expression_heads_probe.py
Depends on: Step 1

**What changes:**
- Validity first, and this probe has more to establish than its siblings because the arm is NOT a
  drop-in.
- **Re-measures the geometry lock live:** asserts every head's `fc1.weight` is `(64, 128 * 86)` and
  that a 150 bp latent therefore cannot be fed directly. Records the arithmetic.
- **Re-measures that each head's `base_model.*` weights are byte-identical to the frozen base**
  (124/124 tensors, max |dev| 0.000e+00 when probed during research). If a head ever diverges, the
  arm stops being "frozen encoder + head" and becomes a fine-tuned-encoder arm, which is a different
  claim and a different arm.
- Measures the two candidate window reductions and the error the cheap one introduces: applying
  `conv1/conv2/pool` to the full 150 latent once and then slicing 86-wide windows before `fc1` is
  NOT identical to a true per-window recompute, because `padding='same'` makes the edges differ. The
  probe computes both on a sample and reports max and median absolute deviation, plus the Spearman
  between them — so the bar can choose from a measurement.
- Records `normalization_params` {mean, std} per head and notes that Spearman is invariant to the
  denormalisation, so it cannot move the verdict — but an absolute predicted value would be wrong
  without it.
- **Records the circularity scope, which is narrower than it first looks:** the `urtecho2018` head
  IS circular on the sigma70 rung (GSE108535) and must never be a comparator there. This arm's
  substrate is GSE144621 peak tiles, and the repo's own withdrawn lead established that no
  genome-wide/tiling Urtecho head exists, so the head is not circular HERE. Recorded as a limit so
  nobody reuses it on the wrong rung.
- Exits non-zero on any failure.

**Test strategy:**
- Unit-test the pure window-slicing arithmetic on synthetic latents with known answers.
- Assert the probe exits non-zero on a synthetic head whose base weights differ from the base.
- The model-touching path skips cleanly when the clone does not resolve.

### Step 9: Freeze the expression-heads bar
Files: scripts/freeze_expression_heads_bar.py, wiki/expression_heads_acceptance_bar.json, tests/test_freeze_expression_heads_bar.py
Depends on: Step 8

**What changes:**
- Freezes the WINDOWING DECISION from Step 8's measurement, before any head output is scored. If the
  cheap full-latent reduction tracks the true per-window recompute within a stated tolerance, it is
  used and the tolerance is recorded; otherwise the true recompute is used and its cost is recorded.
  The choice is derived, not preferred.
- Declares the reduction over the 65 windows (mean of head outputs) as the single PRIMARY, with
  max-over-windows as a pre-declared secondary that is INELIGIBLE as the arm's result.
- **Four heads is a multiplicity of four**, and it is corrected in the bar: each head is judged at
  `alpha/4`, and the arm's metric is the MAX OVER THE DECLARED HEADS — never over heads x reductions.
- Pre-commits the zero-shot framing: these heads were trained on OTHER libraries, so this arm
  measures TRANSFER to the tile rung, and a null result bounds transfer rather than the heads.
- Closed outcome set with reachable == check order, including `INDETERMINATE_GEOMETRY_FAIL` for the
  case where the windowing cannot be validated.

**Test strategy:**
- Assert `verify()` refuses on a mutated probe artifact and on a missing arm.
- Assert the reachable-outcome set EQUALS the check order.
- Assert the Bonferroni divisor equals the number of declared heads, and that the secondary cannot
  enter the verdict function.

### Step 10: Score the expression-heads arm
Files: scripts/expression_heads_score.py, wiki/expression_heads_result_2026-11-02.json, tests/test_expression_heads_score.py
Depends on: Step 9

**What changes:**
- ONE forward pass over the 44,106 tiles computing the shared base latent, then all four heads
  applied to that same latent — the same one-pass-all-rungs design, justified here by the measured
  fact that all four heads share byte-identical base weights.
- Gates first: split equality, then a HEAD-SPECIFIC validity check (each head must produce finite,
  non-constant output over the cohort — a collapsed head would score as a clean null), then the
  shuffled-label null.
- Scores each head's windowed output against the committed atlas per-seed values on byte-identical
  splits, with the bar's corrected interval.
- Checkpointed by batch, PARTIAL-stamped under `--limit`, and the loader refuses a PARTIAL block.

**Test strategy:**
- Drive every `INDETERMINATE_*` branch with a stub head and assert reachability.
- Assert a constant-output head trips the head validity gate rather than scoring ~0.
- Assert the verdict is re-derivable from the artifact alone.
- Assert a head whose output is monotonically transformed scores an IDENTICAL Spearman, which is the
  property that makes the denormalisation irrelevant to the verdict.

### Step 11: Fine-tune feasibility and compute-route probe
Files: scripts/finetune_arm_probe.py, wiki/finetune_arm_feasibility_2026-11-02.json, tests/test_finetune_arm_probe.py
Depends on: Step 1

**What changes:**
- The fine-tune is the long pole and its budget must be MEASURED, not estimated. This probe measures
  training throughput (forward+backward) on one seed for a bounded number of steps, and reports
  seconds/epoch and the extrapolated cost for the bar's seed x epoch grid.
- Measures epochs-to-early-stop on ONE seed using a peak-blocked inner split, so the bar's epoch
  budget comes from an observed convergence curve.
- **Chooses the compute route from the measurement.** Local is CPU-only (`cuda False`, 4 threads);
  the project's established GPU route is a Kaggle T4 kernel, and `scripts/kaggle/` already holds
  seven builders. The probe reports both costs and the bar freezes the route. A local run that would
  exceed a named wall-clock ceiling routes to Kaggle rather than being started and abandoned.
- Records that the encoder is only 1,847,428 parameters, which is why a fine-tune is tractable here
  at all and is NOT comparable to the project's blocked large-LM fine-tunes.
- Asserts the base checkpoint's sha256 still matches the one the positional probe recorded, since
  the clone lives in a volatile temp directory.

**Test strategy:**
- Unit-test the pure extrapolation and route-selection arithmetic on synthetic timings.
- Assert the route selector picks Kaggle when the local extrapolation exceeds the ceiling, and local
  when it does not — both directions, so the rule is not vacuous.
- Assert the probe exits non-zero when the checkpoint hash has moved.

### Step 12: Freeze the fine-tune bar
Files: scripts/freeze_finetune_arm_bar.py, wiki/finetune_arm_acceptance_bar.json, tests/test_freeze_finetune_arm_bar.py
Depends on: Step 11

**What changes:**
- Freezes the seed count, epoch budget, learning rate, early-stopping rule and compute route from
  Step 11's measurements, before any score exists.
- INHERITS the margin 0.02 and the committed atlas comparator; the comparison is the fine-tuned
  encoder's median Spearman against the FROZEN encoder's 0.5076 on byte-identical splits.
- **Pre-commits the two readings that a fine-tune invites and that must not be conflated:** beating
  0.5076 shows the frozen representation was the limit; failing to shows the SUBSTRATE or the
  training budget was, and the bar records that a budget-limited null cannot be read as a
  representation-limited one. The measured epoch curve from Step 11 is what distinguishes them.
- Declares a shuffled-label null on the fine-tune itself — a trained encoder with destroyed labels
  must not score, and this is the arm where that control matters most.
- Closed outcome set with reachable == check order, including `INDETERMINATE_BUDGET_EXHAUSTED` for a
  run that stopped before its frozen epoch budget.

**Test strategy:**
- Assert `verify()` refuses on a mutated probe artifact, a mutated comparator and a missing arm.
- Assert the reachable-outcome set EQUALS the check order.
- Assert the budget-limited-versus-representation-limited distinction is present in the bar text and
  tied to the probe's epoch curve, not left to the reader.

### Step 13: Run and score the fine-tune arm
Files: scripts/finetune_arm_run.py, scripts/kaggle/build_finetune_arm_kernel.py, wiki/finetune_arm_result_2026-11-02.json, tests/test_finetune_arm_run.py
Depends on: Step 12

**What changes:**
- Trains the encoder per seed on byte-identical `leave_peak_out(frac=0.25)` splits, checkpointing
  per seed so a crash costs one seed rather than the run.
- The Kaggle kernel builder pins `machine_shape: NvidiaTeslaT4` explicitly, because
  `enable_gpu: true` alone provisions a P100 that current Kaggle torch cannot run in fp16 — a
  recorded trap. It pushes to a FRESH slug, because a slug holding retrieved output must never be
  re-pushed.
- Emits per-seed scores, the paired bootstrap, the null, and the single verdict from the closed set.
- A capped run is stamped PARTIAL and the scorer refuses it, so a smoke run cannot be read as the
  arm's result.

**Test strategy:**
- Drive the gate chain with a stub trainer and assert every branch is reachable.
- Assert resuming after an interrupted run reproduces the same per-seed scores as an uninterrupted
  one, and that only the unfinished seeds are retrained.
- Assert the PARTIAL stamp is present under `--limit` and that the scorer refuses it.
- Assert the kernel metadata pins the T4 explicitly, and that the slug is not one that already holds
  output.

### Step 14: Record the three arms, render the memos, journal it
Files: soraya_runs/arms-ledger.json, scripts/successor_arms_memo.py, wiki/successor_arms_result_2026-11-02.md, project_state/glm-learned-representation-2026-10-07.md
Depends on: Step 2, Step 3, Step 7, Step 10, Step 13

**What changes:**
- Moves each arm `planned -> measured` by direct dict write with its metric and an evidence string
  naming the verdict, the gates, the limits and the artifact — preserving annotation keys.
- Each arm's metric is the max over ITS OWN declared primaries, never over a ladder and never across
  arms.
- Recomputes `ceiling_verdict` and REPORTS it. If every declared arm is now measured it will read
  `met` again, and that is the same authority call surfacing a second time — the memo states it and
  the run stops rather than recording it in passing.
- Renders one memo covering all three arms FROM their artifacts, ordered by required reader action,
  carrying the observed-before-declared hazard on the as-given arm, the geometry approximation on
  the heads arm, and the budget-versus-representation distinction on the fine-tune.
- Appends one Action-Log row per arm with a class from the 7-value enum, anchored inside the Action
  Log section.
- Re-runs the Step 2 invariant guards and asserts they are STILL green after three status moves,
  which is the whole point of having converted them.

**Test strategy:**
- Assert the ledger round-trips with every annotation key preserved, by semantic diff.
- Assert the memo generator refuses to render any arm on an `INDETERMINATE_*` verdict.
- Assert the Step 2 invariants hold with three more measured arms — the regression the conversion
  exists to prevent.
- Re-run the full suite and diff the failure set BY NAME against the 13-failure baseline.

## Execution Preview

- **Wave 0:** Step 1, Step 4, Step 5 (parallelism 3)
- **Wave 1:** Step 2, Step 3, Step 6, Step 8, Step 11 (parallelism 5)
- **Wave 2:** Step 7, Step 9, Step 12 (parallelism 3)
- **Wave 3:** Step 10, Step 13 (parallelism 2)
- **Wave 4:** Step 14

Total waves 5, max parallelism 5. Critical path 1 -> 8 -> 9 -> 10 -> 14 and 1 -> 11 -> 12 -> 13 -> 14,
both length 5.

**Sequential execution is recommended over parallel worktrees, for the same reason as the preceding
plan.** The stored feature blocks are gitignored and live on `D:`, so a fresh worktree could not run
the parity or geometry gates at all and would silently skip them. Wall-clock is dominated by Step 10
(one forward pass, ~13 min) and Step 13 (the fine-tune, whose budget Step 11 measures).

## Architecture Impact

`sentrux` is present, so the structural gate is eligible. The change adds eight scripts and one
Kaggle kernel builder, with NO new package-level module, so the import-edge count moves by
construction and an aggregate compare cannot distinguish an intended edge from an unintended one.
Step 5 should REDUCE `complex_fn_count`, which is the one metric the last compare flagged
(219 -> 223); this plan's own refactor is the remedy for it. Replacement check, mirroring what the
sibling plans adopted: assert that the set of files importing `positional_readout` is exactly the
intended consumer set, and that `genomewide` gains no new in-package importer.

## Risk Flags

- **The as-given arm is declared AFTER its rung was observed, and that is a real weakness, not a
  formality.** 0.5086 is already visible as a secondary, so its bar is not a blind
  pre-registration. The mitigation is to state it in the bar, predict a TIE rather than a win, and
  keep the single-primary rule so no second rung can be fished. The residual risk is that a reader
  quotes 0.5086 as a positional win; the memo must pre-empt that in its own text.
- **The heads arm's cheap windowing is an approximation and may not be usable.** `padding='same'`
  makes a full-latent reduction differ from a true per-window recompute at the edges. Step 8
  measures it and Step 9 freezes the choice; if the deviation is large the arm costs 65x more and
  that cost lands in Step 10.
- **A fine-tune null is ambiguous between budget and representation**, which is the single most
  over-readable outcome in this plan. Step 11's measured epoch curve is what separates them, and
  Step 12 pre-commits the distinction. Without that curve the arm should not be run.
- **The checkpoint lives in a volatile temp directory** (`.../Temp/pa_probe/PromoterAtlas`). Steps 8
  and 11 re-verify its sha256 against the recorded value; a reboot or cleaner invalidates the whole
  plan's model-touching half and a re-clone is the recovery.
- **Clearing all three arms will flip `ceiling_verdict` to `met` again.** Step 14 surfaces it rather
  than recording it, and the user has now answered this question once — a second `met` is a NEW
  decision about whether the arm set is finally exhausted, not a re-run of the first.
- **Step 5 refactors a script whose output is already published.** The byte-identical re-derivation
  test is the only thing that makes it safe; without it a behaviour-preserving claim is unverified.
- **Multiplicity grows across arms, and the plan does NOT pool it.** Each arm carries its own
  correction; there is no family-wise correction across the three arms, so the probability that at
  least one of three arms clears by chance is higher than any single arm's alpha. Stated here
  because the arms are declared together, which is exactly when that temptation appears.
- **`codex` is absent**, so a pre-execution `/brainstorm` is single-model and must be labelled as
  such rather than presented as a two-engine review.

## Open Questions

- Whether a `met` verdict after all three arms measure should finally be accepted as the criterion's
  terminal, or whether the arm set should keep growing. Deferred to Step 14, which surfaces it with
  the measured numbers in hand rather than deciding it now.
- Whether the fine-tune should run on Kaggle or locally. Deferred to Step 11, which measures both
  and freezes the route from the measurement rather than from preference.
- Whether sentrux should be re-baselined to absorb the new shape. Independent of this plan's
  correctness, still the user's call, and Step 5 may make it unnecessary by reducing the one metric
  that was flagged.

## Verification

- `uv run pytest tests/ -q` with the failure set diffed BY NAME against the 13 at `33b24322`; zero
  arms-guard entries is the pass condition, not a smaller count.
- `uv run python scripts/declare_successor_arms.py` twice, asserting the second run is a no-op.
- `uv run python scripts/freeze_asgiven_arm_bar.py` and the other two freeze scripts exit 0, and each
  exits non-zero on a mutated anchor.
- `uv run python -m scripts.expression_heads_probe` exits 0 and its artifact reports the geometry
  lock arithmetic and the base-weight identity.
- `uv run python -m scripts.positional_arm_score` re-emits an artifact byte-identical to the
  committed one except `wall_seconds` and the RAM readings.
- `git diff --stat` shows the frozen five, `uv.lock`, `genomewide.py` and `expression.py` untouched,
  verified by hash.
- The arms-ledger diff checked SEMANTICALLY — arms added, arms removed, which fields of which arm
  changed.
- `ceiling_verdict` recomputed and reported after Step 1 and again after Step 14, with `best` named
  both times.

## Save-time amendments

Captured at: 2026-11-02
Source: `/save-plan` arguments

**Audit notes only.** `/execute-plan` reads executable work from `## Implementation Steps`
alone; nothing in this section is an instruction. If an amendment ever changes a step's file
list, dependencies or structure, re-validation fires and the plan must go back through
`/technical-plan` before execution.

- Decision answered: declare and run all three successor arms rather than accept `met`; guards converted to invariants rather than re-pinned to counts
<!-- toolkit: check=clean waves=clean gate=fired:open-questions -->
