# Bound The Training Variance That Dominated, And Decide If The Delta Is Resolvable At All

## Lens status
Repo-index: absent — `repo-index-status.sh` returned nothing, so research was direct file reading with no Layer 2 degradation.
Tools: `sentrux` present (plan is non-structural, so the gate will skip), `gh` present, `codex` ABSENT — the pre-exec review will be single-voice again, which for a plan about variance decomposition is a named weakness.
Baseline tests: 12 failed / 9054 passed / 191 skipped, all 12 in four documented absent-data families.
Predecessor: `executed_plans/Measure_The_Generator_Noise_Floor_As_A_Paired_Delta_In_One_Session_Plan/` (as-built, 2026-11-20). This plan executes the follow-on that one explicitly scoped OUT.

## Problem Statement

The paired sweep returned **`DELTA_NOT_ESTABLISHED`** for the carried claim that LoRA fine-tuning
makes the generator's distinguishability worse by +0.0271: mean paired delta +0.0458 against a
spread of 0.2072, 3 of 10 seeds against the committed sign, and at seed 7 the delta was **-0.0142**,
the opposite sign to the committed value at the same seed.

Its pre-registration had recorded, in advance, that the design varies GENERATION only and leaves
**TRAINING variance unmeasured**. That caveat named the dominant term:

- zero-shot at seed 7 reproduced the committed **0.7787 EXACTLY**, so the scoring path is
  deterministic and the instrument is sound
- epoch 3 at seed 7 came out **0.7645** against the committed **0.8058** — a drift of **0.0413**,
  *larger than the 0.0271 effect under test*
- the loss curves differed (**7.43638 / 7.16777 / 6.91858** vs **7.4361 / 7.1655 / 6.9148**) under
  the same method, same data and same generation seed

**Measured cause, not speculation.** There is **no torch seeding anywhere in the kernel** — no
`manual_seed`, no `set_seed` — while the training data order IS deterministic (`for i in range(0,
len(seqs), BATCH)`, no shuffle). So the drift comes from LoRA adapter initialisation plus cuDNN
nondeterminism, and training runs are currently **unlabelled and irreproducible**.

The question this plan answers is not "is the delta positive". It is **whether a delta of this size
is resolvable at all under single-run training**, by separating two variance components that the
previous design confounded into one:

1. **measurement variance** — spread across generation seeds at a FIXED training run
2. **training variance** — spread across training runs at a FIXED generation seed

If (2) alone exceeds the effect, the single-run fluency-delta comparison is **closed as
unmeasurable at this geometry**, which is a durable methodological result for the whole
`glm-generator-finetune` family rather than another refinement of a null.

Non-goals: changing the frozen bar (0.6679 stands); re-reading the surviving null-comparison claim
(the generator is worse than its own seed-matched null at 10 of 10 seeds, mean gap +0.1186 — this
plan does not touch that); forcing bit-deterministic training; raising `n` from 120.

## Codebase Context


See: ./codebase-context.md

## Pre-Change Baseline

- Paired sweep verdict: **`DELTA_NOT_ESTABLISHED`**, mean paired delta **+0.0458**, spread
  **0.2072**, 3 of 10 seeds against the committed sign, delta at seed 7 **-0.0142**
- Trajectory: **`TRAJECTORY_UNRESOLVED`** (step means +0.0016 / +0.0187 / +0.0255, positive counts
  6 / 8 / 9)
- Null comparison, SURVIVING and untouched by this plan: generator worse than its own seed-matched
  null at **10 of 10** seeds, mean gap **+0.1186**
- Committed arms: zero-shot **0.7787**, fine-tuned **0.8058**, both at generation seed 7
- Observed training drift at a fixed generation seed: **0.0413** (one observation, n=1 — which is
  precisely why this plan exists)
- Training variance across repeated runs: **NEVER MEASURED.** No torch seeding exists, so no two
  training runs are even labelled.
- Measured per-unit costs: scoring **138 s**, 3-epoch training **578 s**
- Test suite: 12 failed / 9054 passed / 191 skipped; the 12 are cyp4f2/GeT-RM 7,
  essentiality-ladder 3, README-quickstart 1, pneumo cps DB 1

## Verification Signal

- The emitted body parses for every flag combination, and `build_body` still REFUSES an
  unsubstituted placeholder (the new `__TRAIN_SEEDS__` must be wired into it).
- `scripts/kaggle/glm_finetune_kernel.py` stays byte-identical, sha256 **6bdbd02641abedc5** — it is
  the provenance record of the two committed arms.
- `verdict()`, `EPOCHS`, `LR`, `BATCH`, `DTYPE`, `LORA` and `RECORDED` unchanged by diff.
- Two training runs at the SAME training seed produce the same loss curve to at least 3 decimals
  (labelling works); two runs at DIFFERENT training seeds do not (the seeds are distinct).
- The run returns an epoch-3 score for every (training seed x generation seed) cell with no gaps.
- The reused zero-shot row is VALIDATED, not assumed: one zero-shot cell is re-scored in-run and
  must match the committed value exactly, or the reuse is refused.
- The decomposition reports both components separately and never a single pooled number.
- Frozen five byte-unchanged; `prospective_lock.verify_lock` returns `ok=True`.
- No regression against the 12-failure baseline, reconciled by CAUSE not count.

## Implementation Steps

### Step 1: Pre-register the decomposition reading
Files: wiki/glm_training_variance_bar_2026-11-21.json, tests/test_glm_training_variance_bar.py
Depends on: none

**What changes:**
- A frozen bar, written before any training-variance number exists, naming the decomposition and
  every outcome.
- Decomposition declared explicitly: for epoch-3 score `s(t, g)` at training seed `t` and generation
  seed `g`, and reused zero-shot `z(g)`, the delta is `d(t, g) = s(t, g) - z(g)`.
  **`training_spread(g)`** = spread of `d(t, g)` over `t` at fixed `g`; **`measurement_spread(t)`** =
  spread of `d(t, g)` over `g` at fixed `t`. Headline statistics are the MEDIAN of each across the
  other index, so one outlier cell cannot carry either component.
- Outcomes, ordered, first match wins: `INSTRUMENT_INVALID` (null fails to reproduce at seed 7, any
  cell missing, or the zero-shot reuse check fails); `TRAINING_VARIANCE_DOMINATES` (median
  `training_spread` >= |0.0271|, so a single-run delta cannot resolve the effect — the question is
  CLOSED at this geometry); `MEASUREMENT_VARIANCE_DOMINATES` (median `measurement_spread` > median
  `training_spread` AND median `measurement_spread` >= |0.0271|); `BOTH_SMALL_EFFECT_RESOLVABLE`
  (both medians < |0.0271|).
- A note that `TRAINING_VARIANCE_DOMINATES` is the outcome that CLOSES a question rather than
  refining one, so it must not be read as a failure of the run.
- `claim_status: "DIAGNOSTIC_METHOD_VARIANT"` and the same `not_the_frozen_arm` disclaimer: this
  measures the variance of the frozen method, it does not re-run the frozen arm.
- An explicit statement of what THIS design does not bound in turn: generation-seed count (5, not
  10), epoch coverage (epoch 3 only), and `n` (120).
- A forecast with two-sided reasoning.

**Test strategy:**
- Assert the file parses, carries all four outcomes in a declared order, defines both spreads as
  medians over the opposite index, and states the closing-not-failing note.
- Assert no outcome condition is decided by anything other than a comparison against |0.0271| or
  against the other component — i.e. no newly invented constant.
- Assert the design's own unbounded terms are listed.

### Step 2: Label training runs with a seed, and declare the training-seed set
Files: scripts/kaggle/build_glm_finetune_kernel.py
Depends on: none

**What changes:**
- `TRAIN_SEEDS = (101, 102, 103, 104, 105)` as a builder-level constant, injected via a new
  `__TRAIN_SEEDS__` placeholder wired into `build_body` — the guard will refuse the body otherwise.
  Deliberately disjoint from `SEED_SET` (0-9) so a training seed can never be mistaken for a
  generation seed in a grid cell.
- `step4_lora_shapecheck` gains a keyword-only `train_seed: int | None = None`; when not None it
  calls `torch.manual_seed(train_seed)` and `torch.cuda.manual_seed_all(train_seed)` immediately
  BEFORE `get_peft_model`, which is where the adapter's random `lora_A` init happens. Default None
  leaves the existing path byte-identical in behaviour.
- **Deterministic algorithms are deliberately NOT forced.** `torch.use_deterministic_algorithms(True)`
  would change what is being measured: residual cuDNN nondeterminism IS part of training variance,
  and this plan measures that variance rather than suppressing it. Seeding the init is what makes
  runs LABELLED and DISTINCT, which is all the decomposition needs. Recorded in the kernel comment
  so a later reader does not "fix" it.

**Test strategy:**
- Assert `TRAIN_SEEDS` has no duplicates, has at least 3 entries, and is disjoint from `SEED_SET`.
- Assert `build_body` still refuses an unsubstituted placeholder after `__TRAIN_SEEDS__` is added.
- Assert `torch.manual_seed` appears immediately before `get_peft_model` in the real training path
  and only under the `train_seed is not None` guard.
- Assert the frozen kernel artifact's sha256 is unchanged (module-scoped snapshot, as established).

### Step 3: The training-variance variant
Files: scripts/kaggle/build_glm_finetune_kernel.py
Depends on: Step 2

**What changes:**
- A new builder flag `--training-variance` emitting `scripts/kaggle/glm_training_variance_kernel.py`,
  mutually exclusive with `--diagnostic` and `--paired-sweep` (the existing exclusion check extends
  to three flags).
- For each `t` in `TRAIN_SEEDS`: seed, build the adapter, train once under the frozen method, save
  the per-epoch adapters, merge epoch 3, and score it at every `g` in a 5-seed generation subset
  that MUST contain 7. Each cell is written to the grid immediately, carrying `train_seed`,
  `gen_seed`, `epoch`, `both`, `kmer` and that cell's own `null`.
- **The zero-shot row is REUSED, not re-scored** — it is training-independent (line 447 scores the
  pristine `MODEL_ID`) and 10 cells are already committed. One zero-shot cell at generation seed 7
  is re-scored in-run as a REUSE CHECK and must match 0.7787 exactly; a mismatch aborts rather than
  silently mixing two measurement regimes. This saves 9 scorings, about 21 minutes.
- The loss curve for every training seed is recorded, so same-seed reproducibility and
  different-seed distinctness are both checkable from the artifact.
- Its own reading function; `verdict()` and `paired_reading` are not called.
- Cost at 5 training seeds x 5 generation seeds: 5 trainings (578 s) + 25 scorings (138 s) + 1 reuse
  check + model loads, about **1.9 h**, inside one session.

**Test strategy:**
- Build with `--training-variance`, assert the body AST-parses and contains one training call per
  `TRAIN_SEEDS` entry and one scoring call per (training seed x generation seed) cell.
- Assert the generation subset contains 7 and is a subset of `SEED_SET`.
- Assert the zero-shot reuse check is present and that an exact-match failure aborts.
- Assert the three variant flags are mutually exclusive pairwise.
- Assert `EPOCHS`, `LR`, `BATCH`, `DTYPE`, `LORA` unchanged in the emitted body.

### Step 4: Tests for the variant and the labelling
Files: tests/test_glm_finetune_kernel_builder.py
Depends on: Step 3

**What changes:**
- Add the Step 2 and Step 3 strategies as real tests in the existing module, which already carries
  the module-scoped sha256 no-write invariant on the frozen artifact.
- Extend the existing placeholder-guard non-vacuity test to cover `__TRAIN_SEEDS__`.
- Add a test that the emitted training-variance artifact on disk matches a fresh build of its own
  code, with the base64 payload normalised out (a full-body compare fails for payload reasons, as
  established in the predecessor plan).

**Test strategy:**
- `uv run pytest tests/test_glm_finetune_kernel_builder.py -q` green.
- Each new load-bearing assertion confirmed to FAIL when its subject is reverted, by mutating the
  emitted body and re-running the assertion — not by reading the code.

### Step 5: Run the training-variance grid and commit it verbatim
Files: scripts/kaggle/glm_training_variance_kernel.py, wiki/glm_training_variance_grid_2026-11-21.json
Depends on: Step 1, Step 4

**What changes:**
- Build, push to the NEW slug `glm-training-variance-v1` with `--gpu` and internet ON (the kernel
  pip-installs peft and pulls the model from HuggingFace), poll, pull.
- NEVER re-push: `pa-finetune-arm-v1`/`v2`/`v3`, `pa-finetune-ckpt-probe`/`probe2`/`probe3`,
  `pa-finetune-chain-s1`/`s2`/`s3`/`s4`, `glm-generator-finetune` (the committed arms' provenance),
  or `glm-paired-sweep-v1` (the predecessor grid's provenance).
- Commit the pulled grid verbatim BEFORE any interpretation, as the predecessor did — that
  separation is what let a spoiled in-kernel verdict be corrected without touching the data.
- If the kernel errors with an empty failureMessage and zero output files, record it and STOP rather
  than retrying blind.
- If the session is cut mid-grid, the per-cell records and the saved adapters both survive in the
  collected output, so scoring can resume without retraining.

**Test strategy:**
- Manual verification that the pulled JSON carries one cell per (training seed x generation seed)
  with no gaps, that every loss curve is present, and that the zero-shot reuse check passed.
- Assert the null reproduces 0.6679 at generation seed 7; anchor there ONLY, never per-cell.

### Step 6: Decompose the variance and write the verdict
Files: scripts/glm_training_variance_verdict.py, wiki/glm_training_variance_2026-11-21.md, wiki/glm_training_variance_verdict_2026-11-21.json
Depends on: Step 5

**What changes:**
- A script that reuses `glm_paired_sweep_verdict.py`'s grid helpers, joins the new grid against the
  reused zero-shot row, computes `training_spread(g)` and `measurement_spread(t)`, takes the median
  of each, and applies Step 1's ordered rule with no post-hoc threshold.
- Verify same-seed reproducibility and different-seed distinctness from the recorded loss curves,
  and report them — if the seeds did not actually separate the runs, the decomposition is void and
  the verdict is `INSTRUMENT_INVALID` regardless of the numbers.
- The memo reports the two components SEPARATELY and never pools them into one "variance".
- On `TRAINING_VARIANCE_DOMINATES`: state plainly that the single-run fluency-delta comparison is
  closed at this geometry, that this extends to the whole `glm-generator-finetune` family, and that
  it does NOT touch the surviving null-comparison claim.
- Record limits: 5 training seeds and 5 generation seeds bound both components weakly; epoch 3 only;
  `n=120`; one organism; and the k-mer discriminator's long-range blindness, so a low score would
  remain necessary-not-sufficient.
- Amend Step 1's bar augment-only if the run exposes a defect in it, proving augment-only against
  git HEAD rather than asserting it — as the predecessor did.

**Test strategy:**
- A test asserting the memo's stated verdict is the one Step 1's rule yields from the committed
  grid, so the prose cannot drift from the computation.
- A test asserting the memo reports both components separately and carries the epoch-3-only and
  n=120 limits.

## Execution Preview

5 waves, max parallelism 2, critical path 5 steps (Step 2 -> Step 6).

- Wave 0: Step 1 (pre-registration + its tests), Step 2 (training-seed labelling) — disjoint files
- Wave 1: Step 3 (the variant) — same builder file as Step 2, hence sequenced
- Wave 2: Step 4 (builder tests)
- Wave 3: Step 5 (the run, ~1.9 h attended)
- Wave 4: Step 6 (decomposition + verdict)

Sequential execution is expected again for the predecessor's reasons: gitignored `data/` makes a
worktree's absent-data failure set diverge from the baseline this plan reconciles against, and the
Kaggle SDK lives in the gitignored `.venv` while the run is one host-global session.

## Risk Flags

- **The outcome that closes the question is also the least falsifiable-looking.** If
  `TRAINING_VARIANCE_DOMINATES` fires, the conclusion is "this comparison cannot be made", which is
  unfalsifiable in the wrong hands. Mitigated by the ordered rule being fixed in advance, by
  reporting both components separately, and by the loss-curve reproducibility check — if the
  training seeds did not separate the runs, the decomposition is void and says so.
- **5 training seeds is a weak bound on a spread.** A median over 5 is noisy, and the plan reports
  the median precisely so one outlier cannot carry it, but the honest limit is that this
  distinguishes "larger than the effect" from "smaller than the effect", not a precise variance.
- **Seeding the init may not fully separate runs** if cuDNN nondeterminism dominates the adapter
  init. That is itself a finding, and the loss-curve check surfaces it rather than letting the
  decomposition proceed on unseparated runs.
- **Reusing the zero-shot row mixes two runs' measurements.** Mitigated by an in-run exact-match
  reuse check at generation seed 7 that aborts on mismatch; the predecessor measured that the
  zero-shot path reproduces exactly, so the check should pass, and if it does not the reuse is
  wrong and must not proceed.
- **A third variant flag triples the exclusion surface.** `--diagnostic`, `--paired-sweep` and
  `--training-variance` must be pairwise exclusive; a two-flag check silently admits the third
  combination. Tested pairwise rather than by one `if a and b` line.
- **`uv sync --extra forward` deletes the Kaggle SDK** (declared in `[dev]`). Restore additively
  with `uv pip install kaggle`; a bare `uv sync` would also delete torch.
- **`codex` is absent**, so the pre-exec review is single-voice. For a plan whose core risk is a
  variance decomposition that can be made to say what one wants, that is a real weakness.
- **Pre-existing, non-blocking:** 12 suite failures in four absent-data families, none touched here.
- **Pre-existing and plan-relevant:** `wiki/glm_generator_paired_sweep_2026-11-20.json` carries
  `"verdict": "INSTRUMENT_INVALID"` in its own body, which is a known gate defect, not a statement
  about its data. Anyone joining against that grid must read the corrected verdict artifact instead;
  Step 6 joins against the grid's zero-shot CELLS, not its verdict field.

## Open Questions

- Is the surviving 0.1186 gap to the Markov null the better target? This plan refines the
  understanding of an effect that already failed to resolve, while the robust finding — the
  generator loses to a trivial 3-mer chain at 10 of 10 seeds — is untouched and is arguably the real
  problem. The case for this plan is that it can CLOSE a question; the case against is that closing
  it changes nothing about the 0.1186 gap. Worth a redirect decision before Step 5 spends the session.
- Should `n` rise from 120 to the CLI default 400 in the same run? It attacks measurement variance at
  source rather than characterising it, but it breaks comparability with the reused zero-shot row and
  with both committed arms, which would forfeit the anchors. Named rather than silently scoped out.
- Should epochs 1 and 2 be scored per training seed? It triples the scorings (~3 h) and the
  trajectory is already `TRAJECTORY_UNRESOLVED`, so this plan covers epoch 3 only. If
  `TRAINING_VARIANCE_DOMINATES` fires, the trajectory question is closed along with the delta and the
  extra cost would have been wasted.
- Do the three variance components want a proper nested model rather than two medians of spreads?
  A random-effects decomposition would give variance components with intervals, but it imports a
  modelling choice and an estimator this repo has no precedent for, against a 5x5 grid that cannot
  support it. Medians of spreads are coarse and defensible; naming the alternative so it is not
  rediscovered as novel.

## Verification

1. `uv run pytest tests/test_glm_finetune_kernel_builder.py tests/test_glm_training_variance_bar.py -q`
   green, each new load-bearing test confirmed to fail when its subject is reverted.
2. `git diff` shows `verdict()`, `EPOCHS`, `LR`, `BATCH`, `DTYPE`, `LORA` and `RECORDED` unchanged.
3. `scripts/kaggle/glm_finetune_kernel.py` sha256 still `6bdbd02641abedc5`.
4. Frozen five byte-unchanged; `prospective_lock.verify_lock` returns `ok=True`.
5. Full suite reconciled by CAUSE against the 12-failure baseline, not by count.
6. The pulled grid committed verbatim before any interpretation is written.
7. Every verdict in the memo recomputed by test from the frozen rule and the committed grid.
8. Both variance components reported separately, never pooled into one number.

## Save-time amendments

Captured at: 2026-11-21
Source: `/save-plan` arguments

**Audit-notes-only contract.** This block is provenance for human readers. `/execute-plan` reads
ONLY `## Implementation Steps` for executable work; nothing here is an instruction.

- training-variance decomposition; zero-shot row reused since it is training-independent; redirect question open on whether the 0.1186 null gap is the better target
- The three clauses above are each grounded in the plan body rather than being a gloss on it. The
  decomposition is Step 1's `training_spread(g)` / `measurement_spread(t)` pair. The zero-shot reuse
  rests on a verified fact: line 447 of the emitted kernel scores the pristine `MODEL_ID`, so that
  row cannot depend on training, and the predecessor measured it reproducing 0.7787 exactly -- the
  plan still re-scores one cell in-run as a reuse check that aborts on mismatch. The redirect
  question is the first Open Question, and it is NOT rhetorical: this plan can close a question but
  cannot move the 0.1186 gap, which is the robust finding.
- Unlike the predecessor plan's save, this amendment text was written AFTER the pre-exec review of
  the conversation it summarises, so no clause of it is superseded at save time.

<!-- toolkit: check=clean waves=clean gate=fired:open-questions,reuse-missing -->
