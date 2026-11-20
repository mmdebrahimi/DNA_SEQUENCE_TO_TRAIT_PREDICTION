# Measure The Generator Noise Floor As A Paired Delta, In One Session

## Lens status
Repo-index: absent — `repo-index-status.sh` returned nothing, so research was direct file reading with no Layer 2 degradation.
Tools: `sentrux` present (structural gate eligible), `gh` present, `codex` ABSENT — the pre-exec review of this plan was single-voice, which for a plan whose core risk is self-deception about a statistic is a named weakness, not a clean bill of health.
Feature record: `features/glm-generator-finetune/` exists and governs this family; this plan is conversation-driven and additive to it.
Baseline tests: 12 failed / 9012 passed / 191 skipped, all 12 in four documented absent-data families.
Revision: this is v2. v1 proposed a two-phase gated design; its own pre-exec review refuted three of its load-bearing choices and v1 was never saved. The refutations are recorded under Risk Flags as near-misses.

## Problem Statement

The GLM generator fine-tune arm is a recorded negative. Against a bar frozen before any fine-tuned
number existed, mode-both distinguishability had to fall strictly below the 3-mer Markov null of
**0.6679**. Zero-shot `GENERator-v2-prokaryote-1.2b-base` scored **0.7787**; the LoRA fine-tuned
merge scored **0.8058**. The conclusion carried since is that fine-tuning moved the metric the
WRONG WAY, by **+0.0271**.

**That delta is one draw versus one draw, and its noise has never been measured.** Generation is
seeded — `HFGenerator.__init__(..., seed: int = 0)`, and the falsifier CLI exposes `--seed` with
default **7**, which `run_arm` never passes. So both committed arms are a single sample at seed 7.
A one-sample-versus-one-sample comparison has no error bar.

The existing derived control does not fill the gap. `falsifier.self_control_ceiling` splits the
natural corpus against itself over seeds and is **identical across both arms** (0.5865 mode-both,
0.5447 kmer) because it depends only on the natural corpus and the discriminator seeds. It bounds
*"is this distinguishable at all"* — both arms far exceed it — and says nothing about whether two
different generated sets, drawn from two different weight sets, differ by more than sampling noise.

**The decision quantity is a DIFFERENCE, so the band must come from the variance of that
difference.** `Var(A - B) = Var(A) + Var(B) - 2Cov(A, B)`; one arm's seed spread is not the
uncertainty of `A - B`. The design therefore scores every model at the **same seed set** and judges
the distribution of per-seed paired deltas. Pairing also cancels the shared term: the eval windows,
the discriminator split and the natural corpus are identical within a seed, so a paired delta
isolates the weights.

Two questions, answered by one run rather than two:

1. **Is the +0.0271 delta larger than its own sampling variation?** Measured as the spread and sign
   consistency of `delta_s = score_ft(s) - score_zs(s)` across the seed set.
2. **Does distinguishability rise monotonically across fine-tuning epochs?** The diagnostic that
   separates "the corpus is too small" from "the adaptation is too aggressive". Scored at the same
   shared seeds, so every epoch comparison is paired too.

Non-goals: changing the frozen bar (0.6679 stands); re-running or re-scoring the two committed
arms; presenting any epoch variant as the pre-registered arm; publishing a corrected delta before
the paired measurement exists.

## Codebase Context


See: ./codebase-context.md

## Pre-Change Baseline

Measured from the committed artifacts and the 2026-10-31 run log, not from memory:

- `wiki/glm_generator_finetune_armA_zeroshot_2026-10-31.json` — mode-both **0.7787**, kmer **0.7355**
- `wiki/glm_generator_finetune_armB_finetuned_2026-10-31.json` — mode-both **0.8058**, kmer **0.7559**
- Markov null, both arms: mode-both **0.6679** (the bar), kmer **0.5990**
- Uniform control, both arms: mode-both **0.9641**, kmer **0.9408**
- `control_self_natural_vs_natural_max`: mode-both **0.5865** (0.5348 / 0.5455 / 0.4548 / 0.5053 /
  0.5865), kmer **0.5447** — identical across arms, by construction
- Geometry: `n=120` generated per candidate, `limit_natural=1200`, `length=150`, `context=300`,
  `bins=6`, conditional, fp32, Tesla T4, ONE seed (CLI default 7), `knob_wired: false`
- **Measured cost, from the run log of `emanueleebrahimi/glm-generator-finetune` (status COMPLETE):**
  imports 41 s, split 3 s, model download 27 s, zero-shot scoring **138 s**, 3-epoch LoRA **578 s**
  (3.2 min per epoch), fine-tuned scoring **135 s**; **whole run 944 s**.
- **The 2026-10-31 LoRA adapter still exists** in that run's output
  (`lora_adapter/adapter_model.safetensors`), because `ADAPTER` is under `/kaggle/working`. The
  merged model does not, because it went to `/tmp`.
- Paired generation-sampling variance of the delta: **NEVER MEASURED.** That absence is this plan's
  subject.
- Test suite: 12 failed / 9012 passed / 191 skipped; the 12 are cyp4f2/GeT-RM 7,
  essentiality-ladder 3, README-quickstart 1, pneumo cps DB 1.

## Verification Signal

- The emitted kernel AST-parses (the builder refuses otherwise), and the frozen two-arm kernel body
  is byte-identical when built with no new flag — compared in memory, never by writing.
- `verdict()`, `EPOCHS`, `LR`, `BATCH`, `DTYPE`, `LORA` and `RECORDED` unchanged by diff.
- The run returns a mode-both score for each of 4 models (zero-shot, epoch 1, epoch 2, epoch 3) at
  each seed in the declared set, with no missing cells.
- The seed set INCLUDES 7, so the zero-shot score at seed 7 is comparable to the committed 0.7787
  and the epoch-3 score at seed 7 to the committed 0.8058.
- Per-seed paired deltas are computed and reported individually, not only as a summary.
- The artifact states whether the committed +0.0271 survives its own paired variation, whichever way
  it falls.
- Frozen five (`amr_rules.py`, `calibrated_amr_rules.json`, `mic_tiers.py`,
  `shipped_decoder_surface.py`, `cohort_manifest.py`) byte-unchanged and
  `prospective_lock.verify_lock` still `ok=True`.
- No regression against the 12-failure baseline, reconciled by CAUSE and not by count.

## Implementation Steps

### Step 1: Pre-register the paired reading before any number exists
Files: wiki/glm_generator_paired_noise_bar_2026-11-20.json
Depends on: none

**What changes:**
- A frozen-bar artifact naming, BEFORE the run, every outcome — and keyed on the PAIRED delta, not
  on any single arm's spread.
- The decision quantity is declared explicitly as `delta_s = score_ft_ep3(s) - score_zs(s)`
  evaluated per seed `s`, with the summary statistics being the sign pattern, the mean and the
  spread of those deltas.
- Delta outcomes: `DELTA_ROBUST` when every seed's delta has the same sign as the committed +0.0271
  AND the mean delta exceeds the deltas' own spread (the carried "fine-tuning made it worse" reading
  survives); `DELTA_NOT_ESTABLISHED` when the signs are mixed OR the mean delta is within the
  deltas' own spread (the carried claim is NOT established at one draw per arm — a CORRECTION);
  `INSTRUMENT_INVALID` when the Markov null or the self-control does not reproduce, or any
  model-by-seed cell is missing.
- The rule is deliberately SCALE-FREE: it compares the delta to its own measured variation rather
  than to a guessed threshold. A numeric band is NOT frozen here, because its magnitude is exactly
  what the run measures, and freezing a scale before measuring it is how a bar gets mis-specified.
- Epoch-trajectory outcomes, judged on paired per-seed differences between consecutive epochs:
  `OVERFIT_MONOTONE` (each step up, consistent in sign across seeds — gentler adaptation is the
  lever); `FLAT_NOT_EPOCHS` (no consecutive difference is sign-consistent — epoch count is not the
  cause and the adaptation direction itself is wrong); `HELPS_THEN_DEGRADES` (epoch 1 below
  zero-shot sign-consistently, later epochs above); `TRAJECTORY_UNRESOLVED` (sign-consistent
  nowhere and the deltas' spread exceeds the total epoch-1-to-epoch-3 range).
- `claim_status: "DIAGNOSTIC_METHOD_VARIANT"` plus a `not_the_frozen_arm` note: epochs 1 and 2 are
  NOT the pre-registered method (`method_frozen.epochs = 3`), the bar 0.6679 is NOT restated as
  clearable by them, and no outcome here may be reported as the frozen arm's.
- An anchor rule: if the seed-7 zero-shot score differs from the committed 0.7787, or the seed-7
  epoch-3 score from 0.8058, by more than the deltas' spread, the run is labelled
  `ANCHOR_DRIFT` and no mechanism claim is made.
- A forecast with reasoning, so no result can be retro-fitted.

**Test strategy:**
- A test asserting the artifact parses as JSON and carries all three delta outcomes, all four
  trajectory outcomes, the anchor rule, `claim_status`, and a `frozen_before` timestamp.
- A test asserting the file declares the decision quantity as a per-seed DIFFERENCE and does NOT
  define any band from a single arm's spread.
- A test asserting the file does NOT restate 0.6679 as a target for the variant arms.

### Step 2: Thread an explicit seed through run_arm
Files: scripts/kaggle/build_glm_finetune_kernel.py
Depends on: none

**What changes:**
- `run_arm` gains a keyword-only `seed: int | None = None`; when not None it appends `--seed <n>` to
  the argv. Default None reproduces the current argv EXACTLY, so the frozen two-arm path is
  unchanged byte-for-byte.
- The seed set is a module constant in the builder, not a CLI argument, so the swept seeds are part
  of the committed artifact rather than an invocation detail. It MUST contain 7.

**Test strategy:**
- Assert the no-flag emitted BODY is byte-identical to the committed frozen kernel artifact,
  comparing the in-memory body against the file's contents and writing nothing. Invoking the
  builder's file-writing entry point is forbidden here: its default output path IS that artifact,
  so a build-then-compare test would overwrite what it is checking.
- Assert `--seed` appears in the argv only when a seed is supplied, and never otherwise.
- Assert the declared seed set contains 7 and has no duplicates.

### Step 3: Per-epoch adapter saves and the paired-sweep variant
Files: scripts/kaggle/build_glm_finetune_kernel.py
Depends on: Step 2

**What changes:**
- Inside the existing `for ep in range(EPOCHS)` loop in `step6_finetune`, after the epoch's loss is
  recorded, save the adapter to `ADAPTER + "/ep%d" % (ep + 1)`. Adapters are small and `ADAPTER` is
  under the collected output path, so all three survive the run.
- A new builder flag `--paired-sweep` emitting `scripts/kaggle/glm_paired_sweep_kernel.py`: ONE
  kernel that trains ONCE under the frozen method, then scores FOUR models — zero-shot plus the
  three epoch adapters, each merged in turn — at EVERY seed in the declared set, via `run_arm` with
  an explicit seed.
- Scoring order is model-major within seed-minor, and every model-by-seed cell is written as its own
  record carrying `model`, `epoch` (null for zero-shot) and `seed`, so a truncated session leaves a
  partially-filled grid that is still interpretable rather than an ambiguous average.
- Each epoch merge is released before the next is built, so peak disk holds one 4.8 GB merge.
- The variant calls its OWN reading function and never `verdict()`. Output
  `glm_paired_sweep_result.json` carries `claim_status: "DIAGNOSTIC_METHOD_VARIANT"`, the full grid,
  and the per-seed paired deltas.
- No gate and no second phase: the measured cost is 10 min training plus about 2.3 min per scoring,
  so the whole grid is one session.

**Test strategy:**
- Build with `--paired-sweep` and assert the body AST-parses, saves one adapter per epoch, and
  contains one scoring call per model-by-seed cell.
- Assert `EPOCHS`, `LR`, `BATCH`, `DTYPE` and the `LORA` dict are unchanged in the emitted body.
- Assert `verdict(` appears nowhere in the variant's reading path.
- Assert the no-flag emitted body is still byte-identical to the committed frozen kernel artifact
  after this step, again comparing in memory and writing nothing.

### Step 4: Builder tests, the no-write invariant, and non-vacuity
Files: tests/test_glm_finetune_kernel_builder.py
Depends on: Step 3

**What changes:**
- Add the Step 2 and Step 3 test strategies as real tests in the existing builder test module.
- Add a NON-VACUITY check: deliberately corrupt the emitted body in-test and assert the builder's
  AST guard refuses it, so the guard is proven to bind rather than assumed.
- Add a frozen-constants test that reads the builder source and asserts the five method constants
  and `RECORDED` are byte-unchanged against their committed values.
- Add a NO-WRITE invariant: the module must never call the builder's file-writing entry point,
  because that entry point's default output path IS the committed frozen kernel artifact. Assert it
  by snapshotting that artifact's sha256 before and after the module's tests run.

**Test strategy:**
- `uv run pytest tests/test_glm_finetune_kernel_builder.py -q` green.
- Confirm each new test FAILS when its subject is reverted, so none is vacuous.

### Step 5: Run the paired sweep and commit the raw grid verbatim
Files: scripts/kaggle/glm_paired_sweep_kernel.py, wiki/glm_generator_paired_sweep_2026-11-20.json
Depends on: Step 1, Step 4

**What changes:**
- Build the sweep kernel, push to the NEW slug `glm-paired-sweep-v1` with `--gpu`, poll to
  completion, pull the output.
- NEVER re-push the forbidden slugs: `pa-finetune-arm-v1` / `v2` / `v3`,
  `pa-finetune-ckpt-probe` / `probe2` / `probe3`, `pa-finetune-chain-s1` / `s2` / `s3` / `s4`. Also
  never re-push `glm-generator-finetune`, which holds the two committed arms' provenance.
- Commit the pulled grid verbatim BEFORE any interpretation is written.
- If the kernel errors with an empty failureMessage and zero output files — the documented
  2026-10-31 signature — record that and STOP rather than retrying blind, because that signature is
  indistinguishable from infrastructure failure.
- If the session is cut mid-grid, the per-cell records and the saved epoch adapters both survive in
  the collected output, so scoring can resume without retraining.

**Test strategy:**
- Manual verification that the pulled JSON carries one record per model-by-seed cell with no gaps,
  and that the seed-7 cells for zero-shot and epoch 3 are present.
- Assert the Markov null reproduces 0.6679 and the self-control reproduces 0.5865; if either
  drifts, the run is `INSTRUMENT_INVALID` by Step 1's frozen reading.

### Step 6: Score the paired deltas against the frozen rule and write the verdict
Files: wiki/glm_generator_paired_sweep_2026-11-20.md
Depends on: Step 5

**What changes:**
- Compute per-seed `delta_s`, the sign pattern, the mean and the spread, and apply Step 1's frozen
  rule with no post-hoc threshold.
- Apply the anchor rule against the committed 0.7787 and 0.8058 before any mechanism claim.
- On `DELTA_NOT_ESTABLISHED`: record it as a CORRECTION to a carried claim, state plainly that
  "fine-tuning made it worse" is not established at one draw per arm, and do NOT edit the two
  committed arm artifacts — add a namespace-separate note and leave their measured numbers
  untouched.
- On `DELTA_ROBUST`: record that the carried reading survives, with the measured variation beside it
  so the delta is never again quoted without one.
- Report the epoch trajectory verdict separately from the delta verdict, so a conclusion about
  epochs cannot be read as a conclusion about the bar.
- Record honest limits: one organism, one corpus, ONE training trajectory (seeds vary generation,
  not training, so this bounds generation noise and NOT training-run-to-run noise), upstream windows
  as a promoter PROXY, and a k-mer discriminator blind to long-range organisation so a low score
  remains necessary-not-sufficient.

**Test strategy:**
- A test asserting the memo's stated verdicts are the ones Step 1's rules yield from the recorded
  grid, so the prose cannot drift from the frozen reading.
- A test asserting the memo carries the training-noise limitation, since seed pairing does not
  address it.

## Execution Preview

5 waves, max parallelism 2, critical path 5 steps (Step 2 -> Step 6).

- Wave 0: Step 1 (pre-registration), Step 2 (seed passthrough) — disjoint files
- Wave 1: Step 3 (per-epoch saves + paired-sweep variant) — same builder file as Step 2, sequenced
- Wave 2: Step 4 (builder tests)
- Wave 3: Step 5 (the run)
- Wave 4: Step 6 (verdict)

Wave 3 is a single attended Kaggle session of roughly one hour at ten seeds, derived from the
measured 944 s run: about 10 min training plus about 2.3 min per model-by-seed scoring.

## Risk Flags

- **Near-miss, recorded: v1 of this plan derived the band from one arm's seed spread.** That is the
  documented anti-pattern in `feedback_tolerance_band_from_the_tested_quantity`, whose first
  recorded instance (2026-10-20) was the same generator and the same shape, also caught on review
  before the run. v2's paired design is the fix. Noting it here because the error is silent and will
  look reasonable again next time.
- **Near-miss, recorded: v1 gated an "expensive" second phase that does not exist.** The whole
  committed run was 944 s. A cost premise asserted from the phrase "well under one Kaggle session"
  rather than from the run log produced three unnecessary steps and a gate that would have
  suppressed the trajectory measurement entirely.
- **Seed pairing does NOT bound training variance.** All four scored models come from ONE training
  trajectory, so a `DELTA_ROBUST` verdict establishes that the delta survives GENERATION noise only.
  A second training run at a different torch seed is the test for training noise and is deliberately
  out of scope; the memo must say so rather than implying the delta is fully bounded.
- **Post-hoc deviation from a frozen pre-registration.** Epochs 1 and 2 are not
  `method_frozen.epochs = 3`. Mitigated by a separate bar frozen before any number, a distinct
  `claim_status`, a separate reading function, and never restating 0.6679 as clearable by the
  variant. This remains the plan's largest honesty risk and is why Step 1 precedes everything.
- **The builder's own documented failure mode.** The builder can parse while the emitted kernel does
  not; on 2026-10-31 that cost two sessions and surfaced as ERROR with an empty failureMessage and
  zero outputs. Mitigated by the existing AST guard plus the Step 4 non-vacuity test that proves the
  guard binds. Checking the ARTIFACT not its generator is the recorded lesson.
- **The builder's write path targets the frozen kernel artifact.** Any test that builds with no flag
  and then compares would overwrite its own subject. Addressed as an explicit no-write invariant
  with a sha256 snapshot in Step 4.
- **`uv sync --extra forward` will delete the Kaggle SDK**, which is declared in `[dev]`. Restore
  additively with `uv pip install kaggle`; a bare `uv sync` would also delete torch.
- **CUDA nondeterminism.** `RECORDED` warns that CUDA yields DIFFERENT sampled sequences than CPU
  for this generator. This run is T4-to-T4 against T4 committed numbers, but cuDNN nondeterminism at
  fixed seed is not excluded, which is why the anchors are checked against the measured spread
  rather than for exact equality.
- **`codex` is absent**, so the pre-exec review of this plan was single-voice. For a plan whose core
  risk is a subtly wrong statistic, that is a real weakness; it caught three criticals, but a second
  voice would have been worth having.
- **Pre-existing, non-blocking:** 12 suite failures in four absent-data families. None is touched by
  this plan; reconcile by CAUSE at the end, since a matching count is not a matching set.

## Open Questions

- How many seeds? Five mirrors `self_control_ceiling`'s own `seeds=(0,1,2,3,4)` convention and keeps
  the artifact comparable; ten roughly doubles a ~46-minute scoring budget and materially tightens
  the paired deltas' spread. Recommend ten, since the cost constraint that argued for five was
  measured away. Must include 7 either way.
- Should `--n` rise from 120 to the CLI default 400? It reduces sampling variance at the source
  rather than characterizing it, but it breaks comparability with both committed arms and so
  forfeits the anchors. Running the paired grid at 120 for comparability AND one extra arm at 400 as
  a small-sample probe is affordable; whether that is worth the extra scorings is a scope call.
- If the verdict is `DELTA_NOT_ESTABLISHED`, does the carried claim need correcting beyond this
  plan's memo — in `features/glm-generator-finetune/`, in `CLAUDE.md`, or in the project memory that
  now records 0.7787 -> 0.8058? That is correction-propagation and is the user's call, not mine to
  take silently.
- Is a second training run at a different torch seed worth it, to bound the training variance this
  design explicitly does not? It doubles the training cost (10 min) and is the only way to make a
  delta claim fully bounded. Named rather than silently scoped out.

## Verification

1. `uv run pytest tests/test_glm_finetune_kernel_builder.py -q` green, and each new test confirmed
   to fail when its subject is reverted.
2. `git diff` shows `verdict()`, `EPOCHS`, `LR`, `BATCH`, `DTYPE`, `LORA` and `RECORDED` unchanged.
3. No-flag emitted body byte-identical to the committed `scripts/kaggle/glm_finetune_kernel.py`,
   verified in memory with that artifact's sha256 unchanged on disk.
4. Frozen five byte-unchanged; `prospective_lock.verify_lock` returns `ok=True`.
5. Full suite reconciled by CAUSE against the 12-failure baseline, not by count.
6. The pulled grid committed verbatim before any interpretation is written.
7. Every verdict in the memo recomputed by test from the frozen rule and the recorded grid.

## Save-time amendments

Captured at: 2026-11-20
Source: `/save-plan` arguments

**Audit-notes-only contract.** This block is provenance for human readers. `/execute-plan` reads
ONLY `## Implementation Steps` for executable work; nothing here is an instruction.

- noise floor measured first; epoch trajectory gated on it because both committed arms are one seed-7 draw each
- CORRECTION to the line above, from the pre-exec review: the first and last clauses stand and are
  load-bearing — the noise floor IS measured first, and both committed arms ARE a single seed-7 draw
  each. The middle clause does NOT: the epoch trajectory is NOT gated. The gate was removed because
  the cost premise behind it was false (the whole 2026-10-31 run was 944 s, so there is no expensive
  second phase to gate) and because gating the trajectory on an unresolved delta was backwards — a
  noisy instrument argues for more draws, not fewer measurements.
- The amendment text supplied was written BEFORE that review and describes the superseded v1 design.
  It is preserved verbatim above rather than silently rewritten, and v1 was never saved, so no
  superseded artifact entered `plans-index.md`.
- The band is now derived from the PAIRED per-seed delta, not from one arm's seed spread. v1's
  version of that choice is the documented anti-pattern in
  `feedback_tolerance_band_from_the_tested_quantity`, whose first recorded instance was this same
  generator.

<!-- toolkit: check=clean waves=clean gate=fired:open-questions,reuse-missing -->
