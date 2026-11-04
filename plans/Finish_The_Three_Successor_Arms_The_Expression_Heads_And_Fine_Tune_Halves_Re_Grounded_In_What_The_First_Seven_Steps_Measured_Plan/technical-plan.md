# Finish The Three Successor Arms: The Expression-Heads And Fine-Tune Halves, Re-Grounded In What The First Seven Steps Measured

## Lens status
Dispatch: conversation-driven (no argument). Subject: the seven unexecuted steps of the in-flight plan `plans/Declare_And_Run_All_Three_Successor_Arms_And_Convert_The_Guards_From_State_Pins_To_Invariants_Plan/technical-plan.md`, redesigned against live measurement.
Repo-index: absent — `repo-index-status.sh` returned nothing, so this plan was built by direct file reading, with no Layer 2 summaries consulted.
Reuse signals: no `graphify-out/GRAPH_REPORT.md`, and none of the conventional shared-code directories exist. The survey below is grounded in the sibling arm scripts the first seven steps actually built.
Tools: `sentrux` present, `gh` present, `codex` absent.
Project rules: `project-rules.md` absent, so the `## Project Rules Check` section is omitted by contract.
Model-touching facts in this plan were MEASURED during planning, not inherited.

## DECISIONS NEEDED

**One decision is genuinely yours, and it changes a number inside an acceptance bar before any score exists.**

**D1 — Four expression heads or five?** The plan being superseded says four, and sets the multiplicity correction (the statistical penalty for testing several things at once) to four. Live measurement says there are **five** heads, and all five are identical in the ways that matter: every one has an `fc1.weight` of shape `(64, 11008)`, and every one carries base-encoder weights byte-identical to the frozen base at 124 of 124 tensors, maximum absolute deviation exactly 0.000e+00. The fifth is `promoteratlas-tslpred-kosuri2013.pt`, and the argument for leaving it out is not its shape but its **endpoint**: the other four are `trspred`, transcription prediction, whereas `tslpred` is translation prediction, and it uniquely carries `log_transform: True` in its normalisation parameters. A rung that measures promoter activity is asking a transcription question, so a translation head is arguably predicting a different quantity.

**This plan proceeds on your standing instruction** — "when given a choice we always want to do all the more elaborate version and make sure, and not leave any stones unturned" — so it **declares all five heads and corrects at α/5**, recording the endpoint mismatch as a stated limit on the `tslpred` head rather than as grounds for dropping it. The honest cost of that choice, stated in advance: a five-way correction gives every head a wider confidence interval than a four-way one, so this makes the arm **harder** to pass, not easier. If you would rather scope the arm to the four transcription heads, say so before Step 4 freezes the bar; afterwards it is a frozen-bar edit and a different kind of decision.

**D2 — the `met` verdict will resurface at Step 8, and that is the same authority call you already answered once.** Recording all three arms removes the last `planned` arm again, so `ceiling_verdict` flips from `blocked:untried-arm` back to `met`. Step 8 reports it and stops rather than recording it in passing. No action is needed from you now; it is named here so its reappearance is expected rather than surprising.

## CORRECTIONS

**Two things the superseded plan asserts about the expression-heads arm are refuted by measurement. Both were inherited from research notes, and both would otherwise have been discovered mid-implementation.**

1. **`ExpressionPredictor.from_pretrained` cannot be used at all, so the plan's "`from_pretrained(cls, model_path)` is the loader" is wrong.** Measured: calling `ExpressionPredictor()` under the installed torch 2.7.1+cu118 raises `UnpicklingError` before doing anything else. Its `__init__` runs `torch.load(base_model_path)` with no `weights_only` argument, and torch 2.6 changed that default to `True`; the base checkpoint embeds a `pathlib.PosixPath` because it was saved on Linux, which the restricted unpickler refuses. This is independent of a second, separate defect in the same two lines — the path `trained_weights/base_model/promoteratlas-base.pt` is **relative**, so the class could only ever work with the current directory set to the clone root. The arm must therefore construct the head module itself and load state dictionaries directly. Step 1 exists to own that.

2. **There are five heads, not four.** See D1. The consequence is arithmetic: the multiplicity divisor in the bar changes.

**A third claim was checked and survives, and the check is recorded because the near-miss is instructive.** A state-dictionary element count for the base checkpoint gives **1,847,492**, not the 1,847,428 that `dna_decode/glm/atlas_embed.py:32` pins as `EXPECTED_PARAMS` and enforces as a hard refusal. The repo is right and the discrepancy is not an error: the gap is exactly 64 elements, which is the eight `blocks.N.rotary_attention_block.rotary_emb.inv_freq` tensors at 8 elements each. Those are registered **buffers**, not parameters, and the gate counts `sum(p.numel() for p in model.parameters())`, which excludes buffers. The trap is real and sits directly in the path of Steps 2 and 3, both of which re-verify the checkpoint: a probe that sums the state dictionary gets 1,847,492, trips the refusal, and reports a checkpoint mismatch that does not exist.

## Problem Statement

The in-flight plan declared three successor arms to the ceiling criterion and has completed one of them end to end. The as-given arm is scored (`ASGIVEN_TIES`: median 0.5086 against the incumbent 0.5076, gain −0.0042, 95% confidence interval [−0.0068, +0.0004] containing zero, 4 of 10 seed wins), the eight guards that reddened when the verdict flipped are converted to invariants and green, and the scoring runner is decomposed so a sibling arm can share its gate chain instead of copying it.

What remains is the expensive half: the two arms that need the model itself rather than the committed feature block. The **expression-heads arm** asks whether the five trained prediction heads shipped beside the frozen encoder transfer to this rung at all; it is not a drop-in, because every head is geometry-locked to 86 base pairs while the substrate is 150. The **fine-tune arm** asks whether the frozen representation was the limit; it is the only arm whose cost must be measured before it is allowed to start, because a budget-limited null and a representation-limited null look identical in the artifact and mean opposite things.

The goal is a verdict on each from a bar frozen before any score exists, reached through the same gate chain that produced the first two arms, with the three arms then recorded together.

## Codebase Context

**The substrate and comparator are fixed and must not be recomputed.** The cohort is 44,106 peak tiles at 150 base pairs over 2,433 peaks, selected by `category == "tile" and len(seq) == 150`; `scripts/positional_arm_probe.py:load_arm_cohort()` refuses any other count. The comparator is `data/promoteratlas/local_reconcile_2026-11-01.json` (`verdict: RECONCILED`, `mismatched: []`), whose atlas median is 0.5076 with committed per-seed values and `n_test_per_seed` for seeds 0 through 9 at `frac = 0.25`. Both new bars inherit the margin 0.02 and these per-seed numbers verbatim; neither recomputes them.

**Head facts, all measured during planning** by loading each checkpoint behind a static pickle audit (`dna_decode.glm.atlas_embed.audit_checkpoint`), and by instantiating the clone's class to watch it fail:

| head | `fc1.weight` | base tensors identical | normalisation mean / std |
|---|---|---|---|
| `trspred-hossain2020` | (64, 11008) | 124/124, maxdev 0.000e+00 | −6.497496 / 1.797108 |
| `trspred-lafleur2022` | (64, 11008) | 124/124, maxdev 0.000e+00 | −2.618656 / 0.757282 |
| `trspred-urtecho2018` | (64, 11008) | 124/124, maxdev 0.000e+00 | −1.911309 / 1.378551 |
| `trspred-yu2021` | (64, 11008) | 124/124, maxdev 0.000e+00 | −2.577607 / 1.065576 |
| `tslpred-kosuri2013` | (64, 11008) | 124/124, maxdev 0.000e+00 | 8.036710 / 1.180190, `log_transform: True` |

Each head checkpoint is 11.31 MB, totals 2,814,533 state-dictionary elements of which 967,041 are head-only, and carries exactly `conv1`, `conv2`, `fc1`, `fc2` weights and biases above the base. `11008 = 128 × 86`, so the heads expect an 86-position latent while a 150-position one flattens to 128 × 150 = 19,200; the window count is therefore `150 − 86 + 1 = 65`. The architecture is `conv1 → relu → conv2 → relu → pool → flatten → fc1 → relu → fc2` over the latent, with `conv1`/`conv2` at `kernel_size=8, padding='same'` and `pool = MaxPool1d(7, stride=1, padding=3)` — and `padding='same'` is exactly why the cheap reduction and a true per-window recompute cannot be assumed equal at the edges.

**Loading these checkpoints safely is already solved in-repo, and the obvious alternative is a dead end.** `atlas_embed.audit_checkpoint` disassembles every pickle with `pickletools.genops` **without executing it** and refuses any global outside `SAFE_GLOBAL_PREFIXES`; `load_atlas` then aliases `pathlib.PosixPath = pathlib.WindowsPath` on Windows and loads with `weights_only=False`, which is sound precisely because the static audit already ran. Measured: **all five head checkpoints PASS that audit unchanged**, each naming six globals (`_codecs.encode`, `collections.OrderedDict`, `numpy.core.multiarray.scalar`, `numpy.dtype`, and two torch entries), every one already covered. No allowlist change is needed. The alternative route — `weights_only=True` with `torch.serialization.safe_globals` — was attempted during planning and **cannot be made to work**: torch keys its allowlist on `f"{obj.__module__}.{obj.__name__}"`, the checkpoints were pickled under numpy 1.x as `numpy.core.multiarray.scalar`, and numpy 2.4.4 exposes that object at `numpy._core.multiarray`, so the name the unpickler asks for and the name the allowlist registers can never meet without a stand-in class that lies about its own module. Reuse the existing audit-then-alias pattern rather than rebuilding it.

**The base checkpoint carries its own pretraining recipe, which the fine-tune should anchor to rather than invent.** Its top-level keys are `config`, `epoch`, `model_state_dict`, `optimizer_state_dict`, `val_loss`, and `config` records `learning_rate: 0.001`, `patience: 25`, `batch_size: 1024`, `n_epochs: 2000`, `n_point_masks: 60`, `window_size: 6`, `gradient_accumulation: 1`, `architecture: DNATransformer`, `device: cuda`. Every head checkpoint also carries `val_loss`. Step 3 reads these rather than choosing hyperparameters freely.

**One hazard to carry forward:** `load_atlas` mutates `pathlib.PosixPath` process-globally and never restores it. That is pre-existing and works, but it means a process that loads the base and then the heads inherits the alias, so Step 1 must not depend on `PosixPath` being its real self.

### Reusable-Code Survey

Searched: `scripts/positional_arm_*.py`, `scripts/asgiven_arm_score.py`, `scripts/freeze_asgiven_arm_bar.py`, `scripts/freeze_positional_arm_bar.py`, `dna_decode/glm/atlas_embed.py`, `dna_decode/glm/positional_readout.py`, `scripts/kaggle/` (21 builders), `tests/test_*arm*.py`. No `graphify-out/GRAPH_REPORT.md` and no `src/utils`-style shared directory exists, so there were no injected god-node candidates to start from.

Reuse, rather than rebuild:

- **The gate chain.** `scripts/positional_arm_score.py` exports `GATE_FUNCS`, `gate_split_equality`, `gate_parity`, `gate_shuffled_label_null`, and a `ScoreCtx` dataclass whose `primary_rung` field exists precisely so a sibling arm can point the null at its own primary. `scripts/asgiven_arm_score.py` already consumes it this way and its 20 tests pass, so the pattern is proven rather than hoped for. Steps 6 and 7 import it; neither copies a gate.
- **The statistics.** `paired_bootstrap` (`N_BOOT = 10_000`, `BOOT_SEED = 707`), `order_statistic_median_ci`, `median_of`, `leave_peak_out_indices`, `build_rung`.
- **The bar skeleton.** `freeze_asgiven_arm_bar.py` supplies the shape every bar here follows: a closed `OUTCOMES` tuple, a `CHECK_ORDER`, `reachable_outcomes()` asserted **equal** to the closed set rather than contained in it, a total `classify()`, and a `verify()` that refuses to freeze against drifted anchors. `ceiling_arms_and_criterion(ledger)` is already public in `freeze_positional_arm_bar.py` for exactly this sharing.
- **The feasibility helpers.** `scripts/positional_arm_probe.py` supplies `fits_in_ram` (returns `None` for unknown and never `True`), `free_ram_bytes`, `design_matrix_bytes`, `extrapolate_wall_seconds`, `sha256_of`, `stride_sample`, `load_arm_cohort`.
- **The encoder.** `atlas_embed.load_atlas`, `audit_checkpoint`, `EXPECTED_PARAMS`, `NUC`, `ridge_spearman`, `masked_mean_pool`.
- **The Kaggle pattern.** `scripts/kaggle/build_promoteratlas_tile_arm_kernel.py` is the closest template: it already pins `machine_shape: NvidiaTeslaT4` and already carries the `EXPECTED_PARAMS` refusal and the `PosixPath` guard.

Deliberately **not** reused: `promoter_atlas.models.expression_predictor.ExpressionPredictor` — refuted above.

## Pre-Change Baseline

**Full suite, reused rather than re-run:** 21 failed, 7,514 passed, 180 skipped, measured at commit `b0521d90`, with the failing names recorded in `.claude/execute-plan-state/Declare_And_Run_All_Three_Successor_Arms_And_Convert_The_Guards_From_State_Pins_To_Invariants_Plan.json`. The plan-reference baseline is 13 failed / 7,256 passed at `33b24322`; the eight extra were the guards that reddened when `ceiling_verdict` flipped to `met`.

**Those eight are now cleared, and that is measured, not assumed.** At HEAD `3a85a638`, pytest over the eight arms-related files (`test_arms_ledger_stack_arm.py`, `test_glm_tile_arms_declared.py`, `test_tile_chain_substrate_wall_amendment5.py`, `test_tile_chain_constraint_registered.py`, `test_successor_arms_declared.py`, `test_arms_invariant_guards_non_vacuous.py`, `test_asgiven_arm_score.py`, `test_freeze_asgiven_arm_bar.py`) gives **135 passed in 3.72 s, zero failures**. The expected full-suite failure set at HEAD is therefore the 13-name reference set minus the one tile-chain constraint test that is now green, every remaining member being an absent-data or absent-tool failure (`cyp4f2` 7, essentiality-ladder 3, pneumo 1, README quickstart 1). **That derived figure of 12 must be re-measured at execution time and the 13 names used as the diff reference; do not treat the derived number as a measurement.**

**Arms state:** 16 arms — 12 measured, 3 planned, 1 killed. `ceiling_verdict = blocked:untried-arm`. `best` is `promoteratlas-base-embeddings-frozen-150bp-padded-to-200` at 0.5076, which is 54% of the 0.941 replicate-derived ceiling.

**Architecture gate:** not re-baselined, deliberately, so the final compare stays cumulative across both plans and can show whether the scorer decomposition brought `complex_fn_count` back down from 223.

**Frozen surface:** the five frozen files and `uv.lock` are byte-unchanged and must remain so; `dna_decode/glm/` is not part of that surface, so Step 1 may add a module there.

## Verification Signal

1. **Geometry is re-measured live, not inherited.** The probe asserts every declared head's `fc1.weight` is `(64, 128 × 86)` and exits non-zero otherwise, and asserts 124 of 124 base tensors identical at maxdev 0.000e+00 per head. A head whose base weights have drifted makes the arm a fine-tuned-encoder arm rather than a frozen-encoder one, which is a different claim, so the probe refuses rather than proceeding.
2. **The windowing choice is derived from a measurement, not preferred.** The probe computes the cheap full-latent reduction and a true per-window recompute on a sampled subset and reports maximum deviation, median deviation and the Spearman correlation between them. The bar then selects, and records the tolerance it selected under.
3. **Head output validity is gated before scoring.** Each head must produce finite, non-constant output across the cohort. A collapsed head would otherwise score as a clean null, which is the failure mode that looks most like a finding.
4. **Monotone invariance is asserted.** A head whose output is monotonically transformed must score an identical Spearman, which is the property that makes the denormalisation irrelevant to the verdict. This is a test, not a claim.
5. **The fine-tune budget is measured before the fine-tune is allowed to start.** Seconds per epoch and an observed epochs-to-early-stop curve on one seed, and a route selection exercised in both directions so the rule is not vacuous.
6. **A shuffled-label null runs on the trained encoder itself**, which is the arm where that control matters most: a network with enough capacity to memorise destroyed labels would otherwise look like signal.
7. **Every verdict is re-derivable from its own artifact** plus the frozen bar, with nothing else consulted, and the reachable outcome set is asserted **equal** to the closed set.
8. **The failure set is diffed by name** against the 13-name reference, not by count.
9. **The converted invariant guards stay green** with three more arms measured — the regression the conversion exists to prevent.
10. **The frozen five stay byte-identical** and the prospective lock re-verifies.

## Implementation Steps

### Step 1: A shared head-loading and windowing module
Files: dna_decode/glm/atlas_heads.py, tests/test_atlas_heads.py
Depends on: none

**What changes:**
- Adds `dna_decode/glm/atlas_heads.py` beside `atlas_embed.py`, owning every fact about the five heads so the probe and the scorer cannot drift apart. The sibling arms' defect history is the reason this is a module and not two copies.
- `HEADS` names the five checkpoints with their endpoint (`trspred` or `tslpred`) and the `log_transform` flag, derived from the checkpoint rather than restated.
- `load_head(path)` reuses `atlas_embed.audit_checkpoint` first and only then loads, mirroring `load_atlas` exactly. It returns the state dictionary and the `normalization_params`, and refuses a checkpoint whose `fc1.weight` is not `(64, 11008)` or whose `base_model.*` tensors are not byte-identical to the frozen base.
- **It must not assume `pathlib.PosixPath` is unaliased**, because `load_atlas` rebinds it process-globally and never restores it.
- `EXPECTED_HEAD_ELEMENTS = 2_814_533` and `EXPECTED_HEAD_ONLY = 967_041` are pinned, with a docstring stating plainly that these are **state-dictionary element counts** and that the base's comparable figure is 1,847,492 against `EXPECTED_PARAMS = 1,847,428`, the 64-element gap being eight rotary `inv_freq` buffers. Writing the trap down beside the constant is what stops the next reader re-deriving the wrong number.
- Pure reductions, each refusing non-finite input: `window_starts(length, width)` returning the 65 starts for 150 and 86; `head_forward_full(latent, w)` applying `conv1/relu/conv2/relu/pool` once over the full 150 latent and then slicing 86-wide windows before `fc1` (the cheap route); `head_forward_per_window(latent, w)` slicing the latent first and recomputing the convolutions per window (the true route); and `reduce_windows(values, how)` for `mean` and `max`.
- No file in the frozen five is touched, and `atlas_embed.py` is left unmodified — the measurement above shows its allowlist already suffices.

**Test strategy:**
- Unit-test `window_starts` against the hand-checked answer (65 windows, first start 0, last start 64) and against widths that do not fit, which must refuse rather than truncate.
- Unit-test both forward routes on a synthetic latent built so the correct answer is known in closed form, and assert they **differ at the edges and agree in the interior** — the `padding='same'` consequence. A test asserting they are equal everywhere would be asserting the bug.
- Assert `load_head` refuses a synthetic checkpoint whose base weights differ by one tensor, and one whose `fc1` shape is wrong.
- Assert the pinned element counts reconcile: head-only plus base state-dictionary elements equals the head total.
- Skip the clone-touching tests cleanly when the volatile temp directory does not resolve.

### Step 2: Expression-heads feasibility and circularity probe
Files: scripts/expression_heads_probe.py, wiki/expression_heads_feasibility_2026-11-02.json, tests/test_expression_heads_probe.py
Depends on: Step 1

**What changes:**
- Re-measures, live, every fact the arm rests on, and exits non-zero on any failure: the geometry lock per head, the base byte-identity per head, the base checkpoint's sha256 against `bb413bc8602f05b6…` and its 22,313,814 bytes, and the parameter count through `model.parameters()` rather than the state dictionary.
- Measures the cheap-versus-true windowing deviation on a strided sample of the cohort, reporting maximum absolute deviation, median absolute deviation and the Spearman correlation between the two routes **per head**, so the bar chooses from a number. Reports the wall-clock cost of each route so the choice has a price attached.
- Records each head's `normalization_params` and states that Spearman is invariant under the monotone denormalisation, so it cannot move the verdict — while an absolute predicted expression value would be wrong without it.
- Records the RAM pre-check for the 44,106-row by 65-window design through `fits_in_ram`, which returns `None` for unknown and never `True`.
- **Records the circularity scope, which is narrower than it looks.** The `urtecho2018` head is circular on the sigma70 rung (GSE108535) and must never be a comparator there. This substrate is GSE144621 peak tiles, and the project's own withdrawn lead established that no genome-wide or tiling Urtecho head exists, so it is not circular here. Recorded as a limit so nobody reuses the head on the wrong rung.
- **Records the `tslpred` endpoint mismatch** as a stated limit on that head: it predicts translation, carries `log_transform: True`, and is included under the standing instruction to leave no stones unturned rather than because its endpoint matches the rung.

**Test strategy:**
- Unit-test the pure arithmetic (window counts, deviation summaries, route costing) on synthetic latents with known answers.
- Assert the probe exits non-zero on a synthetic head whose base weights differ from the base, and on a moved checkpoint hash.
- Assert the artifact names all five heads and carries a deviation figure for each, so a head cannot be silently skipped.
- The model-touching path skips cleanly when the clone does not resolve.

### Step 3: Fine-tune feasibility and compute-route probe
Files: scripts/finetune_arm_probe.py, wiki/finetune_arm_feasibility_2026-11-02.json, tests/test_finetune_arm_probe.py
Depends on: none

**What changes:**
- Measures training throughput — forward **and** backward — on one seed for a bounded number of steps, and reports seconds per epoch plus the extrapolated cost of the bar's seed-by-epoch grid. The fine-tune is the long pole and its budget must be measured rather than estimated.
- Measures epochs-to-early-stop on one seed using a peak-blocked inner split, so the bar's epoch budget comes from an observed convergence curve rather than a guess. Without that curve this arm should not run at all, because it is the only thing separating a budget-limited null from a representation-limited one.
- **Anchors the hyperparameters to the checkpoint's own `config`** rather than choosing them freely: `learning_rate: 0.001`, `patience: 25`, `batch_size: 1024`, `n_epochs: 2000`, recorded verbatim from the base checkpoint, with any deviation named and justified in the artifact. A fine-tune whose recipe is unrelated to the pretraining recipe is a weaker test of the same question.
- **Chooses the compute route from the measurement.** Local is CPU-only (`torch.cuda.is_available()` is `False`, 4 threads, measured forward throughput 57.5 sequences per second). The project's established GPU route is a Kaggle T4 kernel and `scripts/kaggle/` holds 21 builders. The probe reports both costs and the bar freezes the route; a local run that would exceed a named wall-clock ceiling routes to Kaggle rather than being started and abandoned.
- Records that the encoder is 1,847,428 parameters by `model.parameters()`, which is why a fine-tune is tractable here and why this is **not** comparable to the project's blocked large-language-model fine-tunes.
- Re-verifies the base checkpoint's sha256, since the clone lives in a volatile temp directory.

**Test strategy:**
- Unit-test the pure extrapolation and route-selection arithmetic on synthetic timings.
- Assert the route selector picks Kaggle when the local extrapolation exceeds the ceiling and local when it does not — **both directions**, so the rule is not vacuous.
- Assert the probe exits non-zero when the checkpoint hash has moved.
- Assert the recorded hyperparameters match the checkpoint's `config` and that a deviation is flagged rather than silently accepted.

### Step 4: Freeze the expression-heads bar
Files: scripts/freeze_expression_heads_bar.py, wiki/expression_heads_acceptance_bar.json, tests/test_freeze_expression_heads_bar.py
Depends on: Step 2

**What changes:**
- Freezes the windowing decision from Step 2's measurement, before any head output is scored. If the cheap route tracks the true recompute within a stated tolerance it is used and the tolerance is recorded; otherwise the true recompute is used and its cost is recorded. The choice is derived, and the artifact says which measurement drove it.
- Declares **mean over the 65 windows** as the single primary per head, with max-over-windows a pre-declared secondary that is **ineligible** as the arm's result — the same shape that kept the positional arm honest when its highest rung was an ineligible secondary.
- **Declares all five heads and corrects at α/5** per decision D1, with the `tslpred` endpoint mismatch recorded as a limit. The arm's metric is the **maximum over the declared heads**, never over heads crossed with reductions.
- Inherits the margin 0.02 and the committed atlas comparator verbatim; `verify()` refuses to freeze against a drifted comparator, an unreconciled one, a missing arm, an arm already carrying a metric, or a mutated probe artifact.
- Pre-commits the zero-shot transfer framing: these heads were trained on other libraries, so a null bounds **transfer to this rung** rather than bounding the heads.
- Closed outcome set with `reachable_outcomes()` asserted equal to the check order, including `INDETERMINATE_GEOMETRY_FAIL` for a windowing that cannot be validated.
- Records its own forecast with probabilities summing to 1, and states that the forecast record in this project is 5 of 9, so the prior is not to be deferred to.

**Test strategy:**
- Assert `verify()` refuses on each mutated anchor and on a missing arm.
- Assert the reachable-outcome set **equals** the check order, and that padding the check order makes the equality fail — the non-vacuity proof.
- Assert the Bonferroni divisor equals the number of declared heads, so changing the head set cannot leave a stale divisor behind.
- Assert the secondary cannot enter the verdict function, by driving it with a secondary that would otherwise win.
- Assert `build_bar` is a pure function of the verified anchors.

### Step 5: Freeze the fine-tune bar
Files: scripts/freeze_finetune_arm_bar.py, wiki/finetune_arm_acceptance_bar.json, tests/test_freeze_finetune_arm_bar.py
Depends on: Step 3

**What changes:**
- Freezes the seed count, epoch budget, learning rate, early-stopping rule and compute route from Step 3's measurements, before any score exists.
- Inherits the margin 0.02 and the committed atlas comparator; the comparison is the fine-tuned encoder's median Spearman against the **frozen** encoder's 0.5076 on byte-identical splits.
- **Pre-commits the two readings that must not be conflated:** beating 0.5076 shows the frozen representation was the limit; failing to shows that the substrate or the training budget was, and a budget-limited null cannot be read as a representation-limited one. The bar ties that distinction to Step 3's measured epoch curve rather than leaving it to the reader.
- Declares a shuffled-label null on the fine-tune itself.
- Closed outcome set with reachable equal to the check order, including `INDETERMINATE_BUDGET_EXHAUSTED` for a run that stopped before its frozen epoch budget.

**Test strategy:**
- Assert `verify()` refuses on a mutated probe artifact, a mutated comparator and a missing arm.
- Assert the reachable-outcome set equals the check order, with the padding non-vacuity proof.
- Assert the budget-versus-representation distinction is present in the bar text **and** references the probe's epoch curve by field name, so it cannot decay into an adjective.

### Step 6: Score the expression-heads arm
Files: scripts/expression_heads_score.py, wiki/expression_heads_result_2026-11-02.json, tests/test_expression_heads_score.py
Depends on: Step 1, Step 4

**What changes:**
- One forward pass over the 44,106 tiles computing the shared base latent, then all five heads applied to that same latent. The one-pass design is justified by the measured fact that all five share byte-identical base weights, not by convenience.
- Imports `GATE_FUNCS` and `ScoreCtx` from `positional_arm_score` with `primary_rung` set to this arm's primary, exactly as the as-given scorer does. Gates run in the bar's order: split equality, then a head-specific validity check (finite and non-constant output per head), then the shuffled-label null.
- Scores each head's windowed output against the committed atlas per-seed values on byte-identical splits with the bar's corrected interval, and takes the arm's metric as the maximum over heads.
- Checkpointed by batch, stamped `PARTIAL` under `--limit`, and the loader refuses a `PARTIAL` block so a smoke run cannot be read as the arm's result.
- Writes per-tile predictions under `data/`, never `wiki/` — a large decimal dump in `wiki/` has previously polluted the contract-number audit's decoy pool and broken two of its guards.

**Test strategy:**
- Drive every `INDETERMINATE_*` branch with a stub head and assert each is reachable.
- Assert a constant-output head trips the validity gate rather than scoring near zero.
- Assert a head whose output is monotonically transformed scores an **identical** Spearman — the property that makes denormalisation irrelevant to the verdict.
- Assert the verdict is re-derivable from the artifact alone.
- Assert the `PARTIAL` stamp is written under `--limit` and that the loader refuses it.

### Step 7: Run and score the fine-tune arm
Files: scripts/finetune_arm_run.py, scripts/kaggle/build_finetune_arm_kernel.py, wiki/finetune_arm_result_2026-11-02.json, tests/test_finetune_arm_run.py
Depends on: Step 5

**What changes:**
- Trains the encoder per seed on byte-identical `leave_peak_out_indices(frac=0.25)` splits, checkpointing per seed so a crash costs one seed rather than the run.
- The Kaggle kernel builder pins `machine_shape: NvidiaTeslaT4` explicitly, because `enable_gpu: true` alone provisions a P100 that current Kaggle torch cannot run in half precision — a recorded trap. It carries the same `EXPECTED_PARAMS` refusal and `PosixPath` guard as the existing tile-arm kernel, and pushes to a **fresh slug**, because Kaggle keeps one output per kernel and a slug holding retrieved output must never be re-pushed.
- Emits per-seed scores, the paired bootstrap, the null, and one verdict from the closed set.
- A capped run is stamped `PARTIAL` and the scorer refuses it.

**Test strategy:**
- Drive the gate chain with a stub trainer and assert every branch is reachable.
- Assert that resuming after an interrupted run reproduces the same per-seed scores as an uninterrupted one, and that only unfinished seeds are retrained.
- Assert the `PARTIAL` stamp under `--limit` and the scorer's refusal of it.
- Assert the kernel metadata pins the T4 explicitly, and that the chosen slug is not one that already holds output.

### Step 8: Record the three arms, render the memo, journal it
Files: soraya_runs/arms-ledger.json, scripts/successor_arms_memo.py, wiki/successor_arms_result_2026-11-02.md, project_state/glm-learned-representation-2026-10-07.md
Depends on: Step 6, Step 7

**What changes:**
- Moves each of the three arms from `planned` to `measured` by **direct dictionary write**, because `arms_ledger.measure_arm` ends in `Arm(**base)` with no field filter and raises `TypeError` on the annotation keys these arms carry — after mutating the dictionary, which is the worse half. Confirmed empirically by test on the positional arm.
- Each arm's metric is the maximum over **its own** declared primaries, never over a ladder and never across arms.
- Recomputes `ceiling_verdict` and reports it. It will read `met` again once no `planned` arm remains; per decision D2 the memo states it and the run **stops** rather than recording it in passing.
- **Expects the `best` pin to fire**, and that is the guard working rather than drift: the as-given arm's 0.5086 displaces the atlas arm's 0.5076, so the leader changes. `best` was deliberately left an exact pin paired with an argmax invariant, with a message instructing a re-pin rather than a relaxation. Re-pin it and say so in the commit.
- Renders one memo covering all three arms **from their artifacts**, ordered by required reader action, carrying the observed-before-declared hazard on the as-given arm, the windowing approximation and the `tslpred` endpoint mismatch on the heads arm, and the budget-versus-representation distinction on the fine-tune. The memo refuses to render any arm on an `INDETERMINATE_*` verdict.
- Appends one Action-Log row per arm with a class from the seven-value enumeration, anchored inside the `<!-- project-state:end:action-log -->` marker.
- Re-runs the converted invariant guards and asserts they are still green after three status moves, which is the whole point of having converted them.

**Test strategy:**
- Assert the ledger round-trips with every annotation key preserved, by semantic diff.
- Assert the memo generator refuses an `INDETERMINATE_*` verdict.
- Assert the invariant guards hold with three more measured arms.
- Re-run the full suite and diff the failure set **by name** against the 13-name reference baseline.

## Execution Preview

- **Wave 0:** Step 1, Step 3 (parallelism 2)
- **Wave 1:** Step 2, Step 5 (parallelism 2)
- **Wave 2:** Step 4, Step 7 (parallelism 2)
- **Wave 3:** Step 6
- **Wave 4:** Step 8

Five waves, maximum parallelism 2, critical path Step 1 → Step 2 → Step 4 → Step 6 → Step 8 at five steps. No two steps in a wave share a file.

**Run this sequentially regardless.** Parallel mode uses git worktrees, and the feature blocks this work reads are gitignored and live on the `D:` drive, so a worktree cannot run the parity or geometry gates and would silently skip them. That is the same reason the in-flight plan ran sequentially, and the conditions have not changed. Step 7 is the long pole and may route to Kaggle, in which case its wall-clock is dominated by queue time outside this session.

## Risk Flags

- **This plan supersedes Steps 8 through 14 of an in-flight plan, and two artifacts will disagree until that is resolved.** The execution state file records the order `[1, 4, 5, 2, 3, 6, 8, 11, 7, 9, 12, 10, 13, 14]` with seven steps complete. Per the save-once rule, prefer editing the original plan in place and setting its index `Status:` field over leaving two live plans covering one body of work; if this is saved as a new plan, mark the original `superseded` with a `Superseded by:` pointer and retire its state file deliberately rather than leaving it resumable. Do not run `/execute-plan` against both.
- **The clone lives in a volatile temporary directory** (`C:/Users/LAPTOP~1/AppData/Local/Temp/pa_probe/PromoterAtlas`). It resolves today and its checkpoint hash matches, but a reboot can remove it. Steps 2, 3, 6 and 7 all re-verify the hash and refuse rather than proceeding; a clone that has vanished blocks the arms rather than silently changing them.
- **`load_atlas` aliases `pathlib.PosixPath` process-globally and never restores it.** Pre-existing and working, but it is load-bearing for every step here, and a future refactor that "cleans it up" without replacing it breaks all four model-touching steps on Windows.
- **The `tslpred` head's endpoint does not match the rung.** Including it is a deliberate choice under decision D1 and is recorded as a limit, not resolved by it.
- **The fine-tune arm can produce an uninterpretable null.** This is the arm most likely to yield a number nobody can read, which is exactly why Step 3 gates Step 7 and why a probe without a convergence curve means the arm should not run.
- **`best` will fire at Step 8.** Expected, not a defect. Treat a re-pin as the correct response and a relaxation as the wrong one.
- **The architecture gate's final compare may report exit 1 advisorily.** The previous compare flagged only `complex_fn_count` (219 → 223) with cycles, god files and coupling flat. Step 1 adds a module and Steps 2 through 8 add scripts, so the count will move again. Advisory only; do not re-baseline without saying so.
- **Do not run the full suite while a torch job holds a model** — pytest dies with an access violation that reads like a test failure. Steps 2, 3, 6 and 7 all hold models.
- `torch` is not in the default dependencies; use `uv sync --extra forward`, and note `[forward]` and `[ml]` are mutually exclusive.

## Open Questions

1. **Four heads or five, and therefore α/4 or α/5?** Carried as decision D1 above. This plan proceeds on all five under your standing instruction; it is answerable up until Step 4 freezes the bar, after which it becomes a frozen-bar edit.
2. **Whether the segmentation head is in scope.** Settled here as no, and recorded rather than asked: `promoteratlas-annotation.pt` is not an expression head, carries no `normalization_params`, names only three pickled globals, and has a different output shape. It is excluded by kind, not by preference.
3. **Where the fine-tune runs.** Deliberately unanswered — Step 3 measures it and Step 5 freezes it. Naming the route now would be the estimate this arm exists to avoid.

## Verification

- Steps 1 through 3 are verified by their own tests plus a non-zero exit on each refusal path, and the two forward routes' interior-agreement/edge-disagreement test is the one that proves the windowing measurement means something.
- Steps 4 and 5 are verified by `verify()` refusing each mutated anchor and by the reachable-equals-closed-set equality with its padding non-vacuity proof.
- Steps 6 and 7 are verified by every `INDETERMINATE_*` branch being driven to reachability with stubs, by the verdict being re-derivable from the artifact alone, and by the monotone-invariance assertion.
- Step 8 is verified by a semantic diff of the ledger showing every annotation key preserved, by the invariant guards staying green across three status moves, and by a by-name diff of the suite's failure set against the 13-name reference.
- Finally: the five frozen files stay byte-identical, the prospective lock re-verifies, and the architecture compare is read as advisory.

## Save-time amendments

Captured at: 2026-11-02
Source: `/save-plan` arguments

These are AUDIT NOTES for human readers, not executable instructions. `/execute-plan` reads only
`## Implementation Steps`; nothing here changes a step's file list, dependencies or structure.

- Supersedes Steps 8-14 of the successor-arms plan
- from_pretrained refuted
- five heads not four so alpha/5
- head loading routed through atlas_embed's static pickle audit
- fine-tune anchored to the checkpoint's own config

<!-- toolkit: check=clean waves=clean gate=fired:open-questions -->
