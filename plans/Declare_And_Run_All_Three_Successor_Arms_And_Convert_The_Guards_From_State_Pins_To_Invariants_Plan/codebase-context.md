**The ledger API, read rather than assumed.** `Arm` is a 7-field dataclass (`name`, `cost_tier`,
`status`, `metric`, `evidence`, `declared_at`, `resolved_at`). `measure_arm` ends in `Arm(**base)`
with **no field filter**, unlike `kill_arm` which filters to `Arm.__dataclass_fields__` — so it
raises `TypeError` on any arm carrying annotation keys, and does so AFTER mutating the dict. That is
why every write in this plan is a direct dict write. `ceiling_verdict` is a pure function of the
sidecar: `planned` non-empty -> `blocked:untried-arm`; else `len(measured) < 2` ->
`blocked:single-arm`; else `met`. `best` is `max(measured, key=metric)`.

**What the eight red guards actually pin.** `test_arms_ledger_stack_arm.py` holds module constants
`MEASURED` (a dict of 11 arm -> metric), `EXPECTED_PLANNED` (`{positional-...-150bp}`) and
`KILLED_ARM`, and asserts exact equality against each. Its own docstring records the standing rule
this plan relies on: *"had re-pinning flipped the verdict the correct action would have been to stop
and escalate instead"* — which is why declaring an arm first makes Step 2 a re-pin rather than a
rewrite. `test_glm_tile_arms_declared.py` reads the JSON **from disk** and never constructs a ledger
to ask it about itself, because three writes once reported success while the file's mtime never
moved; its failing assertion already carries the instruction *"replace this with an assertion about
the settled verdict instead of deleting it"*. `test_tile_chain_substrate_wall_amendment5.py` asserts
the arm left `planned` only by naming its own artifact, and that `blocking()` still names the
criterion.

**The expression heads are geometry-locked, and the lock is arithmetic rather than a convention.**
`ExpressionPredictor.__init__(self, seq_length=86)` computes `fc_input_dim = base_dim * seq_length`,
and all four committed heads have `fc1.weight` of shape `(64, 11008)` where `11008 = 128 x 86`. A
150 bp tile flattens to `128 x 150 = 19200`, so no head can take a 150 bp latent. `forward` runs
`conv1 -> relu -> conv2 -> relu -> pool -> flatten -> fc1 -> relu -> fc2` over the latent, with
`conv1`/`conv2` at `kernel_size=8, padding='same'` and `pool` a `MaxPool1d(7, stride=1, padding=3)`
— all position-local, which is what makes an 86-wide window over the full 150 latent nearly
equivalent to a true per-window recompute, and the `padding='same'` edges are exactly where it is
not. `from_pretrained(cls, model_path)` is the clean loader. Each checkpoint carries
`normalization_params` {mean, std}; Spearman is invariant to that denormalisation, so it cannot move
a verdict, but an absolute predicted value would be wrong without it.

**The heads sit on the FROZEN base, measured not assumed.** Each head checkpoint embeds its own
`base_model.*` tensors, and all four are **byte-identical to `promoteratlas-base.pt`** — 124/124
tensors, max |dev| **0.000e+00**. Head total is 2,814,533 params against the base's 1,847,428, so
the head itself is ~967k params. `__init__` also sets `requires_grad = False` on the base, which
corroborates the identity. This is what makes a heads arm "frozen encoder + trained head" rather
than a fine-tuned-encoder arm, and it is re-checked in Step 8 because a divergence would change the
claim.

**The fine-tune is tractable here precisely because the encoder is small.** 1,847,428 parameters,
not a 1.2B LM, so this is not comparable to the project's blocked large-LM fine-tunes.
`dna_decode/glm/encoder.py::fit_encoder(train, test, *, widths, protocol, seed, device,
shuffle_labels, group_key=tile_peak_key)` already implements the training protocol this arm mirrors,
including early stopping against a **peak-blocked inner split of `train`** with `test` never seen
until the final score, and a `shuffle_labels` validity null that permutes TRAINING targets only.
Local compute is CPU-only (`torch.cuda.is_available()` **False**, 4 threads); forward throughput
measured at **57.5 seq/s**. `kaggle` CLI **2.2.5** resolves and `scripts/kaggle/` already holds seven
kernel builders, so the GPU route is established rather than speculative.

**The stored positional block is reusable for one arm and not the others.**
`D:/dna_decode_cache/promoteratlas_positional/positional_block.npy` is 44,106 x 5,120 float32,
sidecar `COMPLETE`, layout `bins10_mean [0:1280] / bins10_max [1280:2560] / autocorr [2560:5120]`
with the checkpoint sha256 stamped. The as-given arm needs only this block, so it costs minutes. The
heads arm needs the RAW `(B, 128, 150)` latent, which the block does not store, so it needs its own
pass.

**The comparator and split surfaces are unchanged and re-readable.**
`data/promoteratlas/local_reconcile_2026-11-01.json` holds `verdict: RECONCILED`, `mismatched: []`,
the atlas per-seed values (median 0.5076), the GC per-seed values and `n_test_per_seed` for all ten
seeds. `genomewide.leave_peak_out(tiles, frac=0.25, seed=s)` holds out whole PEAKS and membership is
a pure function of the peak id, so disjointness is by construction.
`atlas_embed.ridge_spearman` mirrors `expression.fit_predict_ridge` exactly and already carries the
`shuffle_seed` parameter the nulls use.

**The constraint rail's shape.** `test_the_eleven_pre_existing_constraints_are_still_present` names
11 constraint ids individually — its docstring is right that *"a count alone would pass if one were
swapped"* — and then asserts `len(rows) == 12` against a registry that now holds **14**. The names
are the guard; the trailing total is the only part that drifts on an append.

**The refactor target, measured.** `scripts/positional_arm_score.py::run` is cyclomatic ~58 over 304
lines (AST estimate); the next highest in the same family are `positional_arm_probe.run` at ~24 and
`positional_arm_extract.extract` at ~20. Its 43 existing tests already drive every gate branch, so
they are the regression suite for Step 5 rather than something to be written.

### Reusable-Code Survey

- `dna_decode/glm/positional_readout.py` — `bin_means`, `bin_maxes`, `global_mean_from_bin_means`,
  `channel_autocorrelation`, `orient`, `flatten_bins`. The as-given arm reuses `orient`/`flatten_bins`
  verbatim; nothing new is needed for it.
- `scripts/positional_arm_probe.py` — OWNS the cohort contract (`load_arm_cohort`, which refuses any
  count but 44,106), plus `design_matrix_bytes`, `fits_in_ram`, `free_ram_bytes`,
  `extrapolate_wall_seconds`, `stride_sample`, `parity_report`, `sha256_of`. Every new probe and
  scorer in this plan imports these rather than re-deriving them.
- `scripts/positional_arm_extract.py` — `feature_layout`, `slice_block`, `as_bins`, `rebin`,
  `load_block` (which refuses a PARTIAL or incomplete block). The as-given scorer reuses `load_block`
  unchanged; the heads runner reuses the memmap-plus-sidecar resume pattern.
- `scripts/positional_arm_score.py` — `leave_peak_out_indices`, `paired_bootstrap`,
  `order_statistic_median_ci`, `median_of`, and after Step 5 the extracted per-gate helpers. All
  three new scorers compose these.
- `scripts/freeze_positional_arm_bar.py` — `bonferroni_percentiles`, `derive_parity_tolerance`,
  `reachable_outcomes`, and the freeze-with-`verify()`-refuses-on-drift pattern the three new bars
  copy.
- `dna_decode/glm/atlas_embed.py` — `load_atlas`, `audit_checkpoint`, `EXPECTED_PARAMS`, `NUC`,
  `ridge_spearman`. Reused verbatim; not modified.
- `dna_decode/glm/encoder.py` — `fit_encoder` and its `TrainingProtocol`, the reference the fine-tune
  arm mirrors rather than reinventing.
- `dna_decode/glm/genomewide.py` — `load_peak_tiles`, `leave_peak_out`, `spearman`, `gc_features`.
  READ-ONLY; it is the registered owner of GSE144621 and this plan adds no second reader.
- `scripts/kaggle/build_promoteratlas_tile_arm_kernel.py` — the existing kernel-builder pattern the
  fine-tune's Kaggle route copies, including the `machine_shape: NvidiaTeslaT4` pin.
- `scripts/positional_arm_memo.py` — the render-from-artifact-and-refuse-on-INDETERMINATE pattern the
  combined memo copies.
- Searched and found nothing reusable for the head windowing itself:
  `graphify-out/GRAPH_REPORT.md` (absent), `src/`, `lib/`, `utils/`, `common/`, `helpers/` (none
  exist in this layout), and `.repo-index/` (absent on this host).
