# Resolve The Scoring-Frame Confound Behind The Bacterial-Versus-Human Gap

## Lens status

Execution mode recommendation: **sequential** (see Risk Flags R3 — the substrate lives under gitignored `data/`, so a worktree has no inputs).
gh available: yes · sentrux available: yes (plan is non-structural — gate will skip) · codex: absent
Repo-index: absent — no degradation, direct file reading used throughout
Wiki replay: 3 HIGH-salience entries consulted (tile-substrate-wall 2026-10-21, promcalc regeneration 2026-10-22, zero-shot-vs-supervised 2026-10-23)

## Problem Statement

The eukaryotic arm reported a conv-over-GC gain of **+0.2471** (HepG2) and **+0.3016** (K562) against **+0.0122** on bacterial promoter tiles, and that ratio has been quoted as "20–25× more learnable signal on human regulatory DNA". The arm's own first honest limit states it **cannot** say which of four differences — assay, element length, kingdom, library construction — carries the effect.

A scratch measurement on 2026-11-16 (`wiki/glm_bacterial_scoring_frame_2026-11-16.md`) found a fifth candidate the arm never considered, and it is large. The bacterial label is sharply bimodal: the `active` flag is a **pure threshold on the measurement itself** at 0.5778 with zero overlap either side, **29,193 of 44,106 tiles (66%) sit below it**, and a single binary split of the label recovers **0.8194** of the 0.941 ceiling. So the pooled Spearman is overwhelmingly a *detection* metric. Restricted to the 14,913 active tiles the learned-versus-GC comparison **flips sign unanimously** — a 4-mer ridge loses to GC on 0 of 10 seeds pooled (margin −0.1793) and beats it on 10 of 10 within actives (margin +0.1820) — because GC is a detector and nothing else (pooled +0.3112, within actives **−0.046**).

That scratch result has three gaps that make it suggestive rather than decisive, and this plan closes all three:

1. It used a **4-mer ridge, not the published conv**. The conv is a materially better detector (0.3234 pooled, *above* GC) than the proxy (0.1438, well below), so its active-frame behaviour is genuinely unknown.
2. It never touched the **strongest arm**. The as-given arm (0.5086) and the frozen-embedding arm (0.5076) beat GC by **+0.1964** pooled, so their active-frame margin is the realistic upper bound on this substrate.
3. It could not measure the **human side at all**, which is the quantity that decides whether the 20–25× is a frame artifact.

**Non-goals.** This plan does not re-run any human model, does not retune any bacterial arm, does not edit `scripts/glm_encoder_gate.py::PREREGISTERED`, and does not touch the eukaryotic arm's artifacts or `verdict` field. Re-interpreting a published verdict is an authority call; producing the measurements beside it is not.

## Codebase Context

**`dna_decode/glm/genomewide.py` (559 LOC)** — owns the substrate. `Tile` carries `seq`, `expression`, `peak`, `strand`, `category`, `active`, `rep1`, `rep2`. `load_peak_tiles(cache_dir, *, allow_download)` returns `(list[Tile], LoadStats)` and `ensure_cached` raises `GenomeWideDataError` with the `external_root` report when the configured volume is detached — it refuses to auto-substitute. `leave_peak_out(tiles, *, frac=0.25, seed)` returns `(train, test, n_held_peaks)` and holds out **whole peaks**, so disjointness is by construction. `noise_ceiling(frags)` is the Spearman-Brown step-up over `rep1`/`rep2` and already returns `single_replicate_agreement` + `ceiling_for_a_sequence_model`. Feature builders `gc_features`, `kmer_features(seqs, k)`, `positional_onehot` all present.

**Substrate is LOCAL and verified.** `data/mpra/GSE144621_peak_tile_expression_formatted_std.txt.gz` reproduces every published anchor: `n=44,106` after `category=="tile"` + `len==150`, `14,913` active / `29,193` inactive, single-replicate `0.7945`, ceiling `0.9410`, GC `0.3110` against the published `0.3112`, and **GC seed 0 exactly `0.3000`** — the committed reconcile anchor. The default cache (`D:/dna_decode_cache/mpra`) is detached; pass `data/mpra` explicitly.

**`scripts/glm_encoder_gate.py` (413 LOC)** — the published conv gate, and the active flag is **already plumbed through it**: `_active_labels(test)` exists at L142, `auroc(scores, labels)` with tie handling at L120, and `main()` already prints `{n_active} active / {n_inactive} inactive`. The flag feeds the secondary AUROC framing and **never restricts the Spearman**. `PREREGISTERED` (L65) is frozen and carries `reconcile_anchor = {"arm": "gc", "seed": 0, "expected": 0.3000, "tolerance": 0.01}` plus six verdicts including `INDETERMINATE_RECONCILE_FAILED`. `ridge_gc(train, test)` returns `(spearman, predictions)` and fits `Ridge(alpha=1.0)` on `StandardScaler`-scaled `gc_features` against `np.log1p(expression)`. The established pattern for a new arm is a **flag that leaves the published arms byte-identical** (`--pyramid`, `--extra-widths`).

**`dna_decode/glm/encoder.py` (439 LOC)** — `fit_encoder(train, test, *, widths=DEFAULT_WIDTHS, protocol=REGISTERED_PROTOCOL, seed, device, shuffle_labels, group_key=tile_peak_key)` returns a dict with `spearman`, `predictions`, `n_parameters`, `epochs_ran`, `device`, and the protocol actually used. `DEFAULT_WIDTHS = (6, 20, 75)` = the published 13,217-param arm. `shuffle_labels` permutes **training** targets only — the validity null. Early stopping runs against a peak-blocked inner split via `peak_blocked_val_split`; `test` is never seen until the final score.

**Committed per-tile predictions exist, so the strongest arm needs no retraining.** `data/promoteratlas/asgiven_arm_pertile_2026-11-02.json` (1.83 MB, `rung: bins10_mean_asgiven`) and `data/promoteratlas/per_tile_predictions_2026-11-01.json` (2.04 MB) each hold 10 seeds of `{tile_indices, pred}`. Lengths match the committed split exactly — seed 0 → 11,198, seed 1 → 11,018 — against `wiki/asgiven_arm_result_2026-11-02.json`'s `gates/split_equality/per_seed` (`held_out_peaks: 608` every seed). Committed per-seed Spearmans are available for both the as-given rung (`rung/per_seed`, median **0.5086**) and the atlas comparator (`comparator/atlas_per_seed` = `{0: 0.5461, 1: 0.5499, 2: 0.4983, 3: 0.4986, 4: 0.5075, 5: 0.5163, 6: 0.5078, 7: 0.4933, 8: 0.5016, 9: 0.5210}`, median **0.5076**).

**`dna_decode/glm/lentimpra.py` (444 LOC)** — `DEFAULT_DIR = D:/dna_decode_cache/mpra_human` (detached). `load_cell(cell, data_dir, *, min_join_rate)` joins design to activity and **refuses below the join rate**; `load_activity` returns `{sequence_id: (activity, fold)}` from Table S6; `usable_cells(data_dir)` reports per-cell usability with a reason. `LentiRecord` carries `name`, `seq`, `activity`, `fold`, `cell`. Activity is `log2(RNA/DNA)` and goes negative. **WTC11 refuses** at join rate 0.0000.

**`scripts/fetch_lentimpra.py` already exists and must be reused, not reimplemented.** CLI: `--dest` (default the detached `D:` path, overridable), `--select tables|all` (`tables` is the ~74 MB four-table subset, not the 3.5 GB + 4.33 GB dumps), `--list`, `--yes`. It verifies **every byte against the publisher's Zenodo MD5**, and the MD5 — not the size — is the integrity guard.

**Dependency reality, measured not assumed:** `openpyxl` **3.1.5 is PRESENT**, `torch` 2.7.1+cu118, `sklearn` 1.8.0, `scipy` 1.17.1. `C:` has **33.8 GB free of 511 GB (93.4% full)** — ample for a 74 MB fetch.

### Reusable-Code Survey

Candidates found and adopted rather than rebuilt:
- `dna_decode/glm/genomewide.noise_ceiling` — the replicate/ceiling computation, reused verbatim on tile subsets (`Tile` satisfies its `rep1`/`rep2` duck type; the published 0.9410 comes through it).
- `dna_decode/glm/genomewide.leave_peak_out` — the split, reused unchanged so split discipline is inherited, not reinvented.
- `dna_decode/glm/genomewide.{gc_features, kmer_features}` and `encoder.fit_encoder` — the arms, reused unchanged so the comparison stays architecture-matched to the published +0.0122.
- `scripts/glm_encoder_gate.{ridge_gc, auroc, _active_labels}` — imported, not copied; `ridge_gc`'s docstring records that it was re-implemented once already and that the reconcile gate is what proves fidelity.
- `scripts/fetch_lentimpra.py` — the MD5-verified fetcher, invoked as a subprocess with `--dest` overridden.
- `scripts/asgiven_arm_score.py` + `wiki/asgiven_arm_result_2026-11-02.json` — the committed per-seed values that become Step 2's reconcile anchor.

Searched: `graphify-out/GRAPH_REPORT.md` (absent), `dna_decode/glm/`, `dna_decode/eval/`, `scripts/glm_*`, `scripts/*arm*`, `tests/test_glm_*`.

## Pre-Change Baseline

Published, on the pooled frame (`wiki/glm_encoder_gate_2026-10-08.json`): GC **+0.3112** (2 params), `encoder-w20` **+0.3163** (2,657), `encoder-multi` **+0.3234** (13,217), `encoder-shuffled` **−0.0240**. Paired margin `encoder-multi` − GC = **+0.0125 mean, 6/10 seeds, 95% CI [−0.0099, +0.0349]** (includes zero); frozen verdict `WEAK`. Best arm overall **0.5086**, pooled margin over GC **+0.1964**. Ceiling **0.941**.

Scratch, 2026-11-16, 4-mer-ridge proxy, same splitter and seeds: pooled GC **+0.3112** / 4-mer **+0.1438** / margin **−0.1793** / 0 of 10 seeds; active-only GC **−0.0464** / 4-mer **+0.1320** / margin **+0.1820** / 10 of 10 seeds. Replicate agreement pooled **0.7945** (ceiling 0.9410), active **0.8849** (0.9690), inactive **0.3927** (0.7510). Threshold-only-versus-label Spearman **0.8194** = 87.1% of the ceiling.

Human side (`wiki/eukaryotic_arm_result_{HepG2,K562}_2026-11-01.json`): HepG2 conv **0.4853** / GC **0.2362** / gain **+0.2471** `[+0.2399, +0.2536]`; K562 conv **0.5237** / GC **0.2200** / gain **+0.3016** `[+0.2916, +0.3100]`. **No human replicate ceiling exists** — both artifacts carry the *bacterial* anchor (0.941) as context.

Test suite: **13 failed / 8,717 passed / 204 skipped**, the 13 being the standing absent-data set (`test_claude_md_growth_ceiling` 1, `test_cyp4f2_getrm_concordance` 7, `test_essentiality_report_card_ladder` 3, `test_pneumo_selection_probe` 1, `test_readme_quickstart_output` 1).

## Verification Signal

This plan has worked when all of the following hold:

1. **Step 2 reconciles.** Pooled per-seed Spearman recomputed from the committed per-tile predictions matches the committed per-seed values within 0.001 for all 10 seeds, for whichever arm each file holds — and the arm identity is established **by that reconcile, not by filename**. A mismatch exits non-zero and reports nothing.
2. **Step 3 reconciles on the frame that has a published anchor before changing frames.** GC seed 0 on the **pooled** frame reproduces `0.3000 ± 0.01`, and `encoder-multi` reports 13,217 parameters. Only then is the active frame scored.
3. **Step 3's active-frame validity null is clean** — the shuffled-label arm scores `|spearman| < 0.05` on the active frame. A new frame does not inherit the old frame's null.
4. **Step 4 either measures or refuses by name.** Either the human label's shape is reported with its join rate, or the script exits non-zero naming the specific wall (fetch failure, MD5 mismatch, join rate below threshold). No silent partial.
5. **No new test failures.** `uv run pytest tests/ -q` shows the same 13 absent-data failures and nothing else, plus the new tests passing.
6. **Nothing frozen moved.** The frozen five + `uv.lock` + `wiki/finetune_arm_acceptance_bar.json` byte-identical; `prospective_lock.verify_lock` returns `ok=True` on the 2026-08-31 manifest; `scripts/glm_encoder_gate.py::PREREGISTERED`, `wiki/glm_encoder_gate_2026-10-08.json` and both `wiki/eukaryotic_arm_result_*.json` byte-identical.
7. **Every new guard is proven non-vacuous** — each new test file is shown to fail when its subject is reverted, and the mutation is confirmed to have applied.

## Implementation Steps

### Step 1: Pure active-frame decomposition module
Files: dna_decode/glm/active_frame.py, tests/test_glm_active_frame.py
Depends on: none

**What changes:**
- New module holding every pure function the three measurement steps share, so the frame logic has one definition.
- `threshold_audit(tiles) -> dict` — establishes whether `active` is a clean threshold on `expression`: returns the implied cut, the count on each side, and `overlap_either_side`. **It does not assume a threshold; it measures one.** If overlap is non-zero the audit says so and callers must treat the flag as a label rather than a cut.
- `frame_decomposition(tiles) -> dict` — pooled / active / inactive blocks, each with `n`, `single_replicate_agreement` and `ceiling`, computed by delegating to the existing `noise_ceiling` rather than re-deriving Spearman-Brown.
- `label_bimodality(values) -> dict` — substrate-agnostic shape measures: best single-threshold-versus-label Spearman and the share of a supplied ceiling it recovers, the modal-bin share, and the count of distinct values. Takes raw floats, so Step 4 can call it on human activity where no `active` flag exists.
- `restrict_predictions(tile_indices, preds, tiles, keep) -> tuple[list, list]` — maps stored per-tile predictions onto a tile subset, raising on any index out of range so a mis-aligned prediction file fails loudly instead of scoring a shifted vector.
- Every returned dict carries a `conditions_on_the_outcome: bool` field, true for any block defined by a threshold on the label, so downstream artifacts cannot present a range-restricted number as unconditioned.

**Test strategy:**
- Synthetic tiles with a planted clean threshold → `threshold_audit` recovers the cut and reports zero overlap; with one deliberately mislabelled tile → reports non-zero overlap and does not claim a cut.
- `frame_decomposition` on synthetic replicates with known agreement reproduces the hand-computed Spearman-Brown ceiling.
- `label_bimodality` on a unimodal uniform sample versus a planted two-cluster sample separates them; on an all-identical sample it refuses rather than returning a vacuous cut (the recorded trap: a generic scanner must not return a cut on a degenerate ladder).
- `restrict_predictions` raises on an out-of-range index and on a length mismatch between indices and predictions.
- Real-data test gated on `data/mpra/...` being present: reproduces `44,106 / 14,913 / 29,193`, pooled ceiling `0.9410`, and the `0.5778` cut with zero overlap. Skips with a reason when the file is absent.

### Step 2: Re-score the committed arms on the active frame, with no retraining
Files: scripts/active_frame_rescore.py, tests/test_active_frame_rescore.py
Depends on: Step 1

**What changes:**
- Reads `data/promoteratlas/asgiven_arm_pertile_2026-11-02.json` and `data/promoteratlas/per_tile_predictions_2026-11-01.json`, loads tiles from `data/mpra` via the explicit-cache path, and re-scores each stored prediction vector on the active subset. **No model is trained and the `D:` volume is not needed.**
- **Reconcile gate runs first and is blocking.** For each file and seed, recompute the pooled Spearman from the stored predictions and compare against the committed per-seed values in `wiki/asgiven_arm_result_2026-11-02.json` (`rung/per_seed` for the as-given rung, `comparator/atlas_per_seed` for the atlas arm). Arm identity is **assigned by whichever committed series matches**, and if neither matches within tolerance the script exits non-zero with `INDETERMINATE_ARM_IDENTITY_UNRESOLVED`.
- Also cross-checks each seed's `n_test` against `gates/split_equality/per_seed` so a prediction file from a different split cannot pass on Spearman agreement alone.
- Reports, per arm: pooled per-seed and median (the reconciled values), active-frame per-seed and median, the paired per-seed margin against `ridge_gc` re-scored on the same active test rows, win count, and a paired bootstrap CI on the margin.
- Writes `wiki/active_frame_rescore_<date>.json` under its own schema. Does not write into any existing arm artifact.

**Test strategy:**
- Unit: a synthetic prediction file whose pooled recompute matches a synthetic committed series resolves to that arm; one that matches neither exits with the indeterminate code; one whose `n_test` disagrees fails even when Spearman agrees.
- Unit: the reconcile tolerance is asserted to reject a 0.01 perturbation, so it is not so loose as to be vacuous.
- Real-data test gated on the prediction files: asserts the as-given series reconciles to median `0.5086` and the atlas series to `0.5076`.
- Non-vacuity: perturb one stored prediction vector and confirm the gate fires.

### Step 3: Run the published conv on the active frame, under its own frozen bar
Files: scripts/glm_active_frame_gate.py, tests/test_glm_active_frame_gate.py
Depends on: Step 1

**What changes:**
- A **separate script with its own frozen pre-registration**, importing `ridge_gc`, `auroc`, `_active_labels` from `scripts/glm_encoder_gate.py` and `fit_encoder` / `DEFAULT_WIDTHS` from `dna_decode.glm.encoder`. `scripts/glm_encoder_gate.py::PREREGISTERED` is **not edited** — the active frame breaks its reconcile anchor by construction (GC seed 0 is ≈ −0.04 there, not 0.3000), so sharing the dict would retroactively rewrite a pre-registration and force `INDETERMINATE_RECONCILE_FAILED`.
- **Two-phase reconcile, and the ordering is load-bearing.** Phase 1 scores GC and `encoder-multi` on the **pooled** frame and asserts the published anchor (GC seed 0 `0.3000 ± 0.01`, 13,217 parameters). Only after that passes does Phase 2 score the active frame. Validating on the frame that *has* a published anchor is the only way to prove the data and split are right before the frame changes.
- Arms on the active frame: `gc`, `encoder-multi` (widths `(6, 20, 75)`, the published 13,217-param arm) and `encoder-shuffled` as the active-frame validity null. Registered seeds `range(10)`, matching the published list.
- The frozen bar is stamped `not_a_blind_pre_registration: true` with the reason — the 2026-11-16 scratch measurement motivated it, so the thresholds were set with a directional expectation in view. The repo's own standing rule is that a bar frozen *before* the mechanism analysis is mis-specified; this one is frozen after and says so.
- Pre-committed outcomes: `FRAME_EXPLAINS_THE_GAP` (active-frame margin ≥ +0.15 and ≥ 8 of 10 seeds), `FRAME_CONTRIBUTES` (margin > +0.05, CI excludes zero), `FRAME_DOES_NOT_EXPLAIN` (margin ≤ +0.05 or CI includes zero), `INDETERMINATE_NULL_NOT_CLEAN`, `INDETERMINATE_POOLED_RECONCILE_FAILED`.
- Pins `OMP_NUM_THREADS` / `OPENBLAS_NUM_THREADS` / `MKL_NUM_THREADS` / `NUMEXPR_NUM_THREADS` to 1 and records the per-seed wall time, per the recorded BLAS-oversubscription lesson.
- Writes `wiki/glm_active_frame_gate_<date>.json`.

**Test strategy:**
- Unit on the verdict function: each of the five outcomes is reachable from a synthetic arm dict, and `INDETERMINATE_POOLED_RECONCILE_FAILED` takes precedence over any graded verdict.
- Unit: the script's import of `PREREGISTERED` is read-only — a test asserts `scripts/glm_encoder_gate.py` is byte-identical to `HEAD` after the module is imported, and that the new script defines its own distinct pre-registration dict.
- Unit: a null arm scoring 0.2 forces `INDETERMINATE_NULL_NOT_CLEAN` regardless of the candidate margin.
- Smoke run on a small tile subset with `--seeds 1` to prove the end-to-end path without the full CPU cost.
- Thread-pinning is asserted by reading the environment the script sets, not by timing.

### Step 4: Fetch the human tables locally and measure the human label's shape
Files: scripts/human_frame_shape.py, tests/test_human_frame_shape.py
Depends on: Step 1

**What changes:**
- Invokes the existing `scripts/fetch_lentimpra.py --select tables --dest data/mpra_human --yes` as a subprocess. The fetcher's MD5 verification is the integrity guard and is **not** reimplemented; the `--dest` override exists precisely so the detached `D:` default can be bypassed.
- Loads HepG2 and K562 via `lentimpra.load_cell(cell, data_dir="data/mpra_human")` and reports each cell's join rate. WTC11 is expected to refuse at join rate 0.0000 and that refusal is recorded, not worked around.
- Runs `label_bimodality` from Step 1 on each cell's `activity` values and on the bacterial `expression` values, and prints them side by side. This is the decisive comparison: whether the human label is similarly detection-dominated.
- Reports explicitly whether the deposit carries **replicate-level** columns. If it does, computes a human ceiling via `noise_ceiling`; if it does not, records `human_ceiling: null` with the reason, because the eukaryotic artifacts currently carry the *bacterial* ceiling as context and that substitution must not be silently inherited.
- **Does not invent a human `active` flag.** lentiMPRA publishes no such call, so the comparison is on distributional shape and best-single-threshold recoverable share, never on a threshold this plan chose.
- Writes `wiki/human_frame_shape_<date>.json`.

**Test strategy:**
- Unit: the subprocess invocation is asserted to carry `--select tables` and the overridden `--dest`, so a future edit cannot silently pull the 7.8 GB dumps or write to the detached volume.
- Unit: a simulated fetch failure, a simulated MD5 failure and a join rate below threshold each exit non-zero with distinct codes; none returns a partial shape measurement.
- Unit: `human_ceiling` is `null` with a reason when no replicate columns are found, and a test asserts the artifact never carries the bacterial 0.941 as a human ceiling.
- Real-data test gated on `data/mpra_human/` being present; skips with a reason otherwise.

### Step 5: Write the resolution artifact and make the frame caveat machine-visible
Files: wiki/glm_scoring_frame_resolution.md, wiki/glm_scoring_frame_resolution.json, soraya_runs/arms-ledger.json, tests/test_scoring_frame_artifact.py
Depends on: Step 2, Step 3, Step 4

**What changes:**
- One artifact consolidating Steps 2–4 into an answer to the original question, with the verdict taken from Step 3's frozen bar and the Step 2 upper bound and Step 4 human comparison beside it.
- States plainly which of the eukaryotic arm's four named differences this resolves and which it leaves open, and carries forward the five honest limits from the 2026-11-16 scratch memo that still apply — in particular that the active frame conditions on the outcome and that its absolutes are not comparable to pooled ones.
- Adds **one** annotation key to the `ceiling:` criterion in `soraya_runs/arms-ledger.json` recording that bacterial metrics under this criterion are measured on a detection-dominated frame, pointing at the new artifact. The writing script asserts every pre-existing field unedited and **refuses to write unless the change is augment-only** — the established discipline from the session-limits layering.
- Does **not** edit any arm's `metric`, `status`, `evidence` or `verdict`, and does not touch the eukaryotic arm artifacts.

**Test strategy:**
- Artifact-shape tests: required keys present, every cited `wiki/` path exists, the JSON parses (the recorded trap — an invalid hand-edited JSON artifact was committed and cited once).
- Ledger augment-only test: arm count, key set, and every `metric` / `status` value identical before and after; exactly one new key added at the criterion level.
- A test asserting the artifact does not state a bare ratio without the frame caveat in the same section.
- Non-vacuity: delete the new annotation and confirm exactly the new tests fail while the pre-existing ledger tests pass.

### Step 6: Propagate the result to the reader-facing surfaces
Files: README.md, dna_decode/eval/regime.py, tests/test_regime_boundary.py
Depends on: Step 5

**What changes:**
- Replaces the interim caveat currently in the README's roadmap step 4 with the measured outcome, under whichever of Step 3's pre-committed verdicts fired. The caveat was written as a placeholder on 2026-11-16 and explicitly says the ratio is unresolved; this step resolves it in one direction or the other.
- Updates the `constructed → molecular` row in `dna_decode/eval/regime.py` only if Step 3 returns `FRAME_EXPLAINS_THE_GAP` or `FRAME_CONTRIBUTES`, since that is the row whose cited evidence the bacterial-versus-human contrast informs. If `FRAME_DOES_NOT_EXPLAIN`, the row is left unchanged and the artifact records that the frame hypothesis was tested and rejected.
- Pre-committed so the write-up cannot be shaped by the result: each verdict maps to one README wording and one regime decision, fixed in Step 3's frozen bar before the run.

**Test strategy:**
- README guards already exist and must stay green: `tests/test_advertised_commands.py` + `tests/test_advertised_env_vars.py` (148 tests) resolve every advertised command through its real parser.
- A test asserting the README no longer contains the interim "unresolved" wording once a verdict has been recorded, and that it cites the new artifact.
- `tests/test_regime_boundary.py` must still make all three historical scope-compressions fail loudly; any row edit is checked not to weaken those.

## Execution Preview

| Wave | Steps | Parallelism |
|---|---|---|
| 0 | Step 1 | 1 |
| 1 | Step 2, Step 3, Step 4 | 3 |
| 2 | Step 5 | 1 |
| 3 | Step 6 | 1 |

Total waves: **4**. Max parallelism: **3**. Critical path (toolkit): **Step 1 → Step 2 → Step 5 → Step 6**. The toolkit's path is topological; the **wall-clock** long pole is Step 3 in the same wave — 10 seeds × 3 arms of CPU conv training, though on 14,913 active tiles rather than 44,106, so roughly a third of the published run's cost.

Wave 1 has no intra-wave file overlap: Steps 2, 3 and 4 each own one script and one test file, and all three only read Step 1's module.

## Risk Flags

**R1 — The active frame conditions on the outcome, and no step may hide it.** `active` is a threshold on `expression` itself, so every active-frame number is range-restricted. The defence is structural: Step 1 stamps `conditions_on_the_outcome` on every such block, and the finding rests on the *differential* between GC and the learned arms under the *same* restriction, not on any absolute. Mitigated, not eliminated.

**R2 — `[unverified]` whether the lentiMPRA deposit carries replicate-level data.** If it does not, no human ceiling can be computed and the human-versus-bacterial comparison stays partly unnormalised. Step 4 is required to record `human_ceiling: null` with a reason rather than inherit the bacterial 0.941, which both eukaryotic artifacts currently carry as context.

**R3 — Sequential execution recommended over parallel worktrees.** The substrate (`data/mpra/`) and the committed predictions (`data/promoteratlas/`) sit under gitignored `data/`, so a worktree would start with no inputs and an unreconcilable baseline — the same reasoning that made the previous two plans in this arc sequential.

**R4 — A frozen pre-registration must not be retrofitted.** Step 3 creates a new bar rather than extending `glm_encoder_gate.PREREGISTERED`, and a test pins that file byte-identical to `HEAD`. Extending it would both break its own reconcile anchor and rewrite a commitment after seeing results.

**R5 — Step 3's bar is not a blind pre-registration and says so.** The 2026-11-16 scratch result motivated the thresholds. The artifact carries `not_a_blind_pre_registration: true` with the reason, following the existing field of that name in `wiki/asgiven_arm_result_2026-11-02.json`.

**R6 — BLAS oversubscription will mis-time Step 3 if threads are unpinned.** The recorded instance cost a 2.34× throughput error and a wrong projection. Step 3 pins all four thread variables to 1 and records per-seed wall time.

**R7 — Do not report the 4-mer proxy and the conv as the same thing.** The scratch margin (+0.1820) came from a 4-mer ridge, not the published conv, and the repo has a recorded instance of a zero-shot result being reported as replicating a supervised one. Steps 2, 3 and 5 must label each number with the arm that produced it.

**R8 — Disk.** `C:` is 93.4% full with 33.8 GB free. Step 4 adds ~74 MB and `--select tables` must not drift to `--select all` (7.8 GB); a test pins the flag.

## Open Questions

1. **If Step 4 finds the human label is also strongly bimodal, the 20–25× is largely a frame artifact and the eukaryotic arm's `SUBSTRATE_LIMITED` verdict needs reinterpreting.** The arm's `verdict` field is deliberately untouched by this plan — changing it is a user authority call, not an analysis choice.
2. **Should the as-given arm's active-frame score become a declared arm in `soraya_runs/arms-ledger.json`, or stay a diagnostic?** Declaring it makes it eligible to set `best` under the `ceiling:` criterion, and the criterion's metric is a pooled-frame quantity — the exact non-commensurability the ledger was just annotated against. Recommendation is to keep it a diagnostic; the decision is the user's.
3. **Whether the bacterial substrate should be re-framed going forward** — scoring future arms on actives only, or scoring detection and strength as two separate declared criteria. This changes what every future bacterial number means and is out of scope here.

## Verification

- `uv run pytest tests/ -q` → same 13 absent-data failures, no new ones, plus the new test files passing.
- Each new test file proven non-vacuous by reverting its subject and confirming only that file fails — with the mutation confirmed to have applied, per the recorded `MUTATION-DID-NOT-APPLY` trap.
- Frozen five + `uv.lock` + `wiki/finetune_arm_acceptance_bar.json` byte-identical to `HEAD`; `prospective_lock.verify_lock` on the 2026-08-31 manifest returns `ok=True` with 5 pinned files.
- `scripts/glm_encoder_gate.py`, `wiki/glm_encoder_gate_2026-10-08.json` and both `wiki/eukaryotic_arm_result_*.json` byte-identical to `HEAD`.
- `soraya_runs/arms-ledger.json` diffed before and after: identical arm count, identical key set, no `metric` or `status` changed, exactly one new criterion-level key.
- Step 2's reconcile reproduces as-given median `0.5086` and atlas median `0.5076`; Step 3's pooled phase reproduces GC seed 0 `0.3000 ± 0.01` and 13,217 parameters.
- `uv run python scripts/project_status.py` still exits 0 and reports the registry unchanged at 128 cells.

## Save-time amendments

Captured at: 2026-11-16
Source: `/save-plan` arguments

**Audit-notes-only contract.** This block is provenance for human readers. `/execute-plan` reads
ONLY `## Implementation Steps` for executable work, so nothing here is an executable instruction.
If an amendment changes a Step contract (file lists, dependencies, structure), re-run
`/technical-plan` before `/execute-plan` rather than relying on this section.

- openpyxl is already installed (3.1.5) so no dependency decision is needed
- the committed per-tile predictions under data/promoteratlas/ let the strongest arm be re-scored with no training and no D: drive
- Step 3 gets its own frozen bar rather than extending the published PREREGISTERED

<!-- toolkit: check=clean waves=clean gate=fired:open-questions,unverified,test-strategy-leak -->
