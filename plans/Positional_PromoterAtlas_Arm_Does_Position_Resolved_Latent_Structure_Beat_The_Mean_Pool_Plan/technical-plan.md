# Positional PromoterAtlas Arm: Does Position-Resolved Latent Structure Beat The Mean-Pool 0.5076

## Lens status

Repo-index: absent — direct file reading only, no degradation.
Tools: `sentrux` present (structural gate eligible), `gh` present, `codex` ABSENT (so a second-opinion
`/brainstorm` round will fall back to single-model and must be labelled as such).
`project-rules.md`: absent — the Project Rules Check section is omitted.
`DESIGN.md`: absent, and this plan touches no UI.
Baseline: measured, not asserted — see Pre-Change Baseline.
Substrate: present and loaded during research (44,106 tiles read from `data/mpra`).
Checkpoint: present and the forward signature was exercised on the real model during research.

## Problem Statement

`positional-non-mean-pooled-promoteratlas-150bp` was declared in `soraya_runs/arms-ledger.json` on
2026-11-01 under R6, at the moment its sibling arm ran, and it is now the ONLY `planned` arm on the
tile-chain ceiling criterion (`ceiling_verdict` = `blocked:untried-arm`, 11 measured / 1 planned).
It has no plan. This is that plan.

The measured sibling takes PromoterAtlas's final latent and **masked-mean-pools it over positions**
to 128 features, reaching median Spearman 0.5076 — the best arm on the rung. That arm's own recorded
limit says a mean-pool is the weakest reasonable use of the encoder, because averaging over positions
discards where in the window a signal sits. The question is whether a position-resolved reduction of
the SAME frozen latent beats it by the rung's inherited margin.

**Two scoping corrections the research forced, both of which narrow the arm from how it was
declared.** I wrote that rationale an hour ago and it over-claims on one point.

1. **This arm CANNOT test the −35/−10 spacing claim its rationale asserts.** `Tile` exposes
   `seq, expression, peak, strand, category, active, rep1, rep2` — there is no TSS, no promoter
   register, and no peak-relative offset. Position is therefore relative to an arbitrary 150 bp
   window, not to a −35/−10 frame. What the arm can test is narrower and still worth testing:
   *does a position-resolved reduction beat a position-averaged one.*
2. **Because of (1), the readout that actually matches the stated mechanism is NOT position bins.**
   Spacing is translation-invariant, so an anchor-free reduction — autocorrelation of each latent
   channel across lags — targets spacing directly where bins cannot. It is therefore a declared
   CO-PRIMARY testing a DIFFERENT hypothesis, not a second shot at the same one.

**Non-goals.** No fine-tuning (the clone ships four expression-prediction heads; using them is a
different arm about a fine-tuned encoder, not this one). No new substrate. No change to the frozen
AMR surface, to `genomewide.py`, or to `expression.py`. No attempt to re-read the measured sibling's
0.5076, which stands.

## Codebase Context


See: ./codebase-context.md

## Pre-Change Baseline

Measured, all from committed artifacts or live reads during research:

- Incumbent to beat: atlas mean-pool **median Spearman 0.5076** over 10 seeds, per-seed values above.
- Dumb baseline on the same splits: GC **median 0.3112**.
- Rung ceiling: **0.941**, replicate-derived. Never read a tile-rung number as a fraction of it
  without saying so.
- Inherited margin: **0.02**, taken verbatim from the tile-rung arms so bars stay commensurable.
- `ceiling_verdict`: `blocked:untried-arm`, 11 measured / 1 planned, `best` =
  `promoteratlas-base-embeddings-frozen-150bp-padded-to-200` at 0.5076.
- Suite: **13 failed / 7256 passed / 180 skipped** as of 33b24322. The 13 are pre-existing and are
  absent-data or a known registry-count drift; this plan must not change that set.

## Verification Signal

The plan has succeeded when ALL of these hold, independent of which way the arm's verdict goes:

1. The parity gate passes: the global mean derived from stored equal-size bin means reproduces the
   committed `tile_embeddings_native150_n44106.npy` within float tolerance, AND re-scoring it
   reproduces the committed atlas per-seed values.
2. Split equality holds against `data/promcalc/tile_predictions/` for all ten seeds
   (`n_test_per_seed` identical).
3. The shuffled-label null passes its inherited bar, or the run reports
   `INDETERMINATE_NULL_FAIL` and spends MORE SEEDS rather than re-specifying the bar.
4. A verdict artifact exists with a verdict from the pre-frozen closed set, and the set of reachable
   outcomes EQUALS the check order (not merely is contained in it).
5. The arm moves `planned -> measured` in the arms ledger with its evidence string, and
   `ceiling_verdict` is recomputed and reported — including if it now reads `met`.
6. The suite's failure set is unchanged, diffed by name rather than counted.
7. The frozen five, `uv.lock`, `genomewide.py` and `expression.py` are byte-identical.

## Implementation Steps

### Step 1: Pure positional-readout module
Files: dna_decode/glm/positional_readout.py, tests/test_positional_readout.py
Depends on: none

**What changes:**
- New module of PURE reductions over a `(B, C, L)` array, with no model import and no I/O, so every
  reduction is testable on synthetic tensors whose answer is known by construction.
- `bin_means(z, n_bins)` and `bin_maxes(z, n_bins)` — REFUSE when `L % n_bins != 0` rather than
  truncating or padding a ragged final bin, because a silent ragged bin makes the global mean
  non-recoverable and would break the parity gate quietly.
- `global_mean_from_bin_means(bm)` — the exact inverse for equal bins; this is what the parity gate
  calls, and it exists so the gate tests the STORED features rather than a re-run.
- `channel_autocorrelation(z, lags)` — per-channel Pearson autocorrelation across positions at each
  lag, mean-centred per channel, anchor-free. Returns `(B, C * len(lags))`.
- `orient(features_by_bin, strands)` — reverses the BIN axis for `-` strand rows. Operates on stored
  bin features, so orientation costs no second forward pass.
- Every function refuses a non-finite input rather than propagating NaN, because `spearman` maps NaN
  to 0.0 and would render a broken reduction as a clean null.

**Test strategy:**
- Construct a tensor whose per-bin means are known integers and assert `bin_means` returns them
  exactly; assert `global_mean_from_bin_means` reproduces `z.mean(axis=2)` to float tolerance.
- Assert the ragged refusal fires for `n_bins` not dividing `L` (e.g. 7 into 150).
- Autocorrelation: a pure sinusoid of period p must peak at lag p and trough at lag p/2; white noise
  must give near-zero at all non-zero lags.
- `orient`: applying it twice is the identity; applying it to an all-`+` cohort is a no-op.
- A NaN anywhere in the input raises for every entry point.

### Step 2: Feasibility and parity probe on the real model
Files: scripts/positional_arm_probe.py, wiki/positional_arm_feasibility_2026-11-01.json
Depends on: Step 1

**What changes:**
- Validity-first probe that must pass BEFORE any threshold is frozen, mirroring how the sibling arms
  derived their bars from a feasibility run rather than from judgement.
- Resolves the clone, runs `audit_checkpoint`, and asserts the parameter count equals the module's
  `EXPECTED_PARAMS` (1,847,428) — the checkpoint lives in a volatile temp directory, so its presence
  is a live fact and not an inherited one.
- Asserts the forward contract: `model(X, return_latent=True)[2]` has shape `(B, 128, L)`.
- **Parity on a sample:** for ~512 tiles, extract bin means, derive the global mean, and compare
  against the committed embedding rows for those same tiles. Reports max absolute deviation.
- Measures throughput on that sample and extrapolates the full-pass wall clock.
- Measures free RAM and computes the peak design-matrix bytes for every ladder rung, emitting a
  per-rung `fits_in_ram` flag against a named headroom factor.
- Reports feature counts and the on-disk size of the stored block.
- Exits non-zero on any failure, so a broken probe cannot read as a clean feasibility.

**Test strategy:**
- Unit-test the pure helpers (RAM arithmetic, byte estimates, extrapolation) on synthetic inputs.
- The model-touching path is gated on the clone resolving and skips cleanly when it does not, so the
  file is green on a host without the checkpoint.
- Assert the probe's exit code is non-zero when a synthetic parity deviation exceeds tolerance.

### Step 3: Freeze the acceptance bar
Files: scripts/freeze_positional_arm_bar.py, wiki/positional_arm_acceptance_bar.json
Depends on: Step 2

**What changes:**
- Freezes the bar from the Step 2 feasibility figures, before any score exists.
- INHERITS rather than re-chooses: margin **0.02**; the comparator is the committed atlas per-seed
  values read from `data/promoteratlas/local_reconcile_2026-11-01.json`; the split is
  `leave_peak_out(frac=0.25)` over seeds 0-9; the scorer is `ridge_spearman` unchanged.
- `verify()` re-reads the comparator and the arms ledger and REFUSES to write if the atlas median,
  the GC median or the arm's `planned` status have drifted — a bar anchored to retyped numbers is
  anchored to nothing, and a bar frozen after a run is not a bar.
- **Two co-primaries, testing different hypotheses, with the multiplicity named and corrected:**
  `bins10_mean_oriented` (position-resolved) and `autocorr` (spacing-resolved, anchor-free). Each is
  judged at a Bonferroni-corrected two-primary level rather than at 95%, and the correction is
  written into the bar rather than applied afterwards.
- **Pre-declared secondaries,** judged at the 5-comparison corrected level and explicitly NOT
  eligible to be reported as the arm's result: `bins10_mean_asgiven` (the orientation diagnostic),
  `bins10_max_oriented`, `bins2_mean_oriented`, `bins5_mean_oriented`.
- Closed outcome set, with the reachable set required to EQUAL the check order:
  `POSITIONAL_ADDS` / `POSITIONAL_NO_GAIN` / `POSITIONAL_AMBIGUOUS` /
  `INDETERMINATE_PARITY_FAIL` / `INDETERMINATE_SPLIT_MISMATCH` / `INDETERMINATE_NULL_FAIL`.
- Pre-commitments: the margin will not be re-specified after a number exists; a failing null spends
  MORE SEEDS rather than moving the bar; the measured sibling's 0.5076 is never retro-edited.
- Honest limits fixed in advance, including both scoping corrections from the Problem Statement, the
  orientation convention being a 1-bit inference from a performance difference, and the 689 known
  duplicate sequences being unchanged by this arm.
- Records a forecast with its probabilities and the standing note that my forecast record in this
  project is 3 of 7 and should not be deferred to.

**Test strategy:**
- Assert `verify()` refuses on a mutated comparator, on a mutated arm status, and on a missing arm.
- Assert the reachable-outcome set equals the check order exactly, which is what catches a dead
  validity branch.
- Assert the margin in the written artifact equals the value in the arms-ledger-derived anchors, so
  an inherited number cannot silently become a re-chosen one.

### Step 4: Extraction runner — one forward pass, all rungs
Files: scripts/positional_arm_extract.py
Depends on: Step 1

**What changes:**
- One pass over the 44,106 tiles emitting, per tile, `bins10_mean` (1,280), `bins10_max` (1,280) and
  the autocorrelation block — roughly 5,120 float32 per tile, ~900 MB, written under
  `D:/dna_decode_cache/` per the standing D: rule rather than onto a tight C:.
- Writes a sidecar recording the checkpoint hash, parameter count, geometry, tile count, bin count,
  lag set, dtype and wall clock, so a later scoring run can assert it is reading the features it
  thinks it is.
- Checkpointed by batch so a crash costs one batch, not the pass — the project has a measured
  precedent of a mid-run crash leaving a partial cache that looked complete to its consumer.
- Refuses a ragged cohort, and refuses to overwrite an existing block whose sidecar geometry differs.
- `--limit` for a smoke run, which stamps the sidecar PARTIAL so a partial block can never be
  scored as if it were whole.

**Test strategy:**
- Run `--limit` on a few hundred tiles with a stub encoder and assert the stored block's shape,
  dtype and sidecar fields.
- Assert resuming after an interrupted run reproduces the same bytes as an uninterrupted run over
  the same tiles.
- Assert the PARTIAL stamp is present for a limited run and absent for a full one, and that the
  scoring runner's loader refuses a PARTIAL block.

### Step 5: Scoring runner — gates first, then the ladder
Files: scripts/positional_arm_score.py, wiki/positional_arm_result_2026-11-01.json
Depends on: Step 3, Step 4

**What changes:**
- Executes the frozen bar's check order, VALIDITY FIRST, so a broken run cannot read as a finding.
- Gate A, split equality: recompute `leave_peak_out(frac=0.25)` for seeds 0-9 and assert
  `n_test_per_seed` matches the committed values for all ten; mismatch -> `INDETERMINATE_SPLIT_MISMATCH`
  and stop.
- Gate B, parity: derive the global mean from the stored bin means, score it through `ridge_spearman`
  on those same splits, and assert it reproduces the committed atlas per-seed values; deviation ->
  `INDETERMINATE_PARITY_FAIL` and stop. This is what makes every later number comparable to 0.5076.
- Gate C, shuffled-label null on the primary rung, judged against the inherited bar; failure ->
  `INDETERMINATE_NULL_FAIL` and spend more seeds.
- Then the ladder: each rung scored on byte-identical splits, paired against the committed atlas
  per-seed values, with a paired bootstrap over seeds and the bar's corrected interval.
- Emits per-rung medians, per-seed values, seed-win counts, the paired CI, the RAM actually used, and
  the single verdict from the closed set — plus per-tile predictions for the primary rung written
  under `data/` and NOT into `wiki/`, because a large decimal dump in `wiki/` has previously polluted
  the contract-number audit's decoy pool and broken two of its guards.
- A rung whose design matrix fails the RAM pre-check is SKIPPED and recorded as skipped with its
  byte estimate, never silently dropped.

**Test strategy:**
- Drive the whole gate chain with synthetic features and a stub scorer, asserting each
  `INDETERMINATE_*` branch is reachable and returns the right verdict and exit code.
- Assert a deliberately perturbed stored block trips the parity gate, which proves the gate is not
  vacuous.
- Assert a rung over the RAM budget is recorded as skipped rather than attempted.
- Assert the verdict is a pure function of the measurements plus the frozen anchors, by recomputing
  it from the emitted artifact alone.

### Step 6: Record the arm, render the memo, journal it
Files: soraya_runs/arms-ledger.json, scripts/positional_arm_memo.py, wiki/positional_arm_result_2026-11-01.md, project_state/glm-learned-representation-2026-10-07.md
Depends on: Step 5

**What changes:**
- Moves the arm `planned -> measured` with its metric and an evidence string naming the verdict, the
  gates, the ladder and the limits. The ledger write is a direct dataclass-safe edit because
  `measure_arm` raises on arms carrying off-schema annotation keys, and this arm has them.
- Recomputes `ceiling_verdict` and reports it. If it now reads `met`, that is a SUBSTANTIVE change
  to the criterion and the run stops to surface it rather than recording it in passing.
- Renders the memo FROM the result artifact so no figure is hand-typed, ordered by required reader
  action, and carrying both scoping corrections.
- Appends one Action-Log row with a class from the 7-value enum, anchored inside the Action Log
  section so a candidate row cannot be picked up as the numbering anchor.
- Marks G-A candidate 8 done with its journal marker.
- Re-pins the arms guards if their subject moved, by the documented precedent: re-pin when the
  verdict did not flip, and STOP and escalate if it did.

**Test strategy:**
- Assert the ledger round-trips with the off-schema annotation keys preserved, since losing them is
  the failure mode the direct write exists to avoid.
- Assert the memo generator refuses to render on an `INDETERMINATE_*` verdict, so a withheld verdict
  cannot be dressed as a finding.
- Re-run the arms guards and assert the partition reconciles with the whole.

## Execution Preview

- **Wave 0:** Step 1
- **Wave 1:** Step 2, Step 4 (parallelism 2)
- **Wave 2:** Step 3
- **Wave 3:** Step 5
- **Wave 4:** Step 6

Total waves 5, max parallelism 2. Critical path 1 -> 2 -> 3 -> 5 -> 6.

Wall-clock is dominated by Step 4's single forward pass (~17.5 minutes uncontended at the measured
42 sequences/second, materially longer while the K562 replication holds CPU) and by Step 5's 10
seeds x up to 6 rungs of ridge fits. Step 4 may legitimately run before the bar is frozen because it
produces FEATURES, not scores; no score exists until Step 5.

**Sequential execution is recommended over parallel worktrees.** The stored feature block is
gitignored and lives on D:, so a fresh worktree could not run Step 5 at all and would silently skip
the parity gate — the same reasoning that made the preceding plan sequential.

## Architecture Impact

`sentrux` is present, so the structural gate is eligible. The change adds ONE new leaf module
(`glm/positional_readout.py`) with no inbound dependencies and four scripts, so the import-edge count
moves by construction and an aggregate compare cannot distinguish an intended new edge from an
unintended one. Replacement check, mirroring what the sibling plan adopted after the same problem:
assert that the set of files under `dna_decode/` importing `positional_readout` is exactly the
intended consumer set, and that `genomewide` gains no new in-package importer beyond the three that
already exist.

## Risk Flags

- **The checkpoint lives in a volatile temp directory.** 686 MB under
  `C:/Users/LAPTOP~1/AppData/Local/Temp/pa_probe/PromoterAtlas`, which a reboot or a cleaner can
  remove. Step 2 treats its presence as a live fact and fails loudly; a re-clone is the recovery.
  Mitigation worth taking in Step 2: record the checkpoint hash so a replacement can be verified
  identical rather than assumed.
- **RAM is the binding constraint and it moved during research** — 5.08 GB free with K562 running
  against 10.5 GB recorded in October. The largest rung is a 2,560-feature matrix; the RAM pre-check
  and the skip-and-record behaviour exist because of this, and running Step 5 after the replication
  finishes is the cheap mitigation.
- **The orientation convention of `seq` is unresolved** and is inferred from a performance
  difference between two rungs. That is a 1-bit inference and is flagged as such in the bar; it is
  NOT established biology.
- **Multiplicity is real.** Six rungs against one incumbent is a garden of forking paths if the best
  is reported. Two co-primaries with a named correction and secondaries that are explicitly
  ineligible as the arm's result are the mitigation; the risk is that a reader quotes a secondary
  anyway, which the memo must pre-empt in its own text.
- **A NO_GAIN here is weaker than it looks.** Without a promoter anchor, a null says position bins
  relative to an arbitrary window add nothing — not that positional structure is absent. The bar
  records this so the negative cannot later be over-read as closing positional readouts generally.
- **689 duplicate sequences** are known on this substrate and unchanged by this arm; a
  train/test-straddling duplicate affects this arm exactly as it affects the incumbent, so the
  paired comparison is the right shape and the absolute levels inherit the issue.
- **Clearing the last planned arm may flip `ceiling_verdict` to `met`.** That is a substantive change
  to a criterion, not bookkeeping, which is why Step 6 stops and surfaces it instead of recording it
  quietly.

## Open Questions

- Whether a `met` ceiling verdict should be ACCEPTED as the criterion's terminal, or whether a third
  arm should be declared first. That is an authority call about scope, not a technical one, and Step 6
  surfaces it rather than deciding it.
- Whether the autocorrelation lag set should be dense (1..20) or log-spaced. Deferred to Step 2,
  which measures the feature count and RAM for each and freezes the choice from the measurement
  rather than from preference.
- Whether the four shipped expression-prediction heads in the clone deserve their own arm. Named here
  so it is not lost, explicitly out of scope, and it would be a NEW arm under R6 rather than a
  variant of this one.

## Verification

- `uv run pytest tests/test_positional_readout.py -q` plus the per-step test files, green.
- Full suite diffed by failure NAME against the 13 pre-existing failures at 33b24322, not counted.
- `uv run python -m scripts.positional_arm_probe` exits 0 and its artifact reports parity within
  tolerance.
- `uv run python scripts/freeze_positional_arm_bar.py` exits 0; re-running it on a mutated comparator
  exits non-zero.
- `uv run python -m scripts.positional_arm_score` emits a verdict from the closed set, and the verdict
  is re-derivable from the emitted artifact alone.
- `git diff --stat` shows the frozen five, `uv.lock`, `genomewide.py` and `expression.py` untouched,
  verified by hash rather than by `git status`.
- The arms-ledger diff is checked SEMANTICALLY — arms added, arms removed, and which fields of which
  arm changed — because a `sort_keys` rewrite produces a large line diff that hides what moved.
- `ceiling_verdict` recomputed and reported, with `best` named.

<!-- toolkit: check=clean waves=clean gate=fired:open-questions -->
