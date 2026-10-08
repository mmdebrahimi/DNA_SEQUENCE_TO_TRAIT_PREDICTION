## Lens status
Repo-index: absent (`unconfigured` — update-index CLI not installed) → source-only reading, no degradation.
Baseline tests: `uv run pytest tests/test_glm_*.py -q` → **69 passing** across the four GLM data modules.
sentrux: **absent** → Architecture Gate n/a. codex: present. gh: present. graphify: absent.
Wiki replay: `wiki/decisions-log.md` scanned (32 HIGH-salience entries); the relevant prior is the v0 falsifier
FAIL that still shipped an honest scope-limited artifact with a leakage-safe retrain — consistent with
decision D1 here, not in conflict with it.
Mode: parallel (4 waves, max parallelism 4, a git remote is configured).

## Problem Statement

Build G-A (learned representation) and G-C (condition conditioning) — the two GLM families cleared to start —
with the revisions forced by an adversarial review and three now-settled authority decisions. The v1 plans are
**not** executable as written: the review found that G-C's model class cannot test its own hypothesis, that
G-A's anti-search-artifact protections existed only in prose, and that the MVP bar is satisfiable by a broken
run.

Three fixes are therefore load-bearing rather than cosmetic:

1. **G-C's primary arm must carry `sequence × medium` interactions.** A medium *indicator* appended to a
   linear ridge learns only a global per-medium offset; it cannot modulate sequence weights. The two-model
   comparator fits separate weight vectors per medium and is therefore **strictly more expressive on exactly
   the axis conditioning is meant to exploit** — so the v1 design was rigged against itself and could
   essentially only return `GAIN_IS_DATA_VOLUME` or `NO_GAIN`.
2. **G-A must freeze its training protocol and go multi-seed.** The GC baseline 0.3000 was fixed on seed 0's
   split; evaluating a tunable model on that same split, against a 0.05 margin, is a search surface.
3. **Both MVP bars must require the experiment to have been VALID** (decision D1) — the existing predicates
   (module exists / tests pass / artifact exists / ledger row) are all satisfied by a broken run that writes
   a null-result artifact.

Compute is no longer CPU-bound (free Kaggle GPU available, paid acceptable where it adds value), and per
decision D2 that budget is spent on **rigour, not search**: more seeds and full data, never architecture
hunting.

## Codebase Context


See: ./codebase-context.md

## Pre-Change Baseline

**G-A, on genomic peak tiles** (`wiki/glm_genomewide_oracle_2026-10-07.json`; leave-peak-out, **seed 0 only**,
n_train 32,908 / n_test 11,198, controls and non-150bp dropped):
GC **0.3000** (0.3188 of ceiling) · 4-mer **0.2769** · positional one-hot **0.2691** · 6-mer **0.2022**.
Tile noise ceiling **0.9410** (single-replicate 0.7945). No learned encoder exists.

**G-A grid comparator** (`wiki/glm_expression_oracle_result_2026-10-07.md`; leave-element-out):
GC 0.1414 · 4-mer 0.4137 · 6-mer 0.5872 · one-hot **0.5914**.

**G-C: no per-medium baseline exists yet** — this is the gap Step 4 fills. The only fragment numbers are
LB-only on a 30,000 subsample: GC 0.1023 · 4-mer 0.1313 · 6-mer 0.0919. Fragment ceilings **LB 0.7928 /
M9 0.8014**. Condition effect: noise-matched gap **−0.0614** on 297,868 coordinate-keyed shared fragments,
disattenuated true cross-condition **0.7887**.

**Reference geometry** (`wiki/glm_reference_geometry_2026-10-07.json`, 10 seeds): grid vs tile **1.0000**,
grid vs natural **1.0000**, tile vs natural 0.6664 / **0.8338**; all self-controls at chance (0.494–0.501).

**Tests:** 69 passing across `test_glm_{genomewide,reference_geometry,selection,expression}.py`.

## Verification Signal

**G-A primary:** median encoder Spearman across **≥10 leave-peak-out seeds** minus median GC across **the same
seeds**, required ≥ **0.05**, reported as a fraction of the **0.9410** tile ceiling, with seed dispersion.
**G-A validity (blocking):** the label-shuffled arm's margin over GC must stay below 0.05 across the same
seeds. If it does not, the gate returns `INDETERMINATE_NULL_NOT_CLEAN` and nothing is graded.
**G-A secondary:** the same encoder on the designed grid under leave-element-out — a **consistency check
only**, never evidence of transferable grammar (see Risk Flags for why).

**G-C primary:** Spearman of the interaction model minus the pooled two-model baseline on the **same**
position-blocked split, required ≥ **0.02**, each medium against **its own** ceiling.
**G-C validity (blocking):** the medium-shuffled arm must show no gain, AND the data-matched arm must also
clear the margin. If the shared arm beats the two-model baseline but the data-matched arm does not, the
verdict is `GAIN_IS_DATA_VOLUME`, not a win.

## Implementation Steps

### Step 1: Upgrade both MVP bars with the D1 validity predicates
Files: project_state/glm-learned-representation-2026-10-07.md, project_state/glm-condition-conditioning-2026-10-07.md
Depends on: none

**What changes:**
- Add two predicates to each family's `## MVP Criteria`: `null_clean` (the pre-registered shuffled/permuted
  control arm came out flat) and `registered_protocol_used` (the frozen seed list and training protocol were
  used, not a post-hoc variant).
- Record decision D1 in each ledger's `## Decisions Made` with the reason: the four original predicates are
  all satisfied by a broken run that writes a null-result artifact.
- The bar is frozen BEFORE any run, which is the whole point of sequencing this first.

**Test strategy:**
- A schema check that each ledger still has all 8 required sections and all 6 end markers exactly once
  (reuse the validator already written for these ledgers).
- Assert both new predicate names appear in both ledgers.

### Step 2: Feature-parity drift guard for the two k-mer implementations
Files: tests/test_glm_feature_parity.py
Depends on: none

**What changes:**
- Pin that `falsifier.kmer_features`, `genomewide.kmer_features` and `expression.kmer_feature` return
  identical vectors for the same input and k, across k ∈ {1, 4, 6} and several sequence lengths.
- This is not hypothetical bookkeeping: G-A's gate arms and the reference-geometry triage draw k-mer features
  from *different* modules, and their current agreement is coincidental rather than enforced.

**Test strategy:**
- The test IS the deliverable. Include a non-vacuity assertion that the implementations are genuinely separate
  objects (not the same function re-exported), so the guard cannot pass trivially.

### Step 3: G-A linear-headroom diagnostic (a real stop gate)
Files: scripts/glm_tile_headroom.py, tests/test_glm_tile_headroom.py
Depends on: none

**What changes:**
- Fit ridge on GC, GC+dinucleotide (16), GC+positional-GC-in-bins, and dinucleotide+position — all under
  `leave_peak_out` across the registered seed list, on the full 44,106-tile set.
- Emit `wiki/glm_tile_headroom_<date>.json` with each arm against the 0.9410 ceiling.
- **This step can close G-A in an afternoon.** If composition-only features already reach ~0.30 and adding
  position adds nothing, hypothesis GA-H3 is confirmed and a conv encoder is a predicted-negative.
- Emit an explicit `recommend_user_ratification_before_cnn` boolean so the stop gate is a machine-readable
  output rather than a judgement call buried in prose.

**Test strategy:**
- Unit-test each feature builder for width and length-invariance.
- **Reconcile gate:** assert the GC arm reproduces the committed 0.3000 ± 0.01 on seed 0. If it does not, the
  split or the data differs and the run must refuse rather than report.

### Step 4: G-C shared split, per-medium baselines and per-arm ceilings
Files: scripts/glm_condition_conditioning.py, tests/test_glm_condition.py
Depends on: none

**What changes:**
- Load LB and M9, intersect on the **coordinate key** (the fix already landed in `condition_effect`), and
  build **one** position-blocked split per seed, reused identically by every arm — so arms differ only in the
  model, never in the partition.
- Compute LB and M9 `noise_ceiling` **separately** and pin both; no shared normaliser.
- Fit the two independent per-medium ridges and emit the pooled two-model baseline that the conditioned model
  must beat. Use the **full 297,868 shared fragments**, not a subsample (decision D2).

**Test strategy:**
- Assert every arm receives the identical held-block id set for a given seed.
- Assert each medium's ceiling is derived from that medium's own replicates.
- Refuse if the shared-fragment count falls below a floor, so a degraded intersection cannot be scored.

### Step 5: G-C interaction model plus the three control arms
Files: scripts/glm_condition_conditioning.py, tests/test_glm_condition.py
Depends on: Step 4

**What changes:**
- **`interaction` arm (the C1 fix, now the PRIMARY):** features `[sequence ‖ medium ‖ sequence × medium]`, so
  the model can represent condition-specific sequence response rather than a global offset.
- `indicator` arm: `[sequence ‖ medium]` — **demoted to a diagnostic**, retained precisely to show that an
  indicator alone cannot win.
- `shuffled` arm: the interaction model with the medium label permuted. **Must show no gain.**
- `data_matched` arm: the interaction model trained on a row-count-matched subsample, so a gain cannot be
  attributed to seeing ~2× the rows.
- Freeze the four-branch verdict in the file before any run: `CONDITIONING_HELPS` /
  `GAIN_IS_DATA_VOLUME` / `NO_GAIN` / `INDETERMINATE_NULL_NOT_CLEAN`.

**Test strategy:**
- The verdict function is pure — unit-test all four branches on synthetic arm dicts, including
  `GAIN_IS_DATA_VOLUME`, which exists because it is the likely outcome.
- A planted-interaction synthetic test: construct data where a motif's effect **flips sign** by medium and
  assert the interaction arm beats the indicator arm. Without this the C1 fix is unverified.

### Step 6: G-A encoder module with lazy torch and pinned determinism
Files: dna_decode/glm/encoder.py, tests/test_glm_encoder.py
Depends on: Step 3

**What changes:**
- `SequenceEncoder`: one-hot `(N, 4, L)` → parallel `Conv1d` branches at **kernel widths 6 / 20 / 75** (motif
  / spacer / whole-promoter scale) → ReLU → global max+mean pool per branch → concat → linear head. Small by
  design (~32 channels per branch). Widths produce different output lengths (L − k + 1), so **pooling before
  concatenation is structural, not stylistic**.
- `fit_encoder(train, test, *, widths, protocol, seed) -> dict` mirroring `_fit`'s return shape. `protocol` is
  a frozen dataclass (optimizer, LR, batch size, epoch cap, early-stop rule) — **not loose kwargs**, so a
  post-hoc tweak is visible in a diff.
- **torch is lazy-imported inside `fit_encoder`**, mirroring `generate.py`, so the module's pure helpers and
  their tests stay runnable on a default install where torch is absent.
- Device: CPU is the **determinism reference**; `device="auto"` may use the local GPU or a Kaggle T4 for
  speed, and the artifact records which ran. The determinism claim is scoped to CPU — CUDA convolutions are
  documented as potentially nondeterministic.

**Test strategy:**
- Planted-signal test: a motif at a fixed offset must be learnable to Spearman > 0.8.
- Shuffled-label test: must collapse to ≈ 0.
- CPU determinism: same seed → identical predictions; pinned with `use_deterministic_algorithms`.
- Ragged-input refusal (no silent padding), empty-train refusal.
- `importorskip("torch")` on the training tests only; the pure helpers must pass without torch.

### Step 7: G-A frozen-protocol multi-seed gate with a re-derived baseline
Files: scripts/glm_encoder_gate.py, tests/test_glm_encoder_gate.py
Depends on: Step 6

**What changes:**
- Arms: `gc`, `encoder-w20` (single resolution), `encoder-multi` (6/20/75), `encoder-shuffled` (the null).
- **Registered seed list** (≥10) used by every arm; verdict on the **median**; seed dispersion reported.
- **Re-derive the GC baseline over the same seed list** and report `seed0_gc=0.3000` alongside
  `median_gc_over_registered_seeds` for continuity — the committed single-seed artifact is left untouched and
  pointed at, rather than amended.
- Freeze the verdict branches before any run: `GENERALISES` (margin ≥ 0.05) / `WEAK` (0 < margin < 0.05) /
  `DOES_NOT_GENERALISE` (margin ≤ 0) / `INDETERMINATE_NULL_NOT_CLEAN`.
- Record `registered_protocol_used` and `null_clean` into the artifact so Step 9's MVP predicates read a
  machine-readable field rather than a human claim.
- Runs on tiles (primary) and the grid (consistency check only).

**Test strategy:**
- The verdict function is pure — unit-test all four branches, and assert the null-not-clean branch **refuses
  to grade** rather than returning a band.
- Assert the gate fails closed if the seed list used differs from the registered one.

### Step 8: Execute G-C, commit the artifact, record the verdict
Files: wiki/glm_condition_conditioning_2026-10-08.json, project_state/glm-condition-conditioning-2026-10-07.md
Depends on: Step 1, Step 5

**What changes:**
- Run the four arms across the registered seeds on the full shared-fragment set; commit the artifact whatever
  the verdict.
- Append an Action Log row naming the verdict; update hypotheses GC-H1..GC-H4.
- Carry the honest limits into the artifact: two conditions only; low-ceiling fragment substrate; conditioning
  tested as interactions over k-mer features, not as representation modulation.

**Test strategy:**
- The four MVP predicates (including the two new validity ones) are the test, re-evaluated live.

### Step 9: Execute G-A, commit the artifact, record the verdict
Files: wiki/glm_learned_representation_2026-10-08.json, project_state/glm-learned-representation-2026-10-07.md
Depends on: Step 1, Step 7

**What changes:**
- Run the gate; commit the artifact whatever the verdict.
- Append an Action Log row naming the verdict; resolve GA-H1/GA-H2 and confirm-or-falsify GA-H3.
- If `DOES_NOT_GENERALISE`: write the closed-negative memo and record whether the failure is a **capacity**
  ceiling (which closes G-B) or a merely weak encoder (which only defers it) — that distinction is already a
  pending decision in the G-B ledger and this run is what resolves it.

**Test strategy:**
- The four MVP predicates re-evaluated live (two of them being the new validity predicates), so a run that
  did not use the registered protocol or whose null arm was dirty cannot satisfy the bar.
- Full-suite regression is run from the `## Verification` section rather than named here, so the step's
  contract stays a concrete per-step check.

## Execution Preview

Waves: **4**. Max parallelism **4**.
- **Wave 0:** Step 1 ∥ Step 2 ∥ Step 3 ∥ Step 4
- **Wave 1:** Step 5 (← 4) ∥ Step 6 (← 3)
- **Wave 2:** Step 7 (← 6) ∥ Step 8 (← 1, 5)
- **Wave 3:** Step 9 (← 1, 7)

Critical path: **Step 3 → 6 → 7 → 9** (four steps). G-C completes one wave earlier than G-A, which is the
intended hedge: the family with an already-positive motivating measurement reports first.
Gate classification: all steps `auto` (in-cwd writes/edits, local tests). No money, no `destructive-local`,
no outward-irreversible actions. A Kaggle escalation would remain token/free-tier and is not required by any
step.

## Risk Flags

- **Step 3 may close G-A before Steps 6-9 run.** That is a designed outcome, not a risk to mitigate — but it
  means Steps 6-9 should not begin before Step 3 reports, which is why Step 6 depends on Step 3 rather than
  running in Wave 0 beside it.
- **Framing risk, larger than the architecture risk:** the tile set is 15,171 active / 31,542 inactive, so
  Spearman over the whole set may be substantially ranking noise among inactive tiles. Step 7 additionally
  reports AUROC on the `active` column as a diagnostic. If regression is the wrong frame, that is a finding.
- **The grid arm must not be read as evidence of transferable grammar.** Newly measured: the designed grid is
  separable from genomic sequence at AUROC **1.0000**, and the mechanism is recombination from a small element
  vocabulary (median 13 bp identical block between any two grid sequences vs 3 bp for tiles; per-position
  entropy 1.513 vs 1.987 bits; **no** constant scaffold — that hypothesis was refuted). With ~3×8×8×8×8
  element variants, a model scoring well under leave-element-out can be exploiting the four axes it *did* see.
- **Capacity has already hurt once on this substrate:** 6-mers (4,096 features) scored WORST on tiles
  (0.2022). The encoder is deliberately small, and the shuffled-label null must stay clean.
- **Two k-mer implementations agree only coincidentally** (Step 2 guards it). Until that guard lands, any
  cross-module feature comparison is unverified.
- **torch absent from default deps** — mitigated by the lazy-import convention in Step 6; without it a default
  install would fail at import on the new module.
- **CUDA nondeterminism and tensor layout:** `Conv1d` requires `(N, C, L)` (verified) and CUDA convolutions
  are documented as potentially nondeterministic. CPU is the determinism reference; the device used is
  recorded in the artifact.
- **Local GPU is CC 5.0 / 4 GB** — verified to run `Conv1d` but with no tensor cores. Kaggle T4 is the
  escalation; creds are present but the `kaggle` CLI is **not** on PATH, so the existing `scripts/kaggle/`
  kernel pattern must be followed rather than a CLI assumed.
- **sentrux absent** → no architecture gate for this plan.

## Open Questions

- Whether a G-A negative also closes G-B or merely defers it. Already a pending decision in the G-B ledger;
  Step 9 is what resolves it, and the plan records the distinction rather than pre-judging it.
- Whether the active/inactive classification framing should become a full G-A arm rather than a diagnostic.
  Currently diagnostic; promoting it would change what G-A's bar measures, so it is an acceptance-bar question.
- Whether `GAIN_IS_DATA_VOLUME` on G-C closes conditioning for this substrate or triggers a multi-task
  formulation. Two conditions is the minimum, so a null is weak evidence against conditioning in general.

## Verification

1. `uv run pytest tests/test_glm_feature_parity.py tests/test_glm_tile_headroom.py tests/test_glm_condition.py tests/test_glm_encoder.py tests/test_glm_encoder_gate.py -q`
2. `uv run pytest tests/test_glm_*.py -q` — must not regress the 69 passing today.
3. Live gate runs for both families, then both MVP bars re-evaluated including the two new validity predicates.
4. Frozen-five byte-check and prospective-lock re-verification after any commit (`amr_rules.py`,
   `calibrated_amr_rules.json`, `mic_tiers.py`, `shipped_decoder_surface.py`, `cohort_manifest.py`).

## Save-time amendments

Captured at: 2026-10-08
Source: `/save-plan` arguments

**Audit-notes-only contract:** this block is provenance for human readers. `/execute-plan` reads ONLY
`## Implementation Steps` for executable work, so amendments here are NOT executable instructions. If an
amendment changes step contracts (file lists, dependencies, structure), re-run `/technical-plan` before
`/execute-plan` rather than relying on this section.

- G-A + G-C v2: brainstorm fixes C1/C2/C3 + authority decisions D1/D2/D3 folded in; G-E triage result reframes the grid arm as a consistency check only

<!-- toolkit: check=clean waves=clean gate=fired:open-questions -->
