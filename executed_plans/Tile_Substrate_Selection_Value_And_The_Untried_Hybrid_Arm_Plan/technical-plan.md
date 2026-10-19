# Technical plan — Is the promoter PREDICTION criterion still worth pushing? (tile-substrate selection value + the untried hybrid arm)

## Lens status

Dispatch: conversation-driven (no argument).
Repo-index: absent (`update-index` CLI not installed) — source-only reading, no degradation.
sentrux: present (`~/.local/bin/sentrux` shim → `D:/tools/sentrux/sentrux.exe` v0.5.7). **No `.sentrux/rules.toml` exists in this repo**, so the architecture gate has nothing to enforce; `/execute-plan`'s gate will skip.
`project-rules.md`: absent → `## Project Rules Check` omitted.
`DESIGN.md`: absent; plan touches no UI.
Wiki replay: `wiki/plans-index.md` + `wiki/decisions-log.md` read. One **direct conflict** with a prior recorded decision is surfaced in Risk Flags and Open Questions rather than silently overridden.

## Problem Statement

The GLM portfolio's live `ceiling:` criterion reads *"best achievable **prediction** of bacterial promoter activity from sequence, reported against the replicate-derived 0.941 ceiling."* Nine arms are declared against it; five are measured. The project's own 2026-10-07 measurement (`wiki/glm_oracle_selection_2026-10-07.json`) concluded the opposite priority — verdict `ORACLE_READY_AS_SCORING_HEAD`, 5/5 axes, with the learned scorer realising a cross-axis mean **0.7937** of achievable selection value at k=10 against GC's **0.1328** — and `wiki/alphagenome_lessons_actionable_2026-10-07.md` recorded the consequence as *"the scoring half is adequate and the GENERATOR is the bottleneck."*

That measurement was taken on a **different substrate** (the GSE108535 designed σ70 grid, leave-element-out, learned scorer ρ≈0.59). The arms since declared all run on the **GSE144621 peak-tile** rung, where the best measured scorer is the 346-parameter Promoter Calculator at ρ **0.3968** against a 0.941 ceiling. **Nobody has measured what 0.3968 is worth to a choice on that rung.** Until that number exists, "keep pushing prediction" and "the criterion is the wrong target" are both unevidenced positions.

Two things are therefore missing, and both are cheap because the instruments already exist:

1. **The decision instrument.** Re-run the selection-value simulation on the tile rung with the scorers that actually ship there. If 0.3968 already captures most of the achievable selection value, further prediction work has measured diminishing returns for the generative goal; if it does not, prediction work is justified and we will know why.
2. **The one untried combination.** The Promoter Calculator was measured *against* the conv encoder (+0.0401, 4/4 seeds) and the 55-feature decomposition was measured against the scalar (NO_GAIN, CI includes zero) — but the learned encoder and the biophysical scalar have **never been combined**. The repo's own measured rule is that a naive rank-average lifts only when the modalities are **orthogonal**, and that a modality losing on its own can still lift (GEMME loses to ESM2 and still lifts it on 84% of proteins). That makes this arm worth one run and its precondition worth measuring first.

Non-goals, stated so scope cannot drift:

- **No generator work.** `features/glm-generator-finetune/` exists at status `intake` with a technical plan that is **DEFERRED BY DESIGN** and carries a named unblock trigger (*"the sampling sweep artifact exists AND the 8 GB host is reachable"*). Proposing generator implementation here would be this repo's recorded orientation error #3.
- **No criterion rewrite.** Changing the `ceiling:` wording is an acceptance-bar redefinition and therefore a user authority call. This plan produces the measurement that should inform it, not the edit.
- **No change to any published arm.** The three published encoder arms, the reconcile gate and the frozen AMR surface are untouched; every new file is additive.

## Codebase Context

**The selection instrument already exists and is pure.** `dna_decode/glm/selection.py` ships `simulate_selection(predicted, measured, k=, scorer=, trials=, seed=, percentiles=)` → `SelectionResult`, `permutation_null(measured, k=, trials=, n_perms=, seed=)`, `assert_k1_is_degenerate(result)` and `midrank_percentiles(values)`. `headroom_captured = (picked − random)/(ceiling − random)` IS the "fraction of achievable selection value" quantity; it is `None` at k=1 by construction (0/0), ties break uniformly because `random.sample` returns its subset in random order, and the module is I/O-free so every property is offline-testable. Tested at `tests/test_glm_selection.py`.

**Its existing runner is substrate-bound.** `scripts/glm_oracle_selection.py` wires that instrument to `dna_decode/glm/expression.py` (`load_pairs`, `leave_element_out`, `ELEMENT_COLUMNS`, `fit_predict_ridge`) — the designed-grid loader with a per-element axis split. The tile rung uses a different loader and a different split, so a new runner is required; the instrument itself is reused unchanged.

**The tile rung's API, read from source:** `dna_decode/glm/genomewide.py` exposes `load_peak_tiles(cache_dir=DEFAULT_CACHE, allow_download=True, ...)`, `leave_peak_out(tiles, frac=0.25, seed=0) -> (train, test, held)`, `noise_ceiling(frags) -> dict` (key `ceiling_for_a_sequence_model`) and `spearman(a, b)`. The published filter is `category == "tile" and len(seq) == 150` → 44,106 tiles.

**The join key is POSITION in that filtered list, and it is asserted not assumed.** `scripts/promcalc_decomposition_arm.py` builds `_pos = {id(t): i for i, t in enumerate(tiles)}` and asserts `not set(itr) & set(ite)`, with an in-line note that `leave_peak_out` returns the very objects it was given rather than copies — so a splitter that started copying would mis-join every feature row and still produce a plausible number. Stage A (`scripts/promcalc_export_tiles.py`) exported that same index.

**Targets are asymmetric and must be copied exactly:** training targets are `log1p(t.expression)`, test targets are raw `t.expression`.

**The four tile scorers, as the published arms compute them:**
- `gc` — `scripts/glm_encoder_gate.py::ridge_gc(train, test) -> (spearman, predictions)`. **Already returns predictions.**
- `max_tx` — ridge on the single feature `log1p(max(Tx_rate))` per tile, via `_fit_ridge`.
- `decomposition` — ridge on the 55 features from `tile_features(sites)` / `feature_names()`.
- `conv` — `dna_decode/glm/encoder.py::fit_encoder(train, test, widths=DEFAULT_WIDTHS, protocol=REGISTERED_PROTOCOL, seed=, device=, shuffle_labels=)`, whose returned dict contains `"spearman"` and `"predictions"`. `DEFAULT_WIDTHS = (6, 20, 75)`; `REGISTERED_PROTOCOL` must be used unchanged or the predictions are not the published arm's.

**Neither arm persists per-tile predictions — this is the finding that shapes the plan.** `scripts/glm_encoder_gate.py` computes a second framing from the same fit and then **deliberately discards the vector**: `r["auroc_active"] = auroc(r.pop("predictions"), lab)`, commented *"then DROP the raw predictions -- one float"*. `promcalc_decomposition_arm.py` likewise aggregates to a per-seed Spearman only. So both new measurements require a re-fit, and the economical shape is one shared step that persists an aligned table both consume.

**Substrate availability, verified on disk:** `D:/dna_decode_cache/promcalc/` holds `decomposition.csv.gz` plus resume shards `decomposition.part-01.csv.gz` and `part-02.csv.gz`, and `tiles.tsv`. `load_decomposition(in_dir, smoke=)` reads a GLOB of all of them and refuses coverage < 0.999 — a hardcoded single filename would silently score 88.4%.

**Reconcile anchors exist and are the proof a new consumer reads the same data:** `PREREGISTERED["reconcile"]` in `promcalc_decomposition_arm.py` pins `gc` seed 0 at **0.3000 ± 0.01** and `max_tx` at **0.3968 ± 0.02`, and the verdict function returns `INDETERMINATE_RECONCILE_FAILED` rather than a number when either drifts.

**The combiner also already ships.** `dna_decode/forward/variant_effect.py::rank_average_hybrid(tables: list[dict[str, float]]) -> dict[str, float]` mid-ranks each table over the key intersection, averages ranks, normalises to [0,1], and raises on <2 tables or an empty intersection. It has **no fitted parameters**, so combining two already-out-of-sample prediction vectors introduces no second-level leakage.

### Reusable-Code Survey

Searched: `graphify-out/GRAPH_REPORT.md` (absent), `src/utils`, `src/lib`, `src/common`, `src/helpers`, `utils`, `lib`, `common`, `helpers` (none exist), plus `dna_decode/glm/`, `dna_decode/forward/`, `scripts/glm_*`, `scripts/promcalc_*`, `tests/test_glm_*`.

Candidates adopted, all reused **unchanged**:

- `dna_decode/glm/selection.py` — `simulate_selection`, `permutation_null`, `assert_k1_is_degenerate`, `midrank_percentiles`. The entire selection measurement. No new statistics are written.
- `dna_decode/forward/variant_effect.py::rank_average_hybrid` — the shipped deployable form of the orthogonal-modality finding; exactly the combiner the hybrid arm needs.
- `dna_decode/glm/genomewide.py` — `load_peak_tiles`, `leave_peak_out`, `noise_ceiling`, `spearman`.
- `dna_decode/glm/encoder.py` — `fit_encoder`, `REGISTERED_PROTOCOL`, `DEFAULT_WIDTHS`, `count_parameters`.
- `scripts/glm_encoder_gate.py::ridge_gc` — the GC baseline, already prediction-returning.
- `scripts/promcalc_decomposition_arm.py` — `tile_features`, `feature_names`, `load_decomposition`, `_fit_ridge`, `paired`, `_agg`, and its `PREREGISTERED["reconcile"]` anchors.
- `scripts/glm_oracle_selection.py` — the **verdict pattern** (frozen `PREREGISTERED` dict, noise-aware bar B as `max(0.10, 2·null stdev)`, and the `INDETERMINATE_INSUFFICIENT_*` refusal when the bar is unreachable by construction). Mirrored, not imported, because its bar is per-axis and the tile rung's is per-seed.
- `soraya_runs/arms-ledger.json` via `~/.claude/skills/soraya/scripts/arms_ledger.py` — `ArmsLedger`, `declare_arm`, `measure_arm`, `ceiling_verdict`, `blocking`, and the **explicit `save()`**.

## Pre-Change Baseline

Measured in this session, not recalled:

- **Test suite:** 6,501 tests collected (`uv run pytest tests/ --collect-only -q`, 22.3 s).
- **Arms ledger** (`soraya_runs/arms-ledger.json`, read from disk in a separate process): 9 arms under one `ceiling:` criterion — 5 measured, 4 planned; `ceiling_verdict` = `blocked:untried-arm`.
- **Measured tile-rung Spearmans:** `gc-single-feature` 0.3112 (seed 0 = 0.3000) · `conv-3width-13k-params-150bp` 0.3234 (paired CI vs GC includes zero, 6/10 seeds) · `pyramid-4scale-1.70decades-150bp` 0.3167 · `promoter-calculator-biophysical-346params-fwd_max_tx` **0.3968** · `promcalc-full-decomposition-ridge-150bp` 0.3925 (`NO_GAIN`, paired −0.0040, CI95 [−0.0111, 0.0016], 2/10 seeds).
- **Tile rung:** 44,106 tiles at a fixed 150 bp; replicate-derived ceiling **0.941**. Best measured scorer sits at 42.2% of it.
- **Selection value on the tile rung: NEVER MEASURED.** The only selection-value artifact is `wiki/glm_oracle_selection_2026-10-07.json`, on the designed grid.
- **Conv ⊕ promcalc hybrid: NEVER MEASURED.** No arm in the ledger combines a learned and a biophysical scorer.
- **Per-tile predictions: NOT PERSISTED by any existing arm** (verified in both scripts).
- **Uncommitted tree state:** `scripts/glm_encoder_gate.py` carries +15 additive lines (a sibling session's `--extra-widths` flag). `ridge_gc` is **not** in that diff (verified by `git diff | grep`).

## Verification Signal

The plan has worked when all of these hold:

1. **The shared table reconciles.** Spearman recomputed *from the persisted prediction table* gives `gc` seed 0 within 0.3000 ± 0.01 and `max_tx` within 0.3968 ± 0.02. A drift means the new consumer is not reading the published substrate, and the step must refuse rather than report.
2. **Coverage is complete or refused.** `load_decomposition` coverage ≥ 0.999 over 44,106 tiles, or the run exits non-zero naming the missing count.
3. **The selection harness self-tests.** `assert_k1_is_degenerate` passes at k=1 for every scorer and seed, and `headroom_captured` is `None` there rather than a number.
4. **Every headroom is reported beside its own permutation null** (one null per `(seed, k)`), and no headroom point is quoted without it.
5. **The selection run emits one of the frozen verdicts**, or an explicit `INDETERMINATE_*` when the bar is unreachable — never a headline from a partial run.
6. **Orthogonality is reported before the hybrid verdict** — Spearman between the conv and `max_tx` prediction vectors, per seed, with the stated rule that rank-averaging a near-chance modality is measured to *reduce* signal.
7. **The hybrid arm's outcome reaches the ledger either way.** `measure_arm` is called with the measured number, `save()` is called, and the result is re-read **in a separate process**; a negative is recorded as a measured negative, never deleted.
8. **Suite green and grown:** ≥ 6,501 collected, 0 failures.
9. **Frozen surfaces untouched:** `dna_decode/eval/amr_rules.py` and `dna_decode/data/calibrated_amr_rules.json` byte-identical, and the active prospective lock still verifies.

## Implementation Steps

### Step 1: Declare the hybrid arm before it runs
Files: soraya_runs/arms-ledger.json, tests/test_glm_tile_arms_declared.py
Depends on: none

**What changes:**
- Declare one new arm on the existing `ceiling:` criterion via `ArmsLedger.declare_arm`: `promcalc-plus-conv-rank-average-hybrid-150bp`, `cost_tier="cheap"`, status `planned`. Declaring first is R6: if the hybrid loses, the arm records a measured negative instead of quietly never existing.
- **Call `save()` explicitly.** `declare_arm` and `measure_arm` mutate in memory only; this was measured on 2026-10-14 when three ledger writes reported success and the file's mtime never moved, because every verification called `ceiling_verdict()` on the same in-memory object that had just been mutated.
- Re-read the ledger **in a separate Python process** and confirm 10 arms (5 measured / 5 planned).
- Deliberately do **not** declare the selection-value run as an arm. An arm is a candidate approach to raising the ceiling metric; the selection measurement is a decision instrument about whether that metric matters. Filing it as an arm would make the ledger claim an improvement attempt that was never attempted.

**Test strategy:**
- New `tests/test_glm_tile_arms_declared.py`: load `soraya_runs/arms-ledger.json` **from disk** and assert the hybrid arm is present under the `ceiling:` criterion — a code-owned guard replacing the procedural "remember to save" rule.
- Assert the arm's status is one of `planned` / `measured` (so the test survives Step 5 measuring it) and that no arm was removed: the pre-existing 9 names are all still present.
- Assert `ceiling_verdict()` still reports a `blocked:*` state while any arm is `planned` — pinning that a cheap arm cannot cancel an untried one.

### Step 2: Persist one aligned per-tile prediction table for all four scorers
Files: scripts/glm_tile_predictions.py, tests/test_glm_tile_predictions.py
Depends on: none

**What changes:**
- New script `scripts/glm_tile_predictions.py`. For each registered seed (default 0..9): `load_peak_tiles()` → filter `category == "tile" and len(seq) == 150` → assert 44,106 → `leave_peak_out(tiles, seed=s)` → build the `{id(t): position}` join and assert train/test index disjointness, copying the published pattern verbatim.
- Compute per-test-tile predictions for all four scorers on that one split: `gc` via `ridge_gc`, `max_tx` via `_fit_ridge` on the single `log1p(max Tx_rate)` feature, `decomposition` via `_fit_ridge` on the 55 `tile_features`, `conv` via `fit_encoder(..., protocol=REGISTERED_PROTOCOL, widths=DEFAULT_WIDTHS, seed=s, device=...)`.
- Targets copied exactly: fit on `log1p(expression)`, score against raw `expression`.
- Persist per seed to `D:/dna_decode_cache/promcalc/tile_predictions/seed-NN.json.gz` (D: rule — derived cache, never C:): the test tile positions, the raw measured expression, and one prediction vector per scorer, plus the held-out peak count and the scorer's own reported Spearman.
- **Reconcile gate inside the script:** recompute Spearman *from the persisted vectors* and refuse (non-zero exit, no artifact) unless `gc` seed 0 is 0.3000 ± 0.01 and `max_tx` is 0.3968 ± 0.02. This is what proves the table is the published substrate rather than something that merely resembles it.
- Refuse a partial decomposition sweep: reuse `load_decomposition`'s ≥ 0.999 coverage rule and its resume message rather than re-implementing the check.
- `--seeds N`, `--scorers a,b,c`, `--smoke` and `--device` for cost control; the conv arm is the long pole and must be subsettable without changing the others.
- Write a small sidecar manifest recording the tile count, the filter, the protocol dict, `count_parameters()` for the conv arm, and the shard list — so a later consumer can tell two tables apart.

**Test strategy:**
- New `tests/test_glm_tile_predictions.py`, fully offline on synthetic tiles: the pure helpers (join construction, target transform, table serialisation round-trip, manifest shape) are separated from loading so they are testable with no D: and no model.
- Pin that the reconcile gate **refuses**: feed a deliberately wrong anchor and assert a non-zero exit and that no artifact file is written (non-vacuity — a gate never shown to refuse proves nothing).
- Pin that a train/test index overlap raises rather than scoring.
- Pin the round-trip: predictions written and re-read are index-aligned to the same tiles and recompute the identical Spearman.
- Pin that a missing scorer key is an error, not a silently-dropped column.

### Step 3: Measure what the tile oracle is worth to a CHOICE
Files: scripts/glm_tile_selection.py, tests/test_glm_tile_selection.py
Depends on: Step 2

**What changes:**
- New script `scripts/glm_tile_selection.py` consuming the Step 2 table. Per seed, per scorer, for k in (1, 2, 5, 10, 20): `simulate_selection` on that seed's held-out pool, `assert_k1_is_degenerate` at k=1, and ONE `permutation_null` per `(seed, k)` since a shuffled scorer is a uniform permutation whatever it came from.
- Freeze `PREREGISTERED` **before any number is read**, mirroring `glm_oracle_selection.py`'s bars so the two substrates are directly comparable: bar A = `headroom_captured ≥ 0.50` at k=10 for the best measured scorer; bar B = margin over `gc` ≥ `max(0.10, 2 × null stdev)`; majority of seeds must pass. Both thresholds are ASSERTED, not derived — the ceiling arm that normalises them IS derived, and the artifact must say exactly that.
- Verdicts mirrored: `ORACLE_READY_AS_SCORING_HEAD` / `CRUDE_BASELINE_SUFFICES` / `SELECTION_IS_THE_BOTTLENECK`, plus `INDETERMINATE_INSUFFICIENT_SEEDS` when fewer seeds scored than the bar needs — the refusal that exists because a single-axis run once printed `SELECTION_IS_THE_BOTTLENECK` while reporting that every scored axis passed both bars.
- Emit `wiki/glm_tile_selection_<date>.{md,json}` with per-seed cells, the null band per cell, each scorer's Spearman, the 0.941 ceiling for context, and a comparison row against the designed-grid artifact's cross-axis means (0.7937 learned / 0.1328 gc at k=10).
- **Honest limits, written into the artifact:** the pool is MEASURED peak tiles, so this is an upper bound on a generator's candidates exactly as the designed-grid run was; the tile pool is natural-sequence and peak-enriched rather than fully crossed, so it is a *different* upper bound, not a stricter one; and a per-seed headroom carries scorer-realisation spread that `trials` cannot reduce.

**Test strategy:**
- New `tests/test_glm_tile_selection.py`, offline on a synthetic prediction table.
- Pin the frozen-bar verdict function as a pure function of its inputs across all four branches, including `INDETERMINATE_INSUFFICIENT_SEEDS` — and prove that branch non-vacuous by supplying fewer seeds than the bar needs.
- Pin that a perfect scorer gives `headroom_captured ≈ 1.0` and a shuffled one sits inside its own null band — the instrument's two anchors.
- Pin that k=1 yields `headroom_captured is None` and that a fabricated non-degenerate k=1 cell raises.
- Pin that no verdict is emitted when the table's reconcile manifest is absent.

### Step 4: Measure orthogonality, then run the rank-average hybrid arm
Files: scripts/glm_tile_hybrid_arm.py, tests/test_glm_tile_hybrid_arm.py
Depends on: Step 1, Step 2

**What changes:**
- New script `scripts/glm_tile_hybrid_arm.py` consuming the Step 2 table.
- **Measure the precondition first and report it regardless of outcome:** per-seed Spearman between the `conv` and `max_tx` prediction vectors. The repo's measured rule is that a naive rank-average lifts only for orthogonal modalities, and the PEAR run measured the failure mode — rank-averaging a chance-level structure scorer took ESM2 from 0.352 to 0.204. Conv sits at 0.3234 against GC's 0.3112, so a negative result here is a legitimate and expected-possible outcome, not a bug.
- Combine with `rank_average_hybrid` over `{conv, max_tx}` keyed by tile position as a string, and additionally over `{conv, max_tx, gc}`. No fitting, no calibrator: both inputs are already out-of-sample on the held-out tiles, so the combine adds no second-level leakage.
- Score each combination's Spearman against raw measured expression; compute the per-seed **paired** margin against `max_tx` alone using the existing `paired()` helper (median, win count, CI95) — a difference of medians is not a paired comparison and this repo has a standing lesson on exactly that.
- Freeze the bar to mirror the decomposition arm's, so the two are comparable: `HYBRID_ADDS` iff paired median ≥ +0.02 AND ≥ 8/10 seeds win; `WEAK` if positive but short of either; `NO_GAIN` otherwise; `INDETERMINATE_RECONCILE_FAILED` if the table's anchors drift.
- Add a control that makes a positive interpretable: rank-average `max_tx` with a **shuffled** conv vector. If that also lifts, the gain is an artifact of rank-averaging two vectors rather than of the conv's content.
- Emit `wiki/glm_tile_hybrid_arm_<date>.{md,json}` carrying the orthogonality numbers, all arms, the paired margins with CIs, the control, and the verdict.

**Test strategy:**
- New `tests/test_glm_tile_hybrid_arm.py`, offline on synthetic vectors.
- Pin the headline mechanism in both directions: two genuinely orthogonal synthetic scorers rank-average to a *higher* Spearman than either alone; a near-chance second modality rank-averages to a *lower* one. This is the arm's own hypothesis made falsifiable before real data.
- Pin that the shuffled-conv control does **not** lift, so the control is non-vacuous.
- Pin the frozen verdict function across all four branches.
- Pin that `rank_average_hybrid` is called with ≥2 tables and raises on an empty key intersection (its own contract), and that keys are tile positions rather than list order.

### Step 5: Record both results, measure the arm, and route the findings
Files: soraya_runs/arms-ledger.json, NEXT.md, LESSONS_LEARNED.md
Depends on: Step 3, Step 4

**What changes:**
- Call `measure_arm` for `promcalc-plus-conv-rank-average-hybrid-150bp` with the measured Spearman and an evidence string naming the artifact, the seed count, the paired CI, the orthogonality number and the control's behaviour — written **from the artifact**, never from recall, because under-claiming is as much a trust-surface falsehood as over-claiming.
- **Call `save()` and verify from disk in a separate process.** Re-read `ceiling_verdict()` and report whether it still names a `planned` arm.
- If the hybrid loses, record it as a measured negative with an R7 scope line (the narrowest wording the executed evidence supports: these two scorers, this rung, this split, rank-average combine) and a revival condition that could actually fire — not a decoration.
- Append a short `NEXT.md` entry: the selection verdict, the hybrid verdict, and the one authority fork the selection number informs (whether the `ceiling:` criterion should be rewritten generatively). `NEXT.md` is explicitly transient — prune it, don't grow it.
- Append one concise `LESSONS_LEARNED.md` bullet **only if** a defect was found during execution; otherwise skip rather than manufacture one.
- **Do not add a long `CLAUDE.md` bullet.** The file is at its enforced word ceiling; the findings live in the two `wiki/` artifacts, and whether to compress `CLAUDE.md` to make room is a separate user decision (Open Questions).

**Test strategy:**
- Re-run the Step 1 test, which now reads the arm as `measured` from disk.
- Assert `dna_decode/eval/amr_rules.py` and `dna_decode/data/calibrated_amr_rules.json` are byte-identical via `git diff --quiet` on the **real paths** — `dna_decode/data/amr_rules.py` does not exist, and ad-hoc guards naming it have passed vacuously before.
- Re-verify the active prospective lock.
- Full suite: `uv run pytest tests/ -q`, expect ≥ 6,501 collected and 0 failures. Launch it only when no file is being edited — a prior baseline run was contaminated by concurrent edits in this very tree.

## Execution Preview

Waves (topological, from the `Depends on:` lines):

- **Wave 0 (parallel, 2 steps):** Step 1 (ledger declaration + its disk test), Step 2 (shared prediction table). No file overlap — Step 1 touches `soraya_runs/arms-ledger.json` and its own new test; Step 2 touches two new files.
- **Wave 1 (parallel, 2 steps):** Step 3 (selection value, consumes Step 2), Step 4 (hybrid arm, consumes Steps 1 and 2). No file overlap; two new scripts and two new tests.
- **Wave 2 (1 step):** Step 5 (ledger measurement + routing).

Total waves 3 · max parallelism 2 · critical path Step 2 → Step 3 → Step 5 (the toolkit's deterministic pick; Step 2 → Step 4 → Step 5 is the same length).

Cost shape: Step 2 is the only expensive step, and the conv arm is its long pole — ten `fit_encoder` trainings over ~33,000 tiles each. Steps 3 and 4 are pure arithmetic over a persisted table and run in seconds. Step 2 is deliberately subsettable (`--seeds`, `--scorers`, `--smoke`) so the cheap scorers can land and be verified before the conv arm is committed to.

## Risk Flags

1. **A prior recorded decision conflicts with a claim made in the conversation that produced this plan, and the plan does not resolve it.** `wiki/decisions-log.md` (2026-10-13) records: *"The fragment-rung recommendation flips. At a 0.46 'ceiling' the rung traded 7.3× the data for ~49% of headroom and looked unattractive; at the true 0.79–0.80 it trades for ~85% and is attractive. The wrong number argued against the right move."* This session separately measured that only **3.30%** of the 321,123 fragments exceed 2× the median expression. Those are statements about **different quantities** — a ceiling versus an informative-row fraction — and both can hold at once. Neither supersedes the other, and no step here touches the fragment arm. Carried to Open Questions.
2. **`scripts/glm_encoder_gate.py` is uncommitted and owned by another session.** Step 2 imports `ridge_gc` from it. The +15-line diff is purely additive (`--extra-widths`) and `ridge_gc` is verified absent from it, so importing is safe today — but a sibling's further edit could move it. Mitigation: the Step 2 reconcile gate on `gc` seed 0 = 0.3000 would catch a behavioural change immediately, and refuses rather than reports.
3. **Concurrent agents share this working tree.** A HIGH-salience retrospective records a baseline `pytest` run contaminated by edits made while it was still collecting, and a sibling mutating `dna_decode/glm/encoder.py` between two shell calls in one session. Mitigation: run the full suite only when nothing is being edited; use `TaskStop` by harness task id rather than killing PIDs; verify file identity with `git hash-object`, not file size — a length-preserving constant change is invisible to a size check.
4. **The hybrid arm may legitimately lose, and that must not read as a defect.** Conv at 0.3234 is barely above GC's 0.3112 and its paired CI against GC includes zero, so the orthogonality precondition the whole arm rests on may simply fail. The PEAR replication measured the exact failure mode. The plan handles this by measuring orthogonality first and by recording a negative as a negative (R6/R7) rather than retiring the arm.
5. **`CLAUDE.md` is at its enforced word ceiling** (35,460 of 35,460; 43 bullets totalling 24,677 words are provably stored elsewhere). A new long bullet trips `test_claude_md_growth_ceiling`. Step 5 routes findings to `wiki/` plus a short `NEXT.md` line instead.
6. **The designed-grid 0.7937 may not transfer in either direction.** It was measured with a learned scorer at ρ≈0.59 on a fully-crossed designed grid under leave-element-out; the tile rung is natural, peak-enriched, and its best scorer is biophysical at ρ 0.3968. The plan must not assume the pivot — a tile-rung headroom well below 0.50 would justify continued prediction work, which is the opposite conclusion and an equally useful one.
7. **`sentrux` resolves but has no rules file.** No `.sentrux/rules.toml` exists, so the architecture gate will skip silently. If an architecture baseline is wanted, `/setup-sentrux` must run first; this plan does not depend on it.
8. **Step 2's conv arm is the only real cost and it is CPU-bound here.** If ten seeds prove impractical, the honest degradation is fewer seeds reported as fewer seeds — which Step 3's `INDETERMINATE_INSUFFICIENT_SEEDS` branch and Step 4's 8-of-10 bar both already refuse to paper over.

## Open Questions

1. **Should the `ceiling:` criterion be rewritten generatively?** This is an acceptance-bar redefinition and therefore a user authority call, not an executor task. Step 3 produces the measurement that should decide it; the plan deliberately stops there.
2. **The fragment rung.** Risk Flag 1 names two measurements pointing in different directions about different quantities. The question neither has answered is whether a model trained on the fragment rung actually reaches a useful fraction of its own 0.7928 ceiling, and a reframing is also open: with ~97% of rows at background, a promoter-vs-not classification target may be better posed than a Spearman over mostly-noise. Not in scope here.
3. **Generator work stays parked.** `features/glm-generator-finetune/` is at `intake` with a deliberately deferred plan and a named unblock trigger. Whether to ratify its idea-anchor is a user decision; re-proposing it as new work here would be orientation error #3.
4. **`CLAUDE.md` compression.** 24,677 words across 43 bullets are provably stored elsewhere and the in-repo `wiki/arcs/` extraction pattern is already proven on two bullets. Whether to compress, raise the ceiling, or leave it is the user's call.
5. **Should the selection-value instrument eventually be run on the fragment rung too?** It would need that rung's own ceiling (0.7928) as the normaliser and its own split; cheap once Step 2's table shape generalises, but out of scope until question 2 is settled.

## Verification

1. Run Step 2 with `--smoke` first, then the cheap scorers only (`--scorers gc,max_tx,decomposition`), and confirm the reconcile gate passes on real data before committing the conv arm's compute.
2. Confirm the gate refuses: perturb an anchor in a scratch copy and verify non-zero exit with no artifact written.
3. Run Step 3; confirm `assert_k1_is_degenerate` passes for every scorer/seed, that every headroom is printed beside its null band, and that the emitted verdict is one of the frozen four.
4. Run Step 4; confirm the orthogonality number is reported before the verdict, that the shuffled-conv control does not lift, and that the paired margin carries a CI.
5. Re-read `soraya_runs/arms-ledger.json` in a separate process after Step 5 and confirm the arm's measured value and status are on disk, with an mtime later than the run.
6. `git diff --quiet dna_decode/eval/amr_rules.py dna_decode/data/calibrated_amr_rules.json` — both byte-identical.
7. Re-verify the active prospective lock.
8. `uv run pytest tests/ -q` with nothing being edited: ≥ 6,501 collected, 0 failures.
9. Diff both new `wiki/` artifacts against the claims made in `NEXT.md` — every number in the prose must appear in the JSON it cites, the standing rule `scripts/contract_number_audit.py` enforces elsewhere.

## Save-time amendments

Captured at: 2026-10-14
Source: `/save-plan` arguments

**Audit-notes-only contract.** This block is provenance for human readers. `/execute-plan` reads
ONLY `## Implementation Steps` for executable work — amendments are NOT executable instructions.
If an amendment changed a Step contract (files, dependencies, structure) the toolkit fires on
re-validation and `/technical-plan` must be re-run before `/execute-plan`.

- generator arm removed from scope — features/glm-generator-finetune already carries a deliberately-deferred plan
- fragment-rung conflict with decisions-log 2026-10-13 surfaced unresolved
- selection-value run deliberately NOT declared as a ceiling arm

<!-- toolkit: check=clean waves=clean gate=fired:open-questions,test-strategy-leak -->
