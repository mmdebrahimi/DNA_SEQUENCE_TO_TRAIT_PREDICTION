# Two Unattended Windows: Tonight's 10-12 Hours And The 10-Day Vacation Run

## Lens status
Dispatch: conversation-driven with inline scope ("plan for the next 10-12 hours" + "the longest run you could do" for 2026-11-04 to 2026-11-13).
Repo-index: absent — built by direct file reading, no Layer 2 summaries.
Reuse signals: no `graphify-out/GRAPH_REPORT.md` and no conventional shared-code directory; the survey below is grounded in six existing `run_*_detached.bat` launchers and `scripts/glm_run_lock.py`.
Tools: `sentrux` present, `gh` present, `codex` ABSENT (so no second-opinion review is available for this plan).
Project rules: `project-rules.md` absent, so `## Project Rules Check` is omitted by contract.
Measured during planning: drive capacities, clone size, `PA_REPO` override, resume support, and the host's recorded long-run death history.

## DECISIONS NEEDED

**D1 — One user action gates both windows, and no code can substitute for it: the machine must not sleep.** A detached process cannot run through standby. `powercfg` returned nothing readable for the current scheme during planning, so I am NOT claiming the current setting either way. Before leaving, confirm sleep and hibernate are set to Never on AC for both windows, and that the laptop stays on mains. If that is not acceptable for ten days, say so and the vacation window becomes a Kaggle chain instead, which is strictly worse (see Risk Flags) but survives a sleeping host.

**D2 — The fine-tune's wall-clock ceiling is currently 6 hours, and it is the thing that routes the arm to Kaggle.** `LOCAL_WALL_CEILING_SECONDS = 6 * 3600` is recorded in the probe as *operational, not statistical*: its stated reason is that an overnight run is recoverable while a multi-day local run "tends to be abandoned part-way". A ten-day authorised unattended window falsifies that premise directly. This plan therefore freezes the fine-tune bar against a **window-derived** ceiling instead of a fixed one, which selects the LOCAL route. If you would rather keep the 6-hour constant and accept Kaggle, say so before Step 3 runs; afterwards the bar is frozen and the route cannot move.

**The local route is better here for an operational reason, not a cost one.** Retrieving a Kaggle kernel's output requires someone to run `kaggle kernels output`. I only execute when invoked, and across ten days nothing will invoke me — so a Kaggle run would finish on day one and then sit unread until the 13th, and a 9-hour-per-kernel cap means the full horizon needs a chain that something has to advance. A detached local run writes to disk continuously and needs nobody.

## CORRECTIONS

**The standing "heavy files live on D:" rule is deliberately inverted for the durable clone, and the reason is this host's recorded death history.** `run_tile_conv_detached.bat` documents six deaths of long unattended jobs: four were `STATUS_IN_PAGE_ERROR (0xc0000006)` with torch memory-mapping 4.33 GB of weights off the USB volume, and two produced no event on any of 118 logs and remain unexplained. The single long job that ever reached `rc=0` on this host, at roughly 12 hours, had its data on C:.

Measured now: C: has **30.3 GB free**, D: has **823.8 GB free** and is attached. The clone is **718.8 MB**, of which the weights are only **87.7 MB**. So the durable copy that both windows depend on is ~88 MB plus a small source tree — not a heavy artifact, and 0.3% of C:'s free space. Putting it on C: keeps every byte these runs memory-map on the volume that matched the one survivor, and it leaves the D: rule intact where it was aimed: the bulk caches stay on D:.

This is a favourable-but-not-proven difference, in exactly the register the existing launcher uses. It is not a guarantee, which is why every step below checkpoints and writes a sentinel rather than assuming success.

## Problem Statement

Four steps of the in-flight plan remain — freeze the fine-tune bar, score the expression-heads arm, run and score the fine-tune arm, then record all three arms — and all of them are long or blocked. Two unattended windows are now available: roughly 10-12 hours tonight (2026-11-02 into 2026-11-03), and roughly 10 days from 2026-11-04 to 2026-11-13.

The binding constraint is not compute, it is **attendance**. I execute only when invoked, so anything that needs me mid-flight cannot be scheduled into either window: that rules out a Kaggle chain, rules out harness-tracked background shells (one was already reaped for system memory pressure while the session sat idle), and makes resumability after an unplanned reboot a requirement rather than a nicety.

The goal is to convert the remaining work into two self-contained detached runs that need no attendance, leave a readable state behind, and produce exactly the artifacts the frozen bars require — with the vacation window used for the one measurement that genuinely benefits from days rather than hours.

## Codebase Context

**What the two windows must run, with measured costs.** The expression-heads arm is 4.64 hours over five heads and 44,106 tiles on the faithful per-window route, projected from a 48-tile timing sample; its bar is frozen and its scorer is written and streaming. The fine-tune arm is 13.5 minutes per full-cohort epoch at a measured 40.7 training sequences per second including backward, over 33,080 training tiles, so 10 seeds at 6 epochs is 13.55 hours — already past the current 6-hour ceiling, which is why the route question is live.

**The ten-day window should be ONE long-horizon run, not three runs at three budgets.** Training is sequential and deterministic given a seed, so a run carried to a long horizon with the validation metric recorded at every epoch contains every shorter budget as a prefix: the frozen budget's verdict is read off the same seeds and the same data order, and the longer budgets come free. Three separate runs at E, 2E and 4E epochs would cost 7E epochs per seed for strictly less information. At 13.5 minutes per epoch, 10 seeds at a 60-epoch horizon is roughly 135 hours, which fits inside 240 hours with substantial margin for reboots.

**Why the longer budgets are worth having at all.** The fine-tune bar must pre-commit that a null has two readings — the frozen representation was not the limit, versus the substrate or the budget was — and that only a convergence curve separates them. Today that separation rests on a *subsample* curve. A full-cohort run whose metric is flat from epoch E to epoch 4E closes the budget-limited reading **empirically** instead of by inference. The longer budgets are therefore pre-declared INELIGIBLE diagnostics, exactly as max-over-windows is on the heads arm: they answer "would more budget have helped?" and can never become the arm's result.

**The detached-run convention already exists and has hard-won details that must be copied, not re-derived.** Six `run_*_detached.bat` launchers are in the repo root. Three of their properties are load-bearing:

- Launch is `Start-Process -FilePath <bat> -WorkingDirectory <repo> -WindowStyle Hidden`, which outlives the session. This is NOT a harness background shell and so is not subject to the memory reaper that already killed one job while the session was idle.
- The sentinel redirect must come **before** the echo: `>>"%SENTINEL%" echo EXIT=%RC%`. Written the natural way round, `cmd` reads the digit immediately preceding `>>` as a file handle, redirects stdin, and leaves the sentinel EMPTY. The existing launcher records this as measured on a real run, not theorised.
- `set PYTHONUNBUFFERED=1`, a per-run `.log` appended with `>> "%LOG%" 2>&1`, and a sentinel deleted at start so a stale one cannot be mistaken for this run's result.

**`PA_REPO` is already an environment override**, so pointing both runs at a durable clone needs no code change: `atlas_embed.py:45` reads `os.environ.get("PA_REPO", <temp path>)`.

**Resumability is already most of the way there.** The heads scorer resumes from its sidecar's `n_done` and re-opens the memmap `r+`; a mid-run death costs one batch of 64 tiles. The fine-tune runner must add the per-seed equivalent.

**`scripts/glm_run_lock.py`** provides `acquire`/`release`/`holder`/`pid_alive` and `run_guarded(path, command, label)` with `HELD_EXIT = 9`, which is the right guard against a second launch of the same job over a long window.

### Reusable-Code Survey

Searched: the six `run_*_detached.bat` launchers in the repo root, `scripts/glm_run_lock.py`, `scripts/expression_heads_score.py`, `scripts/finetune_arm_probe.py`, `scripts/freeze_expression_heads_bar.py`, `scripts/positional_arm_probe.py`, `scripts/positional_arm_score.py`, `dna_decode/glm/atlas_embed.py`, `dna_decode/glm/atlas_heads.py`. No `graphify-out/GRAPH_REPORT.md`; no `src/utils`-style directory exists.

Reuse rather than rebuild:

- **The launcher shape** from `run_tile_conv_detached.bat`, including the sentinel-redirect order, the unbuffered env, the stale-sentinel delete and the `Start-Process` usage line in the header comment.
- **`glm_run_lock.run_guarded`** for single-instance protection, with `HELD_EXIT = 9` distinguishing "already running" from a real failure.
- **The streaming extractor and its PARTIAL discipline** in `expression_heads_score.py` — `extract`, `read_sidecar`, `load_scalars`, which already refuses a PARTIAL block.
- **The gate chain and statistics** — `GATE_FUNCS`, `paired_bootstrap` (`N_BOOT = 10_000`, `BOOT_SEED = 707`), `median_of`, `leave_peak_out_indices`.
- **The bar skeleton** from `freeze_expression_heads_bar.py`: closed `OUTCOMES`, `CHECK_ORDER`, `reachable_outcomes()` asserted equal, total `classify()`, a `verify()` that refuses drifted anchors, and the `amendments_before_first_score` record.
- **The probe's pure helpers** — `choose_route`, `project_arm_cost`, `curve_turned`, `deviations_from_recipe`, `seconds_per_epoch`, `free_ram_bytes`, `sha256_of`, `load_arm_cohort`.
- **`atlas_embed.audit_checkpoint`** for verifying the durable clone's weights after the copy.

## Pre-Change Baseline

**Full suite, measured at commit `3a85a638`:** 12 failed / 7,630 passed / 180 skipped in 11m42s. All 12 failures are absent-data or absent-tool cases (cyp4f2 7, essentiality-ladder 3, pneumo 1, README quickstart 1) and are the diff reference by name.

**Work already committed on branch `promoteratlas-exposure-2026-11-01`,** head `ed554966`, nothing pushed: the shared head module (50 tests), the expression-heads feasibility probe (30 tests, artifact written), the fine-tune probe's code (42 pass / 7 skip, measurement in flight), the expression-heads bar (46 tests, frozen then amended), and three review-driven corrections.

**In flight at planning time:** the fine-tune convergence measurement is running detached from this plan; `wiki/finetune_arm_feasibility_2026-11-02.json` is **ABSENT** until it lands. Step 3 below cannot run before it exists, and that is an external prerequisite rather than a plan dependency.

**Arms state:** 16 arms — 12 measured, 3 planned, 1 killed; `ceiling_verdict = blocked:untried-arm`; `best` is the atlas arm at 0.5076, which is 54% of the 0.941 replicate ceiling.

**Host:** C: 30.3 GB free of 475.9; D: 823.8 GB free of 931.5, attached; 15.9 GB RAM with ~6 GB free; CPU-only, 4 torch threads; Kaggle CLI 2.2.5 present.

**Architecture gate:** not re-baselined, deliberately, so the final compare stays cumulative.

## Verification Signal

1. **The durable clone is verified by hash, not by file count.** Every copied weight file re-passes `audit_checkpoint` and the base checkpoint's sha256 matches `bb413bc8602f05b6…` at 22,313,814 bytes. A copy that merely has the right names is not evidence of the right contents.
2. **Both launchers are proven to survive session exit** by a short smoke job that outlives a deliberately-closed shell, with its sentinel read afterwards — not by assuming `Start-Process` does what the comment says.
3. **The sentinel is proven non-empty.** The redirect-order trap is pinned by a test that writes a sentinel both ways and asserts the natural ordering produces an empty file while the reversed one does not.
4. **Single-instance protection is exercised in both directions:** a second launch while the first holds the lock exits `9`, and a launch after release succeeds.
5. **Resume is proven to reproduce, not merely to continue.** The heads run is interrupted and resumed, and the resulting scalars are bit-identical to an uninterrupted run on the same slice; the fine-tune's per-seed resume retrains only unfinished seeds and reproduces the finished ones.
6. **The frozen budget's verdict is read off the long-horizon run's prefix** and is asserted equal to what a run stopped at that budget would have produced, so the diagnostics cost nothing in fidelity.
7. **The longer budgets cannot become the arm's result** — a test drives a diagnostic that would win and asserts the verdict function never sees it.
8. **Every verdict is re-derivable from its own artifact** plus its frozen bar, with the reachable outcome set asserted equal to the closed set.
9. **The failure set is diffed by name** against the 12-name baseline, never by count.
10. **The frozen five and `uv.lock` stay byte-identical**, and the prospective lock re-verifies.

## Implementation Steps

### Step 1: Durable clone on C:, verified by hash
Files: scripts/durable_clone.py, wiki/durable_clone_2026-11-02.json, tests/test_durable_clone.py
Depends on: none

**What changes:**
- Copies the parts of the PromoterAtlas clone both runs actually need — `trained_weights/` (87.7 MB) and `src/promoter_atlas/` — from the volatile temp directory to a durable path under C:, defaulting to `data/pa_clone_durable/` and overridable.
- **Target is C: deliberately, inverting the standing D: rule**, because the host's four `STATUS_IN_PAGE_ERROR` deaths all involved torch memory-mapping weights off the USB volume while the one long survivor kept its data on C:. 88 MB is not a heavy artifact, so the D: rule's purpose — keeping a chronically-full C: clear — is not violated. Both facts are recorded in the artifact so the exception is explicit rather than silent.
- Verifies after copying: the base checkpoint's sha256 against `bb413bc8602f05b6…` and 22,313,814 bytes, every head's sha256 against a freshly computed source digest, and `audit_checkpoint` re-passing on all six checkpoints. A mismatch is a non-zero exit, not a warning.
- Emits the `PA_REPO` value the launchers must export, so the two `.bat` files read it from the artifact rather than hardcoding a path twice.
- Idempotent: a second run re-verifies and copies nothing.

**Test strategy:**
- Unit-test the pure copy-manifest and verification-summary helpers on a synthetic tree with known digests.
- Assert a corrupted destination byte is DETECTED rather than reported as a successful copy.
- Assert a second invocation copies zero files and still re-verifies.
- Clone-touching assertions skip cleanly when the temp directory has already vanished.

### Step 2: The detached launcher harness and its status surface
Files: scripts/detached_launch.py, tests/test_detached_launch.py
Depends on: none

**What changes:**
- A small generator that emits `run_<job>_detached.bat` matching the existing convention exactly: `@echo off`, `cd /d "%~dp0"`, `set PYTHONUNBUFFERED=1`, stale-sentinel delete, append-redirected log, and the sentinel written with **the redirect first** (`>>"%SENTINEL%" echo EXIT=%RC%`), because the natural ordering makes `cmd` treat the digit before `>>` as a file handle and silently produces an empty sentinel.
- Wraps the payload in `glm_run_lock.run_guarded` so a second launch over a long window exits `9` instead of corrupting a half-written block.
- Writes a **heartbeat** line on a fixed cadence (job, pid, unit done, units total, wall seconds, estimated finish) to a status JSON on C:, so on return the state is readable without re-deriving it from a log tail.
- A `status()` reader that classifies each job as `RUNNING` (heartbeat fresh, pid alive), `DEAD_NO_SENTINEL` (pid gone, no sentinel — the unexplained-death signature this host has twice), `FAILED` (sentinel with non-zero), or `COMPLETE`. The `DEAD_NO_SENTINEL` state exists because two of the six recorded deaths left no trace anywhere, and a reader must be able to tell that apart from a clean finish.

**Test strategy:**
- Pin the sentinel-redirect trap: write both orderings and assert the natural one yields an empty file and the reversed one does not.
- Assert the generated `.bat` contains the lock wrapper, the unbuffered env and the stale-sentinel delete, by parsing the emitted text.
- Drive `status()` through all four states with synthetic heartbeats and sentinels, including a fresh heartbeat with a dead pid.
- Launch a 5-second smoke job, close the launching shell, and assert the sentinel appears afterwards — the only way to establish that the process genuinely outlives the session.

### Step 3: Freeze the fine-tune bar against a window-derived ceiling
Files: scripts/freeze_finetune_arm_bar.py, wiki/finetune_arm_acceptance_bar.json, tests/test_freeze_finetune_arm_bar.py
Depends on: none

**What changes:**
- Freezes the seed count, epoch budget, learning rate, early-stopping rule and compute route from the convergence measurement, before any score exists. **External prerequisite:** `wiki/finetune_arm_feasibility_2026-11-02.json` must exist and its verdict must start with `FEASIBLE`; `verify()` refuses otherwise.
- **Replaces the fixed 6-hour ceiling with a window-derived one** and records why: the constant's own stated justification is operational — that a multi-day local run tends to be abandoned part-way — and an authorised ten-day unattended window falsifies that premise. The artifact names the window, the derived ceiling and the resulting route, so a reader can see the route followed from the window rather than from a preference.
- Inherits the margin 0.02 and the committed atlas comparator; the comparison is the fine-tuned encoder's median Spearman against the FROZEN encoder's 0.5076 on byte-identical splits.
- **Inherits the target transform explicitly**: trained on `log1p(expression)`, scored against RAW, matching `ridge_spearman`, whose docstring states any deviation makes the number incomparable with the committed incumbent tables.
- Declares the frozen epoch budget as the arm's result and **pre-declares the longer horizons as INELIGIBLE diagnostics**, with the stated purpose of closing the budget-limited reading empirically.
- Declares a shuffled-label null on the fine-tune itself — and, unlike the heads arm, this one is NOT vacuous: the fine-tune genuinely fits on this substrate, so destroyed labels destroy the fit. The artifact records that the heads arm needed a predictor-side null precisely because a one-feature zero-shot predictor cannot be damaged by relabelling, so the two arms' nulls differ for a measured reason.
- Closed outcome set with reachable equal to the check order, including `INDETERMINATE_BUDGET_EXHAUSTED` for a run that stopped before its frozen budget and `INDETERMINATE_RUN_DIED` for the `DEAD_NO_SENTINEL` case.

**Test strategy:**
- Assert `verify()` refuses a mutated probe artifact, a mutated comparator, a missing arm, an arm already carrying a metric, and a probe whose verdict is INDETERMINATE.
- Assert reachable equals the check order, with the padding and trimming non-vacuity proofs.
- Assert the ceiling is DERIVED from the declared window in both directions — a short window routes to Kaggle, a ten-day window routes local — so the rule is not a disguised constant.
- Assert the budget-versus-representation distinction references the probe's curve field by name and the diagnostics by name.
- Assert the two arms' nulls differ and that the bar says why.

### Step 4: Launch the expression-heads scoring run detached (tonight's window)
Files: run_expression_heads_detached.bat, tests/test_expression_heads_score.py
Depends on: Step 1, Step 2

**What changes:**
- Generates and launches the heads scoring run under the durable `PA_REPO`, writing its block to C: (3.5 MB of scalars — small, and on the volume that matched the survivor), its log and sentinel beside it, and heartbeats into the shared status JSON.
- Runs the **faithful per-window route** the bar froze, at a measured 4.64 hours, which fits tonight's window with roughly a 2× margin.
- The scorer's existing PARTIAL discipline is the safety net: an interrupted run resumes from its sidecar and a capped run is refused as the arm's result.
- Adds the scorer's own test file, which the written code does not yet have: every `INDETERMINATE_*` branch driven to reachability with stub heads, a constant-output head tripping the validity gate rather than scoring near zero, a monotonically-transformed head scoring an identical Spearman, and the verdict re-derivable from the artifact alone.

**Test strategy:**
- Drive the four gates to failure individually with stubs and assert each returns its declared outcome.
- Assert the permutation null's two conditions are both required, by constructing a case that passes the median and fails the magnitude.
- Assert an interrupted extraction resumes and produces scalars bit-identical to an uninterrupted run over the same slice.
- Assert the loader refuses a PARTIAL sidecar and that `--limit` produces one.

### Step 5: The long-horizon fine-tune runner
Files: scripts/finetune_arm_run.py, tests/test_finetune_arm_run.py
Depends on: Step 3

**What changes:**
- Trains the encoder per seed on byte-identical `leave_peak_out_indices(frac=0.25)` splits at the full cohort, on `log1p(expression)`, scoring validation Spearman against RAW expression every epoch.
- **Runs to a long horizon once and records every epoch**, so the frozen budget's result is a PREFIX of the same run rather than a separate run: same seed, same data order, same initialisation. The frozen early-stopping rule is then applied by reading the recorded curve, which is what makes the longer budgets free rather than a 7× cost.
- Horizon is derived from the declared window and the measured throughput, with a hard wall-clock stop so the run cannot overrun the window and leave every seed half-done.
- Checkpoints per seed AND per epoch: a death costs at most one epoch, and a resume retrains only unfinished seeds. Written to C:, each under 8 MB.
- Emits per-seed, per-epoch validation Spearman plus the shuffled-label null at the frozen budget, so the result artifact can be assembled without re-running anything.

**Test strategy:**
- Assert a resumed run reproduces the per-seed curves of an uninterrupted one bit-for-bit on a tiny synthetic cohort, and retrains only the unfinished seeds.
- Assert the frozen-budget verdict read off a long-horizon prefix EQUALS the verdict of a run stopped at that budget — the property the whole design rests on.
- Assert the wall-clock stop fires and stamps the run as bounded rather than complete.
- Assert the longer-horizon diagnostics are marked ineligible and that a diagnostic which would win does not enter the verdict.
- Assert training on raw expression instead of `log1p` is refused, so the C3 correction cannot silently regress.

### Step 6: Launch the vacation run detached
Files: run_finetune_arm_detached.bat
Depends on: Step 1, Step 2, Step 5

**What changes:**
- Generates and launches the long-horizon fine-tune under the durable `PA_REPO`, sized to the 2026-11-04 to 2026-11-13 window with margin: at 13.5 minutes per full-cohort epoch, 10 seeds at a 60-epoch horizon is roughly 135 hours against 240 available, leaving headroom for reboots and for the host running slower than its sampled rate.
- Touches **no external volume**: inputs are the durable clone on C: and the committed cohort; outputs are per-seed checkpoints under 8 MB each. That is the profile of the one long job that survived on this host, and the deliberate avoidance of the profile that died four times.
- Writes heartbeats so that on return the status reader can distinguish COMPLETE from FAILED from the DEAD_NO_SENTINEL signature this host has produced twice with no log evidence anywhere.

**Test strategy:**
- Assert the generated launcher exports the durable `PA_REPO` from Step 1's artifact rather than a hardcoded path.
- Assert it refuses to launch when the durable clone has not been verified.
- Assert a second launch exits 9 while the first holds the lock.

### Step 7: Collect both runs and emit their result artifacts
Files: scripts/collect_unattended_runs.py, wiki/finetune_arm_result_2026-11-02.json, wiki/unattended_runs_status.json, tests/test_collect_unattended_runs.py
Depends on: Step 4, Step 6

**What changes:**
- Reads both runs' checkpoints, sentinels and heartbeats and emits one status artifact plus the fine-tune arm's result artifact, assembled from recorded curves without re-running anything.
- Applies the frozen early-stopping rule to each seed's recorded curve to produce the arm's metric, then reports the longer horizons as the pre-declared ineligible diagnostics.
- **Refuses to emit a verdict for a run in the `DEAD_NO_SENTINEL` state**, returning `INDETERMINATE_RUN_DIED` instead, because a partial curve from an unexplained death is not evidence of convergence or of its absence.
- Runs the paired bootstrap and the gate chain from the frozen bar, so the verdict is re-derivable from the artifact alone.

**Test strategy:**
- Drive every verdict branch from synthetic checkpoint sets, including the died-without-sentinel case.
- Assert the collected frozen-budget metric equals the one the runner would have reported directly.
- Assert the diagnostics appear and are flagged ineligible.
- Assert the status artifact distinguishes all four job states.

### Step 8: Record the three arms, render the memo, journal it
Files: soraya_runs/arms-ledger.json, scripts/successor_arms_memo.py, wiki/successor_arms_result.md, project_state/glm-learned-representation-2026-10-07.md
Depends on: Step 7

**What changes:**
- Moves each of the three arms `planned -> measured` by direct dictionary write, because `measure_arm` ends in `Arm(**base)` with no field filter and raises `TypeError` on the annotation keys these arms carry — after mutating the dictionary, which is the worse half.
- Each arm's metric is the max over ITS OWN declared primaries; never over a ladder, never across arms, and never over an ineligible diagnostic.
- Recomputes `ceiling_verdict` and REPORTS it. With no `planned` arm left it reads `met` again, which is the same authority call surfacing a second time; the memo states it and the run stops rather than recording it in passing.
- **Expects the `best` pin to fire** and treats a re-pin as correct: the as-given arm's 0.5086 displaces the atlas arm's 0.5076. The memo must state that the leader changed by 0.0010 on an arm whose own verdict is TIES with a zero-containing interval, so nobody reads `best` as "the best approach".
- Renders one memo from the artifacts, ordered by required reader action, carrying the observed-before-declared hazard on the as-given arm, the windowing approximation and the translation-endpoint mismatch on the heads arm, and the budget-versus-representation distinction — now settled empirically — on the fine-tune.
- Appends one Action-Log row per arm inside the `<!-- project-state:end:action-log -->` marker, and re-runs the converted invariant guards to assert they are still green after three status moves.

**Test strategy:**
- Assert the ledger round-trips with every annotation key preserved, by semantic diff.
- Assert the memo refuses to render any arm on an `INDETERMINATE_*` verdict.
- Assert the invariant guards hold with three more measured arms.
- Re-run the full suite and diff the failure set BY NAME against the 12-name baseline.

## Execution Preview

- **Wave 0:** Step 1, Step 2, Step 3 (parallelism 3)
- **Wave 1:** Step 4, Step 5 (parallelism 2)
- **Wave 2:** Step 6
- **Wave 3:** Step 7
- **Wave 4:** Step 8

Five waves, maximum parallelism 3, critical path Step 3 → Step 5 → Step 6 → Step 7 → Step 8. No two steps in a wave share a file.

**Run sequentially regardless**, for the same reason as the prior plan: the feature blocks are gitignored and the gates read paths a worktree cannot reach, so parallel mode would silently skip them.

**Window fit.** Tonight: Steps 1, 2, 3 are minutes of work, then Step 4's launch starts a 4.64-hour run that completes inside the 10-12 hours with roughly 2× margin. The vacation window: Step 6's launch starts roughly 135 hours of compute against 240 available. Steps 7 and 8 run on return and need attendance, which is correct — they produce the published verdicts.

## Risk Flags

- **The host has killed six long unattended runs [grounded, from `run_tile_conv_detached.bat`].** Four were `STATUS_IN_PAGE_ERROR (0xc0000006)` while torch memory-mapped 4.33 GB of weights off the USB volume; two left no event on any of 118 logs. Mitigations: nothing in either run touches an external volume, the largest mapped file is 21 MB rather than 4.33 GB, every unit is checkpointed, and the status reader has an explicit `DEAD_NO_SENTINEL` state so an unexplained death is legible rather than invisible. This is favourable-but-not-proven and is not a guarantee.
- **The machine must not sleep, and no code can enforce that.** Carried as D1. A detached process does not run through standby.
- **A Windows Update reboot inside a ten-day window is plausible.** Per-seed and per-epoch checkpointing bounds the loss to one epoch, but nothing restarts the run automatically — on return it may need one resume invocation. A logon-triggered scheduled task would close that gap and is deliberately NOT in scope, because a task that relaunches a training job unattended is a bigger change than the gap it fixes.
- **The harness memory reaper already killed one background job while the session was idle.** Both windows therefore use `Start-Process`-detached launchers rather than harness background shells. If you want extra protection, Claude Code can be started with `CLAUDE_CODE_DISABLE_BG_SHELL_PRESSURE_REAP=1` — but note that only affects harness-tracked jobs and setting it from a shell command has no effect.
- **C: has 30.3 GB free and is chronically full.** Everything this plan writes to C: is small — 88 MB clone, 3.5 MB scalars, under 80 MB of checkpoints — but a third party filling C: during the window would kill both runs. D: has 823.8 GB free and remains the home for bulk caches.
- **The amendment window on both bars closes when their runs start.** The heads bar is already frozen and amended; after Step 4 launches it is fixed. The fine-tune bar is not yet frozen, so D2 is answerable until Step 3 runs and not afterwards.
- **The clone lives in a volatile temp directory until Step 1 copies it.** If a reboot removes it before Step 1 runs, both windows are blocked and the clone must be re-fetched.
- **No second-opinion review is available for this plan** — `codex` is absent and the local model endpoint is unreachable, so the review that found three defects in the previous plan cannot be repeated here. Treat the step designs as single-reviewed.
- **`powercfg` output was unreadable during planning**, so the current sleep setting is `[unverified]` rather than assumed to be Never.
- **Do not run the full suite while a torch job holds a model** — pytest dies with an access violation that reads like a test failure. Both windows hold models for hours, so the suite runs only between windows.

## Open Questions

1. **The sleep setting (D1).** The one gate no code can close. Unverified at planning time.
2. **The ceiling constant versus a window-derived ceiling (D2).** Answerable until Step 3 freezes the bar.
3. **Whether a logon-triggered scheduled task should auto-resume after a reboot.** Deliberately out of scope here: it would make an unattended training job self-relaunching, which is a larger change than the one-invocation gap it removes. Raise it if ten days without a manual resume is unacceptable.
4. **What the frozen epoch budget actually is.** Derived from the convergence measurement that is still in flight, so it is intentionally unanswered rather than guessed.

## Verification

- Steps 1 and 2 are verified by hash checks that detect a corrupted copy, by the sentinel-redirect trap pinned in both orderings, by the four-state status reader driven through every state, and by a smoke job proven to outlive a closed shell.
- Step 3 is verified by `verify()` refusing each mutated anchor, by reachable-equals-declared with its non-vacuity proofs, and by the ceiling being shown to route both ways from the declared window.
- Steps 4 to 6 are verified by interrupt-and-resume reproducing bit-identically, by the PARTIAL refusal, and by the single-instance lock exiting 9 on a second launch.
- Step 7 is verified by the frozen-budget metric read off a prefix equalling a directly-reported one, and by the died-without-sentinel case returning INDETERMINATE rather than a verdict.
- Step 8 is verified by a semantic ledger diff preserving every annotation key, by the invariant guards staying green across three status moves, and by a by-name failure-set diff against the 12-name baseline.
- Throughout: the frozen five and `uv.lock` stay byte-identical, the prospective lock re-verifies, and the architecture compare is advisory only.

## Save-time amendments

Captured at: 2026-11-02
Source: `/save-plan` arguments

These are AUDIT NOTES for human readers, not executable instructions. `/execute-plan` reads
only `## Implementation Steps`; nothing here changes a step's file list, dependencies or
structure. Where an amendment supersedes a FIGURE in the body, the body is left as emitted and
the correction lives here, so the plan stays the artifact that was reviewed.

- Two unattended windows planned: durable clone on C: inverting the D: rule on the host's six-death history; one long-horizon fine-tune whose prefix yields the frozen budget; D1 sleep setting unverified and D2 window-derived ceiling open until Step 3
- SUPERSEDES A FIGURE IN THE BODY. The convergence measurement landed AFTER this plan was emitted and returned `FEASIBLE_CURVE_TURNED`: best epoch 6 at 0.4309 followed by eight non-improving epochs, so Open Question 4 is ANSWERED and the frozen budget is 7 epochs.
- Throughput re-measured at 27.0 train seq/s at batch 128, not the 40.7 at batch 32 the body cites -- a 34% swing on the same host, with the fastest batch size moving. A full-cohort epoch is therefore 20.4 min, not 13.5, and Step 6's illustrative 60-epoch/10-seed horizon is ~204 h rather than ~135 h: 85% of the ten-day window instead of comfortable. On today's rate a safe horizon is nearer 42 epochs per seed (~144 h), still 6x the 7-epoch frozen budget.
- Step 5's design is UNCHANGED by that correction and is what absorbs it: the horizon is already specified as derived at launch from a fresh throughput measurement with a hard wall-clock stop, precisely so a stale rate cannot fix it in advance. The variance turns that stop from prudent into mandatory.
- D2 is now concrete rather than hypothetical: at 23.82 h projected against the recorded 6 h constant the arm routes to Kaggle by 4.0x, while a window-derived ceiling routes it local.

<!-- toolkit: check=clean waves=clean gate=fired:open-questions,unverified -->
