# Run The Fine-Tune Arm As A Bounded Chained Sequence — Technical Plan

> Reach a gradeable verdict on `finetuned-promoteratlas-encoder-150bp` by running it as several self-bounding Kaggle sessions that resume from per-epoch checkpoints, rather than the single session the measured rate has ruled out.

## Lens status
Inputs: conversation (the measured smoke result and the review that followed it).
Degradations: codex not installed and LM Studio unreachable, so that review was single-voice; `repo-index-status.sh` absent on this host, so Step 1 read source directly; `graphify-out/` absent; `project-rules.md` absent, so the Project Rules Check section is omitted; `DESIGN.md` absent and this plan touches no UI; sentrux installed but this change is not structural, so Architecture Impact is omitted.

**This plan SUPERSEDES the unsaved "Harden The Fine-Tune Smoke Into A Faithful Rehearsal" plan emitted earlier in this conversation.** That plan's Problem Statement rested on the claim that the arm cannot start at a session-sized ceiling, which Step 1 research has since falsified; and its central fix — a `run()` rehearsal inside the smoke — is made redundant by a bounded full run, which exercises `run()` end to end while doing real training. Its two surviving items (the `seed_checkpoints` refusal, the session-cap resolution) are carried here as Steps 2 and 4.

## Problem Statement

The smoke succeeded and measured **1276.04 s/epoch** on a Tesla T4 (`sm_75`, pin honoured). That number rules out the plan everything downstream assumed: the bar's full horizon of 10 seeds x 38 epochs projects to **134.69 h**, and even the gradeable prefix — the frozen 7-epoch budget across 10 seeds plus 21 null epochs — costs **32.3 h**. No Kaggle session is that long.

A first reading concluded the arm therefore cannot run on Kaggle at all, because `plan_run` refuses whenever the ceiling leaves fewer than `frozen_epochs` per seed, and it refuses at 8.0 h, 9.0 h, 10.5 h and 12.0 h against this rate. **That conclusion was wrong and the correction is the reason this plan exists.** `run()` calls `plan_run` only when `horizon_epochs is None`; supply the horizon explicitly and planning is bypassed, with the record stamping `"horizon supplied by the caller, not planned from the ceiling"`.

With the horizon supplied, the pieces already in the runner compose into a clean chained run. `advance()` checks `if time.time() >= deadline: bounded = True; break`, so a session stops **gracefully** at its own ceiling and the kernel writes its artifact — no reliance on Kaggle preserving output after a timeout kill. `_write_curve` fires inside the per-epoch loop, so progress persists per **epoch**, not per seed. And `advance()` resumes with `if prior and len(prior["val_spearman_by_epoch"]) >= to_epochs:` — read, never recomputed.

What remains is to bake the horizon so a Kaggle kernel can set it, make the one silent-failure path in checkpoint seeding loud, and run the sequence.

**Scope:** this plan ends when a run artifact exists carrying a gradeable record. Scoring it, registering the arm and recording the result remain Steps 5-8 of `plans/Run_The_Re_Routed_Fine_Tune_Arm_On_Kaggle_T4_And_Score_It_Plan/technical-plan.md`, whose Step 4 this plan replaces.

## Codebase Context

**The bypass, read at the call site.** `scripts/finetune_arm_run.py:505-513`:

```python
if horizon_epochs is None:
    plan = plan_run(ceiling_hours=ceiling, seconds_per_epoch=spe, n_seeds=n_seeds, ...)
    horizon_epochs = plan["horizon_epochs"]
else:
    plan = {"horizon_epochs": int(horizon_epochs), ...,
            "note": "horizon supplied by the caller, not planned from the ceiling"}
```

Supplying the horizon is a supported path that records itself. It is not a loophole, and the `note` is what makes a later reader able to tell the two modes apart.

**The three properties that make chaining safe, each read rather than assumed.** `advance()` breaks on `time.time() >= deadline` and sets `bounded`, so a session self-terminates before Kaggle intervenes. `_write_curve(ckpt_dir, seed, rec, kind=kind)` is called inside the per-epoch loop with `epochs_done: len(curve)`, so a stop at any point leaves every completed epoch on disk. `advance()`'s resume branch reads a prior curve when it already reaches the requested depth. The bar's own resume argument holds by construction because each epoch's shuffle is seeded from `(seed, epoch)` rather than accumulated global RNG.

**The phase order is what makes a bounded stop gradeable.** `run()` orders budget, then nulls, then extension, because the frozen bar makes the shuffled-label null a GATE. A session cut after phases one and two yields a verdict; one cut inside phase one does not, and the scorer's `budget_reached` gate returns `INDETERMINATE_BUDGET_EXHAUSTED` rather than a substantive verdict.

**The cost arithmetic, from the measured rate.** 1276.04 s/epoch. Budget phase = 10 seeds x 7 epochs = 70 epochs = **24.8 h**. Nulls = `N_NULL_SEEDS` (3) x 7 = 21 epochs = **7.4 h**. Gradeable prefix = 91 epochs = **32.3 h**. Full horizon adds 10 x 31 further ineligible epochs for **134.7 h** total. At an 8.0 h ceiling a session completes about **22 epochs**, so the prefix needs roughly **five** sessions.

**The silent-failure path.** `seed_checkpoints()` copies files from `CKPT_IN` whose `p.parent.name == "finetune_checkpoints"`, counts them, and prints only `if n:`. A zero is correct for the default empty `/kaggle/input`; it is a contradiction when `--ckpt-dataset` baked an explicit mount, and it would make a chained session recompute every seed while looking normal. Whether Kaggle's dataset mount preserves that directory name is **[unverified]**.

**Push surface, re-read not inherited.** `scripts/kaggle_push_poll.py` exposes `push <code_py> <slug> [--gpu] [--no-internet] [--dataset <slug>]`, `status`, `pull`; `--gpu` pins `machine_shape = "NvidiaTeslaT4"`; `push` fails closed on `.error` / empty `.ref`. There is no environment field in its metadata, which is why knobs are baked. Kaggle permits **2** concurrent batch GPU sessions.

### Reusable-Code Survey
- `scripts/finetune_arm_run.py::run` — the real entry point, already accepting `horizon_epochs`, `n_seeds`, `n_null_seeds` and `wall_ceiling_hours`; call it, do not fork it.
- `scripts/finetune_arm_run.py::advance` / `_write_curve` / `read_curve` — the resume contract; it needs no change for chaining.
- `scripts/kaggle/build_promoteratlas_finetune_kernel.py` — already carries `--smoke`, `--wall-hours` and `--ckpt-dataset` baking built and tested today; extend it with one more flag rather than a second builder.
- `scripts/kaggle_push_poll.py` — the fail-closed push/status/pull helper, used as-is.
- `tests/test_finetune_kernel_builder.py` — 39 tests including an AST scan of every `bar[...]` chain and the `train_one_seed` / `run` builder-asymmetry pin; extend it.
- `wiki/finetune_arm_kaggle_logs/` — preserved attempt-1 and attempt-2 logs, the evidence any new guard is checked against.

Searched: `scripts/`, `scripts/kaggle/`, `tests/`, `wiki/finetune_arm_kaggle_logs/`, `graphify-out/GRAPH_REPORT.md` (absent), `src/utils`, `src/lib`, `utils`, `lib`, `common`, `helpers` (all absent).

## Pre-Change Baseline

- **Measured rate: 1276.04 s/epoch**, Tesla T4 `sm_75`, `t4_pin_honoured: true`, torch 2.11.0+cu130, from the completed `pa-finetune-arm-v3` smoke.
- **Projection at the bar's horizon: 134.69 h**, `fits_one_session: false` against the 10.5 h ceiling the smoke carried.
- **`plan_run` refuses** at 8.0, 9.0, 10.5, 12.0 and 32.0 h; the lowest passing ceiling is **32.3 h** (horizon 7, 91 epochs).
- **`run()` on Kaggle: never executed.** Zero runs of the real entry point on the real device.
- **No run artifact and no checkpoint dir:** `wiki/finetune_arm_run*.json` is empty; `data/promoteratlas/finetune_checkpoints/` absent.
- **Arm state:** `finetuned-promoteratlas-encoder-150bp` is `planned`; `ceiling_verdict` is `blocked:untried-arm`; comparator median Spearman 0.5076 over 10 seeds at `split_frac` 0.25.
- **Tests:** 39 passing in `tests/test_finetune_kernel_builder.py`; full suite 13 failed / 8564 passed / 204 skipped, failure set recorded in `wiki/finetune_arm_prerun_baseline_2026-11-14.json`.
- **Slugs consumed:** `pa-finetune-arm-v1`, `v2`, `v3` all hold retrieved output and must never be re-pushed.
- **GPU quota consumed so far: [unverified]** — three pushes, two of which died inside two minutes.

## Verification Signal

- A Kaggle session ends with `bounded_by_wall_clock: true` and a written run artifact, proving `run()` self-terminates gracefully rather than being killed.
- A second session seeded from the first reports a **non-zero** seeded-checkpoint count and reaches a higher total epoch count than the first, proving resume works across the dataset boundary.
- A baked checkpoint dataset that seeds zero files refuses loudly instead of recomputing.
- The chain terminates with a record whose `status` is `COMPLETE` or `BOUNDED_BY_WALL_CLOCK`, carrying per-seed `frozen_budget.value` for every seed the bar declares.
- Negative control: if any session dies, no result artifact is written, the arm stays `planned`, and `ceiling_verdict` still reads `blocked:untried-arm`.

## Implementation Steps

### Step 1: Bake the horizon so a Kaggle kernel can bound its own session
Files: scripts/kaggle/build_promoteratlas_finetune_kernel.py
Depends on: none

**What changes:**
- Add a `--horizon-epochs N` build flag that bakes `BAKED_HORIZON_EPOCHS`, following the pattern of the three knobs already baked, with `None` meaning "not baked" so the env read stays live for local dry-runs.
- Pass it through to the kernel's `run(...)` call as `horizon_epochs=HORIZON_EPOCHS`, which bypasses `plan_run` and lets `wall_ceiling_hours` alone decide when the session stops.
- Emit a `horizon_divergence` block into the run artifact recording that planning was bypassed deliberately, quoting `plan_run`'s own refusal text and stating why a bounded run is consistent with the bar: the phase order exists so that a stop after budget and nulls is gradeable.
- Keep the default at `None` so an un-flagged build behaves exactly as today and still refuses — the refusal is correct for anyone who has not chosen to chain.

**Test strategy:**
- Unit: a `--horizon-epochs 38` build bakes the literal and passes it to `run(`; a build without the flag bakes `None` and the `run(` call still carries `horizon_epochs=None`.
- Non-vacuity: assert the two builds differ in exactly the baked constant and the `run(` keyword, and nowhere else.

### Step 2: Make the one silent-failure path in checkpoint seeding loud
Files: scripts/kaggle/build_promoteratlas_finetune_kernel.py
Depends on: Step 1

**What changes:**
- `seed_checkpoints()` raises `SystemExit` when `BAKED_CKPT_IN is not None` and zero files were copied, naming the mount path and saying that a chained session which seeds nothing recomputes every seed while looking normal.
- Print the seeded count unconditionally when a dataset was baked, so the chaining path has a positive signal rather than only an absence.
- Leave the default `/kaggle/input` case silent, because there a zero is the correct answer.

**Test strategy:**
- Unit: the refusal is present and conditioned on `BAKED_CKPT_IN is not None`.
- Non-vacuity: a build without `--ckpt-dataset` does not carry the condition as an unconditional raise, so the default path keeps its correct silence.

### Step 3: Extend the static harness to cover both new behaviours
Files: tests/test_finetune_kernel_builder.py
Depends on: Step 1, Step 2

**What changes:**
- Add tests for every assertion named in Steps 1 and 2, reading call sites and keyword values off the emitted body's **AST** rather than by substring — the discipline earned when a body-wide `build=None` check proved wrong because `run(build=None)` is legitimate while `train_one_seed(build=None)` is fatal.
- Add a test that the baked horizon reaches the `run(` call specifically, not merely the module constants, since a baked value that never reaches its consumer is exactly the class of defect the first two Kaggle failures were.
- Re-assert the existing invariants survive: embed set unchanged across flag combinations, every `bar[...]` chain resolving against the committed bar, the smoke's throwaway-checkpoint guard intact.

**Test strategy:**
- Unit: the file passes with a count above 39, every new assertion represented.
- Non-vacuity: for the horizon pass-through and the zero-copy refusal, reintroduce the defect in a throwaway copy and assert the test fails naming the construct.

### Step 4: Resolve the session cap and the weekly GPU quota
Files: wiki/kaggle_session_limits_2026-11-14.json
Depends on: none

**What changes:**
- Record the GPU session cap and the weekly GPU quota from Kaggle's own documentation and the account settings page, with URL and read date, not from secondary blogs. The repo's pusher attributes 12 h to CPU specifically, and public sources split 12 h against 9 h for GPU.
- Record the quota already consumed, because the gradeable prefix costs 32.3 h and a commonly-cited ~30 h weekly allowance would make the prefix exceed one week's budget — which is a scheduling fact the chain's design depends on.
- If a figure cannot be established, record `UNRESOLVED` with what was tried. An unresolved cap keeps the conservative 8.0 h session ceiling.

**Test strategy:**
- Unit: every number carries a source URL and a read date; an `UNRESOLVED` status carries `what_was_tried`.
- Manual: the recorded quota matches what the account settings page displays.

### Step 5: Verify the checkpoint-dataset round-trip before any session depends on it
Files: wiki/ckpt_dataset_roundtrip_2026-11-14.json
Depends on: Step 3

**What changes:**
- Create a minimal Kaggle dataset containing a directory named `finetune_checkpoints` holding one small file, push a kernel built with `--ckpt-dataset <slug>` and `--smoke`, and record the seeded count the kernel reports.
- Record whether the mount preserved the directory name, which is the assumption `seed_checkpoints()` matches on and which nothing has yet confirmed.
- A zero count here is a **pass of the Step 2 guard and a fail of the mount assumption**; record which of the two occurred, because they call for opposite responses.

**Test strategy:**
- Manual: the kernel log reports a non-zero seeded count, or the Step 2 refusal fires with its named path.
- Unit: assert the artifact distinguishes `mount_preserved_directory_name` from `guard_fired`, so a later reader cannot collapse them.

### Step 6: Run session one as a self-bounding chained session
Files: wiki/finetune_arm_session1_2026-11-14.json
Depends on: Step 3, Step 4

**What changes:**
- Build with `--horizon-epochs 38 --wall-hours 8.0`, push to a FRESH slug with `--gpu`, poll to a terminal state, pull.
- Record `status`, `why_status`, `bounded_by_wall_clock`, the per-seed epoch counts reached, and the total epochs completed against the 91 the gradeable prefix needs.
- Confirm the session **self-terminated** — `bounded_by_wall_clock: true` with a written artifact — rather than being killed, which is the property the whole chain rests on.
- Preserve the pulled log under `wiki/finetune_arm_kaggle_logs/` whatever the outcome.

**Test strategy:**
- Manual: `status` reaches `complete`; the artifact exists; `bounded_by_wall_clock` is true.
- Unit: assert the artifact's recorded epoch total equals the sum of per-seed `epochs_done` read independently from the pulled checkpoints.

### Step 7: Chain further sessions until the prefix is gradeable or a wall is hit
Files: wiki/finetune_arm_chain_2026-11-14.json
Depends on: Step 5, Step 6

**What changes:**
- Upload the pulled checkpoint directory as a Kaggle dataset and push the next session with both `--ckpt-dataset <slug>` and `kaggle_push_poll.py push --dataset <slug>`, since the kernel must know where to look and Kaggle must mount it.
- Assert at each link that the seeded count is non-zero and that total completed epochs **increased**; a session that seeds correctly but adds nothing is a stall, not progress, and must stop the chain.
- Continue until the record reports every seed at the frozen budget and the nulls complete, or until the quota from Step 4 is exhausted.
- Record each link's seeded count, epochs added and cumulative total, so the chain's cost is a measured series rather than an estimate.
- STOP and report rather than continuing if any link seeds zero, adds zero, or the quota binds.

**Test strategy:**
- Unit: assert cumulative epochs is strictly increasing across links and that the final record carries `frozen_budget.value` for every seed the bar declares.
- Manual: confirm the final `status` is one of the two accepted values and that `bounded_by_wall_clock` agrees with it.

## Execution Preview

Wave 0 (2 parallel):  Step 1 — bake the horizon, Step 4 — resolve cap and quota
Wave 1 (1):           Step 2 — make the silent seed path loud
Wave 2 (1):           Step 3 — extend the static harness
Wave 3 (2 parallel):  Step 5 — verify the dataset round-trip, Step 6 — run session one
Wave 4 (1):           Step 7 — chain to a gradeable prefix

Critical path: Step 1 → Step 2 → Step 3 → Step 5 → Step 7 (5 waves, as computed by the toolkit)
Step 6 sits in the same wave as Step 5 and is equally load-bearing for schedule even though the
computed path runs through Step 5.
Max parallelism: 2 agents

Steps 1 and 2 share the builder and are therefore sequenced; Steps 5 and 6 are parallel because they touch different artifacts and consume separate Kaggle sessions, and Kaggle permits two concurrent batch GPU sessions.

## Risk Flags

- Severity: high — **The gradeable prefix costs 32.3 h of GPU, which may exceed one week's entire quota [unverified].** The session cap is not the binding constraint; the weekly allowance is. Step 4 resolves it, and if the quota is ~30 h the chain cannot complete inside one week, which is a scheduling and scope question rather than a technical one.
- Severity: high — **`run()` has still never executed on Kaggle.** Step 6 is its first execution, and a defect in the phase ordering or the artifact write would surface only then. The mitigation is structural rather than preventive: per-epoch checkpoint writes mean even a late crash preserves every completed epoch, so the cost of such a defect is one session's artifact, not one session's compute.
- Severity: medium — **Supplying the horizon bypasses a guard whose refusal text says the run "must not start" [grounded].** That refusal is correct for an unchained run and wrong for a chained one, which is why Step 1 records the divergence in the artifact rather than silently passing a parameter. A reader who sees only the bypass and not the note would reasonably conclude the bar was violated.
- Severity: medium — **The dataset mount directory-name assumption is [unverified]** and `seed_checkpoints()` currently fails by passing. Step 2 makes it loud and Step 5 tests it before Step 7 depends on it; until Step 5 reports, every chaining claim is provisional.
- Severity: medium — **A bounded stop inside the budget phase is not gradeable.** The scorer's `budget_reached` gate returns `INDETERMINATE_BUDGET_EXHAUSTED` until every seed reaches 7 epochs, so sessions one through roughly four produce no verdict by design. This is expected and must not be read as failure, which is why Step 7 tracks cumulative epochs rather than per-session verdicts.
- Severity: medium — **The 8.0 h ceiling is conservative against an unresolved cap** and costs throughput if the true cap is 12 h. Raising it is cheap once Step 4 reports; starting low is the asymmetric choice because a ceiling above the true cap loses the artifact.
- Severity: low — File overlaps: Steps 1 and 2 share the builder and are sequenced; no wave holds two steps touching one file.
- Severity: low — **Three slugs are already consumed.** `pa-finetune-arm-v1` through `v3` all hold retrieved output; every push in this plan needs a fresh slug.
- Severity: low — No code generator writes any committed file; the emitted kernel goes to a scratch path.

## Open Questions

- If Step 4 shows the weekly quota cannot cover 32.3 h, does the arm run across multiple weeks, drop to a smaller `n_seeds`, or get re-routed again? Each changes a frozen pre-registration and is an authority call.
- Should the chain target the gradeable prefix (91 epochs, 32.3 h) and stop, or continue to the full 38-epoch horizon (134.7 h)? The tail is already declared INELIGIBLE by the bar, so the prefix yields the verdict and the tail only closes the budget-limited reading.
- Should checkpoints be uploaded as a dataset after **every** session as insurance, or only when the next link needs them? Uploading always costs storage and time; uploading lazily risks losing a session to an unrelated failure.
- `C:\my_skills\.venv` (232 MB, created accidentally by running `uv` with the cwd inside the skills repo) is still present. Deleting it is a recursive delete and hard-pauses. Confirm or leave it.

## Verification

1. `tests/test_finetune_kernel_builder.py` passes with more than 39 tests, and the horizon pass-through and zero-copy refusal have each been shown to fail when their defect is reintroduced.
2. A `--horizon-epochs 38` build reaches the emitted `run(` call with that value; an un-flagged build still carries `horizon_epochs=None` and still refuses.
3. `wiki/kaggle_session_limits_2026-11-14.json` carries a source URL and read date per number, or `UNRESOLVED` with `what_was_tried`.
4. `wiki/ckpt_dataset_roundtrip_2026-11-14.json` distinguishes `mount_preserved_directory_name` from `guard_fired`.
5. Session one's artifact reports `bounded_by_wall_clock: true` with a written record, proving graceful self-termination.
6. The chain artifact shows strictly increasing cumulative epochs and a non-zero seeded count at every link.
7. Frozen five, `uv.lock` and `wiki/finetune_arm_acceptance_bar.json` byte-identical to `wiki/finetune_arm_prerun_baseline_2026-11-14.json`; full-suite failures a strict subset of its recorded set.
8. If the chain stops short, no result artifact exists, the arm reads `planned`, and `ceiling_verdict` reads `blocked:untried-arm`.

## Save-time amendments

Captured at: 2026-11-15
Source: `/save-plan` arguments

These are AUDIT NOTES ONLY. `/execute-plan` reads executable work exclusively from
`## Implementation Steps`; nothing in this section is an executable instruction. An amendment that
changes a step's file list, dependencies or structure must be folded back in via `/technical-plan`
before execution, and the re-validated toolkit footer below is what reports whether that happened.

- Supersedes the unsaved smoke-hardening plan, whose premise that the arm cannot start at a session-sized ceiling was falsified — plan_run is bypassed when horizon_epochs is supplied; chained bounded sessions adopted, 8.0 h ceiling held conservative against an unresolved GPU session cap

<!-- toolkit: check=clean waves=clean gate=fired:open-questions,unverified,severity-high -->
