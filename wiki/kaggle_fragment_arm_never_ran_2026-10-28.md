# The fragment-pyramid Kaggle run never happened — and the projection defect was not why

**Date:** 2026-10-28
**Status:** correction of a session-held belief; the arm's recorded state was never wrong
**Arm:** `pyramid-4scale-2.22decades-499bp-fragments` (one of three `planned` arms blocking the GLM ceiling)

## The correction

A whole session was spent treating the fragment-pyramid Kaggle result as **computed but unread** —
"`COMPLETE`, 58 minutes, clean exit, artifacts never pulled" — with the work item being *fetch the
artifacts once the HTTP 429 clears*. That belief was false in every part that mattered.

Measured 2026-10-28 through the Kaggle API:

| probe | result |
|---|---|
| `kernels_status("emanueleebrahimi/notebookd7f1ad5a41")` | **`ERROR`** (empty `failure_message`) |
| `kernels_list_files(...)` | **0 files** |
| its `kernel-metadata.json` | `kernel_type: "notebook"`, **`enable_gpu: false`**, `machine_shape: "None"` |
| `current_version_number` | `0` |
| `kernels_list(mine=True)`, 30 kernels | no `fragment`/`pyramid` slug exists; the only 2026-10-28 run is that notebook |

So there were never any artifacts to pull. **The four failed pull attempts and the rate-limit were
chasing files that do not exist**, and the 429 — read at the time as the obstacle between us and a
result — was an obstacle in front of nothing.

The deeper fact is the metadata: `enable_gpu: false`. **The arm has never run on a GPU.** Every
projection, break-even figure and "does it fit a 9 h session" argument in this project assumed a T4.

## What was right, and why that matters

Nothing in the repo asserted the run had completed:

- `arms_ledger` has `pyramid-4scale-2.22decades-499bp-fragments` as **`planned`**, and
  `blocking()` has been naming it as a ceiling blocker the whole time.
- `wiki/decisions-log.md` records *"Step 6 = attended Kaggle real-run DEFERRED per plan"*.

The false belief lived **only** in the Soraya resume checkpoint and in session narrative. The
code-owned arm state was never corrupted by it — which is the R5/R6 rails working exactly as
designed: a narrative drifted, and the thing that gates the ceiling did not move with it.

The operational lesson is about where to look. A status poll is cheap and authoritative; a
remembered status is neither. The remembered one survived a compaction boundary and was then
restated in a checkpoint, which is what gave it the appearance of a measurement.

## The projection defect (real, fixed, and NOT the blocker)

`scripts/kaggle/fragment_pyramid_kernel.py` sized the full run from **one** timing point:

```python
per_fit = smoke_min / max(1, 2 * SMOKE_SEEDS)   # total wall-clock / 6 fits
projected_h = per_fit * (321123 / SMOKE_LIMIT) * (2 * FULL_SEEDS) / 3600.0
```

One timing point cannot separate **fixed** cost (import torch, build the model, parse the frames —
paid once) from **marginal** cost (per row, per fit). Dividing total wall-clock by the fit count
charges the fixed cost to every fit, and the projection then multiplies it by
`scale × 20/6 = 80.3 × 3.33 ≈ **268×**`. Roughly 30 s of setup becomes **+2.2 h** of projected run.
`smoke_min` also held **seconds**, despite the name.

**Fixed** by adding a second, smaller timing point and solving for both unknowns:

```
t(N, S) ≈ fixed + rate · N · (2·S)
rate = (t_smoke − t_cal) / (w_smoke − w_cal)      w = N · 2S  (row-fits)
projected = fixed + rate · w_full                 w_full = 321,123 × 20 = 6,422,460
```

`CAL_LIMIT = 1000, CAL_SEEDS = 1` is the second point. `rate ≤ 0` (the smaller run not faster than
the larger — the signal below timing noise) falls back to the one-point number and **labels itself
`DEGRADED`** rather than projecting from a negative rate. `smoke_min`/`full_min` renamed
`smoke_sec`/`full_sec`, because the unit belongs in the name. Both projections now print side by
side, so the 268× inflation stays visible instead of being silently corrected away.

**Checked against the real recorded measurement** rather than asserted: the CPU smoke was 5.1 min at
N=4,000 / 3 seeds, and the one-point formula reproduces the recorded **22.7 h** exactly; the
two-point solve gives **20.5 h**, recovering `fixed = 30 s` and `rate = 11.5 ms` per row-fit.

**Scope — do not over-read this fix.** It removes 2.22 h of phantom projection and it changes
**no verdict**: on CPU both projections exceed the 8 h bar (correctly — the run genuinely does not
fit), and at any plausible T4 speedup (10–30×) **both** say RUN:

| marginal speedup | two-point | one-point |
|---|---|---|
| 10× | 2.06 h → RUN | 4.28 h → RUN |
| 15× | 1.38 h → RUN | 3.60 h → RUN |
| 30× | 0.69 h → RUN | 2.91 h → RUN |

So the projection defect **was never what stood between this project and a verdict.** It is a real
268× error worth fixing, and it bites near the 8 h boundary, but the blocker was `enable_gpu: false`.
Reporting the fix as the unlock would be the same class of error as the belief it replaces.

## What actually unblocks it

A notebook's accelerator lives in a **settings panel a human must click**, and nothing in a
transcript reveals that it was never clicked. A **script kernel pushed through the API** pins the
accelerator in `kernel-metadata.json`:

```bash
uv run python scripts/kaggle_push_poll.py push \
    scripts/kaggle/fragment_pyramid_kernel.py glm-fragment-pyramid-t4 --gpu
```

Pushed 2026-10-28 and confirmed by reading its metadata back: `kernel_type: "script"`,
`enable_gpu: True`, **`machine_shape: "NvidiaTeslaT4"`**, `enable_internet: True`. It reached
`RUNNING` and stayed there, which is itself informative — both cheap guards (no-CUDA → exit 3;
`verify_pinned_tree()` → exit 3) fire within seconds, so passing them means real CUDA and a valid
pinned tree.

**This is the first time the arm has run on the configuration the analysis assumed.** Its verdict is
not in this memo; the arm stays `planned` until a scored artifact is actually read.

## Traps worth carrying

1. **A re-push destroys the previous version's output.** Kaggle serves output from the latest
   version only, so the window to read artifacts closes the moment anything re-runs that slug. Pull
   before pushing again.
2. **Pin the accelerator in metadata, never in a settings panel.** The panel is invisible to every
   downstream reader of the run.
3. **A remembered status is not a status.** `kernels_status` costs one call.
4. **This environment collapses a doubled backslash inside a quoted heredoc**, so a patch anchored
   on a line containing a backslash escape matches 0× and reports a missing anchor rather than a
   wrong one. Write patch scripts with the Write tool, or keep every anchor backslash-free.
