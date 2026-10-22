# promcalc regeneration — runbook and measured cost correction (2026-10-22)

## DECISION NEEDED

Nothing is blocked right now. One thing will need a call when the run lands, and it is named here so
it is not discovered at the end:

- **If stage C's frozen reconcile gate FAILS, that is a real finding and must not be tuned away.**
  The gate compares the regenerated table against `PREREGISTERED["reconcile"]` — `max_tx_median`
  expected **0.3968 ±0.02** and `gc_seed0` **0.3000 ±0.01**. A failure most likely means
  `promoter_calculator` version drift changed the numbers, in which case the honest outcome is
  `INDETERMINATE_RECONCILE_FAILED`, reported as such. Re-sizing a frozen bar after seeing the result
  is the failure this project has twice had to undo. **Editing `PREREGISTERED` is a user-authority
  call, not an executor one.**

## CORRECTIONS

Two beliefs a reader may be carrying are wrong, and both came from this repo's own records.

### 1. The recorded cost figure is NOT reproducible here. Do not quote ~56 min.

`wiki/promcalc_decomposition_arm_2026-10-13.md` records **0.26 s/tile single-core** and **785
tiles/min on 7 workers**, published as ~74 min for 44,106 tiles. Measured on this host today:

| quantity | measured 2026-10-22 | recorded |
|---|---|---|
| per tile, BLAS pinned to 1 thread | **8.7 core-seconds** | 0.26 s single-core |
| per tile, wall, unpinned single process | ~3.0 s | — |
| constructor cost | **0.000 s** (free) | — |
| 7 workers, unpinned | ~28 tiles/min (proxy; see correction) | 785 tiles/min |
| 7 workers, 1 thread each | **51.9 tiles/min** (~12.7 h) | 785 tiles/min |

785 tiles/min would require roughly **0.5 core-seconds per tile** against the 8.7 measured — a ~17x
gap that no thread setting explains. Either the 2026-10-12 run used a materially faster
implementation, or the figure was derived differently (for instance from a partial-resume window).
**This is left UNRESOLVED rather than reconciled by argument.** What is certain is that the recorded
figure caused a real planning error: the technical plan written this morning budgeted "~1 h" and
called it "the one real cost", on the strength of a number that is off by more than an order of
magnitude.

### 2. The sklearn-isolation rationale was called "refuted". That was too strong.

Stage A's and stage B's docstrings say `promoter_calculator` pulls sklearn 1.9.1 while the project is
pinned at 1.8.0, and a project-env `uv pip install --dry-run` was read as refuting it (it resolved
one package, leaving sklearn at 1.8.0). Measured today in a **fresh** environment,
`promotercalculator==1.2.4` genuinely **does** pull `scikit-learn==1.9.1` (plus numpy 2.4.6, scipy
1.17.1, biopython 1.88). The dry-run stayed at 1.8.0 only because the project's existing pin already
satisfied the resolver. **Both facts are true, and the original rationale was closer to right than my
correction was.** The separate environment is correct for both reasons: the fresh-resolve conflict is
real, *and* `uv.lock` is a frozen reproducibility unit so the dependency cannot be added to it.

## What was run

All three stages, from the restored repo-local MPRA cache. The promcalc directory is
`data/promcalc/`, gitignored via `/data/*` (confirmed with `git check-ignore --no-index -v`).

```bash
# Stage A -- export tiles from the local GSE144621 cache (no network)
uv run python scripts/promcalc_export_tiles.py \
  --cache-dir data/mpra --out-dir data/promcalc --no-download
# -> n_filtered_tiles 44106, n_parsed 46713, replicate_ceiling 0.941
#    Both of Stage A's OWN pinned asserts passed (EXPECTED_TILES, EXPECTED_CEILING),
#    so the re-fetched substrate reproduces the published one exactly.

# a bounded smoke slice, written to *_smoke.* so it cannot collide with the full run
uv run python scripts/promcalc_export_tiles.py \
  --cache-dir data/mpra --out-dir data/promcalc --no-download --limit 40

# Stage B smoke -- the gate that paid for itself (24 s)
"C:/Users/Laptop Portable/tools/promcalc_env/Scripts/python.exe" \
  scripts/promcalc_decompose.py --in-dir data/promcalc --smoke --workers 7

# Stage B full -- detached, resume-safe, thread-pinned
#   Start-Process -FilePath <repo>\run_promcalc_decompose_detached.bat \
#     -WorkingDirectory <repo> -WindowStyle Hidden
#   add -ArgumentList "--resume" on any relaunch
```

### The Stage-B interpreter

`uv venv "C:/Users/Laptop Portable/tools/promcalc_env" --python 3.11` then
`uv pip install promotercalculator` into it. **The documented interpreter
`D:/tools/promcalc_env/Scripts/python.exe` is gone** — it was on the volume that detached, so a
reachable path was needed.

Two naming traps, both already recorded and both still live:

- **install name ≠ import name** — `pip install promotercalculator`, `import promoter_calculator`.
  Installing the import name 404s.
- **`uv venv`'s `python.exe` is a 46 KB trampoline.** Its `pyvenv.cfg` points `home` at the
  uv-managed CPython, and `multiprocessing` spawn therefore launches the **base** interpreter for the
  Pool children, not the venv one. The package still resolves (the venv's `sys.path` travels in the
  spawn preparation data, which the clean smoke run confirms), but **a liveness check that greps for
  the venv path will not match the workers** — this bit the watchdog's first draft.

## The two version breaks, and why the smoke gate is the hero here

The plan's single flagged `[unverified]` risk was that the `promoter_calculator` version used for the
2026-10-12 run is recorded **nowhere**. It materialized immediately, as two independent breaks that
each killed every tile:

1. **`from promoter_calculator import Promoter_Calculator` raises ImportError.** 1.2.4's
   `__init__.py` re-exports only the lowercase `promoter_calculator` *function* from `wrapper`; the
   *class* now lives in the `promoter_calculator.promoter_calculator` submodule.
2. **Each per-TSS value is a `PromoCalcResults` object, not a dict.** It defines `__getitem__`,
   `__setitem__` and `keys()` but **not `get`** — so the worker's five `d.get(field, "nan")` calls
   each raise `AttributeError`. Because `decompose_one` catches broadly, every site turned its whole
   tile into a *failed* tile rather than a visibly wrong one.

**The class itself is UNCHANGED, and that is the load-bearing finding.** 1.2.4's own `wrapper.py`
constructs it with the same `Promoter_Calculator(threads=, verbosity=)` and reads the same
`output()['Forward_Predictions_per_TSS']` key this worker reads. So these are an import-path fix and
a container-shape fix, **not a model swap** — which is what lets stage C's frozen reconcile anchors
remain the right acceptance test rather than a comparison against a different model.

**CORRECTED 2026-10-22, by actually re-running the gate three ways instead of citing the first run.**
"Caught both breaks in 24 s" was true of break 1 ONLY. Measured on 40 tiles at 4 workers, with the
breaks deliberately re-introduced into throwaway copies:

| case | tiles ok | errored | sites/tile | error class | elapsed |
|---|---|---|---|---|---|
| current code | 40 | 0 | 73.0 | -- | 165.8 s |
| break 1 re-introduced | 0 | 40 | -- | `ImportError` | **21.2 s** |
| break 2 re-introduced | 0 | 40 | -- | `AttributeError` | **173.5 s** |

Break 1 fails before any work. **Break 2 fails only AFTER the expensive solve** -- each tile spends its
full ~8.7 core-seconds and then dies formatting the result -- so on the full 44,106 tiles it would have
consumed the ENTIRE ~12.7 h run and returned an empty table. That makes the gate more valuable for
break 2 than the original claim, not less; the stated 24 s figure was simply wrong for that case.

The clean case takes 165.8 s here rather than the 58 s recorded above because the box was running the
real 7-worker job at the time and this used 4 workers. Same work, contended machine.

**The gate is NON-VACUOUS, which is the claim that matters** -- it passes clean code and fails every
tile under either break, with the correct error class each time. A gate that only ever passes proves
nothing.

A third shape difference turned out to be harmless and is pinned so it stays that way:
`hex35_position` comes back as a real list `[25, 31]`, not the string the docstring describes.
`_int_start` stringifies before parsing, so it handles both — 0 of 2,920 rows failed to parse.

## Smoke-gate result after the fix

| check | result |
|---|---|
| rows / tiles | 2,920 / 40 = **exactly 73.0 sites/tile** (the published figure) |
| header | equals the committed `SITE_FIELDS` |
| NaN or blank rows | **0** |
| `hex35_start` parse failures | **0** |
| `spacer_len` range | 15–19 (the canonical σ70 17±2) |
| errors | **0** |

## The thread-pinning finding

The first full launch ran at 28 tiles/min. The diagnosis was measured, and it inverted my own
assumption: the solve is **BLAS-parallel**, so the 3.0 s/tile I had measured in a single process was
quietly using about three cores. With 7 worker *processes* each also spawning its own thread pool, 8
cores were oversubscribed several times over and most of the machine went into context switching.

Pinning `OMP_NUM_THREADS` / `OPENBLAS_NUM_THREADS` / `MKL_NUM_THREADS` / `NUMEXPR_NUM_THREADS` to 1
makes process-level parallelism the only parallelism, which is the level that actually scales here.

**28 → 65 tiles/min, 2.34x** as originally measured — but BOTH arms of that A/B were derived from gzip file size, and the pinned figure is now known wrong (see the correction below). Measured on the real workload — same tiles, same 7
workers, same launcher, only the environment differing — not on a synthetic benchmark. The restart
cost nothing, because `--resume` recovered 495 completed tiles.

**Reusable rule: when a worker pool underperforms its single-process timing, suspect that the
single-process timing was already parallel.** The per-worker rate, not the aggregate, is the
diagnostic.

## Unattended supervision

`scripts/promcalc_watchdog.ps1`, because four long runs died unattended on this host in one week —
three to session teardown, one to a USB paging fault — and the previous full run of *this exact
stage* was killed by teardown at 88.4%. Resume-safety is worth nothing if nothing re-launches.

It distinguishes three states and acts on one: **complete** (coverage reached 44,106 — stop, never
re-launch), **died with retries left** (re-launch with `--resume`), **died with retries exhausted**
(stop and leave the evidence). Two safety properties are deliberate:

- **A retry cap.** An unbounded re-launcher on a job that fails instantly is a hot loop that burns
  the machine all night and buries the real error under thousands of log lines.
- **A coverage check.** A completed run can never be mistaken for a dead one and restarted over its
  own finished output. Undetermined coverage **resumes** rather than declaring completion — failing
  toward a resume is recoverable, declaring a partial run finished is not.

Coverage is read only while the job is **idle**, because decompressing a growing gzip during a live
run is the recorded trap that turned a ~1 h job into a 5.0 h projection.

Both detectors were tested against the live job before being trusted: liveness positive while
running, coverage numeric, and a bogus-directory negative control returning 0 rather than a false
completion. Two launch defects were found and fixed in the process — an indented here-string
terminator (a PowerShell parse error, so the first watchdog silently never started and wrote no log)
and an unquoted `-File` path splitting on the space in "Laptop Portable".

## Honest limits

- **The run was still in flight when this was written.** Every number above is from stage A, the
  smoke gate, and two rate windows on the live stage-B run. **No completion, no coverage of
  44,106/44,106, and no stage-C verdict is claimed here.**
- **The reconcile gate has not run.** Whether the regenerated table reproduces the published
  substrate is the open empirical question, and it is exactly what the gate exists to answer.
- Rate is derived from output-file growth against a bytes-per-tile constant measured on the 40-tile
  smoke (3,654 bytes/tile). Compression ratio can drift with content, so the ETA is approximate.
- **RATE CORRECTED 2026-10-22: 51.9 tiles/min, not 65, and the ETA is ~12.7 h, not ~11 h.** The
  job prints EXACT counts to its own log every 2,000 tiles and I never read them -- I inferred tile
  counts from gzip file size against a bytes-per-tile constant taken from the 40-tile smoke. The log
  says `2000/43611 tiles | 37.9 min` then `4000/43611 tiles | 77.1 min`, i.e. 51.0 tiles/min between
  those marks and 51.9 overall, with the job's own eta at 763 min. **The 5-minute file-size window I
  sampled caught a write/flush burst as the pool spun up, not steady state** -- a short size sample on
  a buffered gzip stream is not a rate measurement.
- **This is the SAME error class as the wiki figure corrected at the top of this document**, committed
  in the same session one commit later: a derived proxy quoted as a measurement when an exact counter
  was already being written to disk. Read the log's progress lines.
- **The 2.34x speedup ratio does NOT rest on those proxy windows** (both arms were proxy-derived, so
  the ratio inherits their lumpiness). It is corroborated independently by the mechanism: 8.7
  core-seconds/tile over 7 workers predicts ~48 tiles/min, and the log-measured 51.9 sits right at
  that hardware efficiency, while unpinned 8-core oversubscription accounts for the remainder. Treat
  the DIRECTION as solid and the exact multiple as approximate.
- 8 workers instead of 7 would be expected to gain roughly another 1.3 h. **Not tried** — the job was
  already running efficiently and a second restart was not worth the progress and attention.
- The 2026-10-12 shard set (88.4% complete) is **not recoverable** — it was on the detached volume.
