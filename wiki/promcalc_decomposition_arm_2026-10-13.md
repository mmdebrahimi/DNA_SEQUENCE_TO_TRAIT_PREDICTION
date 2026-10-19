# The mechanistic model's whole decomposition adds NOTHING over the one scalar reduced from it

**2026-10-13.** Machine-readable twin: `wiki/promcalc_decomposition_arm_2026-10-13.json`. Stages:
`scripts/promcalc_export_tiles.py` (project env) -> `scripts/promcalc_decompose.py` (ISOLATED
`D:/tools/promcalc_env`) -> `scripts/promcalc_decomposition_arm.py` (project env).

## Headline

| arm | median Spearman (10 registered seeds, leave-peak-out) | features |
|---|---|---|
| `gc` (the gate's own `ridge_gc`) | 0.3112 | 1 |
| **`max_tx`** (the 2026-10-12 incumbent, refit here) | **0.3968** | 1 |
| `decomposition` (this arm) | 0.3925 | 55 |
| feature-permutation control | 0.0342 | 55 |
| replicate ceiling | 0.941 | — |

**Verdict `NO_GAIN`.** Paired against the scalar: **−0.0040 at the median, CI95 [−0.0111, +0.0016]
(includes zero), 2/10 seeds.** The Promoter Calculator returns ~73 candidate promoters per 150 bp tile,
each decomposed into nine component free energies plus an initiation rate — ~730 numbers, of which the
incumbent used **one**. Those ~730 carry nothing the strongest-promoter rate does not already carry.

It is not a failure of the features to be informative: the decomposition **beats GC** by +0.0804,
CI95 [0.0698, 0.0926], **10/10 seeds**. It just cannot beat `max(Tx_rate)`.

**What this buys:** the biophysical model is already maximally compressed for this task. The
0.3968 → 0.941 gap is **not** closable by using more of this model's output, so the next lever has to
come from somewhere else. A negative that forecloses a direction cheaply is worth more than another
0.004.

## Both baselines were reproduced INSIDE this pipeline, not quoted

This is the whole reason the comparison means anything — a rich-feature fit judged against a number
carried over from another script confounds the feature set with the pipeline.

| gate | expected | got |
|---|---|---|
| `gc` seed 0 | 0.3000 ± 0.01 (`glm_encoder_gate` anchor) | **0.3000** |
| `max_tx` median | 0.3968 ± 0.02 (`promcalc_vs_learned_2026-10-12`) | **0.3968** |

Refit here on all 10 seeds, `max_tx` beats GC by **+0.0865, CI95 [0.0784, 0.0935], 10/10** — an
independent reproduction of 2026-10-12's biophysical win, and on more seeds than that run had (it
logged GC on 5).

## Honest limits

- **The control passes but is not tightly null** (median +0.0342 against a 0.05 bar; per-seed range
  −0.113 to +0.085). With 55 standardised features, a random coefficient vector still forms a
  combination that correlates with the target by chance. **Harmless for THIS verdict because the
  result is negative** — but a *positive* result near the +0.02 margin would have been within one
  control-median of noise, and should not have been believed without a tighter control.
- **One estimator.** Ridge(alpha=1.0) on 55 standardised features. A non-linear model over the same
  features is untested; the independent 2026-09-30 finding that linear was the *ceiling* and not the
  floor on a different substrate is a reason to expect little, not a measurement here.
- **One rung, one condition.** 150 bp endogenous peak tiles, peak-tile readout. Says nothing about the
  fragment rungs or any designed library.
- **Feature choice is mine.** Three families (winning site's breakdown + geometry, distribution over
  all candidates, candidate count). A different reduction of the same per-site table is not excluded.

## Four defects found by RUNNING this, none by reading it

The per-site table is now **persisted** (`D:/dna_decode_cache/promcalc/decomposition*.csv.gz`,
3,219,801 rows), so these were paid once and future mechanism arms are a cheap refit.

1. **A 20-tile smoke extrapolated the sweep at 2,095 min against a published 74.** `chunksize=16` put
   16 tiles on one worker and nothing on five of seven. *A small smoke measures pool spawn and load
   imbalance, not throughput.*
2. **Then the honest-looking ETA was ALSO wrong — by 5.7x, and my own monitoring caused it.** A
   90-second incremental count gave 138 tiles/min, so 5.0 h. The resume ran at **785 tiles/min**
   (5,132 tiles in 6.5 min), consistent with the published 74 min for the full set. The difference is
   that measuring progress meant decompressing a growing 140 MB gzip repeatedly while 7 workers
   competed for 8 cores. *Observing the sweep slowed the sweep, and the slowed rate was then quoted
   as the sweep's rate.*
3. **A tile with SOME of its rows read as complete.** The killed run left tile 38,973 mid-write, and
   `completed_tiles` keyed on "has at least one row", so a resume skipped it — 44,106/44,106 while one
   tile carried a partial, partly-garbage site list. Fixed by dropping the highest tile_id of any
   truncated shard. **This is the HDF5 embedding-cache trap on a different substrate**
   (`populate(skip_existing=True)` skipped at gene level, so a mid-flush crash left a partial gene set
   that looked like a complete strain): *verify at the unit the consumer aggregates over.*
4. **A `decomposition*` glob also matched `decomposition_smoke.csv.gz`**, merging the 20-tile smoke
   output into the real substrate. Harmless only because the smoke tiles are a prefix of the real
   ones. Shards are now named explicitly.

**A guard premise that was measured FALSE.** A first version refused above "one malformed row per
file", reasoning that truncation can only damage the tail row. The real interruption produced **2,199
copies of a single mangled row** (a float emitted as `999999999993236`) — a cut deflate stream
replays its sliding window. The threshold was arbitrary and the premise wrong; completeness is now
decided **only** by coverage against the tile list, which does not depend on whether a file closed
cleanly. Truncation and missing coverage are different questions, and an earlier version conflated
them so hard that a successful resume would still have refused forever.
