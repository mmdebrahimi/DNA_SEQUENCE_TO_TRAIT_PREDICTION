# Tile-rung selection value — 2026-10-14

**VERDICT: `SELECTION_IS_THE_BOTTLENECK`**

at k=10, 2/10 seeds pass bar A, 10/10 pass bar B | seed 0 max_tx 0.498 vs gc 0.265 (margin needed 0.100; A n B y); seed 1 max_tx 0.480 vs gc 0.338 (margin needed 0.100; A n B y); seed 2 max_tx 0.471 vs gc 0.285 (margin needed 0.106; A n B y); seed 3 max_tx 0.457 vs gc 0.274 (margin needed 0.100; A n B y); seed 4 max_tx 0.513 vs gc 0.285 (margin needed 0.100; A y B y); seed 5 max_tx 0.496 vs gc 0.305 (margin needed 0.100; A n B y); seed 6 max_tx 0.482 vs gc 0.298 (margin needed 0.100; A n B y); seed 7 max_tx 0.457 vs gc 0.299 (margin needed 0.110; A n B y); seed 8 max_tx 0.517 vs gc 0.326 (margin needed 0.100; A y B y); seed 9 max_tx 0.476 vs gc 0.273 (margin needed 0.100; A n B y)

## What this measures

Not a correlation — a CHOICE: propose k candidates, score them, keep the best. `headroom_captured` is the fraction of the ACHIEVABLE selection value the scorer realises at k=10, normalised by the measured best of the same candidate set.

Substrate: GEO GSE144621 peak tiles, 44106 tiles, leave_peak_out(frac=0.25), ceiling 0.941

## Cross-seed mean headroom captured

| scorer | k=1 | k=2 | k=5 | k=10 | k=20 |
|---|---|---|---|---|---|
| `conv` | — | 0.3140 | 0.3378 | 0.3634 | 0.4197 |
| `decomposition` | — | 0.3955 | 0.4360 | 0.4693 | 0.5330 |
| `gc` | — | 0.3090 | 0.2956 | 0.2948 | 0.2952 |
| `max_tx` | — | 0.3880 | 0.4468 | 0.4847 | 0.5366 |

## Against the designed-grid anchor

- GSE108535 designed sigma70 grid, leave-element-out: learned **0.7937**, gc **0.1328** (verdict `ORACLE_READY_AS_SCORING_HEAD`)
- cross-AXIS means over 5 element axes; this run reports cross-SEED means over 10 leave-peak-out splits. Different pooling units -- compare with that in mind.

## Null band (scorer-realisation spread, per seed and k)

A headroom without this band is uninterpretable: `trials` reduces Monte-Carlo error in estimating ONE scorer's headroom and says nothing about how far it sits from zero.

- seed 0 k=10: mean -0.0053, stdev 0.0420, p95 +0.0644
- seed 1 k=10: mean +0.0080, stdev 0.0478, p95 +0.0760
- seed 2 k=10: mean -0.0070, stdev 0.0530, p95 +0.0909

## Honest limits

- UPPER BOUND on the loop's selection value: the candidate pool is MEASURED peak tiles, while a generator's candidates sit off this distribution where no label exists.
- A DIFFERENT upper bound from the designed-grid run's, not a stricter one -- natural peak-enriched sequence versus a fully-crossed design.
- Both bars are ASSERTED, not derived; the per-trial ceiling arm that normalises them IS derived.
- Cross-seed means pool over leave-peak-out splits of ONE experiment; the designed-grid anchor pooled over element axes. Different units.
- Peak tiles are selected from ChIP peaks, so this pool is enriched for real promoters relative to the genome -- selection is easier here than on random sequence.
