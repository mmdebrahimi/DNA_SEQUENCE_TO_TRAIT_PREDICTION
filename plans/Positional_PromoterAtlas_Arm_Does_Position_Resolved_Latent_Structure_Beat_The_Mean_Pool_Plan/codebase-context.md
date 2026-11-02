**The reduction seam, read on the real model during research.** `atlas_embed.masked_mean_pool` calls
`model(X, return_latent=True)`, which returns a 3-tuple whose third element `z` has shape
`(B, 128, w)` — verified live as `(2, 128, 150)`. The incumbent is literally
`z[:, :, :L].mean(dim=2)`. So a positional arm is a different REDUCTION over an identical forward
pass, which is what makes the whole ladder cost one pass.

**Geometry divides evenly.** 150 positions / {2, 5, 10} bins = {75, 30, 15}, all exact, so a
bin-mean ladder is nested and the global mean is EXACTLY recoverable from equal-size bin means. That
recovery is the parity gate: if the derived global mean does not reproduce the committed
`data/promoteratlas/tile_embeddings_native150_n44106.npy`, the forward pass has drifted and no
number may be read.

**The comparator already exists locally and is reconciled.**
`data/promoteratlas/local_reconcile_2026-11-01.json` holds `verdict: RECONCILED`, `mismatched: []`,
the atlas per-seed values (median 0.5076: 0.5461, 0.5499, 0.4983, 0.4986, 0.5075, 0.5163, 0.5078,
0.4933, 0.5016, 0.5210), the GC per-seed values (median 0.3112) and `n_test_per_seed` for all ten
seeds. The paired comparison needs no re-run of the incumbent.

**Split-equality has a committed surface.** `data/promcalc/tile_predictions/` holds `seed-00..09`
plus a manifest declaring `split: leave_peak_out(frac=0.25)`, `n_tiles: 44106`, and scorers
`conv/decomposition/gc/max_tx`. The incumbent arm's bar checks split equality FIRST of all, and this
plan inherits that ordering.

**Scoring protocol is fixed and must not be re-derived.** `atlas_embed.ridge_spearman` mirrors
`expression.fit_predict_ridge` exactly — `StandardScaler` on train, `Ridge(alpha=1.0)`, fit on
`log1p(expression)`, score Spearman against RAW expression — and it already carries a `shuffle_seed`
parameter for the shuffled-label null.

**Substrate facts measured during research, not assumed.** 44,106 tiles at 150 bp over 2,433 peaks;
strand is near-balanced at 21,821 `+` / 22,285 `-`; only 576 of 22,285 minus sequences (2.6%) have
their reverse complement among the plus sequences and ZERO of those share a peak, so plus/minus
tiles are DISTINCT genomic windows and not two orientations of one element. The orientation
convention of `seq` is therefore unresolved by the data and is handled as a declared empirical
variant rather than a guess.

**Resource facts, measured now.** Free RAM is **5.08 GB** while the K562 replication runs, against
10.5 GB recorded in October — and this project already has a registered constraint that declared the
kmer6 rung infeasible at a 9.8 GB design matrix, so a RAM pre-check is house practice rather than
caution. Local throughput for this encoder was measured at 42 sequences/second, so one pass over
44,106 tiles is ~17.5 minutes when uncontended. The clone is 686 MB at
`trained_weights/base_model/promoteratlas-base.pt` and lives in a **volatile temp directory**.

### Reusable-Code Survey

- `dna_decode/glm/atlas_embed.py` — `load_atlas`, `audit_checkpoint`, `ridge_spearman`,
  `masked_mean_pool`. The new readout reuses the loader, the checkpoint audit and the scorer
  verbatim; only the reduction is new.
- `dna_decode/glm/genomewide.py` — `load_peak_tiles`, `leave_peak_out`, `spearman`, `gc_features`.
  READ-ONLY; it is the registered owner of GSE144621 and this plan adds no second reader.
- `dna_decode/glm/expression.py` — `fit_predict_ridge` is the protocol `ridge_spearman` mirrors;
  consulted, not modified.
- `scripts/promoteratlas_local_reconcile.py` — the existing reconcile pattern the Step 5 gate copies.
- `scripts/freeze_exposure_bar.py` + `scripts/freeze_eukaryotic_arm_bar.py` — the freeze-with-verify
  pattern, including re-reading anchors from the ledger and refusing on drift.
- Searched and found nothing reusable for the reduction itself: `graphify-out/GRAPH_REPORT.md`
  (absent), `src/`, `lib/`, `utils/`, `common/`, `helpers/` (none exist in this layout).
