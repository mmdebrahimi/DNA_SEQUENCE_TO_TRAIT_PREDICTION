# The supervised complement's three classes are not equal — PI's arm can never fire

**2026-09-30.** The wiring memo (`wiki/supervised_complement_wired_2026-09-30.md`) measured non-vacuity for
**NNRTI only** and named PI + INSTI as the next cheap step. Done. The answer is not uniform, and the PI
result is a real negative that the shipped doubt layer must not paper over.

## What was measured

For each supported class: partition the isolates by the **deployed** `call_hiv_observed`, keep the
catalog-**SUSCEPTIBLE** ones (that set *is* the blind spot), and ask whether the complement's shipped
`0.5` threshold separates the truly-resistant ones (measured Stanford fold ≥ 3) from the rest.

Genotype extraction is lifted from the complement's **own builder** (`build_hiv_complement_model.build`:
drifted-position detection → `ev.isolate_muts` → `prot[p-1] != aa`), so the tokens scored here are built the
same way as the features the model was trained on.

Both dataset variants are reported because they are **not the same question**: `*_DataSet.Full.txt` is the
set the model was *trained* on (its `n_train` 4222 matches NNRTI exactly — maximally in-sample), while
`*_DataSet.txt` is what the rest of this repo's HIV validation uses.

| class | dataset | blind spot | truly R | base rate | flagged @0.5 | precision | recall | **enrichment** |
|---|---|---|---|---|---|---|---|---|
| NNRTI | validation | 1,108 | 52 | 0.047 | 48 | 0.729 | 0.673 | **15.54×** |
| NNRTI | training | 1,877 | 200 | 0.107 | 129 | 0.806 | 0.520 | **7.57×** |
| INSTI | validation | 369 | 19 | 0.051 | 17 | 0.824 | 0.737 | **15.99×** |
| INSTI | training | 824 | 59 | 0.072 | 43 | 0.884 | 0.644 | **12.34×** |
| **PI** | validation | 614 | **2** | **0.003** | **0** | — | — | **never fires** |
| **PI** | training | 903 | **3** | **0.003** | **0** | — | — | **never fires** |

**The NNRTI validation row reproduces the wiring memo** (1,108 vs 1,109 isolates; 52 R; 48 flagged; 0.729 /
0.673; 15.54× vs 15.55×) from an independently written script — so that memo's headline is confirmed, not
merely restated. The one-isolate difference is a position-range guard, not a substantive disagreement.

## INSTI: non-vacuous, and the strongest of the three

**15.99× enrichment at precision 0.824 / recall 0.737** on the validation set, holding at 12.34× on the
training set. The INSTI arm earns its place in the doubt layer on the same evidence NNRTI does.

## PI: the threshold never fires, and that is NOT a complement defect

**The PI blind spot is almost empty — 2 truly-resistant isolates out of 614 (0.3%).** Max risk over the
whole set is 0.388, below the 0.5 threshold, on both datasets. There is essentially nothing there to rescue.

**Mechanism, and it is the PI catalog's own known trade-off.** PI v0 is **position-based** over 12 major
protease positions, so any substitution at any of them returns R. Protease is 99 aa and heavily mutated in
treatment-experienced isolates, so almost every resistant isolate trips a major position and lands in the
catalog-R bucket: **catalog-R fraction 0.660 (PI) vs 0.503 (INSTI) vs 0.488 (NNRTI)**. A catalog that
over-calls thoroughly leaves a near-empty susceptible bucket by construction.

So PI's catalog looks *complete* here for the wrong reason — it is complete on **sensitivity at the cost of
specificity**, which is precisely what the shipped PI v0.1 work was about (`hiv_pi_v0.1_validation_2026-06-23`:
8/8 PIs improve, the gain on **specificity**, TPV 0.469→0.684). INSTI is position-based too but over a
288 aa integrase with a lower catalog-R fraction, which is why it retains a real blind spot.

## Consequence for the shipped doubt layer — stated, not hidden

- **Do NOT claim the PI complement adds call-time value.** On this evidence its arm is structurally silent;
  a user running a PI drug will effectively never see it. `measured and quiet` is the truthful state for PI,
  and it is the state the signal already reports.
- **No code change is made.** The threshold is the module's own `DEFAULT_THRESHOLD`, and lowering it for PI
  to make the signal fire would be tuning a doubt layer to produce doubt — on 2 positive isolates, which
  cannot support a threshold either way.
- The per-class asymmetry is a property of **each catalog's over-call rate**, not of the complement.

## Honest limits

- **Every number here is IN-DISTRIBUTION** to the Stanford knowledge base — the training rows literally so.
  The deployable claim remains the shipped leave-one-**study**-out blind-spot AUROC — **corrected
  2026-09-30 to NNRTI 0.8102 · INSTI 0.8923 · PI no number**. This line read "(0.81 / 0.89 / 0.89)", but
  only the first was ever measured: `hiv_supervised_deployability.py` hardcodes `DRUG="EFV"` and the NNRTI
  dataset, so the PI and INSTI figures were unsourced literals in a config dict.
  `scripts/hiv_deployability_per_gene.py` measures every gene
  (`wiki/hiv_deployability_per_gene_2026-09-30.json`); INSTI's 0.89 turns out **accurate** (0.8923), while
  **PI is not measurable on this substrate** — its own catalog calls 74.2% of isolates resistant, leaving
  3 R of 910 in its blind spot, which is the same fact this table reports as `THRESHOLD_NEVER_FIRES`.
  The in-sample AUROCs in this run (0.94–0.99) are not the deployable claim and must not be quoted as
  performance.
- R/S here is a **threshold on continuous fold-change** at fold ≥ 3, not categorical AST. PI in particular
  has no per-drug clinical cutoff sourced in-repo (the `ILLUSTRATIVE_FOLD_CUTOFF` caveat), so "2 truly
  resistant" is cutoff-dependent — though the blind spot stays tiny under any nearby choice, because the
  catalog-R fraction is what empties it.
- One drug per class (EFV / LPV / RAL — the trained drug). Other drugs in a class may differ.
- **NRTI and CAI have no complement at all** and remain reported as not-measured.

## Reproduce

```
uv run python scripts/supervised_complement_nonvacuity.py --self-check   # pure guards, no data
uv run python scripts/supervised_complement_nonvacuity.py                # both datasets, all 3 classes
```

Artifact `wiki/supervised_complement_nonvacuity_2026-09-30.json`; 9 tests
`tests/test_supervised_complement_nonvacuity.py`. The script **raises rather than returns a number** when
the catalog partition is degenerate, when any row errors, or when a base rate is zero — the guard that was
missing when two earlier hand-written versions of this scan each produced a confident wrong conclusion.
Exits 2 when the gitignored Stanford data is absent, rather than reporting an empty scan.
