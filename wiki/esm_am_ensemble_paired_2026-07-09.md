# ESM2 + AlphaMissense: the ensemble does not help, and 3B is a regression

**Date:** 2026-07-09 · **Hardware:** Kaggle GPU (Tesla T4, fp16) · **Script:** `scripts/kaggle_esm_am_sweep.py`
**Substrate:** humsavar pathogenic-vs-benign, 7 proteins (MLH1, MSH2, PROC, SLC4A1, MECP2, PINK1, MYOC)

Both GPU runs **reproduced the previously committed 650M and 3B artifacts exactly**, per-protein. ESM
masked-marginal scoring is deterministic, so every delta below is real signal, not seed noise. The
AlphaMissense column is byte-identical across runs, as expected.

## 1. The ensemble does not help (negative result)

`scripts/kaggle_esm_am_sweep.py` added a rank-average ensemble of ESM and AlphaMissense and reported a
positive "ensemble lift vs best-single" at both scales. **That lift was an artifact.** It computed

```
median(ENSEMBLE) - max(median(ESM), median(AM))
```

Those three medians are taken over *different proteins*, so the quantity is not a lift at all. The
correct paired statistic is the per-protein delta `ENS - max(ESM, AM)`:

| model | unpaired "lift" (printed) | **paired median** | **paired mean** | ensemble wins |
|---|---|---|---|---|
| ESM2-650M | +0.005 | **−0.0020** | **−0.0006** | 3 / 7 |
| ESM2-3B   | +0.002 | **−0.0060** | **−0.0104** | 2 / 7 |

At both scales the ensemble **loses to the better single predictor on the majority of proteins**. The
650M ensemble's largest single effect is a 0.025 AUROC *loss* on MLH1. Fixed in `b82aa79`: the verdict
now requires a paired median lift > 0.005 *and* a majority of wins, and prints
`ensemble ~= best single (no paired lift)` on this data.

## 2. Scaling 650M → 3B is a net regression on this panel

| protein | 650M | 3B | delta |
|---|---|---|---|
| MLH1 | 0.715 | 0.671 | −0.044 |
| MSH2 | 0.786 | 0.684 | **−0.102** |
| PROC | 0.559 | 0.538 | −0.021 |
| SLC4A1 | 0.904 | 0.900 | −0.004 |
| MECP2 | 0.937 | 0.938 | +0.001 |
| PINK1 | 0.811 | **0.834** | **+0.023** |
| MYOC | 0.970 | 0.970 | 0.000 |

**Mean paired delta −0.021; 3B is worse on 4/7, better on 2/7, tied on 1.**

The committed headline — *"ESM-3B edges past AlphaMissense; scale trajectory complete"* (median ESM 0.834
vs AM 0.823) — is literally true of the median and **rests entirely on PINK1**, the one protein 3B
improves, which happens to be the median element of both the 650M and the 3B run. The prior artifact's
`caveat` field already flagged that 3B was not uniformly better; this quantifies it. The direction is
consistent with the published ProteinGym curve, which flattens after 650M (650M ~0.47 → 3B ~0.48).

**Reading:** on this panel, model scale past 650M buys nothing, and neither does ensembling with the
supervised predictor. Both were the two "quality levers beyond raw size" the previous commit proposed;
both are now measured, and both are flat-to-negative.

## Caveats

- **n = 7 proteins.** Median AUROC gaps of ~0.01 are inside per-protein noise. The *paired* deltas are
  the trustworthy statistic here, and the 3B regression on MSH2 (−0.102) is far outside that noise.
- AlphaMissense is supervised and has home-field advantage; its training likely overlaps these clinical
  labels. A zero-shot LM matching it is the notable part, not the sign of the median gap.
- Only **rank-average** ensembling was tested. A weighted or stacked ensemble would need labels to fit,
  which is circular against a label-derived benchmark.
- **PROC is ~chance for every predictor** (0.52–0.56) and contributes nothing but variance.

## Provenance

- `wiki/esm_am_ensemble_paired_2026-07-09.json` — machine-readable, incl. per-protein paired lifts.
- Kaggle kernels: `emanueleebrahimi/esm2-650m-vs-alphamissense-ensemble`,
  `emanueleebrahimi/esm2-3b-am-ensemble-run` (private).
- Prior runs reproduced: `wiki/humsavar_am_vs_esm_sweep_650M.json`, `wiki/humsavar_am_vs_esm_sweep_3B.json`.
