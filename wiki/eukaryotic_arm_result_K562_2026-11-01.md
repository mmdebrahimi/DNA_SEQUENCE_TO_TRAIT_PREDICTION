# Eukaryotic arm, K562 cell line -- SUBSTRATE_LIMITED

Artifact: `wiki/eukaryotic_arm_result_K562_2026-11-01.json`  ·  bar: `wiki/eukaryotic_arm_acceptance_bar.json`  ·  date: 2026-11-01

## Verdict: SUBSTRATE_LIMITED

the CI for the human gain lies entirely ABOVE the bacterial gain 0.0122.

Criterion: ceiling: best achievable prediction of bacterial promoter activity from sequence, reported against the replicate-derived 0.941 ceiling

- Human conv median Spearman **0.5237**, GC-only baseline **0.22**
- Gain CI over 10 folds: **[+0.2916, +0.3100]**, mean +0.3016
- Bacterial anchor gain (the quantity it is compared against): **0.0122**

## What this verdict does NOT say

- The artifact records **no** provisional verdict, so there is no would-have-been result to read here.
- The comparison that crosses the substrate boundary is the **within-substrate gain**, not raw Spearman. A within-substrate contrast is the only quantity that crosses a substrate boundary. Raw Spearman does not: the panels have different noise floors and different assays.

## Refused, and what was done instead

### `published_incumbents_head_to_head`

**Why refused:** Table S6 publishes per-fold values for Enformer, MPRALegNet, Sei, MPRAnn and a biochemical feature set and NEVER STATES THE METRIC. Comparing our named Spearman to an unnamed quantity would be an assumption wearing a number's clothes.

**Done instead:** they are recorded verbatim with metric_unstated: true, as context for what this panel supports, and they inform no predicate

## Honest limits

- THE BIGGEST ONE: this is ONE substrate change, not a substrate sweep. A human lentiMPRA panel differs from bacterial tiles in assay, in element length (230 bp vs 150 bp), in kingdom and in how the library was constructed. The arm cannot say WHICH of those differences carries any effect it finds.
- The bacterial anchor 0.0122 is treated as a FIXED POINT; its own uncertainty is not propagated into the CI, so a verdict that lands near the boundary is weaker than the predicate makes it look.
- The conv encoder's hyperparameters were chosen on the BACTERIAL rung. Reusing them on human sequence is the honest apples-to-apples choice for a method-class comparison, but it is NOT a tuned human model and must never be reported as one.
- A constructed panel is the regime this project already knows learned methods work in (constructed -> molecular, Spearman 0.761 on TEM-1). So a SUBSTRATE_LIMITED result CORROBORATES an existing finding rather than discovering something; its value is that it bounds the tile-chain ceiling, not that it is news about constructed panels.
- Nothing here licenses a claim about eukaryotic g->p generally. The project's organism-level eukaryotic embedding arm is a CLOSED NEGATIVE on natural populations (Arabidopsis, de-confounded) and this arm does not reopen it: constructed variation is a different regime, which is the whole reason this arm is scoped to a constructed panel.

## Substrate

- 98332 rows measured of 243778 designed (join rate 1.0)
- sequence: 230 bp -- 230nt sequence (15nt 5' adaptor - 200nt element - 15nt 3' adaptor)
- published folds: 10

## Target reparameterization

- applied: `expression = activity + 3.990495`
- computed over: the whole cell, once (never per fold, which would make folds incomparable)
- why: activity is log2(RNA/DNA) and goes negative; both incumbent fitters apply np.log1p to the target, which is undefined at <= -1, and spearman maps NaN to 0.0 -- so without this the failure would read as a clean null
- evaluation unaffected: Spearman is invariant to monotone transforms of the target; only the fitting scale moves

## Bacterial anchor

- best arm 0.5076, conv 0.3234, GC 0.3112, replicate ceiling 0.941
- fraction of ceiling captured by best: 0.5394
- headroom (best minus GC): 0.1964
- corroborating artifact: `wiki/glm_encoder_gate_result_2026-10-08.md`

## Per fold

| fold | conv | GC | gain | n_train | n_test |
|---|---|---|---|---|---|
| 1 | 0.5106 | 0.2444 | +0.2662 | 88499 | 9833 |
| 2 | 0.5357 | 0.2164 | +0.3193 | 88499 | 9833 |
| 3 | 0.5191 | 0.2138 | +0.3053 | 88499 | 9833 |
| 4 | 0.5175 | 0.2052 | +0.3123 | 88499 | 9833 |
| 5 | 0.5109 | 0.2122 | +0.2987 | 88499 | 9833 |
| 6 | 0.5264 | 0.2231 | +0.3033 | 88499 | 9833 |
| 7 | 0.5333 | 0.2354 | +0.2979 | 88499 | 9833 |
| 8 | 0.5237 | 0.2287 | +0.2950 | 88499 | 9833 |
| 9 | 0.5330 | 0.2112 | +0.3218 | 88499 | 9833 |
| 10 | 0.5166 | 0.2200 | +0.2966 | 88497 | 9835 |

---

Wall time 5469.3 s. Rendered by `scripts/eukaryotic_arm_memo.py` from the artifact named above; every number here is read from it, none recomputed.
