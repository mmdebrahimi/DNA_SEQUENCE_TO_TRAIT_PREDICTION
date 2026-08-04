# FBA essentiality validation: yeast (iMM904) (2026-08-03)

- WT growth **0.2879 /h** (default medium); genes scored **905**; essential prevalence 14.9%
- Gold standard: sgd (1215 experimental essential genes)

| accuracy | MCC | precision | recall | ROC-AUC | PR-AUC |
|---|---|---|---|---|---|
| 0.824 | **0.252** | 0.391 | 0.319 | 0.6301 | 0.3037 |

Confusion (positive=essential): TP 43 FP 67 TN 703 FN 92

**Discrimination: WEAK** (MCC 0.252). Accuracy is flattered by the
imbalanced majority class -- MCC is the honest signal.

## Caveats
- METABOLIC-gene essentiality only; the model's DEFAULT medium.
- Highly class-imbalanced -> MCC (not accuracy) is the discrimination signal.
- In-distribution vs a published knowledge baseline; not an independent-lab claim.
- Essentiality is medium-dependent; a mismatched default medium weakens the metric (a demotion path: set the organism's standard medium + re-score).
