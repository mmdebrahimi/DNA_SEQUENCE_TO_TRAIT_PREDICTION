# FBA essentiality validation: yeast (iMM904) (2026-08-11)

- WT growth **0.9835 /h** (default medium); genes scored **905**; essential prevalence 14.9%
- Gold standard: sgd (1215 experimental essential genes)

| accuracy | MCC | precision | recall | ROC-AUC | PR-AUC |
|---|---|---|---|---|---|
| 0.874 | **0.377** | 0.723 | 0.252 | 0.6237 | 0.3454 |

Confusion (positive=essential): TP 34 FP 13 TN 757 FN 101

**Discrimination: MODERATE** (MCC 0.377). Accuracy is flattered by the
imbalanced majority class -- MCC is the honest signal.

## Caveats
- METABOLIC-gene essentiality only; medium = rich.
- Highly class-imbalanced -> MCC (not accuracy) is the discrimination signal.
- In-distribution vs a published knowledge baseline; not an independent-lab claim.
- Medium: rich (mode=label_matched). Essentiality is medium-dependent -- the SGD labels come from YPD (rich), so a minimal-medium score charges the model for biology; measured effect on yeast/iMM904 was MCC 0.2524 -> 0.3773 with FP 67 -> 13.
