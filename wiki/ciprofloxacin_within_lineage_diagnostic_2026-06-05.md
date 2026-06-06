# ciprofloxacin within-lineage conditional-ranking diagnostic (2026-06-05)

> Conditions on lineage: within each MLST carrying BOTH R and S, fraction of R/S pairs scored
> R>S (≈ lineage-conditioned AUC). If NT's edge persists here, the signal is mechanism, not lineage.
> Source scores: `ciprofloxacin_falsifier_2026-06-05.scores.json`

**Shared lineages used:** 6 · **within-lineage R/S pairs:** 43

| Variant | overall AUROC | within-lineage concordance | perm-null mean | p(≥obs) |
|---|---:|---:|---:|---:|
| NT-XGBoost | 0.914 | 0.605 | 0.506 | 0.365 |
| NT-logreg | 0.863 | 0.628 | 0.505 | 0.326 |
| POINT-XGB | 0.943 | 1.000 | 0.495 | 0.000 |

**NT-best within-lineage concordance 0.605 vs POINT-XGB 1.000 → Δ -0.395**

## Reading
- concordance ~0.5 = no within-lineage discrimination (signal was lineage, not mechanism).
- concordance >>0.5 with low p = the model ranks R above S EVEN within the same lineage = mechanism.
- NT edge persists within lineage if NT concordance > POINT-XGB AND NT p is small. (Δ -0.395, NT p 0.365).
- Small n_pairs ⇒ low power; treat as directional diagnostic, not a gate.