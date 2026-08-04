# FBA cell -> Keio essentiality validation (2026-08-03)

**Claim tested:** the FBA metabolic cell predicts gene-KO essentiality (a cell-level trait) for ANY
iML1515 gene, validated against the free Keio mutant-fitness gold standard.

- Model: **iML1515** (E. coli K-12); WT growth **0.877 /h** (D-Glucose minimal aerobic)
- Label: Keio BW25113 RB-TnSeq mutant fitness (Bernstein 2023); fitness<-2 = essential-on-glucose
- Genes scored (model AND labelled): **1339**  (essential prevalence 7.2% -- highly imbalanced)

## Metrics (assayable gene set)

| accuracy | MCC | precision | recall | F1 | ROC-AUC | PR-AUC |
|---|---|---|---|---|---|---|
| **0.954** | 0.652 | 0.684 | 0.670 | 0.677 | 0.839 | 0.5445 |

Confusion (positive = essential): TP 65 · FP 30 · TN 1212 · FN 32 (n=1339)

**Corroboration:** 101 FBA-essential genes have NO viable mutant in the Keio pool (an independent essentiality signal, excluded from the assayable-set metrics above).

## Caveats (honest scope)
- METABOLIC traits only; glucose M9 aerobic medium.
- Metrics computed on the ASSAYABLE gene set (genes with a measurable mutant); absolutely-essential genes without mutants are corroborated separately.
- Gene essentiality is highly class-imbalanced -> MCC + PR-AUC are more meaningful than ROC-AUC.
- In-distribution vs a published knowledge baseline; not an independent-lab claim.

**Gate:** accuracy 0.954 >= 0.85 -> **PASS**
