# SOTA architectures for bacterial AMR prediction from genome sequence on small cohorts — unsupported claims (V1 invocation)

> Slug: sota-bacterial-amr-prediction-small-cohorts-2026-05-13. Captured 2026-05-13. These rows were rejected during V1 intake.
> Reasons fall into: missing locator, mapping-floor failure, hard-reject banned phrase, low-confidence label.

## Rejected rows

| Row content | Numeric value claimed | Verbatim quote (what intake saw) | Rejection reason | Suggested follow-up |
|---|---|---|---|---|
| Mean token embedding AUROC gain vs summary-token pooling for NT-v2 on DNA sequence classification | 6.8% AUC improvement (median; IQR 3.7-9.6%) | "mean token embedding delivered higher AUROC...in 42 [datasets] for NT-v2" | **mapping-floor failure**: quote contains "42 datasets" but does NOT contain the claimed 6.8% delta. The 6.8% value appears elsewhere in the source (per WebFetch extract), but not in this specific verbatim excerpt. | Re-fetch the source and capture a quote that contains the literal "6.8%" string OR the IQR "3.7-9.6%". Then re-submit. |
| Mean token embedding AUROC gain vs summary-token pooling for DNABERT-2 | 4.0% AUC improvement (median; IQR 2.0-5.5%) | "delivered higher AUROC...in 41 out of the 52 binary sequence classification datasets for DNABERT-2" | **mapping-floor failure**: quote contains "41 of 52" but does NOT contain the claimed 4.0% delta. | Same fix: re-fetch + capture a quote containing "4.0%" or the IQR. |
| TabPFN small-dataset rank vs XGBoost on ≤1250-sample benchmark | 4.88 rank (lower=better) | "TabPFN performs best...with 92.3% accuracy, outperforming the competition" | **mapping-floor failure**: quote contains "92.3% accuracy" (a different metric) but does NOT contain the claimed rank 4.88. | Re-fetch the McElfresh et al. paper and capture a quote that contains "rank 4.88" or the specific rank-table cell. The 92.3%-accuracy quote could anchor a separate row about accuracy instead. |
| Cross-population generalizability gap for cipro AMR prediction (England→Africa) | 87 → 50 (% accuracy delta) | "Random Forest and Light Gradient Boosting Machine were effective for ciprofloxacin (50% accuracy, F1 Score: 0.56)" in Africa-data | **mapping-floor failure (partial)**: quote contains 50% (the Africa side) but does NOT contain the 87% (England side) needed for the delta claim. | Either (a) split into two rows — one for the 87% England baseline (separate quote) and one for the 50% Africa drop, or (b) find a single quote containing both numbers from the same source. |

## Summary

- Total rejected: 4
- Reason breakdown:
  - Missing audit-floor locator(s): 0
  - Mapping-floor failure: 4
  - Hard-reject banned phrase: 0
  - Low confidence: 0

All 4 rejections share a common root cause: during raw-memo writing, the **verbatim quote selected did not contain the numeric value claimed by the row**. The values exist in the source documents (verified during WebFetch extraction), but the wrong quote was captured. This is a recoverable problem — re-fetching the sources and capturing different verbatim excerpts would let these rows pass the mapping floor on resubmission.

**Implication for the user's research question:** the substance of these rows is still useful as directional evidence. The supported memo's Decisions table already incorporates the TabPFN-small-sample finding and the pooling-strategy direction. The rejected rows mostly add quantitative specificity that the supported rows already carry qualitatively. Do not interpret 4 rejections out of 16 as "the research is wrong"; interpret as "the verbatim-quote discipline caught 4 cases where exact numbers weren't quote-supported."
