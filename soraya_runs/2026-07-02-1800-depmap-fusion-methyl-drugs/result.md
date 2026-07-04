# Soraya --until-mvp: Extend DepMap with fusions/methylation; add more drugs

**Verdict: mvp-reached** (all 4 criteria live-met). Run 2026-07-02-1800.

## MVP bar (frozen at run start; execute-mode draft-then-ratify)
- C1 file-exists scripts/depmap_fusion_methylation.py .............. MET
- C2 test-exit-0 pytest tests/test_depmap_fusion_methylation.py .... MET (4 passed)
- C3 file-exists wiki/depmap_fusion_methylation_result_2026-07-02.md MET
- C4 project-state-row (ledger Action Log 304) .................... MET

## What shipped
Two new feature-type modalities for the DepMap decoder, de-confounded by lineage via the promoted
dna_decode.deconfound primitives, on real DepMap/PRISM 19Q4 data (+ freshly fetched CCLE_fusions +
Broad RRBS TSS-1kb methylation, on D:):
- FUSION (binary presence): ALK->crizotinib/alectinib/lorlatinib. Fusion presence is the ONLY modality
  carrying a correct-direction de-confounded signal (mutation/CN/expression of ALK all null). Honestly
  UNDERPOWERED (n+~9; BCR-ABL n=0 in the 569-line PRISM subset).
- METHYLATION (continuous): MGMT->temozolomide (-0.171, canonical) + SLFN11->topotecan (+0.103,
  silencing->resistance); methylation agrees opposite-sign with expression.
The "feature type must match mechanism type" law now spans 5 feature types (mutation / copy-number /
expression / fusion / methylation).

## Wall (named)
Fusion is underpowered in PRISM 19Q4 (BCR-ABL absent). EXTERNAL data wall, CODE-CLOSABLE: a broader
drug-response screen (CTRP or GDSC, which include CML/heme lines) would power the BCR-ABL->imatinib/
dasatinib/nilotinib + ALK cases with no new code (loaders already support it) — only a response-matrix fetch.

## Discipline
Verify-in-batch caught + fixed 3 real-data bugs (duplicate fusion columns, NaN depmap_id, space-padded
'NA' methylation coercion). 4 offline synthetic tests build CCLE-shaped fixtures + run the REAL loaders.
Leak guard 9/9; frozen bacterial/viral/fungal AMR surface byte-unchanged. Commit e6f6326 (pushed).

No-resume: bounded attempts per active session only.
