# Soraya --until-mvp: GDSC/CTRP response matrix to power the fusion cases

**Verdict: mvp-reached** (all 4 criteria live-met). Run 2026-07-02-1900. Commit 403a895.

## MVP bar (execute-mode draft-then-ratify)
- C1 file-exists scripts/gdsc_fusion_decoder.py ......... MET
- C2 test-exit-0 pytest tests/test_gdsc_fusion_decoder.py MET (5 passed)
- C3 file-exists wiki/gdsc_fusion_result_2026-07-02.md .. MET
- C4 project-state-row (ledger 305) .................... MET

## What shipped
Chose GDSC over CTRP (direct xlsx, no redirect; COSMIC IDs the DepMap sample_info bridges). Joined GDSC1+GDSC2
fitted dose-response to the existing CCLE_fusions.csv via sample_info (COSMIC<->DepMap_ID<->lineage), scored
fusion presence vs LN_IC50 de-confounded within lineage. ALL 5 cases POWERED-MATCH:
  imatinib BCR-ABL n=9 t -6.5 / dasatinib n=10 t -14.95 / nilotinib n=9 t -9.09 / ponatinib n=8 t -4.83
  (MWU p~0) / crizotinib ALK n=15 t -4.09 (was -1.44 at n=9 in PRISM).
The BCR-ABL case (n=0 / untestable in PRISM) is now the STRONGEST fusion signal in the project and survives
within-blood-lineage de-confounding. The directive ("turn the directional ALK + untestable BCR-ABL into real
numbers") is fully satisfied. Prior artifact's external-data wall marked CLOSED.

## Notes
Installed openpyxl (free, one-time xlsx->parquet). Data (GDSC1/2 + sample_info) on D:, gitignored.
Frozen bacterial/viral/fungal AMR surface byte-unchanged (leak guard 9/9). No-resume: bounded attempts per
active session only.
