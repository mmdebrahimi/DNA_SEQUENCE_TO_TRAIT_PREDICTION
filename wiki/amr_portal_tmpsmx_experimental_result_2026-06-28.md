# TMP-SMX experimental overlay on the EBI AMR Portal — result (2026-06-28)

Tier-1 of the AMR-Portal unscored-cell triage (`wiki/amr_portal_unscored_triage_2026-06-28.md`): the
5 powered TMP-SMX cells scored with the EXPERIMENTAL `(>=1 acquired sul) AND (>=1 acquired dfr)` overlay
(`dna_decode/data/experimental_drug_rules.py`) on the AMR Portal's free, provenance-disjoint, measured-AST
isolates. Scorer: `scripts/tmp_smx_amr_portal_validate.py` (reuses the Portal scorer's genotype/phenotype/
leakage loaders + the experimental rule). Artifact: `wiki/amr_portal_tmpsmx_experimental_2026-06-28.json`.

**Branding (load-bearing):** EXPERIMENTAL_SCORED / scorer_local / `not_in_shipped_surface=true`. This is a
SEPARATE artifact — it does NOT touch the frozen `amr_portal_independent_*` card (the deployed-surface
validation) nor the frozen AMR surface (this rule is not in it). Namespace-separation by design (the
shared-key silent-overwrite trap).

## Result — 4/5 SCORED (strata-reproduced); Klebsiella honestly INDETERMINATE

| cell | nR / nS | acc | sens | spec | strata reproduced | verdict |
|---|---|---|---|---|---|---|
| Escherichia coli | 2619 / 9269 | 0.926 | 0.727 | 0.983 | ✅ | SCORED |
| Salmonella enterica | 1667 / 24936 | 0.963 | 0.540 | 0.991 | ✅ | SCORED |
| Shigella sonnei | 796 / 343 | 0.874 | 0.837 | 0.959 | ✅ | SCORED |
| Shigella flexneri | 138 / 69 | 0.961 | 0.964 | 0.957 | ✅ | SCORED |
| Klebsiella pneumoniae | 2827 / 1384 | 0.611 | 0.430 | 0.981 | ❌ | **INDETERMINATE** |

Metric = binary measured AST only (the Portal has no MIC → no strict/relaxed tier). The 4-genotype-strata
table + the strata-reproduction gate (sul+dfr is the highest-R stratum AND sul-only R-rate < 0.5) is the
experimental-honesty check, reproducing the Sci234/Oxford `(sul AND dfr)` pattern.

## Reading it (honest)
- **The experimental AND-rule reproduces independently on 4/5 Enterobacterales/Shigella cells** — its
  premise (sul-only and dfr-only strata are NOT resistant; only sul+dfr is) holds on free provenance-disjoint
  measured AST, the strongest cross-cohort confirmation the overlay has (now Sci234 + Oxford + AMR Portal).
- **Specificity is uniformly high (0.957–0.991)** — the AND-rule rarely over-calls (its whole point).
- **Sensitivity varies (0.43–0.96)** — the known, documented fidelity caveat: the rule sees only ACQUIRED
  sul+dfr genes; `folP`/`folA` target-site POINT-mutation TMP-R is AMRFinder-determinant-invisible, so
  R isolates resistant via that route (or via only one acquired family + a point mechanism) are missed.
  This is a true-negative-blind-spot of the determinant approach, not a scoring bug.
- **Klebsiella INDETERMINATE is the honesty gate working:** its strata did NOT reproduce (sul-only R-rate
  too high / sul+dfr not the clean top stratum), so the scorer REFUSES the SCORED label rather than report a
  misleading 0.61-acc number. Klebsiella TMP-SMX has more non-acquired-gene resistance (folA/promoter); the
  AND-rule's premise doesn't hold there → correctly flagged, not force-scored.

## Scope + integrity
- EXPERIMENTAL, not a deployed claim; invisible to the frozen report card by design.
- FROZEN AMR surface (`amr_rules.py` + `calibrated_amr_rules.json`) byte-unchanged; leak guard 9/9.
- 4 pure-logic tests `tests/test_tmp_smx_amr_portal.py` (synthetic strata + gate; no parquet/network).
- Independence: AMR Portal provenance-disjoint (BioSample/ERS/GCA vs CRyPTIC + tuning cohorts); genotype =
  the Portal's own AMRFinderPlus run; phenotype = wet-lab measured AST (non-circular).

## Companion
`wiki/external_validation_tmpsmx_result_2026-06-16.md` (Oxford + Sci234), the prior independent TMP-SMX
validations; this extends the same experimental overlay to the AMR Portal at much larger N.
