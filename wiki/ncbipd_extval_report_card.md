# NCBI-PD no-compute external-validation report card

_Generated 2026-07-21 from 4 organism artifacts. Roll-up of the reusable NCBI-PD external-validation substrate (`scripts/score_ncbipd_extval.py`): a pure metadata join + score of NON-FROZEN organism cells against NCBI Pathogen Detection's OWN published AMRFinderPlus calls, with no-compute lineage-collapse via NCBI-PD SNP clusters._

**Honest tier (no aggregate headline):** provenance-disjoint (frozen accessions excluded) but NOT methodology-independent (same AMRFinderPlus + same cell). RAW sens/spec is clonality-inflated; the **lineage-collapsed** columns (one vote per NCBI-PD SNP cluster) are the honest number. A cell `SCORED_ENDORSED` only if powered (>=5/class), non-degenerate, spec >= 0.85, sens >= 0.5.

**9 of 12 cells SCORED_ENDORSED** across 4 organisms.

| organism | drug | n (R/S) | raw sens/spec | **lineage sens/spec** | discordant | verdict |
|---|---|---|---|---|---|---|
| Campylobacter | gentamicin | 121 (31R/90S) | 0.968/1.000 | **1.000/1.000** | 1.000 | ✅ SCORED_ENDORSED |
| Campylobacter | tetracycline | 121 (66R/55S) | 1.000/0.945 | **1.000/0.933** | 1.000 | ✅ SCORED_ENDORSED |
| Neisseria gonorrhoeae | azithromycin | 156 (110R/46S) | 0.000/1.000 | **0.000/1.000** | 1.000 | DEGENERATE_NOT_ENDORSED |
| Neisseria gonorrhoeae | cefixime | 166 (19R/147S) | 0.789/0.905 | **0.727/0.892** | 2.000 | ✅ SCORED_ENDORSED |
| Neisseria gonorrhoeae | ceftriaxone | 169 (2R/167S) | 0.000/0.946 | **0.000/0.908** | 0.000 | UNDERPOWERED |
| Neisseria gonorrhoeae | ciprofloxacin | 163 (94R/69S) | 0.989/0.986 | **1.000/1.000** | 3.000 | ✅ SCORED_ENDORSED |
| Neisseria gonorrhoeae | penicillin | 31 (14R/17S) | 0.929/0.941 | **0.917/0.933** | 0.000 | ✅ SCORED_ENDORSED |
| Neisseria gonorrhoeae | tetracycline | 60 (34R/26S) | 0.324/1.000 | **0.367/1.000** | 1.000 | SCORED_NOT_ENDORSED |
| Staphylococcus aureus | ciprofloxacin | 82 (22R/60S) | 1.000/0.967 | **1.000/0.976** | 1.000 | ✅ SCORED_ENDORSED |
| Staphylococcus aureus | rifampin | 109 (8R/101S) | 0.875/0.980 | **0.857/0.986** | 0.000 | ✅ SCORED_ENDORSED |
| Streptococcus pneumoniae | erythromycin | 113 (43R/70S) | 0.977/0.971 | **1.000/0.962** | 1.000 | ✅ SCORED_ENDORSED |
| Streptococcus pneumoniae | tetracycline | 113 (10R/103S) | 1.000/1.000 | **1.000/1.000** | 0.000 | ✅ SCORED_ENDORSED |

## Notes
- **Lineage-collapse is no-compute** — NCBI-PD publishes per-isolate SNP clusters (`<PDG>.reference_target.all_isolates.tsv` → `PDS_acc`), collapsed via `clonality.cluster_weighted_confusion` (no Mash/Docker). Every endorsed cell HOLDS at the lineage level → the rules decode mechanism, not clonal structure.
- **DEGENERATE guard**: a cell predicting all-one-class (e.g. gono azithromycin all-S) is never endorsed even at spec/sens 1.0.
- **NON-FROZEN cells**; the frozen decoder surface is byte-unchanged throughout (`verify_lock` OK).
- Per-cell detail: the `wiki/*_ncbipd_extval_*.md` result docs.
