# AMR Portal — INDEPENDENT validation report card

Standing roll-up of the frozen deterministic decoder scored on the EBI AMR Portal (CABBAGE) provenance-disjoint, measured-AST isolates. **21 SCORED_INDEPENDENT** + 4 UNDERPOWERED of 25 cells. NO aggregate headline (per-cell truth only).

## Honesty rails
- **Independent at the ACCESSION level (upper bound).** Disjoint vs CRyPTIC + our tuning cohorts by BioSample/ERS/GCA; BioSample cross-archive resolution would only TIGHTEN it.
- **Genotype = the AMR Portal's own AMRFinderPlus run** (different operator → more independent; AMRFinder-version a named caveat). **Phenotype = wet-lab MIC/disk** (non-circular).
- **Rule applied UNCHANGED**; the frozen surface (`amr_rules.py` + `calibrated_amr_rules.json`) is byte-unchanged. NAMESPACE-SEPARATE from the NCBI-PD / HIV / external-cohort cards.
- **Calibrated registry now independently validated:** `Salmonella|ciprofloxacin` (broad) + `Klebsiella|ciprofloxacin` (qrdr_point + oqxAB-exclusion) — both OPT-IN configs whose provenance asked for an independent cohort; this card IS that cohort. (Promoting them to DEFAULT mutates the sha-pinned frozen file → a deliberate ratify-first freeze-amendment, NOT done here.)

## Cells (per-organism, measured AST, provenance-disjoint)
| Organism | Drug | Tier | Routing | nR / nS | sens (95% CI) | spec (95% CI) | acc |
|---|---|---|---|---|---|---|---|
| Escherichia coli | ceftriaxone | SCORED_INDEPENDENT | drug_rule_default | 1544/7234 | 0.921 [0.906, 0.933] | 0.918 [0.912, 0.924] | 0.919 |
| Escherichia coli | ciprofloxacin | SCORED_INDEPENDENT | drug_rule_default | 3068/12951 | 0.837 [0.823, 0.849] | 0.977 [0.974, 0.979] | 0.950 |
| Escherichia coli | gentamicin | SCORED_INDEPENDENT | drug_rule_default | 1944/13762 | 0.927 [0.915, 0.938] | 0.995 [0.994, 0.996] | 0.987 |
| Escherichia coli | meropenem | SCORED_INDEPENDENT | drug_rule_default | 98/12670 | 0.776 [0.683, 0.847] | 0.988 [0.986, 0.989] | 0.986 |
| Escherichia coli | tetracycline | SCORED_INDEPENDENT | drug_rule_default | 3310/5113 | 0.970 [0.963, 0.975] | 0.987 [0.984, 0.990] | 0.980 |
| Klebsiella pneumoniae | ceftriaxone | SCORED_INDEPENDENT | drug_rule_default | 2496/874 | 0.964 [0.955, 0.970] | 0.905 [0.884, 0.923] | 0.948 |
| Klebsiella pneumoniae | ciprofloxacin | SCORED_INDEPENDENT | calibrated_registry | 2970/1415 | 0.755 [0.739, 0.770] | 0.994 [0.989, 0.997] | 0.832 |
| Klebsiella pneumoniae | gentamicin | SCORED_INDEPENDENT | drug_rule_default | 2026/2904 | 0.928 [0.916, 0.939] | 0.963 [0.955, 0.969] | 0.949 |
| Klebsiella pneumoniae | meropenem | SCORED_INDEPENDENT | drug_rule_default | 1724/2898 | 0.948 [0.936, 0.957] | 0.859 [0.845, 0.871] | 0.892 |
| Klebsiella pneumoniae | tetracycline | SCORED_INDEPENDENT | drug_rule_default | 1059/936 | 0.731 [0.703, 0.757] | 0.980 [0.969, 0.987] | 0.848 |
| Salmonella enterica | ceftriaxone | SCORED_INDEPENDENT | drug_rule_default | 1629/21259 | 0.952 [0.941, 0.961] | 0.998 [0.998, 0.999] | 0.995 |
| Salmonella enterica | ciprofloxacin | SCORED_INDEPENDENT | calibrated_registry | 2434/22538 | 0.907 [0.895, 0.918] | 0.964 [0.962, 0.967] | 0.959 |
| Salmonella enterica | gentamicin | SCORED_INDEPENDENT | drug_rule_default | 1122/24608 | 0.912 [0.894, 0.927] | 0.991 [0.989, 0.992] | 0.987 |
| Salmonella enterica | meropenem | UNDERPOWERED | drug_rule_default | 0/19464 | — — | 1.000 [1.000, 1.000] | 1.000 |
| Salmonella enterica | tetracycline | SCORED_INDEPENDENT | drug_rule_default | 6941/19411 | 0.958 [0.953, 0.962] | 0.992 [0.990, 0.993] | 0.983 |
| Shigella flexneri | ceftriaxone | SCORED_INDEPENDENT | drug_rule_default | 14/373 | 1.000 [0.785, 1.000] | 0.987 [0.969, 0.994] | 0.987 |
| Shigella flexneri | ciprofloxacin | SCORED_INDEPENDENT | drug_rule_default | 50/333 | 0.900 [0.786, 0.957] | 0.997 [0.983, 0.999] | 0.984 |
| Shigella flexneri | gentamicin | UNDERPOWERED | drug_rule_default | 9/373 | 1.000 [0.701, 1.000] | 1.000 [0.990, 1.000] | 1.000 |
| Shigella flexneri | meropenem | UNDERPOWERED | drug_rule_default | 0/375 | — — | 1.000 [0.990, 1.000] | 1.000 |
| Shigella flexneri | tetracycline | SCORED_INDEPENDENT | drug_rule_default | 192/22 | 0.990 [0.963, 0.997] | 0.909 [0.722, 0.975] | 0.981 |
| Shigella sonnei | ceftriaxone | SCORED_INDEPENDENT | drug_rule_default | 426/910 | 0.988 [0.973, 0.995] | 0.995 [0.987, 0.998] | 0.993 |
| Shigella sonnei | ciprofloxacin | SCORED_INDEPENDENT | drug_rule_default | 625/961 | 0.730 [0.693, 0.763] | 0.997 [0.991, 0.999] | 0.892 |
| Shigella sonnei | gentamicin | SCORED_INDEPENDENT | drug_rule_default | 513/967 | 0.472 [0.429, 0.515] | 0.999 [0.994, 1.000] | 0.816 |
| Shigella sonnei | meropenem | UNDERPOWERED | drug_rule_default | 0/883 | — — | 1.000 [0.996, 1.000] | 1.000 |
| Shigella sonnei | tetracycline | SCORED_INDEPENDENT | drug_rule_default | 951/475 | 0.950 [0.934, 0.962] | 0.975 [0.956, 0.985] | 0.958 |

## Provenance
Scores `wiki/amr_portal_independent_scores.json` (`scripts/amr_portal_score_independent.py`, frozen rule). Feasibility `wiki/amr_portal_feasibility_result_2026-06-23.md`; validation memo `wiki/amr_portal_independent_validation_2026-06-23.md`. Rebuild: `uv run python scripts/build_amr_portal_report_card.py`.
