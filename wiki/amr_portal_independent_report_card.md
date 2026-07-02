# AMR Portal — INDEPENDENT validation report card

Standing roll-up of the frozen deterministic decoder scored on the EBI AMR Portal (CABBAGE) provenance-disjoint, measured-AST isolates. **23 SCORED_INDEPENDENT** + 4 UNDERPOWERED of 27 cells. NO aggregate headline (per-cell truth only).

## Honesty rails
- **Independent at the ACCESSION level (upper bound).** Disjoint vs CRyPTIC + our tuning cohorts by BioSample/ERS/GCA; BioSample cross-archive resolution would only TIGHTEN it.
- **Genotype = the AMR Portal's own AMRFinderPlus run** (different operator → more independent; AMRFinder-version a named caveat). **Phenotype = wet-lab MIC/disk** (non-circular).
- **Rule applied UNCHANGED**; the frozen surface (`amr_rules.py` + `calibrated_amr_rules.json`) is byte-unchanged. NAMESPACE-SEPARATE from the NCBI-PD / HIV / external-cohort cards.
- **Calibrated registry now independently validated:** `Salmonella|ciprofloxacin` (broad) + `Klebsiella|ciprofloxacin` (qrdr_point + oqxAB-exclusion) + `Campylobacter|ciprofloxacin` (qrdr_point; added 2026-06-28, C. jejuni acc 0.981 / C. coli 0.995) — OPT-IN configs whose provenance asked for an independent cohort; this card IS that cohort. (Promoting them to DEFAULT mutates the sha-pinned frozen file → a deliberate ratify-first freeze-amendment, NOT done here.)

## Cells (per-organism, measured AST, provenance-disjoint)
| Organism | Drug | Tier | Routing | nR / nS | sens (95% CI) | spec (95% CI) | acc |
|---|---|---|---|---|---|---|---|
| Campylobacter coli | ciprofloxacin | SCORED_INDEPENDENT | calibrated_registry | 437/1331 | 0.982 [0.964, 0.991] | 1.000 [0.997, 1.000] | 0.995 |
| Campylobacter jejuni | ciprofloxacin | SCORED_INDEPENDENT | calibrated_registry | 1533/4358 | 0.943 [0.931, 0.954] | 0.994 [0.991, 0.996] | 0.981 |
| Escherichia coli | ceftriaxone | SCORED_INDEPENDENT | drug_rule_default | 1541/7228 | 0.921 [0.906, 0.933] | 0.918 [0.912, 0.924] | 0.919 |
| Escherichia coli | ciprofloxacin | SCORED_INDEPENDENT | drug_rule_default | 3060/12947 | 0.837 [0.823, 0.850] | 0.977 [0.974, 0.979] | 0.950 |
| Escherichia coli | gentamicin | SCORED_INDEPENDENT | drug_rule_default | 1941/13756 | 0.928 [0.915, 0.939] | 0.995 [0.994, 0.996] | 0.987 |
| Escherichia coli | meropenem | SCORED_INDEPENDENT | drug_rule_default | 97/12664 | 0.773 [0.680, 0.845] | 0.988 [0.986, 0.990] | 0.986 |
| Escherichia coli | tetracycline | SCORED_INDEPENDENT | drug_rule_default | 3308/5112 | 0.970 [0.964, 0.975] | 0.987 [0.984, 0.990] | 0.980 |
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

## Overlay cells (NON-frozen — EXPERIMENTAL / CURATED, NOT in the shipped surface)

NON-frozen rules scored on the SAME AMR-Portal provenance-disjoint measured-AST isolates: the TMP-SMX `(≥1 acquired sul) AND (≥1 acquired dfr)` experimental overlay (`experimental_drug_rules.py`) + curated Tier-3/4 organism rules (`organism_rules/`, e.g. N. gonorrhoeae cipro gyrA-QRDR). **6/7 SCORED**. **Namespace-separate by design:** these are `EXPERIMENTAL_SCORED` / `CURATED_NONFROZEN` / `scorer_local` — NOT counted in the deployed-surface totals above, NOT in `shipped_decoder_surface`, frozen surface byte-unchanged.

- **Gate:** a cell is SCORED only if the 4-genotype strata REPRODUCE the Sci234/Oxford pattern (sul+dfr = highest-R stratum AND sul-only R-rate < 0.5); else INDETERMINATE (honest — the overlay's mechanism doesn't hold there, e.g. Klebsiella).
| Organism | Drug | headline | nR / nS | sens | spec | acc | strata-reproduced |
|---|---|---|---|---|---|---|---|
| Escherichia coli | trimethoprim-sulfamethoxazole | SCORED | 2619/9269 | 0.727 | 0.983 | 0.926 | True |
| Klebsiella pneumoniae | trimethoprim-sulfamethoxazole | INDETERMINATE | 2827/1384 | 0.430 | 0.981 | 0.611 | False |
| Neisseria gonorrhoeae | ciprofloxacin | SCORED | 5618/6406 | 0.943 | 0.990 | 0.968 | True |
| Salmonella enterica | trimethoprim-sulfamethoxazole | SCORED | 1667/24936 | 0.540 | 0.991 | 0.963 | True |
| Shigella flexneri | trimethoprim-sulfamethoxazole | SCORED | 138/69 | 0.964 | 0.957 | 0.961 | True |
| Shigella sonnei | trimethoprim-sulfamethoxazole | SCORED | 796/343 | 0.837 | 0.959 | 0.874 | True |
| Staphylococcus aureus | ciprofloxacin | SCORED | 1563/2095 | 0.981 | 0.884 | 0.926 | True |

## Provenance
Scores `wiki/amr_portal_independent_scores.json` (`scripts/amr_portal_score_independent.py`, frozen rule). Feasibility `wiki/amr_portal_feasibility_result_2026-06-23.md`; validation memo `wiki/amr_portal_independent_validation_2026-06-23.md`. Rebuild: `uv run python scripts/build_amr_portal_report_card.py`.
