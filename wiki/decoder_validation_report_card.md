# Decoder-suite provenance-disjoint validation report card — 2026-06-13

Standing trust surface for the shipped deterministic AMR decoders (Anchor-4). Rows are the DEPLOYED-CLAIM surface (`dna_decode/data/shipped_decoder_surface.py`) unioned with observed scored/census cells. Each cell is the DEPLOYED `call_resistance(organism, drug)` rule scored on a FRESH, leakage-checked, **provenance-disjoint** NCBI-PD cohort (submitters OUTSIDE NARMS/CDC/FDA/GenomeTrakr/PulseNet/USDA).

> **Honest tier (do NOT inflate):** every SCORED cell is an isolate-level provenance-disjoint stress test (different submitter/lab/country). The R classes are **clonally dominated** — the raw-isolate sens/spec is inflated by over-sampled clones, so the **Lineage disclosure** table below reports lineage-effective N + cluster-weighted sens/spec (one vote per lineage) with a Wilson CI. It is **NOT** methodology-independent (most submitters use CLSI broth microdilution) and **NOT** lineage-independent external clinical validation. There is deliberately **no aggregate “X% validated” number** — read the grid cell by cell.

## State legend

| state | meaning |
|---|---|
| `SCORED` | Stage-2 provdisjoint run exists — acc/sens/spec shown |
| `POWERED_UNSCORED` | censused ≥ 20/class both classes; not yet scored |
| `UNDERPOWERED` | censused < 20/class (surveillance-dominated organism) |
| `ABSTAINS_BY_DESIGN` | registry EXPRESSION_FLOOR — rule refuses what it can't decode |
| `NOT_CENSUSED` | bacterial + census-able; no census yet |
| `LABEL_CONFOUNDED` | phenotype label is an unreliable surrogate (oxacillin AST vs mecA) |
| `NO_FREE_PHENOTYPE_SOURCE` | fungal/antiviral/antimalarial — no free isolate-level AST (structural non-cell) |

## State counts

| state | cells |
|---|---|
| `SCORED` | 10 |
| `UNDERPOWERED` | 3 |
| `ABSTAINS_BY_DESIGN` | 2 |
| `LABEL_CONFOUNDED` | 1 |
| `NO_FREE_PHENOTYPE_SOURCE` | 11 |

## Cells

| organism | drug | state | acc | sens | spec | n | detail |
|---|---|---|---|---|---|---|---|
| acinetobacter | meropenem | `ABSTAINS_BY_DESIGN` | — | — | — | — | registry verdict EXPRESSION_FLOOR (broad@1) — rule refuses expression-driven R it cannot decode |
| campylobacter | ciprofloxacin | `SCORED` | 1.0 | 1.0 | 1.0 | 40 | TP20 FP0 TN20 FN0 |
| candida_auris | caspofungin | `NO_FREE_PHENOTYPE_SOURCE` | — | — | — | — | fungal_fks1; no free isolate-level AST source (structural non-cell) |
| candida_auris | fluconazole | `NO_FREE_PHENOTYPE_SOURCE` | — | — | — | — | fungal_erg11; no free isolate-level AST source (structural non-cell) |
| candida_auris | micafungin | `NO_FREE_PHENOTYPE_SOURCE` | — | — | — | — | fungal_fks1; no free isolate-level AST source (structural non-cell) |
| candida_auris | voriconazole | `NO_FREE_PHENOTYPE_SOURCE` | — | — | — | — | fungal_erg11; no free isolate-level AST source (structural non-cell) |
| escherichia_coli_shigella | ceftriaxone | `SCORED` | 0.967 | 0.967 | 0.967 | 60 | TP29 FP1 TN29 FN1 |
| escherichia_coli_shigella | ciprofloxacin | `SCORED` | 0.817 | 0.933 | 0.7 | 60 | TP28 FP9 TN21 FN2 |
| escherichia_coli_shigella | gentamicin | `SCORED` | 0.95 | 0.9 | 1.0 | 60 | TP27 FP0 TN30 FN3 |
| escherichia_coli_shigella | tetracycline | `SCORED` | 0.933 | 0.933 | 0.933 | 60 | TP28 FP2 TN28 FN2 |
| influenza_a | oseltamivir | `NO_FREE_PHENOTYPE_SOURCE` | — | — | — | — | influenza_na; no free isolate-level AST source (structural non-cell) |
| influenza_a | peramivir | `NO_FREE_PHENOTYPE_SOURCE` | — | — | — | — | influenza_na; no free isolate-level AST source (structural non-cell) |
| influenza_a | zanamivir | `NO_FREE_PHENOTYPE_SOURCE` | — | — | — | — | influenza_na; no free isolate-level AST source (structural non-cell) |
| klebsiella | ceftriaxone | `SCORED` | 0.95 | 1.0 | 0.9 | 60 | TP30 FP3 TN27 FN0 |
| klebsiella | ciprofloxacin | `SCORED` | 0.967 | 0.967 | 0.967 | 60 | TP29 FP1 TN29 FN1 |
| klebsiella | gentamicin | `SCORED` | 0.933 | 0.933 | 0.933 | 60 | TP28 FP2 TN28 FN2 |
| klebsiella | meropenem | `SCORED` | 0.683 | 0.467 | 0.9 | 60 | TP14 FP3 TN27 FN16 |
| klebsiella | tetracycline | `SCORED` | 0.883 | 0.8 | 0.967 | 60 | TP24 FP1 TN29 FN6 |
| plasmodium_falciparum | artemisinin | `NO_FREE_PHENOTYPE_SOURCE` | — | — | — | — | pf_kelch13; no free isolate-level AST source (structural non-cell) |
| plasmodium_falciparum | artesunate | `NO_FREE_PHENOTYPE_SOURCE` | — | — | — | — | pf_kelch13; no free isolate-level AST source (structural non-cell) |
| plasmodium_falciparum | chloroquine | `NO_FREE_PHENOTYPE_SOURCE` | — | — | — | — | pf_pfcrt; no free isolate-level AST source (structural non-cell) |
| plasmodium_falciparum | dihydroartemisinin | `NO_FREE_PHENOTYPE_SOURCE` | — | — | — | — | pf_kelch13; no free isolate-level AST source (structural non-cell) |
| pseudomonas_aeruginosa | meropenem | `ABSTAINS_BY_DESIGN` | — | — | — | — | registry verdict EXPRESSION_FLOOR (broad@3) — rule refuses expression-driven R it cannot decode |
| salmonella | ciprofloxacin | `UNDERPOWERED` | — | — | — | — | censused 4R/87S provenance-disjoint (< MIN/class) — surveillance-dominated |
| salmonella | gentamicin | `UNDERPOWERED` | — | — | — | — | censused 5R/86S provenance-disjoint (< MIN/class) — surveillance-dominated |
| salmonella | tetracycline | `UNDERPOWERED` | — | — | — | — | censused 4R/81S provenance-disjoint (< MIN/class) — surveillance-dominated |
| staphylococcus_aureus | oxacillin | `LABEL_CONFOUNDED` | — | — | — | — | phenotype LABEL is an unreliable surrogate (oxacillin AST vs mecA; cefoxitin is the CLSI surrogate) |

## Lineage disclosure (clonality-corrected)

Raw sens/spec counts one vote per ISOLATE; clones inflate it. Below: lineage-effective N (greedy-representative Mash clustering — chaining-resistant, NOT single-linkage) + cluster-weighted sens/spec (one vote per same-label lineage; mixed-label clones are DISCORDANT, never majority-voted) with a 95% Wilson CI. Weighted N is tiny — the CI is the point. Weighted metrics shown at Mash 0.005 (conservative); the JSON carries 0.001 too.

| organism | drug | raw N | eff lineages R/S @.001 | eff lineages R/S @.005 | wtd sens [95% CI] (n) | wtd spec [95% CI] (n) | discordant | grade |
|---|---|---|---|---|---|---|---|---|
| campylobacter | ciprofloxacin | 40 | 16/17 | 15/14 | 1.0 [0.796–1.0] (n=15) | 1.0 [0.785–1.0] (n=14) | 0 | moderate (>=15 effective lineages) |
| escherichia_coli_shigella | ceftriaxone | 60 | 23/25 | 11/17 | 1.0 [0.741–1.0] (n=11) | 1.0 [0.816–1.0] (n=17) | 3 | limited (8-14 effective lineages) |
| escherichia_coli_shigella | ciprofloxacin | 60 | 14/27 | 4/21 | 0.5 [0.15–0.85] (n=4) | 0.8 [0.584–0.919] (n=20) | 1 | scarce (3-7 effective lineages) |
| escherichia_coli_shigella | gentamicin | 60 | 15/23 | 5/9 | 0.6 [0.231–0.882] (n=5) | 1.0 [0.701–1.0] (n=9) | 4 | scarce (3-7 effective lineages) |
| escherichia_coli_shigella | tetracycline | 60 | 23/26 | 17/19 | 0.882 [0.657–0.967] (n=17) | 1.0 [0.832–1.0] (n=19) | 5 | moderate (>=15 effective lineages) |
| klebsiella | ceftriaxone | 60 | 21/30 | 16/21 | 1.0 [0.806–1.0] (n=16) | 0.95 [0.764–0.991] (n=20) | 2 | moderate (>=15 effective lineages) |
| klebsiella | ciprofloxacin | 60 | 9/23 | 2/18 | 0.5 [0.095–0.905] (n=2) | 1.0 [0.824–1.0] (n=18) | 1 | clonal (<3 effective lineages) |
| klebsiella | gentamicin | 60 | 16/13 | 11/7 | 1.0 [0.741–1.0] (n=11) | 0.857 [0.487–0.974] (n=7) | 2 | limited (8-14 effective lineages) |
| klebsiella | meropenem | 60 | 14/23 | 6/21 | 1.0 [0.61–1.0] (n=6) | 0.952 [0.773–0.992] (n=21) | 4 | scarce (3-7 effective lineages) |
| klebsiella | tetracycline | 60 | 24/28 | 19/27 | 0.842 [0.624–0.945] (n=19) | 0.963 [0.817–0.993] (n=27) | 1 | moderate (>=15 effective lineages) |

## Provenance

- Row set: `dna_decode/data/shipped_decoder_surface.py` (deployed-claim surface) ∪ observed cells.
- SCORED cells: `wiki/provenance_disjoint_validation_*.json` (Stage-2 `provenance_disjoint_validate.py`).
- Powering: `wiki/provdisjoint_census_results.json` (Stage-1 `ncbi_pd_provenance_census.py`).
- ABSTAINS: `dna_decode/data/calibrated_amr_rules.json` (EXPRESSION_FLOOR verdicts).
- Lineage disclosure: `wiki/provdisjoint_lineage_metrics.json` (`scripts/compute_lineage_metrics.py`).
- Rebuild: `.venv/Scripts/python.exe scripts/build_validation_report_card.py` (read-only roll-up; re-run as cells land).
