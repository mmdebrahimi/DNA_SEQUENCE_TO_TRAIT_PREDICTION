# Decoder-suite provenance-disjoint validation report card — 2026-06-10

Standing trust surface for the shipped deterministic AMR decoders (Anchor-4). Each cell is the DEPLOYED `call_resistance(organism, drug)` rule scored on a FRESH, leakage-checked, **provenance-disjoint** NCBI-PD cohort (submitters OUTSIDE NARMS/CDC/FDA/GenomeTrakr/PulseNet/USDA).

> **Honest tier (do NOT inflate):** every SCORED cell is provenance-disjoint (different submitter/lab/country), **NOT** methodology-independent (most submitters use CLSI broth microdilution) and **NOT** external clinical validation. There is deliberately **no aggregate “X% validated” number** — read the grid cell by cell.

## State legend

| state | meaning |
|---|---|
| `SCORED` | Stage-2 provdisjoint run exists — acc/sens/spec shown |
| `POWERED_UNSCORED` | censused ≥ 20/class both classes; not yet scored |
| `UNDERPOWERED` | censused < 20/class (surveillance-dominated organism) |
| `ABSTAINS_BY_DESIGN` | registry EXPRESSION_FLOOR — rule refuses what it can't decode |
| `NOT_CENSUSED` | bacterial + census-able; no census yet |
| `NO_FREE_PHENOTYPE_SOURCE` | fungal/antiviral/antimalarial — no free isolate-level AST (structural non-cell) |

## State counts

| state | cells |
|---|---|
| `SCORED` | 3 |
| `POWERED_UNSCORED` | 3 |
| `UNDERPOWERED` | 1 |
| `ABSTAINS_BY_DESIGN` | 2 |
| `NO_FREE_PHENOTYPE_SOURCE` | 6 |

## Cells

| organism | drug | state | acc | sens | spec | n | detail |
|---|---|---|---|---|---|---|---|
| acinetobacter | meropenem | `ABSTAINS_BY_DESIGN` | — | — | — | — | registry verdict EXPRESSION_FLOOR (broad@1) — rule refuses expression-driven R it cannot decode |
| campylobacter | ciprofloxacin | `SCORED` | 1.0 | 1.0 | 1.0 | 40 | TP20 FP0 TN20 FN0 |
| escherichia_coli_shigella | ciprofloxacin | `SCORED` | 0.817 | 0.933 | 0.7 | 60 | TP28 FP9 TN21 FN2 |
| klebsiella | ceftriaxone | `POWERED_UNSCORED` | — | — | — | — | censused 505R/410S provenance-disjoint (>=MIN/class); not yet scored |
| klebsiella | ciprofloxacin | `SCORED` | 0.967 | 0.967 | 0.967 | 60 | TP29 FP1 TN29 FN1 |
| klebsiella | gentamicin | `POWERED_UNSCORED` | — | — | — | — | censused 317R/339S provenance-disjoint (>=MIN/class); not yet scored |
| klebsiella | tetracycline | `POWERED_UNSCORED` | — | — | — | — | censused 182R/132S provenance-disjoint (>=MIN/class); not yet scored |
| pseudomonas_aeruginosa | meropenem | `ABSTAINS_BY_DESIGN` | — | — | — | — | registry verdict EXPRESSION_FLOOR (broad@3) — rule refuses expression-driven R it cannot decode |
| salmonella | ciprofloxacin | `UNDERPOWERED` | — | — | — | — | censused 4R/87S provenance-disjoint (< MIN/class) — surveillance-dominated |
| candida_auris | fluconazole | `NO_FREE_PHENOTYPE_SOURCE` | — | — | — | — | fungal ERG11 BLAST target-site; no free isolate-level AST source on NCBI-PD (structural non-cell) |
| candida_auris | voriconazole | `NO_FREE_PHENOTYPE_SOURCE` | — | — | — | — | fungal ERG11 BLAST target-site; no free isolate-level AST source on NCBI-PD (structural non-cell) |
| candida_auris | caspofungin | `NO_FREE_PHENOTYPE_SOURCE` | — | — | — | — | fungal FKS1 BLAST target-site; no free isolate-level AST source on NCBI-PD (structural non-cell) |
| candida_auris | micafungin | `NO_FREE_PHENOTYPE_SOURCE` | — | — | — | — | fungal FKS1 BLAST target-site; no free isolate-level AST source on NCBI-PD (structural non-cell) |
| influenza_a | oseltamivir | `NO_FREE_PHENOTYPE_SOURCE` | — | — | — | — | influenza NA-inhibitor target-site; no free isolate-level AST source on NCBI-PD (structural non-cell) |
| plasmodium_falciparum | artemisinin | `NO_FREE_PHENOTYPE_SOURCE` | — | — | — | — | P. falciparum kelch13 target-site; no free isolate-level AST source on NCBI-PD (structural non-cell) |

## Provenance

- SCORED cells: `wiki/provenance_disjoint_validation_*.json` (Stage-2 `provenance_disjoint_validate.py`).
- Powering: `wiki/provdisjoint_census_results.json` (Stage-1 `ncbi_pd_provenance_census.py`).
- ABSTAINS: `dna_decode/data/calibrated_amr_rules.json` (EXPRESSION_FLOOR verdicts).
- Rebuild: `.venv/Scripts/python.exe scripts/build_validation_report_card.py` (read-only roll-up; re-run as cells land).
