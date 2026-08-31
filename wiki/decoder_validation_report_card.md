# Decoder-suite provenance-disjoint validation report card — 2026-08-31

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

`blind.` = determinant-invisible fraction (of the scored measured-R, the fraction the cell calls non-R = FN/(TP+FN) = 1−sens) — the honest 'how much resistance this cell structurally misses'; DESCRIPTIVE, not an endorsement input. The truly-invisible vs rule-limited split is in `wiki/determinant_blindness_atlas.md` (NCBI-PD cells).

| organism | drug | state | acc | sens | spec | n | blind. | detail |
|---|---|---|---|---|---|---|---|---|
| acinetobacter | meropenem | `ABSTAINS_BY_DESIGN` | — | — | — | — | — | registry verdict EXPRESSION_FLOOR (broad@1) — rule refuses expression-driven R it cannot decode |
| campylobacter | ciprofloxacin | `SCORED` | 1.0 | 1.0 | 1.0 | 40 | 0.0 | TP20 FP0 TN20 FN0 |
| candida_auris | caspofungin | `NO_FREE_PHENOTYPE_SOURCE` | — | — | — | — | — | fungal_fks1; no free isolate-level AST source (structural non-cell) |
| candida_auris | fluconazole | `NO_FREE_PHENOTYPE_SOURCE` | — | — | — | — | — | fungal_erg11; no free isolate-level AST source (structural non-cell) |
| candida_auris | micafungin | `NO_FREE_PHENOTYPE_SOURCE` | — | — | — | — | — | fungal_fks1; no free isolate-level AST source (structural non-cell) |
| candida_auris | voriconazole | `NO_FREE_PHENOTYPE_SOURCE` | — | — | — | — | — | fungal_erg11; no free isolate-level AST source (structural non-cell) |
| escherichia_coli_shigella | ceftriaxone | `SCORED` | 0.967 | 0.967 | 0.967 | 60 | 0.033 | TP29 FP1 TN29 FN1 |
| escherichia_coli_shigella | ciprofloxacin | `SCORED` | 0.817 | 0.933 | 0.7 | 60 | 0.067 | TP28 FP9 TN21 FN2 |
| escherichia_coli_shigella | gentamicin | `SCORED` | 0.95 | 0.9 | 1.0 | 60 | 0.1 | TP27 FP0 TN30 FN3 |
| escherichia_coli_shigella | tetracycline | `SCORED` | 0.933 | 0.933 | 0.933 | 60 | 0.067 | TP28 FP2 TN28 FN2 |
| influenza_a | oseltamivir | `NO_FREE_PHENOTYPE_SOURCE` | — | — | — | — | — | influenza_na; no free isolate-level AST source (structural non-cell) |
| influenza_a | peramivir | `NO_FREE_PHENOTYPE_SOURCE` | — | — | — | — | — | influenza_na; no free isolate-level AST source (structural non-cell) |
| influenza_a | zanamivir | `NO_FREE_PHENOTYPE_SOURCE` | — | — | — | — | — | influenza_na; no free isolate-level AST source (structural non-cell) |
| klebsiella | ceftriaxone | `SCORED` | 0.95 | 1.0 | 0.9 | 60 | 0.0 | TP30 FP3 TN27 FN0 |
| klebsiella | ciprofloxacin | `SCORED` | 0.967 | 0.967 | 0.967 | 60 | 0.033 | TP29 FP1 TN29 FN1 |
| klebsiella | gentamicin | `SCORED` | 0.933 | 0.933 | 0.933 | 60 | 0.067 | TP28 FP2 TN28 FN2 |
| klebsiella | meropenem | `SCORED` | 0.683 | 0.467 | 0.9 | 60 | 0.533 | TP14 FP3 TN27 FN16 |
| klebsiella | tetracycline | `SCORED` | 0.883 | 0.8 | 0.967 | 60 | 0.2 | TP24 FP1 TN29 FN6 |
| plasmodium_falciparum | artemisinin | `NO_FREE_PHENOTYPE_SOURCE` | — | — | — | — | — | pf_kelch13; no free isolate-level AST source (structural non-cell) |
| plasmodium_falciparum | artesunate | `NO_FREE_PHENOTYPE_SOURCE` | — | — | — | — | — | pf_kelch13; no free isolate-level AST source (structural non-cell) |
| plasmodium_falciparum | chloroquine | `NO_FREE_PHENOTYPE_SOURCE` | — | — | — | — | — | pf_pfcrt; no free isolate-level AST source (structural non-cell) |
| plasmodium_falciparum | dihydroartemisinin | `NO_FREE_PHENOTYPE_SOURCE` | — | — | — | — | — | pf_kelch13; no free isolate-level AST source (structural non-cell) |
| pseudomonas_aeruginosa | meropenem | `ABSTAINS_BY_DESIGN` | — | — | — | — | — | registry verdict EXPRESSION_FLOOR (broad@3) — rule refuses expression-driven R it cannot decode |
| salmonella | ciprofloxacin | `UNDERPOWERED` | — | — | — | — | — | censused 4R/100S provenance-disjoint (< MIN/class) — surveillance-dominated |
| salmonella | gentamicin | `UNDERPOWERED` | — | — | — | — | — | censused 5R/86S provenance-disjoint (< MIN/class) — surveillance-dominated |
| salmonella | tetracycline | `UNDERPOWERED` | — | — | — | — | — | censused 5R/81S provenance-disjoint (< MIN/class) — surveillance-dominated |
| staphylococcus_aureus | oxacillin | `LABEL_CONFOUNDED` | — | — | — | — | — | phenotype LABEL is an unreliable surrogate (oxacillin AST vs mecA; cefoxitin is the CLSI surrogate) |

## Source-concentration disclosure (how many sources back each SCORED number)

The lineage table above corrects for clonal domination WITHIN a cohort. This asks the question one level up: how many independent SOURCES does the cohort draw on at all? Every cell here is provenance-DISJOINT from the tuning data — that was the design goal and it was met. Provenance-DIVERSE is a different property, was never claimed, and is what this measures.

WHY IT MATTERS, measured: `escherichia_coli_shigella x gentamicin` is 95% one BioProject and contains ZERO carriers of the `rmt` determinant family. It reports sens 0.893; two source-diverse measurements of the same cell with the same frozen rule report 0.429 and 0.523. A cohort with no carriers of a determinant family cannot detect a rule blind to that family.

The error is NOT directional: a 97%-single-source cipro cell reads PESSIMISTIC (spec 0.700 vs 0.988 on an 8-BioProject set). A single-site estimate is an estimate of that site. **These rows change no metric and no cell state.**

| organism | drug | N | BioProjects | largest share | dominant | unknown provenance |
|---|---|---|---|---|---|---|
| campylobacter | ciprofloxacin | 40 | 1 | 100%  **SINGLE-SOURCE** | PRJNA560409 | 0 |
| escherichia_coli_shigella | ceftriaxone | 60 | 4 | 55% | PRJNA278886 | 0 |
| escherichia_coli_shigella | ciprofloxacin | 60 | 2 | 97%  **SINGLE-SOURCE** | PRJNA278886 | 0 |
| escherichia_coli_shigella | gentamicin | 60 | 4 | 95%  **SINGLE-SOURCE** | PRJNA278886 | 0 |
| escherichia_coli_shigella | tetracycline | 60 | 6 | 53% | PRJNA662792 | 0 |
| klebsiella | ceftriaxone | 60 | 2 | 50% | PRJNA278886 | 0 |
| klebsiella | ciprofloxacin | 60 | 36 | 18% | PRJNA529587 | 0 |
| klebsiella | gentamicin | 60 | 32 | 48% | PRJNA278886 | 0 |
| klebsiella | meropenem | 60 | 15 | 37% | PRJNA504784 | 0 |
| klebsiella | tetracycline | 60 | 2 | 50% | PRJNA278886 | 0 |

**3 of 10** cells rest on ONE BioProject holding ≥80% of the cohort.


## Catalog-completeness disclosure (L2 doubt — can the RULE even represent it?)

The three tables above all ask *how good is the evidence for this cell's number*. This asks a different question: **does the deployed rule fail to represent a determinant family that is present in the data at all?** That is the failure which produced the gentamicin `rmt` blind spot, and no amount of better cohort evidence would have surfaced it.

Rows are **drug-level** — the screen runs across the whole cached determinant index, NOT this cell's cohort. `STRONG` means a family the rule cannot represent whose labelled carriers are uniformly resistant, **after a family-wise correction** over the families screened for that drug. The correction is load-bearing: the raw purity signature fires on 5 families and 4 are coincidences (cipro `qnrA1` at 4R/0S is p=0.030 against ~125 families screened).

**These rows change no metric and no cell state.**

| organism | drug | families rule can't represent | raw signature | **STRONG** | families |
|---|---|---|---|---|---|
| acinetobacter | meropenem | 317 | 0 | 0 | — |
| campylobacter | ciprofloxacin | 125 | 2 | 0 | — |
| escherichia_coli_shigella | ceftriaxone | 216 | 2 | 0 | — |
| escherichia_coli_shigella | ciprofloxacin | 125 | 2 | 0 | — |
| escherichia_coli_shigella | gentamicin | 131 | 1 | 1  **KNOWN GAP** | `rmtE1` |
| escherichia_coli_shigella | tetracycline | 89 | 0 | 0 | — |
| klebsiella | ceftriaxone | 216 | 2 | 0 | — |
| klebsiella | ciprofloxacin | 125 | 2 | 0 | — |
| klebsiella | gentamicin | 131 | 1 | 1  **KNOWN GAP** | `rmtE1` |
| klebsiella | meropenem | 317 | 0 | 0 | — |
| klebsiella | tetracycline | 89 | 0 | 0 | — |
| pseudomonas_aeruginosa | meropenem | 317 | 0 | 0 | — |
| salmonella | ciprofloxacin | 125 | 2 | 0 | — |
| salmonella | gentamicin | 131 | 1 | 1  **KNOWN GAP** | `rmtE1` |
| salmonella | tetracycline | 89 | 0 | 0 | — |

Across 5 screened drugs, **1** determinant family survives the correction. **Honest limit:** this project has exactly ONE independently confirmed completeness gap, so recovering it is a single case and **not a rate** — it bounds nothing about gaps never confirmed.


## Prospective-lock disclosure (temporal — leakage-free BY CONSTRUCTION)

A SEPARATE arm from the provenance-disjoint numbers above, not a replacement for them. Every isolate here became public STRICTLY AFTER the decoder was frozen and sha256-pinned (`wiki/prospective_lock_manifest_2026-06-22.json`), so the decoder cannot have been tuned to it — the leakage argument is temporal, not statistical. `verify_lock` re-hashes the live decoder on every scoring run and hard-fails on drift.

HONEST SCOPE: N is small and ACCRUES over time; this is a temporal stress test, NOT lineage-independent clinical validation, and these rows are NOT clonality-corrected (the lineage table above applies to the provdisjoint cohorts only).

| organism | drug | lock date | N (R/S) | acc | sens | spec | abstain | powering | as of |
|---|---|---|---|---|---|---|---|---|---|
| escherichia_coli_shigella | ciprofloxacin | 2026-06-13 | 61 (24R/37S) | 0.967 | 0.917 | 1.000 | 0 | POWERED | 2026-08-24 |
| escherichia_coli_shigella | gentamicin | 2026-06-13 | 62 (49R/13S) | 0.532 | 0.429  **REGRESSION** | 0.923 | 0 | POWERED | 2026-08-24 |

A LOW prospective sens with HIGH spec means the rule under-calls — it is missing determinants, not mislabelling. Diagnose the false negatives' features before reading it as decay; see `wiki/prospective_lock_first_accrual_2026-08-24.md`, where exactly that diagnosis located a real catalog gap rather than drift.


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

## Curated layer value-over-naive-baseline

The deployed `call_resistance` rule vs NAIVE AMRFinder use ('any drug-class determinant → R', no subclass/point/threshold refinement) on the SAME labels, balanced accuracy. The curated layer must BEAT naive tool use on INDEPENDENT data, else the number only proves the tool works (the validate-wrapper-vs-underlying-tool rail). Reconciled cells only.

| surface | organism | drug | frozen balacc | naive balacc | Δ | verdict |
|---|---|---|---|---|---|---|
| Oxford ext. MIC | Escherichia_coli_Shigella | ciprofloxacin | 0.949 | 0.7665 | 0.1825 | CURATED_LAYER_ADDS_VALUE |
| Oxford ext. MIC | Escherichia_coli_Shigella | ceftriaxone | 0.827 | 0.5015 | 0.3255 | CURATED_LAYER_ADDS_VALUE |
| Oxford ext. MIC | Escherichia_coli_Shigella | gentamicin | 0.9585 | 0.714 | 0.2445 | CURATED_LAYER_ADDS_VALUE |
| provdisjoint | Campylobacter | ciprofloxacin | 1.0 | 1.0 | 0.0 | NAIVE_TIES_CURATED |
| provdisjoint | Escherichia_coli_Shigella | ciprofloxacin | 0.8165 | 0.65 | 0.1665 | CURATED_LAYER_ADDS_VALUE |
| provdisjoint | Escherichia_coli_Shigella | ceftriaxone | 0.967 | 0.7335 | 0.2335 | CURATED_LAYER_ADDS_VALUE |
| provdisjoint | Escherichia_coli_Shigella | gentamicin | 0.95 | 0.6335 | 0.3165 | CURATED_LAYER_ADDS_VALUE |
| provdisjoint | Escherichia_coli_Shigella | tetracycline | 0.933 | 0.85 | 0.083 | CURATED_LAYER_ADDS_VALUE |
| provdisjoint | Klebsiella | ciprofloxacin | 0.967 | 0.7 | 0.267 | CURATED_LAYER_ADDS_VALUE |
| provdisjoint | Klebsiella | meropenem | 0.6835 | 0.5 | 0.1835 | CURATED_LAYER_ADDS_VALUE |

_9 cells the curated layer adds value, 1 ties, 0 naive-beats. Sources: `wiki/external_validation_oxford_naive_comparator_*.json` + `wiki/provdisjoint_naive_comparator_*.json`; full synthesis `wiki/curated_vs_naive_value_add_synthesis_2026-06-27.md`._

## Provenance

- Row set: `dna_decode/data/shipped_decoder_surface.py` (deployed-claim surface) ∪ observed cells.
- SCORED cells: `wiki/provenance_disjoint_validation_*.json` (Stage-2 `provenance_disjoint_validate.py`).
- Powering: `wiki/provdisjoint_census_results.json` (Stage-1 `ncbi_pd_provenance_census.py`).
- ABSTAINS: `dna_decode/data/calibrated_amr_rules.json` (EXPRESSION_FLOOR verdicts).
- Lineage disclosure: `wiki/provdisjoint_lineage_metrics.json` (`scripts/compute_lineage_metrics.py`).
- Rebuild: `.venv/Scripts/python.exe scripts/build_validation_report_card.py` (read-only roll-up; re-run as cells land).
