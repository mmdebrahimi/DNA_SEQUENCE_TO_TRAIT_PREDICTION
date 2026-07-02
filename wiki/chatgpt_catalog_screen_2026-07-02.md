# ChatGPT genome->phenotype catalog — screened through the 8-gate scorecard (2026-07-02)

The catalog is USEFUL, and its own central warning (label leakage / evidence circularity) VALIDATES this project's hard-won discipline. But routed through the scorecard, most of it is (a) USER-authority-gated acquisition, (b) circular-label catalogs needing independent validation, or (c) human-variant FEATURE sources — NOT free executor-actionable phenotype substrates. Exactly ONE new source clears all 8 gates like the yeast substrate: **DepMap/CCLE**.

## Label-substrate sources (scorecard-scored)

| candidate | verdict | decoder | depth | fails | unknowns | usage class |
|---|---|---|---|---|---|---|
| DepMap / CCLE | **PASS** | learned-niche | 1000 | — | — | A — EXECUTOR-PILOTABLE NOW (yeast analog) |
| GDC / TCGA (open tier) | **VERIFY** | neither | 10000 | — | G3_sampling_independent | A' — pilotable, verify G3 confounds first |
| All of Us / UK Biobank / FinnGen / dbGaP / EGA | **REJECT** | neither | 500000 | G1_accessible | G3_sampling_independent | C — USER-AUTHORITY acquisition (freeze fwd-path #1) |
| CIViC | **REJECT** | deterministic | catalog | G2_non_circular | G5_provenance_separable | B — deterministic catalog + independent validation (NOT a label) |
| ClinVar / ClinGen | **REJECT** | deterministic | catalog | G2_non_circular | G5_provenance_separable | B — deterministic catalog + temporal-holdout validation (NOT a label) |

## Non-label sources (categorized, NOT scored — they aren't genotype→phenotype label substrates)

- **D — Feature sources** (regulatory/expression features for a HUMAN variant-effect model, not organismal-phenotype labels): GTEx, ENCODE, SCREEN, eQTLGen, 4D Nucleome, FANTOM5, Factorbook, Ensembl VEP, CADD, SpliceAI. Free; relevant only if the project pivots to human noncoding variant-effect prediction.
- **E — Benchmark / anti-circularity harnesses**: GIAB, CAGI, precisionFDA. Align with the project's leakage discipline; adopt as an evaluation harness if going human.
- **Commercial / gated**: HGMD, COSMIC, OncoKB, Tempus, 23andMe. License/consent-gated; USER-authority.

## Notes per scored source

- **DepMap / CCLE** — WEB-VERIFIED free CSVs (depmap.org): PRISM 4,518 compounds x ~578 lines, CRISPRGeneEffect, CCLE mutations/expression for >1,000 lines, joined by depmap_id. The CANCER ANALOG of the yeast win. De-confound by lineage/tissue-of-origin (like yeast clades). EXECUTOR-PILOTABLE NOW. Sources: https://depmap.org/portal/download/, https://doi.org/10.6084/m9.figshare.9393293.v4
- **GDC / TCGA (open tier)** — Open tier has somatic MAF + clinical/survival. G3 UNKNOWN: tumor purity/heterogeneity + site/batch confounds. Controlled tier (germline/BAM) is USER-gated. Sources: https://gdc.cancer.gov/
- **All of Us / UK Biobank / FinnGen / dbGaP / EGA** — G1 FAIL by design: CONTROLLED ACCESS (researcher agreement / DAC approval / institutional affiliation, possibly cost). These are NOT executor-fetchable -- they are the USER-AUTHORITY acquisition path (reproducibility-freeze forward-path #1). Real measured phenotypes (G2 pass). G3 UNKNOWN: EHR ascertainment/site confounds. Highest supervised-training value IF acquired. Sources: https://www.researchallofus.org/, https://www.ukbiobank.ac.uk/
- **CIViC** — CC0 (free). Same circularity shape as ClinVar (curated interpretations). REJECT as label; deterministic-catalog use with independent validation only. Sources: https://civicdb.org/
- **ClinVar / ClinGen** — G2 FAIL: WEB-VERIFIED type-1 circularity -- assertions derived from the in-silico predictors a decoder competes with (ACMG PP3/BP4). REJECT as a LEARNED label. Usable only as a DETERMINISTIC catalog (WHO-TB analog) validated against a TEMPORAL-HOLDOUT / independent measured label. Sources: https://www.ncbi.nlm.nih.gov/clinvar/, https://pubmed.ncbi.nlm.nih.gov/36413997/

## Headline
**DepMap/CCLE is the one new free, measured, joinable, deep, non-circular substrate** — the cancer analog of the yeast win, pilotable exactly like it. Everything else is user-gated acquisition (biobanks), circular-as-label catalogs (ClinVar/CIViC — deterministic-catalog use only, with independent validation), or human-variant feature sources. The catalog's circularity warning is the project's own wall, independently confirmed.
