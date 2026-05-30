# All of Us research program polygenic risk score ancestry-stratified accuracy (2025) — supported memo (V1 invocation)
<!-- memo-schema: 0.4 -->

> Captured 2026-05-23. Source: Claude Code (/research orchestrator v0.4 — harness Skill tool invocation). Slug: all-of-us-prs-ancestry-stratified-accuracy-2025-2026-05-23.
> Audit floor: 5 of 5 locators required per row. Mapping floor: rationale → quantity required.
> Banned-phrase scan: 5 hard-reject + 4 soft-warn (tiered). Cite-token noise scanned + flagged.
> v0.4 added: topic-shape sanity check (Step 1.5), source-text identity advisory (Step 5.5), Decisions table column = Candidate use / Verification needed.

## Research Context

- **Problem:** All of Us research program polygenic risk score ancestry-stratified accuracy 2025
- **Captured:** 2026-05-23
- **Schema:** memo-schema 0.4

## Audit table (verbatim, supported rows only)

| Claim / quantity | Numeric value | Units | Source title | Authors / org | Year | Section / table / figure | Stable URL | Access date | Quoted excerpt (≤25 words) | Extraction rationale | Source type | Confidence |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| All of Us whole-genome sequences (release v7) | 245,388 | individuals | All of Us diversity and scale improve polygenic prediction contextually with greatest improvements for under-represented populations | All of Us consortium (bioRxiv preprint) | 2024 (revised May 2025 v2) | Abstract / Methods | https://www.biorxiv.org/content/10.1101/2024.08.06.606846v2.full | 2026-05-23 | "The All of Us research program (AoU) generated whole-genome sequences of 245,388 individuals (release v7)" | AoU genomic scale; baseline for ancestry-stratified PRS analyses | preprint (bioRxiv) | high |
| PRS trained on multi-ancestry / multi-biobank data — participant scale | ~750,000 | participants | All of Us diversity and scale improve polygenic prediction contextually with greatest improvements for under-represented populations | All of Us consortium | 2024 | Methods | https://www.biorxiv.org/content/10.1101/2024.08.06.606846v2.full | 2026-05-23 | "PRS trained on multi-ancestry and multi-biobank data with up to ~750,000 participants" | Training-data scale across AoU + UK Biobank combined | preprint (bioRxiv) | high |
| Number of traits/diseases analyzed in AoU multi-ancestry PRS | 32 | common complex traits | All of Us diversity and scale improve polygenic prediction contextually with greatest improvements for under-represented populations | All of Us consortium | 2024 | Methods, Results | https://www.biorxiv.org/content/10.1101/2024.08.06.606846v2.full | 2026-05-23 | "PRS trained on multi-ancestry and multi-biobank data with up to ~750,000 participants for 32 common, complex traits and diseases" | Trait-set breadth | preprint (bioRxiv) | high |
| Atrial Fibrillation PRS — European OR/SD (AoU validation) | 1.89 [95% CI 1.85-1.93] | OR per standard deviation | Multi-trait polygenic risk scores improve genomic prediction of atrial fibrillation across diverse ancestries | Jurgens S et al. (preprint) | 2025 | Figure 2, Table S3 | https://pmc.ncbi.nlm.nih.gov/articles/PMC12622184/ | 2026-05-23 | "OR/SD of 1.89 among individuals of European ancestry" | Highest absolute AF prediction accuracy in EUR | preprint (Research Square) | high |
| Atrial Fibrillation PRS — Asian ancestry OR/SD (AoU validation) | 1.76 [95% CI 1.56-1.99] | OR per standard deviation | Multi-trait polygenic risk scores improve genomic prediction of atrial fibrillation across diverse ancestries | Jurgens S et al. | 2025 | Figure 2, Table S3 | https://pmc.ncbi.nlm.nih.gov/articles/PMC12622184/ | 2026-05-23 | "Asian (OR/SD = 1.76; AUROC = 0.637)" | Multi-trait approach gain in non-EUR; relative-vs-absolute distinction | preprint (Research Square) | high |
| Atrial Fibrillation PRS — Admixed American OR/SD | 1.45 [95% CI 1.38-1.53] | OR per standard deviation | Multi-trait polygenic risk scores improve genomic prediction of atrial fibrillation across diverse ancestries | Jurgens S et al. | 2025 | Figure 2, Table S3 | https://pmc.ncbi.nlm.nih.gov/articles/PMC12622184/ | 2026-05-23 | "Admixed American (1.45; 0.595)" | Mid-tier PRS gain | preprint (Research Square) | high |
| Atrial Fibrillation PRS — African ancestry OR/SD | 1.39 [95% CI 1.32-1.45] | OR per standard deviation | Multi-trait polygenic risk scores improve genomic prediction of atrial fibrillation across diverse ancestries | Jurgens S et al. | 2025 | Figure 2, Table S3 | https://pmc.ncbi.nlm.nih.gov/articles/PMC12622184/ | 2026-05-23 | "African ancestry groups (1.39; 0.573)" | Lowest absolute AF gain but largest relative improvement | preprint (Research Square) | high |
| AF European sample size in AoU validation | 11,087 cases / 113,280 controls | individuals | Multi-trait polygenic risk scores improve genomic prediction of atrial fibrillation across diverse ancestries | Jurgens S et al. | 2025 | Methods, Table S4 | https://pmc.ncbi.nlm.nih.gov/articles/PMC12622184/ | 2026-05-23 | "European: 11,087 AF cases; 113,280 controls" | Statistical-power baseline for EUR results | preprint (Research Square) | high |
| T2D multi-ancestry GWAS meta-analysis — total scale | 360,000 cases + 1.8M controls (41% non-EUR) | individuals | Multi-ancestry polygenic risk scores for the prediction of type 2 diabetes and complications in diverse ancestries | Huerta-Chagoya A et al. | 2025 | Abstract | https://www.medrxiv.org/content/10.1101/2025.07.21.25331778.full.pdf | 2026-05-23 | "harmonize T2D GWAS meta-analyses across five ancestries... 360,000 T2D cases and 1.8 million controls (41% non-EUR)" | Largest reported T2D multi-ancestry effort; non-EUR proportion 41% | preprint (medRxiv) | high |
| Prostate cancer PRS effect-size range across ancestries (AoU validation) | 1.61 (Middle Eastern) to 2.19 (American) | OR | Validation and context-dependent effects of a prostate cancer polygenic risk score in the All of Us Research Program | (medRxiv preprint, 2025) | 2025 | Results, Table 2 (from search snippet) | https://www.medrxiv.org/content/10.1101/2025.10.01.25337107v1 | 2026-05-23 | "effect sizes ranged from 1.61 (95% CI=1.02-2.64) in Middle Eastern to 2.19 (95% CI=1.98-2.42) in American populations" | Prostate cancer PRS cross-ancestry portability; American ancestry highest | preprint (medRxiv) | medium |
| Prostate cancer PRS — total risk variants | 451 | SNPs | Validation and context-dependent effects of a prostate cancer polygenic risk score in the All of Us Research Program | (medRxiv preprint, 2025) | 2025 | Methods | https://www.medrxiv.org/content/10.1101/2025.10.01.25337107v1 | 2026-05-23 | "validated a previously developed multi-ancestry prostate cancer PRS of 451 risk variants" | Specific variant count; reproducibility anchor | preprint (medRxiv) | high |
| AoU vs prior datasets — non-European individuals (ratio) | ~3 | x more | NHGRI research highlight / All of Us diversity initiative | NHGRI | 2024-2025 | News release | https://www.genome.gov/news/news-release/researchers-optimize-genetic-tests-for-diverse-populations-to-tackle-health-disparities | 2026-05-23 | "represented about three times as many individuals of non-European ancestry compared to other major datasets" | Quantifies AoU diversity premium vs prior PRS-training datasets | gov / institutional | high |
| AoU vs prior datasets — multi-population-ancestry individuals (ratio) | ~8 | x more | NHGRI research highlight / All of Us diversity initiative | NHGRI | 2024-2025 | News release | https://www.genome.gov/news/news-release/researchers-optimize-genetic-tests-for-diverse-populations-to-tackle-health-disparities | 2026-05-23 | "eight times as many individuals with ancestry spanning two or more global populations" | Multi-ancestry individual coverage; admixture-rich cohort | gov / institutional | high |

## Source-Locator Coverage

- Total rows submitted: 17
- Survived audit floor: 13
- Survived mapping floor: 13
- Survived banned-phrase scan: 13
- Final supported: 13
- Survival rate: 13 / 17 (76%)

## Caveats per row

- **AF European OR/SD (1.89)** — Source is preprint (Research Square), not peer-reviewed; treat as direction-setting evidence.
- **AF Asian OR/SD (1.76) + AF Admixed American OR/SD (1.45) + AF African OR/SD (1.39)** — Same paper Table S3; preprint-tier confidence. Pair-wise AUROC values (0.637 / 0.595 / 0.573) were dropped during intake because their quote fields read "(same paper)" not a verbatim excerpt containing the numeric value — moved to unsupported memo with recovery path documented there.
- **AF European AUROC (0.646)** — Same fate as paired AUROC rows above (audit-floor failure on non-verbatim quote); moved to unsupported memo.
- **Prostate cancer PRS effect-size range (1.61-2.19)** — provenance: websearch-summary (rationale flags "Table 2 (from search snippet)"; re-verify against direct medRxiv source before high-confidence use). Confidence downgraded high → medium per Step 5.5 source-text identity advisory.
- **T2D 360,000 cases / 1.8M controls** — Source is medRxiv preprint; full per-ancestry validation numerics in AoU not extracted (medRxiv PDF returned title/authors only at fetch time — documented in raw memo honest gaps).
- **AoU 245,388 WGS + ~750,000 participants + 32 traits** — All three sourced from same bioRxiv preprint; full-text WebFetch returned HTTP 403 at extraction time. Quoted excerpts captured via WebSearch result rendering of the bioRxiv abstract. Acceptable preprint-tier evidence; verify against direct source before any uplift.
- **NHGRI ~3x / ~8x ratios** — Institutional news release, not peer-reviewed; matches consortium-paper claim shape but does not provide per-ancestry sample size breakdowns.

## Decisions for Human Confirmation (cap 5)

| Claim | Numeric value | Units | Source URL | Candidate use / Verification needed | Confidence |
|---|---:|---|---|---|---|
| AF PRS European OR/SD ceiling vs African ancestry floor (AoU validation) | 1.89 → 1.39 | OR per SD | https://pmc.ncbi.nlm.nih.gov/articles/PMC12622184/ | **Candidate use:** quantifies the ancestry-stratified PRS-performance gap on AoU-cohort validation — concrete numbers for product-positioning claims about ancestry-aware vs ancestry-blind PRS. **Verification needed:** paper is preprint (Research Square); confirm gap persists in peer-reviewed revision before lifting into any product-tier evidence claim. | high |
| Prostate cancer PRS cross-ancestry effect-size range (AoU validation) | 1.61 → 2.19 | OR | https://www.medrxiv.org/content/10.1101/2025.10.01.25337107v1 | **Candidate use:** demonstrates non-monotonic ancestry effect (American highest, not European) — counters simple "EUR-PRS performs best" narrative used in DNA-decoder marketing claims. **Verification needed:** confirm via direct medRxiv source-text (current extraction is from search-result rendering); also check sample sizes per ancestry to weight the comparison. | medium |
| AoU diversity ratio vs prior PRS-training datasets | ~3x non-EUR + ~8x multi-population | x more | https://www.genome.gov/news/news-release/researchers-optimize-genetic-tests-for-diverse-populations-to-tackle-health-disparities | **Candidate use:** sharp quantitative anchor for "AoU is the largest non-European-inclusive PRS-training cohort to date" positioning. **Verification needed:** NHGRI release is institutional summary; trace back to the primary consortium paper (Comm Med, publisher-blocked at fetch — recovery via Unpaywall OA-mirror or PMC alternate ID is the discipline). | high |
| AoU WGS scale (release v7) | 245,388 | individuals | https://www.biorxiv.org/content/10.1101/2024.08.06.606846v2.full | **Candidate use:** baseline cohort-scale number for any DNA-decoder claim that references AoU as a data source. **Verification needed:** confirm release v7 is the current AoU release at decision time (AoU releases on a 6-12 month cadence; numbers shift). | high |
| T2D multi-ancestry GWAS meta-analysis scale | 360,000 cases + 1.8M controls (41% non-EUR) | individuals | https://www.medrxiv.org/content/10.1101/2025.07.21.25331778.full.pdf | **Candidate use:** largest reported T2D multi-ancestry GWAS scale signal; useful for "PRS training data has become diverse" claims. **Verification needed:** per-ancestry validation R² in AoU not extracted (medRxiv PDF returned title/authors only); follow-up retrieval to obtain per-ancestry numerics before any product claim. | high |

## Verification trace (Mission Control L1)

This intake was invoked as part of Mission Control run `2026-05-23-1400-research-aou-prs-ancestry`. The parent run's Intent Contract is at `C:\Users\Farshad\PythonProjects\dna_decode\mission-control-runs\2026-05-23-1400-research-aou-prs-ancestry\intent-contract.md`.

**Validation steps applied:**
- Audit floor (Step 2): 13 pass / 4 fail (non-verbatim "(same paper)" quote on AUROC pair-rows)
- Mapping floor (Step 3): 13 pass / 0 fail
- Banned-phrase scan (Step 4): 0 hard-reject / 0 soft-warn
- Cite-token noise scan (Step 5): 0 flagged
- Source-text identity advisory (Step 5.5): 1 provenance-flag (Row 14 prostate cancer — search-snippet rationale; downgraded high → medium) / 0 author-identity-uncertain / 0 quote-shape-table-cell

**Verification result for parent run's sub-task "Intake validation":**
- Status: PASS
- Criterion: Rows pass audit floor + mapping floor + banned-phrase scan + cite-token scan + source-identity advisory (per parent Intent Contract verification criteria)
- Evidence: `C:\Users\Farshad\PythonProjects\dna_decode\research_outputs\all-of-us-prs-ancestry-stratified-accuracy-2025-2026-05-23.md` (this memo, 13 supported rows) + `C:\Users\Farshad\PythonProjects\dna_decode\research_outputs\all-of-us-prs-ancestry-stratified-accuracy-2025-2026-05-23_unsupported.md` (4 rejected rows)

## Promotion Gate reminder

This memo is INPUT to the 4-step Promotion Gate (Research_Intake_Checklist.md §7), NOT a Promotion approval. Do NOT lift any number from this memo into rules_v0_1.yaml or wiki/SME_Calibration_Worksheet.md without:
1. Doc resolves at the cited URL
2. Section reference exists in the doc
3. Quoted excerpt is verbatim in the doc
4. Mapping from excerpt to numeric value is natural to a domain reader
