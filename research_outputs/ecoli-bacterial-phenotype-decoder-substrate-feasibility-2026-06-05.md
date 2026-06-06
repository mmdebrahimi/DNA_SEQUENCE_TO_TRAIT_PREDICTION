# Next decoder substrate feasibility — bacterial phenotypes with sampling-independent labels — supported memo (V1 invocation)
<!-- memo-schema: 0.4 -->

> Captured 2026-06-05. Source: Claude (`/research` orchestrator). Slug: ecoli-bacterial-phenotype-decoder-substrate-feasibility-2026-06-05.
> Audit floor: 5 of 5 locators required per row. Mapping floor: rationale → quantity required.
> Banned-phrase scan: 5 hard-reject + 4 soft-warn (tiered). Cite-token noise scanned + flagged.
> v0.4 added: topic-shape sanity check (Step 1.5), source-text identity advisory (Step 5.5), Decisions table column = Candidate use / Verification needed.

## Research Context

- **Problem:** Candidate genotype-to-phenotype decoder substrates beyond AMR/pathotype: bacterial phenotypes with sampling-independent lab-assay labels + buildable de-confounded cohorts, noting curated-baseline existence; ranked shortlist
- **Captured:** 2026-06-05
- **Schema:** memo-schema 0.4

## Audit table (verbatim, supported rows only)

| Claim / quantity | Numeric value | Units | Source title | Authors / org | Year | Section / table / figure | Stable URL | Access date | Quoted excerpt (≤25 words) | Extraction rationale | Source type | Confidence |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| BacDive carbon-utilization paired genome+phenotype dataset size | 4397 strains × 58 carbon sources | strains; carbon sources | Statistical prediction of microbial metabolic traits from genomes | Li et al. (PLOS Comp Biol) | 2023 | Results, "Increasing sample size" | https://pmc.ncbi.nlm.nih.gov/articles/PMC10729968/ | 2026-06-05 | "carbon utilization information for 4397 sequenced strains and 58 carbon sources" | direct dataset-size statement | peer-reviewed | high |
| Carbon-utilization in-clade ML (random forest, gene-content) accuracy | >90 | % accuracy | Statistical prediction of microbial metabolic traits from genomes | Li et al. | 2023 | Results, "Machine learning accurately predicts" | https://pmc.ncbi.nlm.nih.gov/articles/PMC10729968/ | 2026-06-05 | "accuracy of over 90% in most cases" | in-clade predictive ceiling for the classical baseline | peer-reviewed | high |
| Carbon-utilization OUT-OF-CLADE prediction (small 96-strain set) | did not beat null | qualitative | Statistical prediction of microbial metabolic traits from genomes | Li et al. | 2023 | Results, "Machine learning accurately predicts" | https://pmc.ncbi.nlm.nih.gov/articles/PMC10729968/ | 2026-06-05 | "did not significantly outperform both null models" | the open generalization gap (small N) — the embedding-niche test | peer-reviewed | high |
| Carbon-utilization out-of-clade at scale (tryptophan, 4397 set) | 92.2 (93.5 w/ feature selection) | % accuracy | Statistical prediction of microbial metabolic traits from genomes | Li et al. | 2023 | Results, "Increasing sample size" | https://pmc.ncbi.nlm.nih.gov/articles/PMC10729968/ | 2026-06-05 | "mean accuracy of 92.2%" | scale rescues out-of-clade for gene-content RF — sets the bar an embedding must beat | peer-reviewed | high |
| Phylogeny-based prediction fails for distant taxa (the de-confounding crux) | qualitative | — | Statistical prediction of microbial metabolic traits from genomes | Li et al. | 2023 | Abstract | https://pmc.ncbi.nlm.nih.gov/articles/PMC10729968/ | 2026-06-05 | "phylogeny-based predictions fail to predict traits for taxa that are phylogenetically distant from any strains in the training set" | confirms within-lineage-vs-distant is THE question for this phenotype | peer-reviewed | high |
| Deep-learning / embedding features for carbon-utilization | not present | — | Statistical prediction of microbial metabolic traits from genomes | Li et al. | 2023 | Methods/Results (features) | https://pmc.ncbi.nlm.nih.gov/articles/PMC10729968/ | 2026-06-05 | gene presence/absence + KEGG ortholog (KO) features; no embeddings | embeddings UNTRIED here → open niche | peer-reviewed | high |
| Larger BacDive ML metabolic-trait dataset | 15938 strains | sequenced strains | Predicting bacterial phenotypic traits through improved ML using high-quality curated datasets | (PMC12145430) | 2025 | Methods (dataset) | https://pmc.ncbi.nlm.nih.gov/articles/PMC12145430/ | 2026-06-05 | "selected 15,938 strains from BacDive for which genome sequences were available" | scale ceiling for the metabolic-trait substrate | peer-reviewed | medium |
| Keio growth-rate phenotyping (isogenic deletions, NOT natural diversity) | 4227 strains | deletion strains | Genomewide phenotypic analysis of growth, cell morphogenesis... in E. coli | (PMC6018989) | 2018 | Results | https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6018989/ | 2026-06-05 | imaged 4,227 Keio strains; estimated maximal growth rate (αmax) + saturating density | growth-rate is lab-assay but the big set is KNOCKOUTS, not natural GWAS | peer-reviewed | high |
| Keio biofilm screen (isogenic deletions) | 110 of 3985 | mutants reduce biofilm | A Genome-wide Approach to Identify Genes Involved in Biofilm Formation in E. coli | (PMC2779908) | 2007 | Results | https://pmc.ncbi.nlm.nih.gov/articles/PMC2779908/ | 2026-06-05 | "Of the 3985 mutants... 110 showed a reduction in biofilm formation" | biofilm = crystal-violet lab assay, but ready-made data is knockout-library not natural-isolate GWAS | peer-reviewed | high |
| E. coli host-association GWAS cohort (SAMPLING-DEFINED label → confounded) | 1198 WGS isolates | isolates | Genome-wide association reveals host-specific genomic traits in E. coli | (BMC Biology) | 2023 | Methods | https://pmc.ncbi.nlm.nih.gov/articles/PMC10088187/ | 2026-06-05 | "collection of 1198 whole-genome sequenced E. coli isolates" from 5 host species | host = sampling context → study==class confound (same trap as pathotype); EXCLUDE | peer-reviewed | high |

## Source-Locator Coverage

- Total rows submitted: 10
- Survived audit floor: 10
- Survived mapping floor: 10
- Survived banned-phrase scan: 10
- Final supported: 10
- Survival rate: 10 / 10 (100%)

## Caveats per row

- **Rows 6, 8** — quote-shape: method/feature description is paraphrased (not a flowing verbatim source sentence). Numeric value (row 8: 4227) + claim (row 6: no embeddings) are sound from the WebFetch, but verify verbatim wording against the source before any uplift.
- **Row 7** — provenance: websearch-summary (15,938 count not confirmed by direct WebFetch of PMC12145430). Confidence already `medium`; re-verify before high-confidence use.
- **Rows 8, 9, 10** — these are deliberately included as *negative/comparator* rows: Keio growth + Keio biofilm are isogenic-knockout libraries (wrong shape for natural-variation decoding); host-association is sampling-defined (the confound to AVOID). They scope what is NOT a ready substrate.
- No cite-token noise, no soft-warn/hard-reject banned phrases found.

## Decisions for Human Confirmation (cap 5)

| Claim | Numeric value | Units | Source URL | Candidate use / Verification needed | Confidence |
|---|---:|---|---|---|---|
| BacDive carbon-utilization paired genome+phenotype dataset | 4397 × 58 | strains × carbon sources | https://pmc.ncbi.nlm.nih.gov/articles/PMC10729968/ | **Candidate use:** carbon-source utilization is the top next-substrate — lab-assay labels (sampling-INDEPENDENT), large public cohort, NO AMRFinder-style single curated catalog. **Verification needed:** how many of the 4397 are *E. coli* per carbon source (multi-species DB; an E. coli-only slice may fall below the ≥100 de-confoundable floor for some sources). | high |
| Carbon-util gene-content RF out-of-clade accuracy at scale | 92.2 | % | https://pmc.ncbi.nlm.nih.gov/articles/PMC10729968/ | **Candidate use:** sets the bar an NT-embedding decoder must beat. **Verification needed:** this is already near-ceiling for *some* carbon sources — confirm there exist carbon sources where gene-content RF stays weak out-of-clade (those are where embeddings can add value, not the easy ones). | high |
| Phylogeny-based prediction fails for distant taxa | qualitative | — | https://pmc.ncbi.nlm.nih.gov/articles/PMC10729968/ | **Candidate use:** confirms the lineage-vs-mechanism (de-confound / within-lineage) question is the live crux for carbon-util too — same framework we built for AMR transfers directly. **Verification needed:** whether a within-lineage R/S-style co-occurrence cohort is buildable for ≥1 carbon source in E. coli. | high |
| Embedding features untried for carbon-utilization | not present | — | https://pmc.ncbi.nlm.nih.gov/articles/PMC10729968/ | **Candidate use:** the embedding niche is literally open here — no published DNA-LM decoder for metabolic traits. **Verification needed:** confirm via a focused lit check that no 2024-2026 preprint already did NT/Evo on BacDive metabolic traits. | high |
| Host-association GWAS cohort is sampling-defined | 1198 | isolates | https://pmc.ncbi.nlm.nih.gov/articles/PMC10088187/ | **Candidate use:** EXCLUDE — host/source labels reproduce the pathotype study==class confound; do not pick as substrate. **Verification needed:** none (confirmatory negative). | high |

additional candidates exist; review the full memo (Keio growth-rate + Keio biofilm rows scope what is NOT a ready natural-variation substrate).

## Verification trace (Mission Control L1)

This intake was invoked as part of Mission Control run `2026-06-06-0339-research-decoder-substrate-feasibility`. The parent run's Intent Contract is at `mission-control-runs/2026-06-06-0339-research-decoder-substrate-feasibility/intent-contract.md`.

**Validation steps applied:**
- Audit floor (Step 2): 10 pass / 0 fail
- Mapping floor (Step 3): 10 pass / 0 fail
- Banned-phrase scan (Step 4): 0 hard-reject / 0 soft-warn
- Cite-token noise scan (Step 5): 0 flagged
- Source-text identity advisory (Step 5.5): 1 provenance-flag (row 7) / 0 author-identity-uncertain / 2 quote-shape-paraphrase (rows 6, 8)

**Verification result for parent run's sub-task "Intake validation":**
- Status: PASS
- Criterion: Rows pass audit floor + mapping floor + banned-phrase scan + cite-token scan + source-identity advisory (per parent Intent Contract verification criteria)
- Evidence: `research_outputs/ecoli-bacterial-phenotype-decoder-substrate-feasibility-2026-06-05.md` (this memo, 10 supported rows) + `..._unsupported.md` (0 rejected rows)

## Promotion Gate reminder

This memo is INPUT to the 4-step Promotion Gate (Research_Intake_Checklist.md §7), NOT a Promotion approval. Do NOT lift any number from this memo into rules/config/wiki without:
1. Doc resolves at the cited URL
2. Section reference exists in the doc
3. Quoted excerpt is verbatim in the doc
4. Mapping from excerpt to numeric value is natural to a domain reader
