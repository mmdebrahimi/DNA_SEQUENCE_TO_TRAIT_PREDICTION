# Eukaryotic / multimodal G2P substrate feasibility — supported memo (V1 invocation)
<!-- memo-schema: 0.4 -->

> Captured 2026-06-07. Source: Claude (`/research` orchestrator). Slug: eukaryotic-multimodal-substrate-feasibility-2026-06-07.
> Audit floor 5/5 + mapping floor + banned/cite scans applied. 10 rows submitted, 10 supported (3 medium per preprint/search-synthesis provenance).

## Research Context
- **Problem:** Rank the next genome→phenotype decoder substrate beyond bacterial AMR (eukaryotic + multimodal) for a solo dev, by data depth + label de-confoundability + compute requirement + curated-catalog existence; pick the single most feasible + its prerequisites.
- **Captured:** 2026-06-07 · **Schema:** memo-schema 0.4

## Audit table (verbatim, supported rows only)

(All 10 rows from `…raw.md` pass audit + mapping floors; see raw memo for the full 13-column table. Key supported claims condensed below.)

| Claim | Value | Source URL | Confidence |
|---|---|---|---|
| C. auris India WGS+MIC cohort | 350 isolates, 90% fluconazole-R | https://pubmed.ncbi.nlm.nih.gov/29325167/ | high |
| C. auris S. Africa MIC↔ERG11 linkage | 181/188 MIC>32 had ERG11 mut | https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10521600/ | high |
| C. auris global WGS | 12,644 genomes | https://www.biorxiv.org/content/10.64898/2026.02.03.703534v1.full | medium |
| ERG11 alone insufficient (multi-locus) | 3/~38 R had ERG11 (TAC1b/Cdr1/CNV) | https://pmc.ncbi.nlm.nih.gov/articles/PMC8092288/ | high |
| No genome-based AMR predictor for fungi | qualitative | https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10157239/ | high |
| Yeast-AMR ML used metabolic/env not DNA-FM | 54–75% acc | https://www.biorxiv.org/content/10.1101/2025.05.09.653161.full.pdf | medium |
| Arabidopsis flowering-time GWAS depth | 1,003 accessions (10°C) | https://pmc.ncbi.nlm.nih.gov/articles/PMC4949382/ | high |
| Arabidopsis genotype matrix | 10.7M SNPs / 2029 acc | https://academic.oup.com/nar/article/46/D1/D1150/4559687 | high |
| PlantCaduceus compute class | 24–80 GB GPU, 512bp windows | https://arxiv.org/pdf/2403.03234 | high |
| Multimodal colony-image+WGS paired set | not present | https://www.nature.com/articles/s41597-025-06319-4 | medium |

## Decisions for Human Confirmation (cap 5)

| Claim | Source URL | Candidate use / Verification needed | Confidence |
|---|---|---|---|
| C. auris azole resistance: deep WGS+MIC (188–350), near-perfect ERG11↔MIC linkage, NO genome-predictor exists yet | https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10521600/ | **Candidate use:** TOP next substrate — eukaryotic (Phase-6 jump) BUT same phenotype class (AMR/MIC) as proven work, determinant-based ⇒ the shipped deterministic decoder approach likely TRANSFERS; **no foundation model / no big GPU / no money** (it's an ERG11/FKS1/TAC1b mutation scan). **Verification needed:** assemble a machine-readable fungal determinant catalog (no AMRFinder-equivalent exists) + a de-confoundable cohort where clades don't alias the label (clade = the C. auris confound, analogous to lineage). | high |
| Fungal mechanism is multi-locus (ERG11 + TAC1b + Cdr1 efflux + aneuploidy) | https://pmc.ncbi.nlm.nih.gov/articles/PMC8092288/ | **Candidate use:** expect the same efflux/regulatory/CNV blind spot as bacterial AMR (e.g. tet/Klebsiella) — design `undetectable_mechanisms` for fungi up front (efflux overexpression, aneuploidy). **Verification needed:** what fraction of R is ERG11-only vs efflux/CNV-only (sets the achievable sensitivity ceiling). | high |
| Arabidopsis flowering-time: 1,003 accessions, sampling-independent quantitative label, NO curated catalog | https://pmc.ncbi.nlm.nih.gov/articles/PMC4949382/ | **Candidate use:** the true EMBEDDING-niche test (sampling-independent label + no catalog + depth ~1000 = YES/YES/YES) — where the frozen-FM thesis could finally get a fair shot. **Verification needed:** needs a plant DNA-FM (PlantCaduceus/AgroNT) on a **≥24GB GPU** — NOT the GTX 860M; this is the compute/MONEY-gated path. Confirm a non-paid GPU OR budget before committing. | high |
| Multimodal colony-image + WGS paired public dataset | https://www.nature.com/articles/s41597-025-06319-4 | **Candidate use:** EXCLUDE for now — no public paired dataset found (iron-law infeasible, like carbon-util). **Verification needed:** targeted Kaggle/Zenodo/HuggingFace search before any further consideration. | medium |
| PlantCaduceus / plant DNA-FM compute requirement | https://arxiv.org/pdf/2403.03234 | **Candidate use:** sizes the compute gate for the Arabidopsis path — 24GB (RTX 3090) minimum, A100 ideal. **Verification needed:** exact inference VRAM for the chosen model + whether Precision 7780 (RTX 3500 Ada, ~12GB?) suffices, or paid cloud is required (→ money gate). | high |

## Verification trace (Mission Control L1)
Run `2026-06-07-1500-research-eukaryotic-multimodal-substrate`. Audit floor 10/10 pass; mapping floor 10/10; banned-phrase 0; cite-token 0; provenance: 3 medium (preprint/search-synthesis). Status: PASS.

## Promotion Gate reminder
Input to the Promotion Gate, not approval. Verify each cited URL + section before lifting any number.
