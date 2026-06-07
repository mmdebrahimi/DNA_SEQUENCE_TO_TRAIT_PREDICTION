# Eukaryotic / multimodal G2P substrate feasibility (raw research, 2026-06-07)

> Captured 2026-06-07. Source: Claude Code (`/research` orchestrator). Topic: ranked feasibility of the next genome→phenotype substrate beyond bacterial AMR (eukaryotic + multimodal), for a solo dev — by data depth, label de-confoundability, compute requirement, and curated-catalog existence. Slug: eukaryotic-multimodal-substrate-feasibility-2026-06-07.
> Web search + one-pass synthesis (~4 calls). "Thorough Google" depth.

## Audit table (verbatim, all candidate rows)

| Claim / quantity | Numeric value | Units | Source title | Authors / org | Year | Section | Stable URL | Access date | Quoted excerpt (≤25 words) | Extraction rationale | Source type | Confidence |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| C. auris India multicentre WGS+MIC cohort | 350 | isolates | A multicentre study of antifungal susceptibility among 350 C. auris isolates | (J Antimicrob Chemother) | 2018 | Abstract | https://pubmed.ncbi.nlm.nih.gov/29325167/ | 2026-06-07 | "isolates (n = 350) from 10 hospitals in India ... 90% of C. auris were fluconazole resistant" | depth + MIC for fungal-AMR substrate | peer-reviewed | high |
| C. auris S. Africa WGS: MIC↔ERG11 linkage | 181 of 188 | isolates MIC>32 with ERG11 mut | C. auris outbreak neonatal unit South Africa | (Emerg Infect Dis) | 2023 | Results | https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10521600/ | 2026-06-07 | "all 181/188 isolates that had a fluconazole MIC >32 µg/mL had ERG11 mutations" | near-perfect genotype↔MIC determinant linkage | peer-reviewed | high |
| C. auris global WGS collection | 12644 | genomes | Global genomic epidemiology of Candida auris | (bioRxiv) | 2026 | Abstract | https://www.biorxiv.org/content/10.64898/2026.02.03.703534v1.full | 2026-06-07 | "analysis of 12,644 whole genome sequences from 1997-2024" | scale ceiling for fungal-AMR substrate | preprint | medium |
| ERG11 alone insufficient (clade II) | 3 of ~38 | R isolates with known ERG11 mut | Genome-wide analysis experimentally evolved C. auris | (PMC8092288 / clade II study) | 2021 | Results | https://pmc.ncbi.nlm.nih.gov/articles/PMC8092288/ | 2026-06-07 | "62.3% of 61 isolates ... only 3 isolates harbored a known azole-conferring mutation in ERG11" | mechanism is multi-locus (TAC1b/Cdr1/CNV) — blind-spot like efflux | peer-reviewed | high |
| Genome-based AMR predictors absent for fungi | qualitative | — | Challenges of genome-based identification of antifungal resistance | (Frontiers Microbiol) | 2023 | Intro | https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10157239/ | 2026-06-07 | "genome-based predictors of antimicrobial resistance (AMR) are not available for fungal pathogens" | NO AMRFinder-equivalent for fungi → catalog must be hand-built | peer-reviewed | high |
| Yeast-resistance ML used metabolic/env, NOT DNA embeddings | 54–75 | % accuracy | ML identifies signatures of antifungal drug resistance in Saccharomycotina | (bioRxiv) | 2025 | Results | https://www.biorxiv.org/content/10.1101/2025.05.09.653161.full.pdf | 2026-06-07 | "532 yeast species ... eight antifungal drugs ... accuracies of 54-75%" | DNA-FM untried for fungal AMR; current SOTA is feature-based | preprint | medium |
| Arabidopsis flowering-time GWAS depth | 1003 | accessions (10°C) | 1,135 Genomes Reveal Global Pattern of Polymorphism in A. thaliana | (Cell / 1001 Genomes Consortium) | 2016 | Results | https://pmc.ncbi.nlm.nih.gov/articles/PMC4949382/ | 2026-06-07 | "flowering time GWAS ... four replicates in 1,003 (10°C) and 971 (16°C) lines" | deep sampling-independent quantitative phenotype | peer-reviewed | high |
| Arabidopsis genotype matrix scale | 10709466 | SNPs (2029 accessions) | AraGWAS Catalog | (Nucleic Acids Research) | 2018 | Methods | https://academic.oup.com/nar/article/46/D1/D1150/4559687 | 2026-06-07 | "SNP-Matrix for 2029 accessions on 10,709,466 segregating markers" | public, GWAS-ready; embedding-niche candidate (no curated catalog) | peer-reviewed | high |
| PlantCaduceus compute class | 3090 / A100 | GPU (24–80 GB) | Caduceus: Bi-Directional Equivariant Long-Range DNA Sequence Modeling | (Schiff et al., arXiv) | 2024 | Methods | https://arxiv.org/pdf/2403.03234 | 2026-06-07 | "3090, A5000, A6000, V100, and A100 GPUs"; "512 base pair windows" | plant DNA-FM needs ≥24GB GPU → NOT the GTX 860M (4GB) | peer-reviewed | high |
| Multimodal colony-image + WGS paired public dataset | not present | — | (web search synthesis, 2026-06-07) | — | 2026 | — | https://www.nature.com/articles/s41597-025-06319-4 | 2026-06-07 | "did not find a single dedicated public dataset that pairs colony images with WGS" | substrate-infeasible (no paired data) — iron-law kill | search-synthesis | medium |

## Highest-confidence rows (top 5)
1. Row 2 — C. auris S. Africa 181/188 MIC>32 ↔ ERG11: near-perfect genotype↔phenotype determinant linkage (the fungal-AMR transfer signal).
2. Row 5 — no genome-based fungal AMR predictor exists: NO AMRFinder-equivalent → catalog hand-built, but the niche is open.
3. Row 7 — Arabidopsis 1,003-accession flowering-time GWAS: deep, sampling-independent, no curated catalog → the embedding-niche candidate.
4. Row 4 — ERG11 alone insufficient (TAC1b/Cdr1/CNV): the fungal mechanism has the same efflux/regulatory blind spot as bacterial AMR.
5. Row 9 — PlantCaduceus needs ≥24GB GPU: the compute gate for the Arabidopsis/embedding path.

## Low-confidence rows
- Row 3 (12,644) + Row 6 (yeast ML) + Row 10 (multimodal absence): medium — preprint / search-synthesis, not direct-fetched.

## Honest gaps
- **Human T2D PRS (UK Biobank / All of Us):** not assessed in depth — both are access-gated (application/approval, not freely downloadable) + PRS has established baselines; lower solo-feasibility, deprioritized this pass.
- **Multimodal colony-image+WGS:** no public paired dataset found (Kaggle/Zenodo/HuggingFace not exhaustively searched). Treat as infeasible until a paired set is located.
- **Fungal determinant catalog completeness:** ERG11/FKS1/TAC1b mutation lists are documented across papers but NOT consolidated into one machine-readable catalog — building it is the main prerequisite for a deterministic fungal-AMR decoder.
- VRAM exact figures for PlantCaduceus inference not pinned (GPU class only).
