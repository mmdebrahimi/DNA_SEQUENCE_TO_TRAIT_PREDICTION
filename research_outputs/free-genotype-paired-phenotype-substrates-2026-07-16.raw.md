# Free genotype-paired phenotype substrates for a deterministic decoder — beyond AMR (raw research, 2026-07-16)

> Captured 2026-07-16. Source: Claude Code (`/research` orchestrator, executed by hand per Soraya). Topic: "Do free, independent, genotype-paired phenotype substrates exist in unexploited modalities (visible traits / human quantitative traits / plant + yeast common-garden traits) that clear the dna_decode label wall AND fit the deterministic curated-catalog paradigm — path B of the VOI menu". Slug: free-genotype-paired-phenotype-substrates-2026-07-16.
> Web search via WebSearch + WebFetch, single-session skeptical extraction. "Thorough Google" depth.

## Audit table (verbatim, all candidate rows)

| Claim / quantity | Numeric value | Units | Source title | Authors / org | Year | Section | Stable URL | Access date | Quoted excerpt (≤25 words) | Extraction rationale | Source type | Confidence |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| PGS Catalog scale | 5,022 scores / 656 traits | scores/traits | The Polygenic Score Catalog (systematic-eval paper) | PGS Catalog / Lambert et al. | 2025 | abstract | https://www.pgscatalog.org/ | 2026-07-16 | "the database has included 5,022 polygenic risk scores associated with 656 different traits" | free deterministic-weight predictors across 656 traits | peer-reviewed/db | medium |
| PGS Catalog license | qualitative | — | PGS Catalog — About | PGS Catalog | 2026 | About/Terms | https://www.pgscatalog.org/about/ | 2026-07-16 | "used in accordance with any licensing restrictions set by the authors" | open DB, per-score author licensing (not universally open) | database | high |
| PGS scoring-file contents | qualitative | — | PGS Catalog — About | PGS Catalog | 2026 | About | https://www.pgscatalog.org/about/ | 2026-07-16 | "weights/effect sizes, effect allele, genome build" | each score = a deterministic apply-from-VCF formula (same shape as HIrisPlex-S) | database | high |
| PGS portability loss (Martin) | 1.6× / 1.7× / 2.5× / 4.9× lower | rel. accuracy | Polygenic scores in biomedical research | Martin et al. (via review PMC9391275) | 2019 | — | https://pmc.ncbi.nlm.nih.gov/articles/PMC9391275/ | 2026-07-16 | "1.6-fold lower in Hispanic/Latino... 2.5-fold lower in East Asians, and 4.9-fold lower in Africans" | PGS accuracy decays with genetic distance = ancestry/lineage confound | peer-reviewed | medium |
| PGS predictive-power drop | 24% / 37% / 51% of European | % power | Dissecting Predictive Accuracy of Polygenic Indexes across Ancestries | PMC12622163 | 2025 | results | https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12622163/ | 2026-07-16 | "lowest in African (24%), followed by East Asian (37%) and South Asian (51%)" | the polygenic space is intrinsically lineage-confounded | peer-reviewed | medium |
| AraPheno Arabidopsis | 1135 accessions | accessions | AraPheno — 1001 Genomes flowering study | Seren/Korte et al. | 2017 | study/12 | https://arapheno.1001genomes.org/study/12/ | 2026-07-16 | "study of 1135 Arabidopsis accessions from the 1001 Genomes Consortium" | free common-garden flowering phenotype; CSV/PLINK/JSON | peer-reviewed/db | medium |
| Yeast 1011 phenome | 1011 genomes / 223 traits | strains/traits | Predicting natural variation in the yeast phenotypic landscape | PMC12583546 / bioRxiv | 2024 | dataset | https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12583546/ | 2026-07-16 | "223 traits measured across 1011 genome-sequenced Saccharomyces cerevisiae strains" | free deep lab-measured growth/life-cycle phenome | peer-reviewed | medium |
| Y1000+ yeast | 1154 strains / 24 conditions | strains/conditions | Y1000+ Project | Wisc./Hittinger | 2026 | — | https://y1000plus.wei.wisc.edu/ | 2026-07-16 | "high-throughput quantitative phenotyping in 24 growth conditions across... 1,154 yeast strains" | free cross-species quantitative phenotyping | database | low |
| HIrisPlex-S pigmentation | 41 SNPs, eye AUC ~0.9 | SNPs/AUC | HIrisPlex-S developmental validation | Walsh et al., FSI:Genetics | 2018 | — | https://www.sciencedirect.com/science/article/abs/pii/S1872497318302205 | 2026-07-16 | "17 skin colour predictive SNPs and the previous HIrisPlex assay for 24 eye and hair colour predictive SNPs" | curated causal-locus catalog + free coefficients = the WINNING paradigm | peer-reviewed | medium |
| openSNP validation cohort | 6401 genotypes / 668 phenotypes | files/phenotypes | openSNP | Greshake et al. | 2026 | — | https://opensnp.org/ | 2026-07-16 | "6,401 genotype files across 668 phenotypes" incl. eye/hair colour | free independent genotype+observed-phenotype cohort (self-reported) | database | medium |

## Highest-confidence rows (top 5)

1. Row 2 (PGS license, direct fetch) — the free-polygenic space is real and open-ish.
2. Row 3 (PGS scoring-file contents, direct fetch) — each PGS is a deterministic apply-from-VCF formula, same shape as the project's cells.
3. Row 5 (PGS ancestry power-drop) — the polygenic space carries the project's own 0-for-5 lineage confound baked in.
4. Row 9 (HIrisPlex-S) — the one free substrate that fits the winning curated-catalog paradigm.
5. Row 10 (openSNP) — the free independent validation cohort for the pigmentation cell.

## Low-confidence rows

- Row 8 (Y1000+ 1154/24) — from search summary, not fetched; directional only.

## Honest gaps

- **No fetched live count** for PGS Catalog scores/traits (About page omits the number; the 5,022/656 is from a dated 2025 publication via search summary → medium).
- **PGS per-score licensing not enumerated** — "open DB" is true at the catalog level; individual-score reuse depends on author terms (not audited per-score here).
- **Arabidopsis flowering deterministic-locus decoding UNTESTED** by this project — the embedding test was a closed NEGATIVE (2026-06-12); whether a curated FLC/FRI/FT-locus deterministic decoder works is not established either way (an open question, not a found substrate).
- **Yeast growth-trait curated causal-locus catalogs** — largely ABSENT (growth is organism-polygenic); the drug-resistance subset MAY have known loci (analogous to AMR) but was not enumerated here.
