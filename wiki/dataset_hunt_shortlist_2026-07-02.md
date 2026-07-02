# Dataset hunt — ranked shortlist (F3, 2026-07-02)

Scored by `scripts/dataset_candidate_scorecard.py` (the F1 8-gate rubric). PASS = all gates pass + a viable decoder paradigm; VERIFY = no fail but gaps to confirm; REJECT = a gate fails or no paradigm. Gate legend: **G1_accessible** free + fetchable, no DUA/paywall; **G2_non_circular** phenotype measured independently, not derived from a tool the decoder competes with; **G3_sampling_independent** phenotype not confounded with sampling context (site/source/study/date); **G4_unit_joinable** per-individual genotype joined to per-individual phenotype (not aggregate); **G5_provenance_separable** a leakage-free split exists (temporal/accession/cohort); **G6_depth_or_catalog** >=~100 same-organism units (learned) OR a curated determinant catalog exists; **G7_genotype_fetchable** actual sequence/variants per unit are downloadable; **G8_label_not_censored** quantitative labels are tierable (not all operator-censored at one bound).

| rank | candidate | creature | verdict | decoder | depth | fails | unknowns |
|---|---|---|---|---|---|---|---|
| 1 | 1002 Yeast Genomes (Peter 2018) | S. cerevisiae (yeast) | **PASS** | learned-niche | 1011 | — | — |
| 2 | DGRP2 (Drosophila Genetic Reference Panel) | D. melanogaster (fruit fly) | **PASS** | learned-niche | 205 | — | — |
| 3 | Rice 3000 Genomes + IRRI phenotypes | Oryza sativa (rice) | **VERIFY** | neither | 3000 | — | G3_sampling_independent |
| 4 | CaeNDR (Caenorhabditis Natural Diversity Resource) | C. elegans (+briggsae/tropicalis) worm | **VERIFY** | learned-niche | 300 | — | G4_unit_joinable |
| 5 | ClinVar | Homo sapiens | **VERIFY** | deterministic | catalog | — | G2_non_circular, G5_provenance_separable |
| 6 | Mouse Phenome Database / Collaborative Cross | Mus musculus (mouse) | **VERIFY** | learned-niche | 100 | — | G4_unit_joinable, G6_depth_or_catalog, G7_genotype_fetchable |
| 7 | Arabidopsis flowering-time (1001G + AraPheno) | A. thaliana (plant) | **REJECT** | learned-niche | 1003 | G5_provenance_separable | — |

## Notes + sources per candidate

- **1002 Yeast Genomes (Peter 2018)** (S. cerevisiae (yeast)) — growth/fitness across ~35 lab conditions (quantitative). VCF matrix + phenotype matrix both free (1002genomes.u-strasbg.fr; ENA PRJEB13017). Phenotype = controlled lab growth across conditions (common-garden -> sampling-independent). Deepest single-species substrate found. Sources: http://1002genomes.u-strasbg.fr/, https://www.nature.com/articles/s41586-018-0030-5, https://www.ebi.ac.uk/ena/browser/view/PRJEB13017
- **DGRP2 (Drosophila Genetic Reference Panel)** (D. melanogaster (fruit fly)) — 31 harmonized quantitative organismal phenotypes (12 studies; DGRPool aggregates more). 205 inbred lines, full genotypes public (dgrp2.gnets.ncsu.edu). Single source population (Raleigh NC) -> G5 met by line-level holdout, not geography. DGRPool = harmonized phenotypes. Sources: http://dgrp2.gnets.ncsu.edu/, https://www.nature.com/articles/nature10811, https://elifesciences.org/reviewed-preprints/88981
- **Rice 3000 Genomes + IRRI phenotypes** (Oryza sativa (rice)) — agronomic quantitative traits. G3 UNKNOWN: field/environment confound is the classic plant-GWAS risk (multi-site trials); need controlled/BLUP phenotypes. FROM-KNOWLEDGE, verify. Sources: https://snp-seek.irri.org/
- **CaeNDR (Caenorhabditis Natural Diversity Resource)** (C. elegans (+briggsae/tropicalis) worm) — quantitative wild-isolate traits (drug/toxin response etc.). MIT-licensed, AWS Open Data, per-strain VCF. G4 UNKNOWN: confirm a HOSTED per-strain trait corpus exists (the GWAS tool is BYO-phenotype; Andersen-lab published trait sets likely qualify). Sources: https://caendr.org/, https://registry.opendata.aws/caendr/, https://pmc.ncbi.nlm.nih.gov/articles/PMC10767927/
- **ClinVar** (Homo sapiens) — variant -> disease/pathogenicity (curated). Deterministic-path (curated catalog). G2 UNKNOWN: clinical assertions can be predictor-derived (ACMG uses in-silico tools) -> circularity risk, the exact wall the HIV cell had to dodge. Sources: https://www.ncbi.nlm.nih.gov/clinvar/
- **Mouse Phenome Database / Collaborative Cross** (Mus musculus (mouse)) — large measured phenome across strains. G4/G7 UNKNOWN: per-strain genotype<->phenotype join granularity + fetchable per-strain genotypes need confirming. FROM-KNOWLEDGE, verify. Sources: https://phenome.jax.org/
- **Arabidopsis flowering-time (1001G + AraPheno)** (A. thaliana (plant)) — flowering time (quantitative). CLOSED NEGATIVE for embeddings (G2 2026-06-12: embedding learned population structure, not the causal signal). G5 marked FAIL for the LEARNED path (structure not separable from signal). Deterministic Mendelian sub-traits only; do NOT re-run embeddings here. Sources: https://arapheno.1001genomes.org/

## Headline
Top substrate: **1002 Yeast Genomes (Peter 2018)** (PASS, learned-niche, depth 1011). It is the F4 pilot-fetch target.

Honest scope: the top 3 are WEB-VERIFIED (F2 2026-07-02); the rest are FROM-KNOWLEDGE first-pass reads with UNKNOWN on their risk gate — a follow-up sweep must confirm before they rank as PASS.
