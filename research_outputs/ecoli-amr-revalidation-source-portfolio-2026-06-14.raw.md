# Independent E. coli AMR re-validation source portfolio (raw research, 2026-06-14)

> Captured 2026-06-14. Source: Claude Code (`/research` orchestrator, breadth run #2). Topic: "breadth of independent external E. coli WGS + AST cohorts with public genomes to re-validate the frozen v0.5.0 AMR decoder (multi-source validation portfolio)". Slug: ecoli-amr-revalidation-source-portfolio-2026-06-14.
> Web search via WebSearch (~6 calls). Goal = MANY independent cohorts (breadth), each provenance-disjoint from the decoder's US-NCBI-PD tuning. Criteria reused from acquirable-ecoli-phenotype-label-sources-2026-06-14: sampling-independent / non-circular / ≥100-150 same-organism w/ downloadable assemblies / provenance-disjoint. For BINARY re-validation, categorical S/I/R AST is sufficient; MEASURED-MIC cohorts additionally allow breakpoint-controlled labels.

## Audit table (verbatim, all candidate rows)

| Claim / quantity | Numeric value | Units | Source title | Authors / org | Year | Section | Stable URL | Access date | Quoted excerpt (≤25 words) | Source type | Confidence |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Oxford cohort: WGS + measured MIC (re-validation gold) | 2875 | isolates | AMR genes ↔ MIC in E. coli | Lipworth et al. | 2025 | Methods | https://www.thelancet.com/journals/lanmic/article/PIIS2666-5247(25)00111-9/fulltext | 2026-06-14 | "2875 isolates with linked whole genome sequencing data and MIC-level phenotyping" (genomes open at ENA PRJNA604975) | peer-reviewed | medium |
| Spain PROBAC: E. coli blood, broth-microdilution MIC, 16 drugs | 224 | isolates | WGS of E. coli bacteraemia (sepsis/septic shock) in Spain | Maldonado et al. (Lancet Microbe) | 2024 | Methods + Data availability | https://www.thelancet.com/journals/lanmic/article/PIIS2666-5247(23)00369-5/fulltext | 2026-06-14 | "Susceptibility testing was performed by broth microdilution ... deposited in the European Nucleotide Archive under the Bioproject PRJEB62601" | peer-reviewed | high |
| Spain PROBAC: drugs with MIC incl cipro/cef/gent | 16 | drugs | (same) | Maldonado et al. | 2024 | Methods | https://pubmed.ncbi.nlm.nih.gov/38547882/ | 2026-06-14 | "fosfomycin, amoxicillin–clavulanic acid, piperacillin–tazobactam, cefoxitin, cefotaxime, ceftriaxone, cefepime, cefuroxime ... ciprofloxacin, amikacin, gentamicin, trimethoprim–sulfamethoxazole" | peer-reviewed | high |
| 234-isolate European blood cohort: WGS + ISO broth-microdilution MIC, 11 drugs | 234 | isolates | Genotypic vs phenotypic resistance in 234 E. coli | (Sci Rep) | 2023 | Methods | https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9829913/ | 2026-06-14 | "234 E. coli isolates from positive blood cultures ... microdilution for 11 clinically relevant antibiotics ... ISO 20776-1 standard broth microdilution method as recommended by EUCAST" | peer-reviewed | high |
| Denmark "One day in Denmark": clinical E. coli WGS, public | 699 | isolates | One day in Denmark: WGS of E. coli from clinical settings | (J Antimicrob Chemother) | 2025 | Data availability | https://academic.oup.com/jac/article/80/4/1011/7989499 | 2026-06-14 | "Raw sequence data have been submitted to the European Nucleotide Archive ... under study accession no. PRJEB37711" (91.4% urine) | peer-reviewed | high |
| Netherlands bacteraemia E. coli WGS, public | 281 | isolates | ESBL vs non-ESBL E. coli bacteraemia in the Netherlands 2014-2016 | (BMC) | 2020 | Data availability | https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6959556/ | 2026-06-14 | "All generated raw reads were submitted to the European Nucleotide Archive ... under the study accession number PRJEB35000" (212 non-ESBL + 69 ESBL) | peer-reviewed | high |
| Norway NORM: nationwide E. coli BSI WGS | 22512 | BSI sampled | Emergence/dissemination of AMR in E. coli BSI in Norway 2002-17 | Gladstone et al. (Lancet Microbe) | 2021 | Methods | https://www.thelancet.com/journals/lanmic/article/PIIS2666-5247(21)00031-8/fulltext | 2026-06-14 | "isolates from 22,512 E coli bloodstream infections included in the Norwegian surveillance programme on resistant microbes (NORM) from 2002 to 2017" | peer-reviewed | medium |
| Jakarta BSI E. coli WGS, public BioProject | 22 | isolates | E. coli bloodstream infection WGS, Cipto Mangunkusumo Hospital Jakarta | (Indonesia) | 2020 | Data availability | https://pmc.ncbi.nlm.nih.gov/articles/PMC8529152/ | 2026-06-14 | "All FASTA raw files are available from BioProject accession number PRJNA596854 and Sequence Read Archive accession numbers SRR10761126–SRR10761147" | peer-reviewed | medium |
| Comparator: Scotland E. coli BSI WGS (prior search) | 162 | isolates | (Oxfordshire comparator set) | — | — | (citation) | https://www.ebi.ac.uk/ena/browser/view/PRJEB12513 | 2026-06-14 | "162 E. coli BSI isolates from Scotland (PRJEB12513)" | peer-reviewed | low |

## Highest-confidence rows (top 5)

1. Spain PROBAC (row 2-3) — 224 blood, EUCAST broth-microdilution MIC, 16 drugs, ENA PRJEB62601 CONFIRMED. Strongest *measured-MIC* non-UK addition.
2. 234-isolate European cohort (row 4) — ISO 20776-1 broth-microdilution MIC, 11 drugs, explicit genotype-vs-phenotype. Measured MIC.
3. Denmark PRJEB37711 (row 5) — 699 clinical (urine-heavy), public WGS. Large categorical-AST validation set.
4. Netherlands PRJEB35000 (row 6) — 281 blood, public WGS. Categorical-AST validation set.
5. Oxford PRJNA604975 (row 1) — the re-validation gold (measured MIC, already the top candidate from run #1).

## Low-confidence rows

- Norway NORM (row 7): huge, national, non-US, but the public WGS accession + whether isolate-level AST is MIC vs categorical needs the paper's data-availability statement (not pinned this pass).
- Jakarta PRJNA596854 (row 8): genomes public + accession confirmed, but small (22) and AST method unconfirmed.
- Scotland PRJEB12513 (row 9): surfaced as a comparator accession; size + AST detail not directly verified.

## Honest gaps

- **AST method per cohort** is confirmed measured-MIC only for Spain + the 234-set + Oxford. Denmark/Netherlands/Norway likely report categorical S/I/R (sufficient for BINARY decoder re-validation, but not breakpoint-controlled MIC). Confirm per cohort before treating as MIC sources.
- **Per-cohort downloadable-assembly vs raw-reads-only**: most deposit raw reads (ENA study accessions); assembling adds a step vs the decoder's assembly-accession fetch path. Confirm whether assemblies or only reads are deposited.
- **Sweden (Svarm) + a single combined NethMap/Svarm MIC+WGS resource** not found as one artifact — national surveillance reports don't link isolate-level MIC to public genomes in one place.
- **Leakage check**: each candidate's accessions must pass the project's `cohort_manifest` exact-identity gate vs the decoder's tuning/validation cohorts before use (almost certainly disjoint — all non-US — but must be verified, not assumed).
