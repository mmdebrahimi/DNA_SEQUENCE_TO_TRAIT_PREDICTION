# C. auris WGS + fluconazole-MIC cohort sources (for Gate G1) — raw research, 2026-06-07

> Captured 2026-06-07. Source: Claude Code (`/research` orchestrator). Topic: "Candida auris clinical
> isolates per-isolate fluconazole MIC with NCBI genome assembly accession and clade from published
> whole-genome-sequencing supplementary tables (South Africa outbreak and India cohorts)". Slug:
> cauris-wgs-mic-cohort-sources-2026-06-07.
> Web search via Claude Code WebSearch + WebFetch — single-session sourcing reconnaissance for the EP-7
> Gate-G1 cohort. GOAL was a buildable per-isolate [genome_accession, fluconazole_MIC, clade] table.

## Headline finding (load-bearing for G1 scoping)

A *buildable* G1 cohort needs three things joined per isolate: (a) a downloadable **assembled** genome,
(b) a fluconazole **MIC**, (c) a **clade**. The published S.Africa cohorts provide (b)+(c) richly but the
genomes are deposited as **SRA raw reads, not assemblies** — so G1 is **compute-gated** (de-novo assembly
of ~90-188 fungal genomes), NOT the no-compute supplementary-fetch the EP-7 plan assumed. This refines the
EP-7 substrate assumption.

## Audit table (candidate sources)

| Claim / quantity | Numeric value | Units | Source title | Authors / org | Year | Section | Stable URL | Access date | Quoted excerpt (<=25 words) | Extraction rationale | Source type | Confidence |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| S.Africa neonatal outbreak isolates WGS-sequenced | 188 | isolates | Candida auris Clinical Isolates ... Neonatal Unit ... South Africa | Govender et al. (NICD) / Emerg Infect Dis | 2023 | Results | https://pmc.ncbi.nlm.nih.gov/articles/PMC10521600/ | 2026-06-07 | "188 (66%) were processed for whole-genome sequencing" | the survey-cited 188-isolate cohort | peer-reviewed | high |
| Isolates with fluconazole MIC>32 that had ERG11 mutations | 181/188 | isolates | (same PMC10521600) | (same) | 2023 | Results/Table 3 | https://pmc.ncbi.nlm.nih.gov/articles/PMC10521600/ | 2026-06-07 | "all 181/188 isolates that had a fluconazole MIC >32 had ERG11 mutations" | confirms near-perfect ERG11<->MIC linkage | peer-reviewed | high |
| PMC10521600 per-isolate accession+MIC table availability | qualitative | - | (same PMC10521600) | (same) | 2023 | Appendix | https://pmc.ncbi.nlm.nih.gov/articles/PMC10521600/ | 2026-06-07 | "Appendix ... 23-0181-Techapp-s1.pdf" | only AGGREGATE tables in main text; no BioProject for the 188; appendix is PDF | peer-reviewed | high |
| Bloodstream-infection C. auris WGS isolates | 92 | isolates | In Vitro Antifungal Resistance ... Bloodstream Infections, South Africa | van Schalkwyk et al. / AAC | 2021 | Results | https://pmc.ncbi.nlm.nih.gov/articles/PMC8370198/ | 2026-06-07 | "the 92 C. auris isolates ... BioProject accession number PRJNA737309" | the ONLY cohort with a stated BioProject | peer-reviewed | high |
| PRJNA737309 deposited as assemblies vs raw reads | 0 assemblies / 115 SRA | records | NCBI Assembly + SRA (eutils) | NCBI | 2026 | eutils esearch | https://www.ncbi.nlm.nih.gov/bioproject/PRJNA737309 | 2026-06-07 | (eutils: assembly Count=0; sra Count=115) | THE blocker: caller needs assembled FASTA; reads need de-novo assembly (compute) | database | high |
| PMC8370198 per-isolate MIC supplementary | qualitative | - | (same PMC8370198) | (same) | 2021 | Table S1 | https://pmc.ncbi.nlm.nih.gov/articles/instance/8370198/bin/aac.00517-21-s0001.pdf | 2026-06-07 | "Table S1: ... fluconazole MICs for clade III isolates with ... mutations" | per-isolate MIC exists but in a PDF; clade-III subset | peer-reviewed | medium |
| National clade-distribution WGS genomes | 115 | genomes | Clade distribution of C. auris in South Africa ... WGS | Maphanga et al. | 2021 | Results | https://pmc.ncbi.nlm.nih.gov/articles/PMC8253216/ | 2026-06-07 | "WGS SNP analysis of 115 South African C. auris genomes" | another candidate cohort; clade-focused (MIC less central) | peer-reviewed | medium |
| Clade-specific ERG11 azole-R substitutions | qualitative | - | Simultaneous Emergence ... 3 Continents | Lockhart et al. / Clin Infect Dis | 2017 | Results | https://pmc.ncbi.nlm.nih.gov/articles/PMC5215215/ | 2026-06-07 | "F126T with South Africa, Y132F with Venezuela, Y132F or K143R with India and Pakistan" | confirms the catalog's clade<->mutation map | peer-reviewed | high |

## Highest-confidence rows (top 5)
1. Row 1 (188 WGS) + Row 2 (181/188 MIC>32 w/ ERG11) — the survey's headline linkage, primary source confirmed.
2. Row 4 (PRJNA737309) — the only cohort with a stated, downloadable BioProject.
3. Row 5 (PRJNA737309 = 0 assemblies / 115 SRA) — eutils-verified; THE structural blocker for no-compute G1.
4. Row 8 (Lockhart clade<->mutation map) — independently confirms the fungal_amr.py catalog numbering.

## Honest gaps
- **No single accessible table joins [assembled genome accession + fluconazole MIC + clade] per isolate.**
  PMC10521600 = aggregate-only + no BioProject for the 188; PMC8370198 = BioProject (raw reads) + MIC in a
  PDF; the accession<->MIC join is manual either way.
- **India cohort (350):** not separately resolved this pass (search centered on the S.Africa outbreak that
  the substrate survey headlined). A dedicated India-cohort search is a follow-up if the S.Africa path is
  abandoned.
- **Assembled + MIC-labeled C. auris genomes:** not found. The 841 NCBI assemblies (ledger row 2) do NOT
  carry fluconazole MIC in BioSample metadata (verified earlier this session).

## G1 path options (consequence of the above)
1. **Compute path:** download PRJNA737309 SRA reads -> de-novo assemble (~30 min/genome) -> parse Table S1
   PDF for per-isolate MIC -> fuzzy-join isolate IDs -> run `build_fungal_cohort.py`. Modest compute, real
   manual curation. No longer "no-compute".
2. **User-curated path:** user supplies a per-isolate [isolate_id, assembly_accession, fluconazole_mic,
   clade] TSV (e.g. from an authoritative supplementary they can extract) -> one-command G1.
3. **Assembled-genome path (cheapest if it exists):** find any C. auris study depositing GCA assemblies
   WITH per-isolate MIC. Not found this pass; low probability given MIC isn't in NCBI metadata.
