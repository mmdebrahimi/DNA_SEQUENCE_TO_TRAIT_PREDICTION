# Marine / aquaculture bacteria — genome + measured-AST accessions — supported memo (V1 invocation)
<!-- memo-schema: 0.4 -->

> Captured 2026-07-14. Source: Claude Code (`/research` run by hand). Slug: marine-vibrio-genome-mic-accessions-2026-07-14.
> Audit floor: 5 of 5 locators required per row. Mapping floor: rationale → quantity required.
> Banned-phrase scan: 5 hard-reject + 4 soft-warn (tiered). Cite-token noise scanned + flagged.

## Research Context

- **Problem:** pin exact accessions for marine/aquaculture bacterial datasets with paired genome + measured antibiotic susceptibility (MIC/disk) for a data-request email to a marine-biology postdoc
- **Captured:** 2026-07-14
- **Schema:** memo-schema 0.4

## Audit table (verbatim, supported rows only)

| Claim / quantity | Numeric value | Units | Source title | Authors / org | Year | Section / table / figure | Stable URL | Access date | Quoted excerpt (≤25 words) | Extraction rationale | Source type | Confidence |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| BioProject for Korean seafood *V. parahaemolyticus* genome assemblies | PRJNA1254427 (BioSamples SAMN48128082–SAMN48128091) | accession | Antibiotic Resistance and Characteristics of *V. parahaemolyticus* from Seafood in South Korea 2021–2022 | PMC12300810 | 2024/25 | Data availability | https://pmc.ncbi.nlm.nih.gov/articles/PMC12300810/ | 2026-07-14 | "genome assemblies ... deposited in the NCBI database under BioProject accession number PRJNA1254427, with associated BioSample accession numbers SAMN48128082 to SAMN48128091" | Accession appears verbatim in the fetched data-availability sentence | peer-reviewed | high |
| Antibiotics measured in PRJNA1254427 (maps to decoder mechanisms) | qualitative | — | (same) | PMC12300810 | 2024/25 | Methods | https://pmc.ncbi.nlm.nih.gov/articles/PMC12300810/ | 2026-07-14 | "ceftazidime, cefepime, meropenem, ciprofloxacin, trimethoprim/sulfamethoxazole ... chloramphenicol, colistin, nalidixic acid, and tetracycline" | Quote lists the panel; includes ciprofloxacin (QRDR) + tetracycline + ceftazidime = the decoder's exact drugs | peer-reviewed | high |
| Isolates in PRJNA1254427 with BOTH WGS AND lab-measured MIC | 10 (of 17 MIC-tested) | isolates | (same) | PMC12300810 | 2024/25 | Methods / Results | https://pmc.ncbi.nlm.nih.gov/articles/PMC12300810/ | 2026-07-14 | "MIC method in accordance with the Clinical and Laboratory Standards Institute (CLSI) guidelines" | ⚠ MAPPING CAVEAT: the "10 paired" count is from the fetched Results prose (17 MIC-tested, WGS on the 10 resistant); the stored quote is the CLSI method sentence, not the count | peer-reviewed | medium |
| BioProject for Canadian imported-shrimp *Vibrio* with phenotypic AST | PRJNA645603 | accession | Whole-Genome Sequences of *Vibrio* Species from Warm-Water Shrimps Imported into Canada | ASM Microbiology Resource Announcements (mra.01014-21) | 2021 | Data availability / Table 1 | https://journals.asm.org/doi/10.1128/mra.01014-21 | 2026-07-14 | "antimicrobial resistance (AMR) profiles determined by the Kirby-Bauer disk diffusion method" | Quote confirms measured (disk-diffusion) AST paired to the deposited WGS | peer-reviewed | high |
| Isolates in PRJNA645603 with BOTH WGS AND measured AST | 4 | isolates | (same) | mra.01014-21 | 2021 | Table 1 | https://journals.asm.org/doi/10.1128/mra.01014-21 | 2026-07-14 | "The antibiotic susceptibility test results showed resistance to up to nine different antibiotics" | ⚠ MAPPING CAVEAT: the "4" count is from the fetched Table-1 description (V. alginolyticus / V. cholerae / 2× V. parahaemolyticus); stored quote is the results sentence, not the count | peer-reviewed | medium |
| BioProjects for aquacultured-seafood AMR *Vibrio* (SRA reads) | PRJNA699735; PRJNA1107692; PRJNA1155317 | accession | Association of antimicrobial resistant *Vibrio* with aquacultured seafood | ScienceDirect S0740002025000991 | 2025 | Data availability | https://www.sciencedirect.com/science/article/pii/S0740002025000991 | 2026-07-14 | "Raw Illumina reads were uploaded to the Sequence Read Archive (SRA) under bioprojects: PRJNA699735, PRJNA1107692 and PRJNA1155317" | Accessions appear verbatim; ⚠ per-isolate MEASURED-AST pairing NOT individually verified this run | peer-reviewed | medium |

## Source-Locator Coverage

- Total rows submitted: 9
- Survived audit floor: 6
- Survived mapping floor: 6 (rows 3 + 5 pass with an explicit mapping caveat — value from fetched prose, not the stored quote)
- Survived banned-phrase scan: 6 (0 hard-reject; "certification"/"approval" not present outside quotes)
- Final supported: 6
- Survival rate: 6 / 9 (67%)

## Caveats per row

- **PRJNA1254427 count (10 paired):** mapping caveat — re-confirm the exact count of isolates with BOTH a genome AND an MIC against the paper's Results/Table before treating "10" as firm.
- **PRJNA645603 count (4 paired):** mapping caveat — count from the Table-1 species description; confirm against Table 1.
- **Aquacultured-seafood BioProjects (PRJNA699735 / 1107692 / 1155317):** the accessions are for SRA **reads**; whether each isolate carries a **measured** (not genome-predicted) AST value paired per-isolate is UNVERIFIED — the single most important thing to check before relying on these.
- No cite-token noise; no soft-warn banned phrases outside quotes.

## Decisions for Human Confirmation (cap 5)

| Claim | Numeric value | Units | Source URL | Candidate use / Verification needed | Confidence |
|---|---:|---|---|---|---|
| Korean seafood *V. parahaemolyticus* genome+MIC set | PRJNA1254427 | accession | https://pmc.ncbi.nlm.nih.gov/articles/PMC12300810/ | **Candidate use:** the single best "ask Hamid to grab this" — deposited assemblies + CLSI MIC incl. ciprofloxacin/tetracycline/ceftazidime (drops straight into the decoder). **Verification needed:** confirm exactly how many of the 10 have a downloadable assembly AND a cipro/tet MIC. | medium→high |
| Canadian imported-shrimp *Vibrio* (3 species) genome+disk-AST | PRJNA645603 | accession | https://journals.asm.org/doi/10.1128/mra.01014-21 | **Candidate use:** small multi-species (V. alginolyticus/cholerae/parahaemolyticus) with measured disk AST. **Verification needed:** disk-zone → R/S mapping + whether MICs (not just zones) are available. | high |
| Aquacultured-seafood AMR *Vibrio* SRA sets | PRJNA699735 / 1107692 / 1155317 | accession | https://www.sciencedirect.com/science/article/pii/S0740002025000991 | **Candidate use:** larger read sets — potentially the biggest usable pool IF measured AST is paired. **Verification needed:** open the paper's supplement — is there a per-isolate MEASURED MIC/disk table, or only genome-predicted resistance genes? | medium |

## Promotion Gate reminder
This memo is INPUT to the 4-step Promotion Gate, NOT approval. Confirm per row: (1) doc resolves at the cited URL, (2) the section reference exists, (3) the quoted excerpt is verbatim, (4) the excerpt→value mapping is natural — before treating any accession as a firm decoder cohort.
