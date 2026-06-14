# Acquirable E. coli phenotype-label sources — supported memo (V1 invocation)
<!-- memo-schema: 0.4 -->

> Captured 2026-06-14. Source: Claude (`/research` orchestrator → `/research-intake`). Slug: acquirable-ecoli-phenotype-label-sources-2026-06-14.
> Audit floor: 5 of 5 locators per row. Mapping floor: rationale → quantity. Banned-phrase scan: clean. Cite-token noise: none.
> v0.4 source-text identity advisory APPLIED: all rows derived from WebSearch result summaries (direct WebFetch of medRxiv/Lancet was HTTP-403 blocked this session) → provenance caveat on every row; Oxford-cohort rows downgraded high→medium pending direct-source verification.

## Research Context

- **Problem:** Ranked shortlist of acquirable E. coli phenotype label sources scored on 4 non-negotiable criteria (sampling-independent, non-circular, organism-depth ≥100-150 with downloadable assemblies, provenance-disjoint-feasible) + acquisition path + N-after-filters, to reopen the frozen v0.5.0 AMR decoder beyond the saturated public-label AMR grid.
- **Captured:** 2026-06-14
- **Schema:** memo-schema 0.4

## Audit table (verbatim, supported rows only)

| Claim / quantity | Numeric value | Units | Source title | Authors / org | Year | Section | Stable URL | Access date | Quoted excerpt (≤25 words) | Source type | Confidence |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Oxford E. coli cohort: isolates with linked WGS + measured MIC | 2875 | isolates | Estimating the association of AMR genes with MIC in E. coli | Lipworth et al. (Oxford / Lancet Microbe) | 2025 | Methods | https://www.thelancet.com/journals/lanmic/article/PIIS2666-5247(25)00111-9/fulltext | 2026-06-14 | "2875 isolates with linked whole genome sequencing data and MIC-level phenotyping for at least one antibiotic" | peer-reviewed | medium |
| Oxford cohort: blood-culture isolates, OUH UK, 2013-2018 | 2410 | isolates | (same) | Lipworth et al. | 2025 | Methods | https://www.medrxiv.org/content/10.1101/2024.05.15.24307162v1 | 2026-06-14 | "2410 were isolated from blood cultures ... Oxford University Hospitals NHS Foundation Trust (OUH) between 2013-2018" | peer-reviewed | medium |
| Oxford cohort: urine isolates, OUH, 2020 | 465 | isolates | (same) | Lipworth et al. | 2025 | Methods | https://www.medrxiv.org/content/10.1101/2024.05.15.24307162v1 | 2026-06-14 | "465 isolates ... all urine cultures sent to the OUH laboratory which were positive for E. coli ... 2020" | peer-reviewed | medium |
| Oxford cohort: antibiotics with measured MIC | 8 | drugs | (same) | Lipworth et al. | 2025 | Methods | https://www.thelancet.com/journals/lanmic/article/PIIS2666-5247(25)00111-9/fulltext | 2026-06-14 | "ampicillin, amoxicillin-clavulanate, ceftriaxone, cefuroxime, ciprofloxacin, gentamicin, piperacillin-tazobactam, and trimethoprim-sulfamethoxazole" | peer-reviewed | medium |
| Oxford cohort: MIC is measured graded value (interval regression) | qualitative | — | (same) | Lipworth et al. | 2025 | Methods | https://www.thelancet.com/journals/lanmic/article/PIIS2666-5247(25)00111-9/fulltext | 2026-06-14 | "multivariable interval regression models to estimate the change in MIC ... predict MIC and binary resistance/susceptibility" | peer-reviewed | medium |
| ecoref panel: E. coli strains with growth phenotyping | 894 | strains | Phenotype inference in an E. coli strain panel | Galardini et al. (eLife) | 2017 | Results | https://elifesciences.org/articles/31035 | 2026-06-14 | "a large genetic reference panel of 894 E. coli strains able to capture the genetic and phenotypic diversity of the species" | peer-reviewed | medium |
| ecoref: growth measured across conditions | 214 | conditions | (same) | Galardini et al. | 2017 | Results | https://elifesciences.org/articles/31035 | 2026-06-14 | "measured growth across 214 conditions, as well as obtained the genomic sequences for the majority of the strains (N = 696)" | peer-reviewed | medium |
| ecoref: strains with public genome sequence | 696 | genomes | (same) | Galardini et al. | 2017 | Results | https://elifesciences.org/articles/31035 | 2026-06-14 | "obtained the genomic sequences for the majority of the strains (N = 696)" | peer-reviewed | medium |
| von Mentzer ETEC: preassembled genomes virulence-profiled | 1083 | genomes | ETECFinder database construction | Joffré / CGE et al. | 2023 | Methods | https://journals.asm.org/doi/10.1128/jcm.00570-23 | 2026-06-14 | "1,083 preassembled genomes (sequenced by Astrid von Mentzer) ... 890 ETEC matches ... 193 non-ETEC matches" | peer-reviewed | medium |
| von Mentzer ETEC: confirmed public BioProject subset | 439 | genomes | (same) | Joffré / CGE et al. | 2023 | Methods | https://journals.asm.org/doi/10.1128/jcm.00570-23 | 2026-06-14 | "the two BioProjects PRJNA421191 with 305 and PRJNA416134 with 134 sequences" | peer-reviewed | medium |
| Pfizer ATLAS (Vivli): measured MICs, NO genome assemblies | 6500000 | MICs | Seeking patterns of antibiotic resistance in ATLAS | Catalán et al. / Vivli-Pfizer (Nat Commun) | 2022 | Abstract | https://www.nature.com/articles/s41467-022-30635-7 | 2026-06-14 | "ATLAS holds 6.5M minimal inhibitory concentrations (MICs) for 3,919 pathogen-antibiotic pairs ... 633k patients ... 2004 and 2017" | peer-reviewed | medium |
| NARMS: real broth-microdilution MIC + NCBI-linked genomes | qualitative | — | NARMS Technical Report | FDA/CDC/USDA NARMS | 2023 | Methods | https://www.fda.gov/media/164290/download | 2026-06-14 | "NARMS uses standardized broth microdilution for MIC determination ... uploads isolate data to NCBI" | regulatory | high |

## Source-Locator Coverage

- Total rows submitted: 13
- Survived audit floor: 12
- Survived mapping floor: 12
- Survived banned-phrase scan: 12
- Final supported: 12 (1 rejected: Oxford-accession row, quote = "VALUE NOT PRESENT" + low confidence)
- Survival rate: 12 / 13 (92%)

## Caveats per row

- **ALL rows — provenance: websearch-summary.** Quotes are the search engine's rendering of the source text, not a direct WebFetch read (medRxiv + Lancet Microbe both returned HTTP 403 this session). Re-verify each quote against the direct source before any high-confidence use. This is why all Oxford/ecoref/von-Mentzer/ATLAS rows are capped at `medium`.
- **Oxford rows** — the cohort's exact public-deposit accession is the load-bearing unknown (see unsupported memo). Criterion-3 ("downloadable assemblies") is UNCONFIRMED until the data-availability statement is read directly.
- **ATLAS row** — the claim is well-supported, but the source is a DISQUALIFIER: ATLAS provides genotype markers "only when tested," NOT whole-genome assemblies → fails criterion 3 (no assemblies at scale).
- **NARMS row** — well-supported, but NARMS IS the surveillance ecosystem the project excludes for provenance-disjointness → fails criterion 4 by definition (it is the in-cohort/tuning source, not a disjoint validation source).
- **von Mentzer rows** — ETEC toxin/pathotype labels are partly molecularly defined; closer to the study==class / sampling-defined trap than a pure lab measurement. Medium-value; needs the per-strain label-provenance check the project already applied to pathotype.

## Scorecard — the actual deliverable (4-criteria pass/fail per candidate)

| Candidate | (1) sampling-indep | (2) non-circular | (3) ≥100-150 + downloadable assemblies | (4) provenance-disjoint-feasible | Clears ALL 4? | Acquisition path | N after filters (est.) | Biggest risk |
|---|---|---|---|---|---|---|---|---|
| **Oxford 2875 WGS+MIC cohort** | ✅ measured clinical MIC | ✅ wet-lab, not genome-derived | ⚠️ N huge; **assembly-deposit UNCONFIRMED** | ✅ UK OUH, outside US ecosystem | **PROVISIONAL YES** (gated on assembly deposit) | free download IF deposited; else data-request to OUH/Oxford | up to ~2875 (per drug ≥100s) if public; ~0 if controlled-access | reads may be controlled-access (clinical) → drops to MTA |
| **ecoref / eLife panel (894 / 696 genomes / 214 growth conditions)** | ✅ wet-lab growth assay | ✅ phenotype is growth, not genome-derived | ✅ 696 public genomes, ecoref repo | ⚠️ diverse natural isolates; per-label ≥20/class + clonality UNCHECKED | **PROVISIONAL YES** for a NEW trait | free download (ecoref repository) | ~696 genomes; per-growth-condition class balance TBD | a continuous growth phenotype may not be cleanly decodable / may track lineage (the embedding-0-for-4 risk) |
| von Mentzer ETEC (≥439 confirmed public, ~1083 total) | ⚠️ toxin type partly molecular | ⚠️ partly molecular label | ✅ ≥439 public genomes | ⚠️ pathotype = sampling-adjacent | NO (criteria 1/2 weak) | free download (PRJNA421191 + PRJNA416134) | ~439-1083 | re-introduces the pathotype label-circularity trap already closed |
| Pfizer ATLAS / Vivli (6.5M MIC) | ✅ measured MIC | ✅ measured | ❌ NO genome assemblies | n/a | **NO** | Vivli application | ~0 (no assemblies) | structural: phenotype without genomes |
| NARMS raw MIC | ✅ measured BMD | ✅ measured | ✅ NCBI-linked | ❌ IS the excluded ecosystem | **NO** | free (FDA/NCBI) | ~0 disjoint | it is the tuning source, not disjoint |

## Decisions for Human Confirmation (cap 5)

| Claim | Numeric value | Units | Source URL | Candidate use / Verification needed | Confidence |
|---|---:|---|---|---|---|
| Oxford E. coli cohort with linked WGS + measured clinical MIC, UK hospitals | 2875 | isolates | https://www.thelancet.com/journals/lanmic/article/PIIS2666-5247(25)00111-9/fulltext | **Candidate use:** the strongest RE-VALIDATION substrate for the frozen decoder — measured MIC on cipro/cef/gent (+5 more), UK provenance outside the US surveillance ecosystem, large N. **Verification needed:** READ THE DATA-AVAILABILITY STATEMENT — confirm genome assemblies/reads are PUBLICLY downloadable (ENA/BioProject) vs controlled-access; if controlled, this becomes an MTA/data-request (user action). | medium |
| ecoref panel: growth phenotyping across conditions + public genomes | 696 | genomes | https://elifesciences.org/articles/31035 | **Candidate use:** a genuinely NEW decodable trait (growth/fitness) at depth with free public genomes — opens a non-AMR axis. **Verification needed:** pick ONE growth condition with ≥20/class balance, check it survives Mash-clonality + isn't lineage-confounded (the 0-for-4 embedding risk), confirm ecoref genome accessions resolve. | medium |
| von Mentzer ETEC reference genomes, toxin-typed | 439 | genomes | https://journals.asm.org/doi/10.1128/jcm.00570-23 | **Candidate use:** a pathotype reference set with public genomes. **Verification needed:** confirm the toxin label is WET-LAB (not molecular gene-call) per strain — if molecular, it re-enters the closed pathotype-circularity trap and should be dropped. | medium |
| Pfizer ATLAS measured MIC (DISQUALIFIED) | 6500000 | MICs | https://www.nature.com/articles/s41467-022-30635-7 | **Candidate use:** none for this decoder. **Verification needed:** none — confirmed it provides genotype markers only, no assemblies → fails criterion 3. Recorded so it is not re-investigated. | medium |
| NARMS raw broth-microdilution MIC (DISQUALIFIED for disjoint) | qualitative | — | https://www.fda.gov/media/164290/download | **Candidate use:** none as a DISJOINT source — it is the surveillance ecosystem already excluded. **Verification needed:** none — fails criterion 4 by definition. Recorded so it is not re-investigated. | high |

## Promotion Gate reminder

This memo is INPUT to the 4-step Promotion Gate, NOT approval. The single highest-leverage verification (do this first): **read the Oxford cohort's data-availability statement** at the Lancet Microbe / medRxiv source to confirm whether the 2875 assemblies are publicly downloadable. That one fact decides whether the top candidate is a free-download re-validation (executor-doable) or an MTA/data-request (user action). Do NOT treat any source here as cleared until its quote is verified against the direct source (all rows are websearch-summary provenance).
