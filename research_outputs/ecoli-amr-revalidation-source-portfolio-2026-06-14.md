# Independent E. coli AMR re-validation source PORTFOLIO — supported memo (V1 invocation)
<!-- memo-schema: 0.4 -->

> Captured 2026-06-14. Source: Claude (`/research` breadth run #2 → intake). Slug: ecoli-amr-revalidation-source-portfolio-2026-06-14.
> Audit floor 5/5 + mapping floor + banned-phrase scan: clean. Provenance: all rows websearch-summary (direct fetch of Lancet/medRxiv 403-blocked) → verify quotes against direct source before any uplift.
> Goal = BREADTH (user-confirmed): a multi-source portfolio of INDEPENDENT external cohorts to re-validate the frozen v0.5.0 AMR decoder, not a single pilot.

## Research Context

- **Problem:** Build a breadth portfolio of independent external E. coli WGS + AST cohorts (public genomes, provenance-disjoint from the decoder's US-NCBI-PD tuning) to re-validate the frozen v0.5.0 AMR decoder across many labs/countries — a "validated across N independent sources" trust claim.
- **Captured:** 2026-06-14
- **Schema:** memo-schema 0.4

## The portfolio (ranked; all clear sampling-independent + non-circular + provenance-disjoint)

| # | Cohort | N (E. coli) | AST type | Genomes (accession) | Country | Tier | Biggest risk |
|---|---|---|---|---|---|---|---|
| 1 | **Oxford** (Lipworth) | 2875 | **measured MIC**, 8 drugs | ENA **PRJNA604975** (open, confirmed) | UK | **MIC — gold** | MIC table location (supplement vs extract) |
| 2 | **Spain PROBAC** (Maldonado 2024) | 224 | **measured MIC** (EUCAST broth microdilution), 16 drugs incl cipro/cef/gent | ENA **PRJEB62601** (open, confirmed) | Spain | **MIC** | sepsis-only sampling; per-isolate MIC table availability |
| 3 | **234-isolate European** (Sci Rep 2023) | 234 | **measured MIC** (ISO 20776-1 BMD), 11 drugs | accession TBC (PMC9829913 data stmt) | Europe | **MIC** | accession + country not pinned this pass |
| 4 | **Denmark "One day in Denmark"** | 699 | categorical AST (likely) | ENA **PRJEB37711** (open, confirmed) | Denmark | binary | urine-heavy (91%); AST is S/I/R not MIC (binary-validation only) |
| 5 | **Netherlands bacteraemia** | 281 | categorical AST (likely) | ENA **PRJEB35000** (open, confirmed) | Netherlands | binary | ESBL-enriched selection; categorical AST |
| 6 | **Norway NORM BSI** (Gladstone 2021) | ~thousands (22,512 sampled) | NORM surveillance AST | accession TBC | Norway | binary (large) | national-surveillance program (still non-US/disjoint); accession + MIC-vs-categorical unconfirmed |
| 7 | **Jakarta BSI** | 22 | unconfirmed | ENA **PRJNA596854** (open, confirmed) | Indonesia | binary (small) | small N; AST method unconfirmed |

(Comparator surfaced, not scored: Scotland E. coli BSI `PRJEB12513`, ~162 isolates.)

## Source-Locator Coverage

- Total candidate rows: 9 (7 cohorts + 2 comparator/context)
- Supported (audit+mapping+banned-phrase pass, medium/high conf): 8
- Unsupported: 1 (Scotland comparator — low confidence, size/AST unverified)
- Survival: 8/9 (89%)

## Caveats per row

- **ALL rows: websearch-summary provenance** (Lancet/medRxiv direct-fetch 403'd). Re-verify each accession + AST method against the direct data-availability statement before use.
- **MIC tier (Oxford/Spain/234)** = measured MIC → supports both binary R/S re-validation AND breakpoint-controlled labels. **Binary tier (Denmark/Netherlands/Norway/Jakarta)** = categorical S/I/R → supports binary decoder re-validation only.
- **Norway NORM** is a national surveillance program — but it is the *Norwegian* one, NOT the US NARMS/CDC/GenomeTrakr ecosystem the decoder excludes, so it is still provenance-disjoint from the decoder's tuning. Confirm no accession overlap via `cohort_manifest`.
- Most cohorts deposit **raw reads** (ENA study accessions) → assembling adds a step vs the decoder's assembly-accession fetch path.

## Decisions for Human Confirmation (cap 5)

| Claim | Source URL | Candidate use / Verification needed | Confidence |
|---|---|---|---|
| Spain PROBAC: 224 blood, EUCAST broth-microdilution MIC (16 drugs), ENA PRJEB62601 | https://www.thelancet.com/journals/lanmic/article/PIIS2666-5247(23)00369-5/fulltext | **Candidate use:** 2nd MIC re-validation cohort (non-UK, Spain) — strengthens "validated across countries". **Verification needed:** confirm per-isolate MIC table is downloadable; confirm PRJEB62601 holds assemblies or reads. | high |
| 234-isolate European cohort: ISO broth-microdilution MIC (11 drugs), genotype-vs-phenotype | https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9829913/ | **Candidate use:** 3rd MIC cohort; already framed as genotype↔phenotype. **Verification needed:** pin the ENA/SRA accession + country from the paper's data statement. | high |
| Denmark PRJEB37711: 699 clinical E. coli, public WGS | https://academic.oup.com/jac/article/80/4/1011/7989499 | **Candidate use:** large categorical-AST binary-validation cohort (urine-dominant — also tests UTI provenance). **Verification needed:** confirm AST fields (S/I/R per drug) are in the metadata + map to cipro/cef/gent. | high |
| Netherlands PRJEB35000: 281 blood E. coli, public WGS | https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6959556/ | **Candidate use:** categorical-AST binary-validation cohort; note ESBL-enriched (cef-skewed). **Verification needed:** AST availability per isolate; correct for ESBL-selection in cef cell. | high |
| Norway NORM nationwide BSI WGS | https://www.thelancet.com/journals/lanmic/article/PIIS2666-5247(21)00031-8/fulltext | **Candidate use:** largest non-US national cohort — high-N binary validation. **Verification needed:** find the public WGS accession + whether isolate-level AST (MIC or categorical) is released; confirm disjoint from tuning. | medium |

*(additional candidate exists: Jakarta PRJNA596854 — small, AST unconfirmed; lower priority.)*

## Portfolio synthesis (the breadth answer)

For the user's stated goal (breadth = multi-source validation portfolio), the decoder can be re-validated across **at least 6 independent cohorts spanning UK, Spain, Denmark, Netherlands, Norway, Indonesia** — none from the US surveillance ecosystem it was tuned on. Recommended execution order by VOI:
1. **Oxford (MIC, UK)** + **Spain PROBAC (MIC, ES)** + **234-European (MIC)** = the measured-MIC core (breakpoint-controlled, highest-quality labels).
2. **Denmark + Netherlands** (categorical AST, confirmed open accessions) = cheap high-N binary breadth.
3. **Norway NORM** (largest, needs accession confirmation) + **Jakarta** (small) = extend the geographic spread.

A decoder that holds across this set = a genuine "validated across N independent labs/countries" trust claim — exactly the breadth goal. Each cohort is `auto`-executor-doable (open ENA download + local AMRFinder via Docker + the shipped `call_resistance`), no money.

## Promotion Gate reminder

INPUT, not approval. Per-cohort, the first verification is the same: confirm (a) the AST/MIC table is obtainable, (b) the accession holds assemblies or reads, (c) it passes the `cohort_manifest` leakage gate. Do not treat any accession as cleared until its data-availability statement is read directly (all rows are websearch-summary provenance).
