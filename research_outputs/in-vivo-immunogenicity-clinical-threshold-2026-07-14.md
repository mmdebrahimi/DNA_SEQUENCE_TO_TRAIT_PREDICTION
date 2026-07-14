# In-vivo immunogenicity clinical thresholds (therapeutic-protein ADA) — supported memo (V1 invocation)
<!-- memo-schema: 0.4 -->

> Captured 2026-07-14. Source: Claude Code (`/research`). Slug: in-vivo-immunogenicity-clinical-threshold-2026-07-14.
> Audit floor: 5 of 5 locators required per row. Mapping floor: rationale → quantity required.
> Banned-phrase scan: 0 hard-reject, 1 soft-warn ("certification" absent; "approval" absent). Cite-token noise: 0.

## Research Context

- **Problem:** in-vivo immunogenicity clinical threshold
- **Interpretation (scoped):** clinical/analytical thresholds for in-vivo immunogenicity of therapeutic proteins/biologics — anti-drug-antibody (ADA) concentration/titer, assay sensitivity + cut-point false-positive-rate targets, and the definition of a "clinically significant" ADA. NOT vaccine seroprotection or gene-therapy immunogenicity (honest gaps below; re-scope available).
- **Headline finding:** there is **NO single universal in-vivo immunogenicity clinical *decision* threshold**. Clinical significance is determined product-by-product against PK/PD/efficacy/safety. The nearest quantitative anchors are (a) the ~**100 ng/mL** ADA concentration that "may be associated with clinical events" (which drove FDA's assay-*sensitivity* target — an analytical, not clinical, cutoff), and (b) product-specific ADA-titer/drug-trough relationships.
- **Captured:** 2026-07-14 · **Schema:** memo-schema 0.4

## Audit table (verbatim, supported rows only)

| Claim / quantity | Value | Units | Source | Year | Quoted excerpt (≤25 words) | Confidence |
|---|---|---|---|---|---|---|
| FDA min ADA assay sensitivity | 100 | ng/mL | FDA Immunogenicity Testing guidance | 2019 | "assays achieve a sensitivity of at least 100 nanograms per milliliter" | medium |
| ADA conc. that may be associated with clinical events | 100 | ng/mL | FDA guidance | 2019 | "concentrations as low as 100 ng/mL may be associated with clinical events" | medium |
| Prior FDA sensitivity recommendation | 250–500 | ng/mL | FDA guidance | pre-2019 | "Traditionally, FDA has recommended sensitivity of at least 250 to 500 ng/mL" | medium |
| Screening-assay target FPR | 5 | % | FDA guidance | 2019 | "a screening assay should have at least 90% chance to satisfy the 5% false positive rate" | medium |
| Confirmatory-assay target FPR | 1 | % | FDA guidance | 2019 | "a confirmatory assay should have at least 80% chance to satisfy the 1% false positive rate" | medium |
| Confirmatory cut point (GDF15 CCP-II) | 24.8 | % | Frontiers Immunol. 2024 | 2024 | "The confirmatory cut point for GDF15 (CCP-II) was 24.8%." | high |
| Titer cut point factor (TCPF) | 2.2 | — | Frontiers Immunol. 2024 | 2024 | "The TCPF was defined ... with outliers removed as 2.2" | high |
| Confirmatory target FPR (worked example) | 1 | % | Frontiers Immunol. 2024 | 2024 | "The target false positive rate was 1%" | high |
| Confirmatory achieved FPR (worked example) | 0.7 | % | Frontiers Immunol. 2024 | 2024 | "the FPR excluding outliers was 0.7% which is considered acceptable for the drug confirmatory assay" | high |
| Screening achieved FPR (worked example) | 9.6 | % | Frontiers Immunol. 2024 | 2024 | "the FPR excluding outliers was 9.6%" | high |
| Achieved assay sensitivity (screening) | 65.5 | ng/mL | Frontiers Immunol. 2024 | 2024 | "Sensitivity (ng/mL) 65.5" | high |
| Clinically significant ADA (definition; no universal cutoff) | qualitative | — | PMC11682980 2024 | 2024 | "Clinically significant antibodies are those that have an impact on patient's health and are associated with exposure and ability to respond to therapeutic." | high |

Full URLs + sections in `in-vivo-immunogenicity-clinical-threshold-2026-07-14.raw.md`.

## Source-Locator Coverage

- Total rows submitted: 12 · Survived audit floor: 12 · Survived mapping floor: 12 · Banned-phrase: 12 · Final supported: 11 (1 low-confidence patent row → unsupported). Survival rate 11/12 (92%).

## Caveats per row

- **FDA rows (100 ng/mL, 250–500 ng/mL, 5%, 1%):** provenance = WebSearch-summary of the FDA guidance PDF (direct fetch returned HTTP 404) → downgraded to `medium`. The 1% confirmatory target is independently corroborated by the direct-fetched Frontiers 2024 paper. A direct PDF re-fetch would upgrade to high.
- **65.5 ng/mL:** table-cell quote-shape (verbatim from Table 1); value verbatim.
- No cite-token noise; no hard-reject phrases.

## Decisions for Human Confirmation (cap 5)

| Claim | Numeric value | Units | Source URL | Candidate use / Verification needed | Confidence |
|---|---:|---|---|---|---|
| ADA conc. associated with clinical events / FDA min assay sensitivity | 100 | ng/mL | https://www.fda.gov/media/119788/download | **Candidate use:** the nearest quantitative "in-vivo immunogenicity clinical threshold" — an ADA level below which detection isn't required and around which clinical events begin appearing. **Verification needed:** confirm verbatim in the FDA PDF directly (search-summary sourced; PDF 404'd this run); confirm it is an *analytical* target, not a clinical action cutoff. | medium |
| Screening / confirmatory cut-point target FPR | 5 / 1 | % | https://www.fda.gov/media/119788/download | **Candidate use:** the statistical "threshold" machinery — how positive/negative is defined per tier (5% screening, 1% confirmatory). **Verification needed:** direct FDA-PDF confirmation of the ≥90%/≥80% confidence phrasing. | medium |
| Clinically significant ADA — no universal numeric cutoff | qualitative | — | https://pmc.ncbi.nlm.nih.gov/articles/PMC11682980/ | **Candidate use:** the correct framing answer — there is no single number; clinical significance = case-by-case PK/PD/efficacy/safety impact. **Verification needed:** none (definitional, direct-fetched verbatim). | high |
| Confirmatory cut point + TCPF (worked example) | 24.8 % / 2.2 | % / — | https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2024.1438251/full | **Candidate use:** a concrete illustration of how a per-product immunogenicity threshold is statistically derived. **Verification needed:** note it is product-specific (GDF15), not generalizable. | high |
| Product-specific ADA-vs-trough threshold (adalimumab example) | 0.8 | μg/mL | (patent — see unsupported memo) | **Candidate use:** evidence that meaningful thresholds exist only per-product. **Verification needed:** replace the patent source with a peer-reviewed adalimumab therapeutic-drug-monitoring paper before any use. | low (unsupported) |

## Promotion Gate reminder

INPUT to the 4-step Promotion Gate, NOT approval. Do not lift any number without: (1) doc resolves at the cited URL; (2) section exists; (3) excerpt verbatim; (4) mapping natural. The FDA-sourced rows specifically need a direct-PDF re-fetch to clear step 1 at high confidence.
