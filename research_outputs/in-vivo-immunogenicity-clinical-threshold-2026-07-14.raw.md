# In-vivo immunogenicity clinical thresholds for therapeutic proteins — ADA concentration/titer, assay sensitivity + cut-point FPRs (raw research, 2026-07-14)

> Captured 2026-07-14. Source: Claude Code (`/research` orchestrator). Topic: "in-vivo immunogenicity clinical threshold". Slug: in-vivo-immunogenicity-clinical-threshold-2026-07-14.
> Web search via Claude Code WebSearch + WebFetch — single-session, one-pass synthesis. Scoped to therapeutic-protein/biologic ADA immunogenicity (see scope note in the intent contract).
> KEY FRAMING (surfaces in the supported memo): there is NO single universal in-vivo immunogenicity *clinical decision* threshold — clinical significance is determined product-by-product. The closest quantitative anchors are (a) the ~100 ng/mL ADA concentration that "may be associated with clinical events" (an assay-sensitivity benchmark, not a clinical cutoff) and (b) product-specific ADA-titer/drug-trough relationships.

## Audit table (verbatim, all candidate rows)

| Claim / quantity | Numeric value | Units | Source title | Authors / org | Year | Section / table / figure | Stable URL | Access date | Quoted excerpt (≤25 words) | Extraction rationale | Source type | Confidence |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| FDA recommended minimum ADA assay sensitivity | 100 | ng/mL | Immunogenicity Testing of Therapeutic Protein Products — Developing and Validating Assays for ADA Detection | FDA / CDER | 2019 | Sensitivity | https://www.fda.gov/media/119788/download | 2026-07-14 | "recommending that screening and confirmatory IgG and IgM ADA assays achieve a sensitivity of at least 100 nanograms per milliliter" | FDA's current analytical sensitivity target | Regulatory guidance | medium |
| ADA concentration that may be associated with clinical events | 100 | ng/mL | (same FDA guidance) | FDA / CDER | 2019 | Sensitivity rationale | https://www.fda.gov/media/119788/download | 2026-07-14 | "recent data suggest that concentrations as low as 100 ng/mL may be associated with clinical events" | closest thing to an in-vivo clinical concentration anchor | Regulatory guidance | medium |
| Prior FDA sensitivity recommendation (superseded) | 250–500 | ng/mL | (same FDA guidance) | FDA / CDER | pre-2019 | Sensitivity | https://www.fda.gov/media/119788/download | 2026-07-14 | "Traditionally, FDA has recommended sensitivity of at least 250 to 500 ng/mL" | historical benchmark before the 100 ng/mL lowering | Regulatory guidance | medium |
| Screening-assay target false-positive rate | 5 | % | (same FDA guidance) | FDA / CDER | 2019 | Cut point determination | https://www.fda.gov/media/119788/download | 2026-07-14 | "a screening assay should have at least 90% chance to satisfy the 5% false positive rate" | screening-tier cut-point FPR target | Regulatory guidance | medium |
| Confirmatory-assay target false-positive rate | 1 | % | (same FDA guidance) | FDA / CDER | 2019 | Cut point determination | https://www.fda.gov/media/119788/download | 2026-07-14 | "a confirmatory assay should have at least 80% chance to satisfy the 1% false positive rate" | confirmatory-tier cut-point FPR target | Regulatory guidance | medium |
| Confirmatory cut point (GDF15 CCP-II, worked example) | 24.8 | % | Immunogenicity assessment strategy for a chemically modified therapeutic protein in clinical development | Frontiers in Immunology (Rondon et al. style team) | 2024 | Results, ADA validation results | https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2024.1438251/full | 2026-07-14 | "The confirmatory cut point for GDF15 (CCP-II) was 24.8%." | study-specific confirmatory cut point | Peer-reviewed | high |
| Titer cut point factor (TCPF) | 2.2 | dimensionless | (same Frontiers paper) | Frontiers in Immunology | 2024 | Results / Table 1 | https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2024.1438251/full | 2026-07-14 | "The TCPF was defined with a robust-parametric approach in the log-transformed data set with outliers removed as 2.2" | factor above screening CP used to report titer | Peer-reviewed | high |
| Confirmatory assay target FPR (direct, worked example) | 1 | % | (same Frontiers paper) | Frontiers in Immunology | 2024 | Methods, confirmatory cut points | https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2024.1438251/full | 2026-07-14 | "The target false positive rate was 1%" | corroborates the FDA 1% confirmatory target | Peer-reviewed | high |
| Confirmatory assay achieved FPR (worked example) | 0.7 | % | (same Frontiers paper) | Frontiers in Immunology | 2024 | Results, ADA validation results | https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2024.1438251/full | 2026-07-14 | "the FPR excluding outliers was 0.7% which is considered acceptable for the drug confirmatory assay" | achieved confirmatory FPR vs 1% target | Peer-reviewed | high |
| Screening assay achieved FPR (worked example) | 9.6 | % | (same Frontiers paper) | Frontiers in Immunology | 2024 | Results, ADA validation results | https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2024.1438251/full | 2026-07-14 | "the FPR excluding outliers was 9.6%" | achieved screening FPR (5% target, conservative accept) | Peer-reviewed | high |
| ADA assay sensitivity achieved (screening, worked example) | 65.5 | ng/mL | (same Frontiers paper) | Frontiers in Immunology | 2024 | Table 1 | https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2024.1438251/full | 2026-07-14 | "Sensitivity (ng/mL) 65.5" | achieved sensitivity better than the FDA 100 ng/mL target | Peer-reviewed | high |
| Clinically significant ADA (definition — no universal numeric cutoff) | qualitative | — | What are clinically significant anti-drug antibodies and why is it important to identify them | (PMC11682980 authors) | 2024 | "What are clinically significant ADAs" | https://pmc.ncbi.nlm.nih.gov/articles/PMC11682980/ | 2026-07-14 | "Clinically significant antibodies are those that have an impact on patient's health and are associated with exposure and ability to respond to therapeutic." | clinical significance is case-by-case, not a fixed number | Peer-reviewed | high |
| Adalimumab free-drug diagnostic cut-off (drug-specific) | 0.8 | μg/mL | Method to optimize the treatment of patients with biological drugs (patent) | (USPTO patent) | n/a | diagnostic cut-off | https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/10215762 | 2026-07-14 | "there is never a positive titer of anti-adalimumab antibodies above 0.8 μg/ml of free adalimumab" | product-specific ADA-vs-trough threshold | Patent / secondary | low |

## Highest-confidence rows (top 5)

1. Row 6 (GDF15 confirmatory cut point 24.8%) — direct WebFetch verbatim, peer-reviewed worked example.
2. Row 7 (TCPF 2.2) — direct WebFetch verbatim, Table 1.
3. Row 8 (confirmatory achieved FPR 0.7%) — direct WebFetch verbatim.
4. Row 11 (clinically-significant-ADA definition) — direct WebFetch verbatim; anchors the "no universal number" finding.
5. Row 10 (achieved sensitivity 65.5 ng/mL) — direct WebFetch verbatim (table cell — quote-shape caveat).

## Low-confidence rows

- Row 12 (adalimumab 0.8 μg/mL) — patent source + WebSearch-summary provenance; drug-specific. Routed to unsupported.
- Rows 1–5 (FDA canonical 100 ng/mL / 5% / 1% / 250–500 ng/mL) — the FDA PDF direct-fetch returned HTTP 404; values captured from WebSearch summary of the FDA guidance + corroborated by the Frontiers paper (which cites the same 1% target directly). Marked medium (provenance: search-summary) pending a direct PDF re-fetch.

## Honest gaps

- **There is NO single universal in-vivo immunogenicity *clinical decision* threshold.** The 100 ng/mL figure is an *assay sensitivity* benchmark ("concentrations as low as 100 ng/mL may be associated with clinical events"), NOT a clinical action cutoff. FDA + the clinical-significance review both state clinical relevance is determined product-by-product against PK/PD/efficacy/safety.
- **Vaccine seroprotection titers** (e.g., anti-HBs ≥10 mIU/mL, rabies ≥0.5 IU/mL, influenza HAI ≥1:40) — a different meaning of "immunogenicity clinical threshold" — NOT covered under this therapeutic-protein/ADA scoping. Re-scope available if that was the intent.
- **Gene-therapy / AAV immunogenicity** (pre-existing NAb titer eligibility cutoffs) — NOT covered; product- and vector-specific.
- **Neutralizing-antibody clinical titer thresholds** — no universal numeric; product-specific.
- FDA PDF direct fetch 404'd → the canonical FDA rows rest on a search-summary of that guidance (medium confidence); a direct re-fetch would upgrade them.
