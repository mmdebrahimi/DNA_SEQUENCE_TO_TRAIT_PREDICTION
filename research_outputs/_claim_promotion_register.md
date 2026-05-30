# Claim-Promotion Register — DNA decoder research outputs
<!-- register-schema: 0.1 -->

> Generated 2026-05-23. Source memos scanned: 3 (AMR small-cohorts 2026-05-13; PRS cross-ancestry 2026-05-22; AoU PRS ancestry-stratified 2026-05-23).
> Filters the 15 follow-up-queue candidates into 5 promotion tiers. The classification governs which claims are OK to use in which contexts.
> NOT a Promotion Gate output. Each row still requires the per-memo 4-step Promotion Gate before lift into rules / wiki / external copy.

## Tiers

| Tier | Definition | Allowed surfaces |
|---|---|---|
| **T1 — Internal direction-setting only** | Direction-setting research signal; treat as hypothesis. Preprint-only, search-snippet-provenance, author-identity-uncertain, or single-source. | Claude conversation context; internal Slack / brainstorm; not committed to any PM doc, sales deck, or code. |
| **T2 — PM-doc usable** | Direction-setting + traceable to a verifiable source; safe inside internal PM verdicts / plans with explicit "preprint-tier evidence" or "single-source" caveats. | PM verdicts; planning docs; technical plans; CLAUDE.md context notes. NOT external. |
| **T3 — Sales-deck usable** | Peer-reviewed or institutional-summary source, verbatim-quote verified, no provenance flags. | Internal sales decks; pilot proposals; investor decks (with citation). NOT public marketing copy. |
| **T4 — Public-claim usable** | T3 standard PLUS: primary source directly verified, cross-ancestry implications stated honestly, no overclaim risk. | Public marketing; website; press releases; press interviews; FDA communications. |
| **T5 — Rules / code usable** | T4 standard PLUS: numeric values empirically reproduced on dna_decode's own cohort, OR cited as field-wide-known threshold with peer-reviewed primary. | `rules_v0_1.yaml`; `wiki/SME_Calibration_Worksheet.md`; production code defaults. |

## Classification

| # | Claim | Source memo | Locator | Tier | Justification |
|---|---|---|---|---|---|
| 1 | AF PRS European OR/SD 1.89 → African 1.39 (AoU validation) | aou-prs-2026-05-23 | https://pmc.ncbi.nlm.nih.gov/articles/PMC12622184/ | **T2** | Direct WebFetch verified; preprint (Research Square), not yet peer-reviewed. PM-doc usable with "preprint-tier" caveat. Promote to T3 when peer-review status flips. |
| 2 | Prostate cancer PRS effect-size 1.61 (Middle Eastern) → 2.19 (American) | aou-prs-2026-05-23 | https://www.medrxiv.org/content/10.1101/2025.10.01.25337107v1 | **T1** | Provenance: websearch-summary (downgraded). medRxiv URL returned HTTP 403 on verification attempt 2026-05-23. NOT PM-doc-usable until direct source-text confirms. |
| 3 | AoU diversity ratios: ~3x non-EUR + ~8x multi-population vs prior datasets | aou-prs-2026-05-23 | https://www.genome.gov/news/news-release/researchers-optimize-genetic-tests-for-diverse-populations-to-tackle-health-disparities | **T3** | Verbatim sentences VERIFIED against NHGRI institutional summary (2026-05-23). NHGRI cites Lennon NJ et al., Nature Medicine 2024 (https://www.nature.com/articles/s41591-024-02796-z) — NOT Comm Med as the raw memo guessed. Sales-deck usable with NHGRI citation. Promote to T4 when primary Nature Med paper directly verified (publisher-blocked at fetch 2026-05-23). |
| 4 | AoU WGS release v7 = 245,388 individuals | aou-prs-2026-05-23 | https://www.biorxiv.org/content/10.1101/2024.08.06.606846v2.full | **T2** | bioRxiv preprint; WebFetch returned HTTP 403; quote captured via WebSearch result rendering. Will move to T3 when AoU PRS Workbench docs or PMC alternate confirms cohort scale. AoU releases on 6-12 month cadence — verify current release at time of claim. |
| 5 | T2D multi-ancestry GWAS: 360K cases + 1.8M controls (41% non-EUR) | aou-prs-2026-05-23 | https://www.medrxiv.org/content/10.1101/2025.07.21.25331778.full.pdf | **T1** | medRxiv PDF returned title/authors only; per-ancestry validation numerics not extracted. Not actionable until per-ancestry R² obtained. |
| 6 | European GD-decile PGS accuracy decrease 14% (Ding 2023) | prs-cross-ancestry-2026-05-22 | https://pmc.ncbi.nlm.nih.gov/articles/PMC10284707/ | **T3** | Ding et al. Nature 2023, peer-reviewed, PMC-resolvable. Verbatim verified via PMC direct fetch in source run. Sales-deck usable with Nature citation. Promote to T4 if/when ancestry-claim-governance section becomes external. |
| 7 | Cross-ancestry overlap: closest HL decile ≈ furthest EA decile (both 0.71 average r̂i²) | prs-cross-ancestry-2026-05-22 | https://pmc.ncbi.nlm.nih.gov/articles/PMC10284707/ | **T3** | Same Ding 2023 paper; same provenance discipline as row 6. |
| 8 | Average GD→accuracy R = −0.95 across 84 traits | prs-cross-ancestry-2026-05-22 | https://pmc.ncbi.nlm.nih.gov/articles/PMC10284707/ | **T3** | Same paper; abstract-level finding. Strong trait-wide consistency signal. |
| 9 | Concordant-SNP R² improvement, African ancestry +0.0039 (p=6.63e-18) | prs-cross-ancestry-2026-05-22 | https://pmc.ncbi.nlm.nih.gov/articles/PMC12820924/ | **T2** | Momin 2026 Genetic Epidemiology; peer-reviewed but specific numeric value would benefit from independent replication before sales-deck use. |
| 10 | Method recommendation: GBLUP + PRS-CSx best for polygenic traits | prs-cross-ancestry-2026-05-22 | https://pmc.ncbi.nlm.nih.gov/articles/PMC12820924/ | **T2** | Qualitative claim, peer-reviewed; "best" depends on metric. PM-doc usable as direction-setting for method selection. |
| 11 | Random Forest standard classifier on frozen DNA-FM embeddings | amr-small-cohorts-2026-05-13 | https://pmc.ncbi.nlm.nih.gov/articles/PMC12663285/ | **T1** | Author-identity-uncertain flag in source memo. Not actionable until author + title confirmed. |
| 12 | XGBoost rank ~8.30 vs TabPFN ~4.88 on small datasets (≤1250 samples) | amr-small-cohorts-2026-05-13 | https://arxiv.org/abs/2305.02997 | **T2** | McElfresh NeurIPS 2024, peer-reviewed; specific to dna_decode's classifier-head decision. Move to T5 if/when reproduced on a dna_decode cohort. |
| 13 | GBDT performance positively correlated with samples-to-features ratio (NN-class opposite) | amr-small-cohorts-2026-05-13 | https://arxiv.org/abs/2305.02997 | **T2** | Same McElfresh paper. Direction-setting for dna_decode architecture flag in CLAUDE.md. |
| 14 | E. coli cipro + AMRFinderPlus + k-mer + XGBoost = >90% accuracy on 256-genome cohort | amr-small-cohorts-2026-05-13 | https://pubmed.ncbi.nlm.nih.gov/39320197/ | **T2** | Talamantes-Becerra 2024, peer-reviewed; sets the classical-baseline floor. Move to T5 once dna_decode reproduces or beats this on a comparable cohort. |
| 15 | For ciprofloxacin in E. coli, SNP tables may outperform gene presence-absence | amr-small-cohorts-2026-05-13 | https://pmc.ncbi.nlm.nih.gov/articles/PMC11684616/ | **T2** | Orcales PLOS Comput Biol tutorial, peer-reviewed; direction-setting for feature engineering in dna_decode. |

## Tier rollup

| Tier | Count | Notes |
|---|---|---|
| T1 (internal-only) | 3 | Rows 2, 5, 11 — all blocked on a single verification step each |
| T2 (PM-doc usable) | 8 | Rows 1, 4, 9, 10, 12, 13, 14, 15 |
| T3 (sales-deck usable) | 4 | Rows 3, 6, 7, 8 — all peer-reviewed Nature/Nature Med + verbatim verified |
| T4 (public-claim usable) | 0 | Promotion gate: primary-source direct verification + cross-ancestry honesty + no overclaim |
| T5 (rules / code usable) | 0 | Promotion gate: empirical reproduction on dna_decode cohort |

## Promotion paths (what would move a row up a tier)

- **Row 1 T2 → T3:** AF PRS paper peer-review publication (currently Research Square preprint). Watch for journal acceptance.
- **Row 2 T1 → T2:** direct medRxiv source-text verification of "1.61 → 2.19 effect sizes" sentence. medRxiv access via alternative client (curl + spoofed UA, or Sci-Hub equivalent if academic-fair-use permits).
- **Row 3 T3 → T4:** direct verification of Lennon et al. Nature Medicine 2024 primary paper (publisher-blocked at 2026-05-23 fetch; PMC alternate may exist or wait for institutional access).
- **Row 4 T2 → T3:** AoU PRS Workbench docs OR PMC alternate ID for the bioRxiv preprint.
- **Row 6/7/8 T3 → T4:** ancestry-claim-governance section in DNA decoder verdict becomes external copy (sales deck or website).
- **Row 14 T2 → T5:** dna_decode reproduces the >90% accuracy baseline on a comparable 150-250 strain cohort.

## What this register does NOT do

- Does NOT promote anything autonomously. Each row's tier is a CLASSIFICATION, not a Promotion Gate decision.
- Does NOT replace the per-memo Promotion Gate (4-step doc-resolves / section-exists / quote-verbatim / mapping-natural).
- Does NOT cover Bombardier-internal (athena-*) memos — those use a different evidence regime (SME-only, no public verification).

## Update discipline

- Append a new row when `/research-followup` adds a candidate to `_followup_queue.md`.
- Bump a row's tier when its promotion-path condition is met (record the verification evidence in the row's Justification column).
- Demote a row if its source is retracted, paywalled, or refuted (record the demotion event).
- Stale-flag any T1 row that has not moved in 60 days (likely a dead lead — drop or document why retained).
