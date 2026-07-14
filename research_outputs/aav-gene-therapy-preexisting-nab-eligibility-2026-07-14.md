# AAV gene-therapy pre-existing neutralizing-antibody eligibility cutoffs — supported memo (V1 invocation)
<!-- memo-schema: 0.4 -->

> Captured 2026-07-14. Source: Claude Code (`/research`). Slug: aav-gene-therapy-preexisting-nab-eligibility-2026-07-14.
> Audit floor: 5 of 5 locators per row. Banned-phrase scan: 0 hard-reject. Cite-token noise: 0.

## Research Context

- **Problem:** in-vivo immunogenicity clinical threshold → re-scope 2: gene-therapy / AAV pre-existing-NAb eligibility.
- **Headline finding:** **NO single universal AAV NAb eligibility threshold** — this MIRRORS the parent therapeutic-protein-ADA finding (not the vaccine one). Cutoffs are program-, vector-, dose-, and assay-specific; NAb assays are explicitly unstandardized. Approved products split between a **qualitative binary** (Roctavian AAV5: Detected/Not Detected — no titer in the label) and an **assay-dependent titer** (Zolgensma AAV9 ~1:50 at one lab, ≥1:25 "elevated" at another). The field is actively debating that binary exclusion is too conservative (high-dose AAV5 hem-B transduces at high titers).
- **Captured:** 2026-07-14 · **Schema:** memo-schema 0.4

## Audit table (supported rows only)

| Program / context | Threshold | Vector | Shape | Source | Confidence |
|---|---|---|---|---|---|
| **Roctavian** (hemophilia A) | anti-AAV5 **Detected → ineligible** (no numeric titer) | AAV5 | qualitative binary (FDA companion Dx) | FDA label / ARUP AAV5 DetectCDx | medium |
| Roctavian pivotal exclusion rate | 26/181 (**14%**) ineligible for pre-existing anti-AAV5 | AAV5 | real-world | Blood Advances 2024 | medium |
| **Zolgensma** (SMA) | **>1:50** positive (one lab); **≥1:25** elevated (another) — assay-dependent, NOT official label | AAV9 | titer (lab-specific) | Kruzik 2022 (PMC8933338) | medium |
| **Etranacogene** (hemophilia B, HOPE-B) | comparable FIX activity **up to titer ~678 (<700)** — permissive | AAV5 | titer (high vector dose) | NEJM HOPE-B 2023 | medium |
| Etranacogene lead-in cohort | median NAb titer **58** (range **9–3,440**) | AAV5 | real distribution | PMC12441693 2024 | **high** |
| Etranacogene lead-in seroprevalence | **47.8%** NAb-positive | AAV5 | real cohort | PMC12441693 2024 | **high** |
| Historic AAV2 (early hem B) | **>1:5** exclusion; ~1:2–1:11 already neutralized | AAV2 | titer (very low) | Kruzik 2022 | medium |
| AAV9 Pompe (mouse-model-informed) | propose exclude **>1:100**; ~**1:1000** blocks transduction | AAV9 | proposed | PMC10976115 2024 | medium |
| Serotype seroprevalence (context) | AAV5 ~**40%** (lowest); AAV8 ~40–58%; AAV1/2 up to **70%** | — | population | Global seroprev. 2024 | medium |

Full URLs + verbatim quotes in `aav-gene-therapy-preexisting-nab-eligibility-2026-07-14.raw.md`.

## Source-Locator Coverage

- Rows submitted 10 · audit floor 10 · mapping floor 10 · banned-phrase 10 · final supported 9 (1 low → unsupported: "official Zolgensma label = 1:50", unconfirmed). Survival 9/10 (90%).

## Caveats per row

- **Etranacogene median 58 / seroprevalence 47.8%** — HIGH (direct-fetch verbatim). The **<700 (678)** permissive-efficacy threshold is the HOPE-B efficacy analysis (medium via review synthesis) — the natural-history paper confirmed the cohort distribution but returned VALUE NOT PRESENT for the <700 efficacy-comparability figure.
- **Zolgensma 1:50** — one laboratory's positivity threshold, NOT a standardized/label cutoff; a second lab uses ≥1:25 "elevated". Assay-dependent.
- **No universal number** — NAb assays are explicitly unstandardized (a more-sensitive luciferase reporter reclassified GFP-assay-negative patients as positive); cross-trial cutoff comparison is unreliable.

## Decisions for Human Confirmation (cap 5)

| Claim | Numeric value | Units | Source URL | Candidate use / Verification needed | Confidence |
|---|---:|---|---|---|---|
| No universal AAV NAb eligibility threshold — program/vector/dose/assay-specific | qualitative | — | https://pmc.ncbi.nlm.nih.gov/articles/PMC8933338/ | **Candidate use:** the framing answer if "immunogenicity clinical threshold" meant gene therapy — mirrors the ADA "no universal number" finding, not the vaccine finding. **Verification needed:** none (well-established). | high |
| Roctavian eligibility is BINARY (Detected → ineligible), no titer | qualitative | — | https://www.fda.gov/media/169937/download | **Candidate use:** the one approved product with a clean regulatory eligibility rule — but it is qualitative, not a threshold. **Verification needed:** confirm verbatim in the current FDA label. | medium |
| Etranacogene real cohort NAb distribution | median 58 (9–3,440) | titer | https://pmc.ncbi.nlm.nih.gov/articles/PMC12441693/ | **Candidate use:** the best-evidenced real titer distribution; direct-fetch verbatim. **Verification needed:** none. | high |
| Etranacogene permissive-efficacy threshold | ~678 (<700) | titer | https://www.nejm.org/doi/full/10.1056/NEJMoa2211644 | **Candidate use:** the most permissive approved-program threshold (high vector dose overcomes NAbs). **Verification needed:** confirm the 678 figure in the NEJM HOPE-B efficacy analysis directly. | medium |
| AAV9 Pompe proposed exclusion | >1:100 | titer | https://pmc.ncbi.nlm.nih.gov/articles/PMC10976115/ | **Candidate use:** a concrete preclinical-informed titer proposal. **Verification needed:** note it is a proposal from a mouse model, not an approved criterion. | medium |

## Promotion Gate reminder

INPUT to the 4-step Promotion Gate, NOT approval. The two direct-fetch rows (median 58, seroprevalence 47.8%) are high-confidence; every titer *cutoff* is program/assay-specific and must be read against that specific trial's protocol + validated assay — cross-program transfer of a numeric cutoff is a value-fidelity error here.
