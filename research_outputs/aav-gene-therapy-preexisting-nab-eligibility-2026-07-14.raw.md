# AAV gene-therapy pre-existing neutralizing-antibody eligibility cutoffs (raw research, 2026-07-14)

> Captured 2026-07-14. Source: Claude Code (`/research` orchestrator). Topic: "gene-therapy / AAV pre-existing neutralizing-antibody eligibility cutoffs" (re-scope 2 of the immunogenicity-threshold topic). Slug: aav-gene-therapy-preexisting-nab-eligibility-2026-07-14.
> KEY FRAMING: like therapeutic-protein ADA (and UNLIKE vaccines), AAV gene-therapy eligibility has NO single universal numeric NAb threshold — cutoffs are program-, vector-, dose-, and assay-specific. Approved products split between a QUALITATIVE binary (Roctavian: Detected/Not Detected) and ASSAY-DEPENDENT titers (Zolgensma ~1:50 at one lab); the field is actively debating whether binary exclusion is too conservative (high-dose AAV5 hem-B transduces at high titers).

## Audit table (verbatim, all candidate rows)

| Claim / quantity | Numeric value | Units | Source title | Authors / org | Year | Section / table / figure | Stable URL | Access date | Quoted excerpt (≤25 words) | Extraction rationale | Source type | Confidence |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Roctavian (AAV5, hemophilia A) eligibility = qualitative, not a numeric titer | qualitative | — | ROCTAVIAN Prescribing Information / AAV5 DetectCDx (ARUP) | FDA / BioMarin / ARUP | 2023 | Indications / companion Dx | https://www.fda.gov/media/169937/download | 2026-07-14 | "anti-AAV5 antibody positive (result of Detected) are not eligible ... antibody negative (Not Detected) are eligible for treatment with ROCTAVIAN" | approved product uses binary detected/not-detected, not a titer | Regulatory guidance | medium |
| Roctavian pivotal: ineligible for pre-existing anti-AAV5 | 26/181 (14%) | fraction | Roctavian gene therapy for hemophilia A (review) | Blood Advances | 2024 | Trial screening | https://ashpublications.org/bloodadvances/article/8/19/5179/516942 | 2026-07-14 | "of 181 potential participants, 26 (14%) were ineligible because of preexisting anti-AAV5 antibodies" | real-world exclusion rate | Peer-reviewed | medium |
| Zolgensma (AAV9, SMA) positivity threshold (assay-dependent, NOT official label) | >1:50 (one lab); ≥1:25 elevated (another) | titer | Testing preexisting antibodies prior to AAV gene transfer therapy | Kruzik et al. (Mol Ther M&C) | 2022 | Assay variability | https://pmc.ncbi.nlm.nih.gov/articles/PMC8933338/ | 2026-07-14 | "samples >1:50 dilution are reported as positive, whereas samples ≥1:25 dilution are reported as elevated at the other laboratory" | commonly cited but lab-specific, not a standardized label cutoff | Peer-reviewed | medium |
| Etranacogene dezaparvovec (AAV5, hem B, HOPE-B) permissive threshold | <700 (~678) | NAb titer | HOPE-B durability analyses (via review synthesis) | CSL/uniQure | 2023 | efficacy-titer relationship | https://www.nejm.org/doi/full/10.1056/NEJMoa2211644 | 2026-07-14 | "FIX activity was similar in participants without and with pre-existing NAbs to AAV5 up to a titer of 678" | most permissive approved-program threshold (high vector dose) | Peer-reviewed | medium |
| Etranacogene lead-in: median pre-existing AAV5 NAb titer | 58 (range 9–3,440) | NAb titer | Natural history of preexisting AAV5 antibodies (etranacogene phase 3 lead-in) | (CSL/uniQure authors) | 2024 | Abstract | https://pmc.ncbi.nlm.nih.gov/articles/PMC12441693/ | 2026-07-14 | "48% (32/67) of enrolled participants had detectable NAbs (NAb+) with a median titer of 58 (range: 9–3,440)" | real-cohort titer distribution | Peer-reviewed | high |
| Etranacogene lead-in: AAV5 NAb seroprevalence | 47.8 | % | (same natural-history paper) | 2024 | Results | https://pmc.ncbi.nlm.nih.gov/articles/PMC12441693/ | 2026-07-14 | "The seroprevalence of NAb+ participants at screening was 47.8%" | population seropositive fraction | Peer-reviewed | high |
| Historic AAV2 (early hem B) — low titers block transduction; common exclusion | >1:5 (exclusion); ~1:2–1:11 neutralized | NAb titer | Testing preexisting antibodies (review) | Kruzik et al. | 2022 | Historical | https://pmc.ncbi.nlm.nih.gov/articles/PMC8933338/ | 2026-07-14 | "low titers (titer < 1:2 to 1:11 ...) of preexisting anti-AAV2 neutralizing antibody completely neutralized doses of AAV2" | very low titers historically disqualifying | Peer-reviewed | medium |
| AAV9 (Pompe, mouse-model-informed) proposed exclusion threshold | >1:100 (exclude); ~1:1000 blocks | NAb titer | Using an In Vivo Mouse Model to Determine Exclusion Criteria of Anti-AAV9 NAb | (PMC10976115 authors) | 2024 | Conclusions | https://pmc.ncbi.nlm.nih.gov/articles/PMC10976115/ | 2026-07-14 | "propose to preclude patients with NAb titers > 1:100"; "high-level NAb, a titer about 1:1000, ... completely blocked transduction" | preclinical-informed exclusion proposal | Peer-reviewed | medium |
| AAV serotype seroprevalence (context for eligibility) | AAV5 ~40% (lowest); AAV8 ~40–58%; AAV1/2 up to 70% | % | Global seroprevalence of NAbs against AAV serotypes | (Mol Ther M&C) | 2024 | Results | https://www.sciencedirect.com/science/article/pii/S2329050124000895 | 2026-07-14 | "NAb activity directed against AAV5 was the least prevalent ... prevalence at the 1:2 dilution was 60.8% for AAV1 and 40.8% for AAV5" | why AAV5 chosen; population-level context | Peer-reviewed | medium |
| Official Zolgensma label numeric NAb cutoff (1:50) | 1:50 | titer | (unconfirmed as official label — see unsupported) | — | — | — | — | 2026-07-14 | source could not confirm 1:50 as the official Zolgensma label criterion (it is one lab's positivity threshold) | routed to unsupported | Secondary | low |

## Highest-confidence rows (top 5)

1. Row 5 (etranacogene median NAb titer 58, range 9–3,440) — direct WebFetch verbatim; HIGH.
2. Row 6 (etranacogene seroprevalence 47.8%) — direct WebFetch verbatim; HIGH.
3. Row 1 (Roctavian binary Detected/Not-Detected eligibility) — FDA-approved companion Dx design; the key qualitative finding.
4. Row 4 (etranacogene permissive up to titer 678) — the most permissive approved-program threshold.
5. Row 8 (AAV9 Pompe proposed >1:100) — preclinical-informed concrete titer proposal.

## Low-confidence rows

- Row 10 (official Zolgensma label 1:50) — the searches explicitly could NOT confirm 1:50 as the official label criterion (it is one laboratory's positivity threshold). Routed to unsupported.

## Honest gaps

- **No single universal AAV NAb eligibility threshold exists** — this MIRRORS the parent-topic finding (therapeutic-protein ADA) more than the vaccine finding: cutoffs are program-, vector-, dose-, and assay-specific, and NAb assays are explicitly NOT standardized (a more-sensitive luciferase reporter reclassified GFP-assay-negative patients as positive).
- Approved-product eligibility splits into two shapes: **qualitative binary** (Roctavian AAV5: Detected/Not Detected) vs **assay-dependent titer** (Zolgensma AAV9 ~1:50 at one lab). No numeric threshold appears in Roctavian's label at all.
- The field is actively debating that binary exclusion is TOO conservative: high-dose AAV5 hem-B (etranacogene, 2×10¹³ vg/kg) transduces at titers that would exclude a lower-dose program — a nonlinear dose×titer threshold.
- NEJM HOPE-B (<700 / 678) carried at medium via review-synthesis; the direct natural-history paper confirmed the cohort distribution (median 58) but not the <700 efficacy-comparability threshold (VALUE NOT PRESENT there — it is the HOPE-B efficacy analysis).
