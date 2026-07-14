# Vaccine seroprotection titer thresholds — supported memo (V1 invocation)
<!-- memo-schema: 0.4 -->

> Captured 2026-07-14. Source: Claude Code (`/research`). Slug: vaccine-seroprotection-titer-thresholds-2026-07-14.
> Audit floor: 5 of 5 locators per row. Banned-phrase scan: 0 hard-reject. Cite-token noise: 0.

## Research Context

- **Problem:** in-vivo immunogenicity clinical threshold → re-scope 1: vaccine seroprotection titers.
- **Headline finding:** UNLIKE therapeutic-protein ADA (no universal clinical threshold), vaccine immunogenicity **HAS established, widely-used per-pathogen numeric "correlates of protection."** Each is pathogen-specific + assay-specific + mostly derived decades ago from natural-immunity / challenge / efficacy data, and several are **population/setting-dependent** (a correlate ≠ an absolute individual guarantee).
- **Captured:** 2026-07-14 · **Schema:** memo-schema 0.4

## Audit table (supported rows only — the correlate-of-protection reference table)

| Pathogen / antigen | Threshold | Units | Assay | Source | Confidence |
|---|---|---|---|---|---|
| **Hepatitis B** (anti-HBs) | ≥10 | mIU/mL | CLIA/ELISA | CDC MMWR RR-62-10 (2013) | medium |
| **Tetanus** (anti-tetanus) | ≥0.1 | IU/mL | ELISA (conservative) | trial/ECDC convention | medium |
| **Diphtheria** (anti-diphtheria) | ≥0.1 (ELISA); ≥0.01 (neutralization) | IU/mL | ELISA / Vero-cell neut. | PMC7563378 | medium |
| **Rabies** (RVNA) | ≥0.5 | IU/mL | RFFIT | WHO (verify direct) | low |
| **Influenza** (adults) | ≥1:40 | HAI titer | hemagglutination inhibition | Ng et al. 2015 (PMC4498268) | medium |
| **Measles** (against infection) | ≥120 | mIU/mL | plaque-reduction neutralization | Bolotin et al. J Infect Dis 2020 | **high** |
| **Hib** (anti-PRP) | ≥0.15 (short); ≥1.0 (long) | μg/mL | ELISA | Hib conjugate convention | medium |
| **Pneumococcus** (infants, WHO aggregate) | ≥0.35 | μg/mL | anti-capsular IgG ELISA | Andrews/Siber 2007 | medium |
| **Meningococcus** (human complement) | ≥1:4 | hSBA titer | serum bactericidal (human C′) | Goldschneider / Borrow 2003 | medium |
| **Meningococcus** (rabbit complement, MenC) | ≥1:8 | rSBA titer | serum bactericidal (rabbit C′) | Borrow 2003 | medium |

Full URLs + verbatim quotes in `vaccine-seroprotection-titer-thresholds-2026-07-14.raw.md`.

## Source-Locator Coverage

- Rows submitted 12 · audit floor 12 · mapping floor 12 · banned-phrase 12 · final supported 10 (2 low → unsupported: measles ≥200 severe-disease [unconfirmed by primary], polio ≥1:8 [not fetched]). Survival 10/12 (83%).

## Caveats per row

- **Measles 120 mIU/mL** — HIGH (direct-fetch verbatim). The commonly-cited **≥200 mIU/mL for severe-disease** protection was NOT confirmed by the fetched systematic review (VALUE NOT PRESENT) → routed to unsupported.
- **Influenza 1:40** — ~50% protection *in adults*; children need ~1:110 (1:40 ≈ only 22%); household exposure lowers it. A population/setting-dependent correlate, not an individual guarantee.
- **Rabies 0.5 IU/mL** — WHO-attributed but not directly fetched (WHO-confirm search hit an AUP content filter) → low; verify against the WHO rabies position paper.
- **Meningococcus / pneumococcus** — assay-source dependent (hSBA vs rSBA) and (pneumo) a population aggregate, not a per-serotype individual cutoff.
- CDC HepB MMWR direct fetch returned 403 → HepB at medium via search-summary.

## Decisions for Human Confirmation (cap 5)

| Claim | Numeric value | Units | Source URL | Candidate use / Verification needed | Confidence |
|---|---:|---|---|---|---|
| Vaccine correlates of protection EXIST per-pathogen (unlike therapeutic-protein ADA) | qualitative | — | https://academic.oup.com/jid/article/221/10/1576/5610904 | **Candidate use:** the framing answer if "immunogenicity clinical threshold" meant vaccines — established numeric correlates exist, pathogen-by-pathogen. **Verification needed:** none (well-established). | high |
| Measles PRN correlate of protection | 120 | mIU/mL | https://academic.oup.com/jid/article/221/10/1576/5610904 | **Candidate use:** the single best-evidenced correlate; direct-fetch verbatim. **Verification needed:** none. | high |
| Influenza HAI 50%-protection benchmark (adults) | 1:40 | titer | https://pmc.ncbi.nlm.nih.gov/articles/PMC4498268/ | **Candidate use:** the classic regulatory correlate. **Verification needed:** note it is adult-specific (~1:110 in children; less in household exposure). | medium |
| Hepatitis B seroprotection | 10 | mIU/mL | https://www.cdc.gov/mmwr/preview/mmwrhtml/rr6210a1.htm | **Candidate use:** the operational "vaccine responder" cutoff. **Verification needed:** direct CDC-page re-fetch (403'd this run; canonical + multiply-corroborated). | medium |
| Rabies adequate RVNA response | 0.5 | IU/mL | https://www.who.int/teams/immunization-vaccines-and-biologicals/diseases/rabies | **Candidate use:** WHO adequacy threshold. **Verification needed:** direct WHO-position-paper confirmation (not fetched this run — AUP filter). | low |

## Promotion Gate reminder

INPUT to the 4-step Promotion Gate, NOT approval. The medium/low rows (HepB, rabies, tetanus, diphtheria, Hib, pneumo, meningo) rest on search-summary or convention and need a direct authoritative re-fetch to clear step 1 at high confidence; measles is the only direct-fetch-verbatim high-confidence row.
