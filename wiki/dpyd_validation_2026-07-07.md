# DPYD v0.1 validation — 1000G/gnomAD AF-corroboration (2026-07-07)

_KNOWLEDGE_BASELINE AF-corroboration (like VKORC1/SLCO1B1) — confirms the caller's 4 actionable variants are real + correctly-positioned at CPIC-expected population frequencies; NOT an independent per-sample diplotype concordance number. NOT clinical._

**Verdict: AF_CORROBORATED** — 4/4 actionable variants within the CPIC/gnomAD-expected EUR frequency band.

- **GeT-RM concordance:** EXTERNAL_WALL: GeT-RM DPYD consensus (Pratt 2016 CDC) is a paper-supplement table, NOT in the CYP-only ursaPGx benchmark → needs manual curation for a clean automated concordance; deferred
- **Deployment:** decoded on 5 real PGP-UK humans (all *1/*1 NM — consistent with these low freqs)

| allele | rsid | function | EUR AF (1000G) | expected band | verdict |
|---|---|---|---|---|---|
| *2A | rs3918290 | no_function | 0.00497 | 0.002-0.015 | IN_BAND |
| *13 | rs55886062 | no_function | 0.00099 | 0.0-0.005 | IN_BAND |
| c.2846A>T | rs67376798 | decreased | 0.00696 | 0.002-0.015 | IN_BAND |
| HapB3 | rs75017182 | decreased | 0.02386 | 0.01-0.06 | IN_BAND |

_AFs: Ensembl 1000GENOMES:phase_3:EUR (2026-07-07). Expected bands: CPIC DPYD allele-functionality (Amstutz 2018) + gnomAD. NOT a clinical tool._
