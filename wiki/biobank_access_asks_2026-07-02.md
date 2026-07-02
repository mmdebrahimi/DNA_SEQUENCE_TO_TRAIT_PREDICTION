# Biobank / controlled-access acquisition asks (Track C, 2026-07-02)

Drafted for the USER to action. These are the catalog's highest supervised-training-value sources but they
FAIL gate G1 (controlled access) — they are the reproducibility-freeze **forward-path #1** ("acquire a
non-public label source"), which is a USER-authority decision (identity + institutional affiliation + any
fee), NOT an executor task. Soraya can draft/prepare; **the user submits + authorizes.** Ordered by
friction (lowest first).

## 1. FinnGen summary statistics — FREE, NO application (do first, zero friction)
- **What:** genome-wide association summary statistics across ~2,400+ disease endpoints (founder population).
- **Access:** public download NOW — no DAC, no fee. `https://www.finngen.fi/en/access_results`.
- **Executor-usable:** YES (summary stats are open). Individual-level data is controlled (skip unless needed).
- **Use:** locus-to-trait priors, endpoint definitions, ancestry-transfer robustness checks. A free immediate
  add; does NOT need the user's authority for the summary tier.

## 2. All of Us Researcher Workbench — Registered Tier (lowest-friction controlled access)
- **What:** >535k short-read WGS + EHR + surveys + physical measures + Fitbit + RNA-seq/proteomics.
- **Ask (USER actions):** (a) create a Researcher Workbench account with a verified institutional/eRA/login.gov
  identity; (b) complete the Responsible Conduct of Research training; (c) sign the Data User Code of Conduct.
  Registered Tier = no genomic-individual cost; Controlled Tier (genomics) needs an added agreement.
- **Friction:** MEDIUM — requires a US institutional affiliation OR the "citizen scientist" path; training ~2h.
- **Draft research statement (user pastes):** *"Genotype→phenotype prediction of common-disease risk +
  biomarker levels from WGS, with explicit ancestry-shift + leakage auditing; deterministic-rule baselines
  before any learned model."* (Aligns with the project's honesty discipline.)

## 3. UK Biobank — Approved Researcher Application (highest value, has a fee)
- **What:** ~500k participants, array+exome+genome, imaging, biomarkers, linked health records.
- **Ask (USER actions):** (a) register on the Access Management System (AMS); (b) an institution must be
  registered; (c) submit an application with a research abstract; (d) **pay the access fee** (tiered by
  data scope — MONEY GATE: the user authorizes the spend).
- **Friction:** HIGH — institutional + fee + ~weeks approval. Highest common-disease value.
- **MONEY:** the access fee is a real purchase → **user-authorized spend** (Soraya does not commit it).

## 4. dbGaP / EGA — per-study controlled access ("dataset multipliers", do AFTER a first model)
- **What:** thousands of disease-specific study cohorts (dbGaP ~3.3k studies / EGA broad EU catalogue).
- **Ask (USER actions):** (a) an eRA Commons (dbGaP) / DAC (EGA) identity + a signed institutional signing
  official; (b) a Data Access Request per study with a research use statement. No blanket fee, but per-study
  approval.
- **Friction:** HIGH per study; treat as a multiplier AFTER All of Us / UKB proves a pipeline. The ingestion
  layer should expect heterogeneous submitter schemas + consent-group partitioning.

## Recommended sequence for the user
1. **FinnGen summary stats** — grab now, free, no authority needed (Soraya can fetch on request).
2. **All of Us Registered Tier** — the lowest-friction controlled unlock; user completes account + training.
3. **UK Biobank** — only if the user authorizes the access fee (money gate).
4. **dbGaP/EGA** — later, per-study, after a first model exists.

## What Soraya can do without user authority
- Fetch + integrate FinnGen summary stats (free/open) into a locus-priors layer.
- Draft the All of Us / UK Biobank research abstracts + application field content (above) for the user to paste.
- Build the ingestion layer (heterogeneous-schema-tolerant) so it's ready when access lands.

## What requires the user (authority boundary)
- Creating the accounts / signing the Data User Code of Conduct / institutional signing official.
- Authorizing the UK Biobank access **fee** (money gate — inviolable).
- Submitting the applications under the user's identity.
