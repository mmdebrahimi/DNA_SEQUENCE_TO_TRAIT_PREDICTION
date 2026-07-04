# UK Biobank application prep — the chosen acquisition path (2026-07-04)

**User directive 2026-07-04:** "Go straight to UK Biobank — the individual/international-friendly one."
Chosen over All of Us after the DURA check (York University not registered → needs a slow institutional
agreement; UWO registered but user is alumnus-only, no active affiliation). **UK Biobank needs NO
institution-level agreement** — the individual applies as PI. This is target #2 in
`wiki/label_acquisition_anchor_2026-07-04.md`, now promoted to the active path.

**Authority boundary (unchanged):** Soraya drafts + prepares. The USER registers, proves identity, signs the
MTA, and authorizes the access fee (money gate — inviolable). Nothing here spends money.

---

## 1. Eligibility verdict — YES, via the York professorship

UK Biobank runs an **open-access** model: any *bona fide researcher*, any country not under UK/US/EU
sanctions (Canada qualifies), academia/charity/public/commercial. Health-related, public-interest research.
**No DURA / no institution-wide agreement.** BUT registration as a bona fide researcher requires three things
your **York University** faculty role supplies:

| Requirement | Your vehicle |
|---|---|
| Personal **institute email** | your `@yorku.ca` faculty email |
| **Proof of researcher status** — CV/résumé OR a link to an institute profile page | your York faculty profile page (or CV) |
| **Authorised signatory** for the Material Transfer Agreement (MTA) | York's Research/Legal/Tech-Transfer office co-signs (per-project MTA, not a DURA) |

The field mismatch (comparative literature vs. genomics) does **not** disqualify you — UKB gates on
*bona fide researcher* + *health-related, public-interest project*, not department. Represent it truthfully:
a genuine independent research project you are PI on.

**One real dependency on York:** the MTA needs an *authorised signatory* (an institutional official). Confirm
early that York's research office will co-sign a per-project MTA for a faculty-led project — this is standard
and light, but it is the one non-you step.

## 2. The AMS process (steps in order)

1. **Register** in the UK Biobank **Access Management System (AMS)** — free. Upload CV **or** institute
   profile-page link + register with the `@yorku.ca` email. UKB aims to verify within **10 working days**.
2. **Create the application** in AMS — you become **Applicant PI** and select the **data tier** (fields you
   want; see §4).
3. **Approval → access-fee request raised** + MTA sent to you **and York's authorised signatory**.
4. **Pay the access fee** (money — your authorization) + return the **signed MTA**.
5. **Data released** on the **UKB Research Analysis Platform (UKB-RAP)** — cloud analysis environment.
   Full pipeline averages **~15 weeks** submission→release.

## 3. Fees (money gate — user authorizes)

Cost-recovery, varies by **data fields + volume + academic vs. commercial** (you = academic, the cheaper
tier). **Registration is free**; the fee is raised only *after* approval — so you can register + apply + see
the exact quote **before** any spend. Notes:
- **£1,000 platform credits/project** auto-apply — BUT the programme is stated "until **Summer 2026**"
  (i.e. likely expiring now — confirm at application).
- **Student tier £500** (data-only) — applies only if a *student* under you is the applicant; not you as
  faculty PI.
- Full WGS/exome fields cost more than array+imputed. §4's phased wishlist starts cheap (array+imputed +
  a few biomarker/record fields) so the first project is low-fee.

## 4. Drafted project — grounded in the decoder's validated human cells

The north-star fit: UK Biobank is the **free independent LAB-MEASURED label** the human decoder needs. The
current human cells are openSNP **self-report** (PILOT tier); UKB's genotype × measured-phenotype at n≈500k
upgrades them to real external validation — the human analogue of the HIV Stanford / TB AMR-Portal
independent numbers. It also gives a substrate for the one unrun *fair learned-decoder* test (deep,
lab-measured trait with no curated determinant catalog).

### Lay summary (draft — for the AMS "lay summary" field)
> We are testing whether a transparent, rules-based computer method can read a person's DNA and correctly
> predict specific, well-understood traits and disease risks — without using an opaque "black-box" model.
> Our method already works on infections (predicting antibiotic resistance from a microbe's DNA) and on a
> few human traits using self-reported data. UK Biobank's large set of genetic data paired with measured
> health information (for example, blood cholesterol) lets us check, on hundreds of thousands of people,
> how accurate and how honest-about-its-uncertainty the method is. The goal is genetic interpretation that
> a clinician or patient can actually understand and trust, and a fair benchmark of simple rules against
> more complex machine-learning approaches. This is health-related and in the public interest: interpretable,
> auditable variant interpretation for clinical genomics.

### Scientific abstract (draft — for the AMS project description)
> **Background.** We have built a deterministic, interpretable genotype→phenotype decoder that emits a
> calibrated call + explicit abstention + provenance, validated across bacterial/viral/fungal antimicrobial
> resistance (curated-determinant catalogs) and a set of human single-locus cells (ClinVar Mendelian variants;
> single-SNP traits). Human cells are currently validated only on self-report data.
> **Aim.** Independently validate the deterministic human decoder against UK Biobank genotype × lab-measured
> phenotype, and benchmark it against learned (embedding/PRS) baselines on a deep, measured trait lacking a
> curated determinant catalog.
> **Cells (phase 1, low-cost fields).**
> (a) **Familial hypercholesterolemia** — ClinVar pathogenic/likely-pathogenic *LDLR / APOB / PCSK9* variants
> (from exome/imputed genotypes) → **measured LDL-cholesterol** + lipid-lowering medication + ICD-coded
> hypercholesterolemia. A clean quantitative deterministic-rule × measured-phenotype test at scale.
> (b) **Lactase persistence** — *LCT* rs4988235 → self-reported milk/dairy intake + lactose-intolerance codes
> (replicates our openSNP lactase WIN on measured intake).
> (c) **ABO / imputation ancestry-limit** — validate the committed LD-imputation map (rs8176719 from rs657152)
> across UK Biobank's genetic-ancestry strata, quantifying the ancestry-dependent purity limit we already
> found (EUR tag ~0.97 vs AFR ~0.33).
> **Methods.** Frozen deterministic rules applied unchanged (no tuning on UKB); calibration (sens/spec,
> Wilson CIs), ancestry-stratified performance, abstention rate, and a clonality/relatedness-aware confusion.
> Learned baselines (PRS, protein-language-model scores) scored on the identical held-out split.
> **Public benefit.** Transparent, auditable variant interpretation; an honest benchmark of deterministic
> rules vs. black-box models on independent measured labels.

### Data-field wishlist (phase 1 — start cheap)
Resolve exact Field IDs in the UKB **Data Showcase** at application time; representative set:
- **Genomics:** Genotyping array + **imputed genotypes** (cheap; covers rs4988235, rs8176719, rs657152) →
  add **exome sequences** later for the ClinVar LDLR/APOB/PCSK9 rare variants (higher-cost, phase 2).
- **Biomarkers:** LDL-direct, total cholesterol (for the FH cell).
- **Health records:** ICD-10 diagnoses (hospital in-patient / first-occurrence) for hypercholesterolemia +
  lactose-intolerance.
- **Diet/lifestyle:** milk-type-used (for the lactase cell).
- **Ancestry:** genetic principal components / genetic ethnic grouping (for the ancestry-stratified analysis).
- **Medications:** lipid-lowering / statin use.

## 5. The single next action (yours)
**Register in AMS** with your `@yorku.ca` email + your York faculty profile-page link (or CV) — free, ~10
working days to verify. In parallel, **email York's Research/Tech-Transfer office** to confirm they'll act as
authorised signatory on a per-project MTA. Once registered, paste the AMS application fields back here and
I'll fit the drafted lay summary + abstract + wishlist into the exact form, and build the UKB-RAP ingestion
scaffold so the pipeline is ready the instant data lands.

**What Soraya can still do now (no authority, reversible):** build the heterogeneous-schema-tolerant UKB-RAP
ingestion + deterministic-rule scoring scaffold against the drafted field wishlist, so it's ready pre-access.
