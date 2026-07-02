# All of Us Researcher Workbench — application draft (2026-07-02)

**Draft-then-ratify:** Soraya wrote this content; **you submit it under your own identity.** Everything below
the "Account setup" section is paste-ready text for the Researcher Workbench "Research Purpose" fields. The
account-creation steps are yours to do (identity verification + training + Code of Conduct = your authority).
Registered Tier is **free** (no fee); the Controlled Tier — needed for individual-level genomics — adds a
separate Data Use Agreement but still no cost.

---

## Part A — Account setup (YOUR steps; not something I can do)
1. Go to `https://www.researchallofus.org/register/` → "Register".
2. **Identity verification** — via login.gov **or** ID.me (government-ID verification). This is required and
   is yours to complete.
3. **Affiliation** — attach an institution if you have one; if not, All of Us supports unaffiliated
   researchers (select the appropriate independent/citizen-scientist path in the form).
4. **Responsible Conduct of Research training** (~1–2 h, in-workbench).
5. **Sign the Data User Code of Conduct.**
6. Once in the Workbench, create a **Workspace** and fill the "Research Purpose" fields with Part B below.
7. For genomics (WGS): request **Controlled Tier** access (extra DUA, no fee).

---

## Part B — Research Purpose (paste-ready)

### Workspace name
`Interpretable genotype-to-phenotype decoding with leakage & ancestry-shift auditing`

### Primary purpose of your project (check all that apply)
- [x] **Research on a specific disease / phenotype**
- [x] **Methods / tool development** (interpretable prediction methods)
- [ ] Ancestry / population research (secondary — used only for de-confounding/robustness, see below)
- [ ] Commercial purpose  ·  [ ] Educational  ·  [ ] For-profit

### Plain-language summary (for the public Research Hub)
We are building an **interpretable "DNA decoder"**: given a person's genetic data, it predicts measurable
health-related traits (e.g. quantitative biomarkers and common-disease risk) **and names the specific
gene(s) responsible**, rather than acting as a black box. Our prior work established, across several
non-human organisms and human cancer cell lines, a clear principle: such a decoder reliably identifies the
true causal gene **when the genetic feature it uses matches how that gene actually works** (point mutations,
gene copy number, or expression). We will use All of Us's deep, measured phenotypes to test whether this
interpretable, mechanism-matched approach transfers to human traits — with **explicit auditing for the two
ways this kind of work most often goes wrong: hidden ancestry/population structure and data leakage.**

### Scientific approach
1. Select quantitative, well-measured phenotypes (e.g. lipid/metabolic biomarkers, blood indices) with clear
   biological candidate loci — traits where a *mechanistic* prediction is testable, not just a black-box fit.
2. Build interpretable, mechanism-matched genotype features (variant-level for point-mutation traits; copy
   number / expression proxies where the mechanism is dosage-driven), following our validated
   "feature-type-matches-mechanism-type" methodology.
3. **De-confound rigorously:** report within-ancestry-group (not just global) performance, use continuous
   relatedness/PCA-based structure controls, and permutation nulls — the same discipline we apply throughout
   our work — so any reported signal is mechanism, not population structure.
4. **Leakage discipline:** treat published pathogenicity scores / catalog labels as *features to audit*, never
   as ground truth; use temporal / cohort holdouts to keep evaluation honest.
5. Benchmark against a domain-knowledge baseline; report calibration and subgroup performance.

### Anticipated findings
- Whether an interpretable, mechanism-matched decoder attributes human quantitative traits to their canonical
  loci **after de-confounding** (the honest test), and where it does *not* — an equally informative negative.
- Ancestry-stratified performance and calibration, to characterize (not paper over) generalization limits.

### Underrepresented-populations / demographics
We will report **ancestry-stratified** results specifically to characterize and improve generalization across
groups — All of Us's diversity is central to auditing whether a genotype→phenotype method is learning biology
vs population structure. We do **not** target or single out any group for a sensitive/stigmatizing trait.

### Commercial purpose
**No.** This is methods/research; the resulting decoder is an interpretable research tool. (If commercial use
is ever considered, we will re-attest and follow All of Us policy.)

### Sensitive / stigmatizing-trait review
Our planned phenotypes are quantitative biomarkers and common physiological traits — not stigmatizing
categories. If any candidate phenotype could be sensitive, we will route it through the All of Us Resource
Access Board review before use.

### Research ethics attestation
We attest to using the data only for the stated research purpose, honoring the Data User Code of Conduct, not
attempting re-identification, and reporting results in aggregate.

---

## Part C — Honest scope + boundary (Soraya)
- **I drafted; you submit.** Creating the account, verifying identity, completing training, and signing the
  Code of Conduct are yours (authority + identity). I cannot and did not do those.
- **No money required** for Registered Tier or the Controlled-Tier genomics DUA — no fee gate here (unlike UK
  Biobank, which does charge). If any step ever shows a cost, that's a separate decision for you.
- **Honest framing of value:** All of Us is the highest-value *free-to-access* deep human phenotype substrate
  (WGS + EHR + biomarkers) and is the credible unlock for extending the interpretable decoder to human traits
  — the reproducibility-freeze forward-path #1 ("acquire a non-public / access-controlled measured-label
  source"). This draft makes the *work* ready; the *access* is your decision to action.
- **Next after access lands:** I can build the ingestion + phenotype-harmonization layer (heterogeneous
  schema-tolerant) so it's ready the moment your Workbench access is approved — no waiting.
