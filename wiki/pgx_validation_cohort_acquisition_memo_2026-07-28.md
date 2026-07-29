# Deep-dive: acquiring a clinical/WGS validation cohort for the human-PGx decoder (2026-07-28)

**Question posed:** acquire a real clinical/WGS validation cohort to push human-PGx precision past the
1000G-30x-panel limits.

**Headline verdict (R2 reframe):** the "acquire a cohort" need **splits into 3 tiers, and the two cheapest
are FREE + mostly already-reachable — we are UNDER-USING free data, not blocked on acquisition.** A genuinely
*independent measured-phenotype* cohort (the HIV-PhenoSense analog) largely **does not exist for PGx in
free individual-level form**; that tier is a controlled-access / money / institutional **authority fork**,
not an executor task.

---

## The gap, precisely (two different needs)

1. **Validation COMPLETENESS/breadth** — the panel-limited sites (TPMT `*6`/`*12`/`*40` etc. not genotyped
   in the NYGC 30x *phased VCF*) + more genes/samples than the subset we currently score.
2. **Independent EVIDENCE TIER** — real patients with a **measured wet-lab metabolizer phenotype** (probe-drug
   PK / EHR outcome), to upgrade our phenotype mapping from "faithful-to-CPIC-guideline" → "independent-measured"
   (what HIV PhenoSense did for the AMR track).

These have very different answers.

---

## Tier 1 — FREE, reachable NOW, closes need #1 (RECOMMENDED first; $0, no authority)

**1a. Recover the panel-limited sites from the 1000G CRAMs we already fetch.** The missing sites are a
**documented artifact**, not a data wall: the Star Allele Search paper states *"several star-allele-defining
variants were present in the Phase 3 10x dataset but absent from the NYGC phased 30x VCF files"* — exactly our
TPMT `*6`/`*12`/`*40`. The **1000G 30x CRAMs are full WGS**, already resolvable (`scripts/resolve_1000g_cram.py`)
and already read-level-callable (`scripts/cyp2d6_pileup_gen.py`, `samtools mpileup`). Generalizing that pileup
caller to the sentinel sites recovers the validation the phased VCF filters out. **FREE, already-tooled.**

**1b. Adopt the FULL free GeT-RM truth (we use only a subset).** The **GeT-RM Consolidated PGx Table** =
**363 samples × 34 PGx genes/loci** (incl. CYP2C9/2C19/2D6/3A4/3A5, TPMT, NUDT15, DPYD…), free Excel at
`cdc.gov/lab-quality/php/get-rm/`; plus **Star Allele Search** = star-allele calls for **all 3,202 1000G
biospecimens**, free at `coriell.org/StarAllele/Search`. We currently score a small slice (the ursaPGx TSV +
3 consensus TSVs). Ingesting the full table = **far more genes + samples of validation truth, FREE.**

> Tier-1 caveat: GeT-RM/Star-Allele truth is **consensus star-allele** (tool-derived), not a measured
> phenotype — so it validates CALLING + the CPIC phenotype MAPPING stays "faithful-to-guideline". It does
> NOT provide need #2. But it is the honest, free way to make the calling validation broad + complete.

## Tier 2 — Independent MEASURED phenotype (the HIV-analog): mostly REACHABILITY-BLOCKED (honest negative)

Measured probe-drug PK datasets pairing metabolic ratio (MR) with genotype **exist**:
- **Estonian Biobank recall study** (npj Genomic Medicine 2025) — omeprazole→CYP2C19, metoprolol→CYP2D6 MRs
  by diplotype (individual data is **Estonian-Biobank-controlled**; publication is aggregate/stratified).
- **Chinese dextromethorphan study** (PMC5411458, n=235) — **verified 2026-07-28: group-means only, NO
  individual genotype↔MR download** (only a primer table supplement).
- **Escitalopram TDM** (n=5067) / antidepressant meta-analysis (13 studies) — measured serum MRs + genotype,
  but again **aggregate / TDM-service-controlled**, not individual-level open data.

**Finding:** unlike HIV PhenoSense (a free, individual-level, isolate-level wet-lab label), PGx measured-
phenotype cohorts publish **aggregate** data; the **individual-level** genotype+MR lives in **controlled
biobanks**. So the "free independent measured-phenotype win" is **not readily available** for PGx — this tier
collapses into Tier 3 (controlled access). (A targeted figshare/Zenodo hunt could still surface an exception;
low-probability, worth a single cheap probe if the user wants to try.)

## Tier 3 — Controlled / paid clinical cohorts (AUTHORITY / MONEY forks — user's decision)

Real patients, WGS + EHR-derived phenotypes — the genuine independent clinical cohort — but all gated:

| cohort | content | access | cost | fit |
|---|---|---|---|---|
| **All of Us** (NIH) | WGS + EHR + star-allele + metabolizer phenotype (NM/IM/PM + activity score) | **free platform**, but the **institution must sign a DURA** (VUMC; ~**months**) + individual credentialing (ID.me/Login.gov + training) for the **Controlled Tier** (WGS) | free access; **GCP compute** beyond the $300 initial credit = money | **best free-access clinical cohort**; measured-ish (EHR) phenotypes at scale |
| **UK Biobank** | WGS + deep phenotypes | application + approval | **application fee (money)**; noted externally-walled 2026-07-05 | strong but paid + walled |
| **eMERGE** (dbGaP) | WGS/array + EHR + PGx | dbGaP **controlled** (IRB + DUC + data-access committee) | free-ish but IRB/committee | strong; heavier application |
| **GeT-RM Coriell DNA** | the reference cell lines themselves | purchasable | **money (per-sample)** | only needed to generate NEW wet-lab data; the consensus is already free (Tier 1b) |

**Concrete path for All of Us:** the institutional DURA needs a **signing institution**. Earlier this session
you noted your **sister is a tenured professor at York University** — York could be the signing institution,
making All of Us the cleanest Tier-3 route. Still a **months-long institutional agreement + your authority
decision + GCP compute cost**, not an executor task.

---

## Recommendation (ranked)

1. **DO Tier 1 now (free, $0, mine to execute):** generalize the CRAM read-level caller to recover the
   panel-limited sentinel sites (closes the exact gap that motivated this) + ingest the full GeT-RM
   Consolidated Table + Star Allele Search for broad calling validation. This is what "acquire a validation
   cohort" actually resolves to for us today — the free data we under-use.
2. **Accept the Tier-2 honest negative:** a free individual-level measured-phenotype PGx cohort ≈ doesn't
   exist; don't spend chasing it (optional single figshare/Zenodo probe only).
3. **Tier 3 = your authority/money fork:** if you want a genuine independent CLINICAL cohort, **All of Us via
   a York-University DURA** is the best free-access option (months + credentialing + GCP compute). UK Biobank
   (fee) + eMERGE/dbGaP (IRB) are heavier. I can DRAFT the DURA-interest note / data-access rationale, but
   the institutional agreement + any spend is yours.

## Sources
- Star Allele Search (30x-VCF-drops-defining-variants caveat; 3,202 biospecimens): https://pmc.ncbi.nlm.nih.gov/articles/PMC10811916/
- GeT-RM Consolidated PGx Table (363 samples / 34 genes): https://pmc.ncbi.nlm.nih.gov/articles/PMC12103986/ · https://www.cdc.gov/lab-quality/php/get-rm/index.html
- Cyrius (1000G CRAM WGS CYP2D6 calling): https://www.nature.com/articles/s41397-020-00205-5
- All of Us access (DURA/credentialing/$300 GCP): https://support.researchallofus.org/hc/en-us/articles/35013049400468-Access-and-DURA-Support-Questions · https://www.researchallofus.org/data-tools/workbench/
- Estonian Biobank recall (measured MRs): https://www.nature.com/articles/s41525-025-00549-6
- Chinese dextromethorphan (aggregate-only, verified): https://pmc.ncbi.nlm.nih.gov/articles/PMC5411458/
