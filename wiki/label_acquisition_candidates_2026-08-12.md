# Label acquisition candidates — verified, gate-screened (2026-08-12)

Every URL below was **fetched or API-queried during this session**, not recalled. Sizes, licences and
accessions are what the servers actually returned. Each candidate is screened against the project's own
8 rejection gates (`wiki/negative_results_map_2026-06-13.md`).

> **Disk (Care check):** `C:` has **5.7 GB free** of 254 GB. `D:` has **4,541 GB free**.
> **Download everything below to `D:`.** The top candidate alone is 2.3 GB compressed.

---

## 1. Fitness Browser — RB-TnSeq compendium ⭐ TOP PICK

**7,552 genome-wide fitness experiments across 46 bacteria + 2 archaea.**

| | |
|---|---|
| URL | <https://figshare.com/articles/dataset/February_2024_release_of_the_Fitness_Browser_RB-TnSeq_data_for_diverse_bacteria_and_archaea/25236931> |
| DOI | `10.6084/m9.figshare.25236931` (published 2024-02-16) |
| **Licence** | **CC BY 4.0** — verified via figshare API, genuinely reusable with attribution |
| Size | `feba.db.gz` **2,302 MB** (this is the one that matters — per-GENE fitness) · 48 × `db.StrainFitness.*.gz` (1.8–104 MB each, per-STRAIN, optional) · `aaseqs.gz` 43 MB |
| Total | 3,722 MB |

**Download just `feba.db.gz`** unless you want per-strain granularity. Per-gene fitness lives in `feba.db`;
the `StrainFitness` files are a finer cut we don't need yet.

**Why this is the top pick — it attacks three live gaps at once:**

1. **It is two-sided by construction.** Continuous fitness values, not "essential / not". This directly
   fixes the one-sided benchmark that currently blocks carbon-source gap-fill validation
   (`fba_carbon_growth_validation`: 21 positives, recall 1.0, **zero negatives**).
2. **It is multi-condition at scale.** The conditional-essentiality cell shipped this week scores 4 media ×
   1,075 genes. This is thousands of conditions — the same metric with ~1000× the resolution.
3. **It is multi-organism.** 48 organisms. The FBA cell's cross-organism weakness currently rests on
   **yeast alone** (every other organism is MODEL_WALLED or LABEL_WALLED).

**The ingestion path is already partly proven** — `db.StrainFitness.Keio.gz` in this release is the same
E. coli Keio fitness source `scripts/fba_carbon_growth_validate.py` already consumes a slice of. We have
been using one file from a compendium of 48.

**Gate screen: CLEARS ALL 8.** G1 clean (barcode-sequencing wet-lab measurement, not gene-call-derived) ·
G3 clean (assay reading, not sampling context) · G6 clean (continuous, no breakpoint censoring) ·
G4/G8 n/a (not surveillance isolate-population data).

> **Known data caveat, carry it:** the authors withdrew the original **sucrose** and **D-mannitol** stocks
> after publication (bad stock solutions — BW25113 grew on the old sucrose stock but not a fresh one).
> Post-Sept-2021 data uses fresh stocks. The Feb-2024 release is after that fix, but if you ever pull the
> 2017 or 2020 snapshots, filter those two carbon sources out.

**Note:** the live site `fit.genomics.lbl.gov` returns **HTTP 403 to scripts** (both the CGI and the root,
even with a browser User-Agent). That is a bot-block, not a licence wall — a real browser gets through,
which is exactly why this is a manual download. The figshare mirror is scriptable if you prefer.

### Direct URLs + checksums (verify your download)

| file | size | md5 | direct URL |
|---|---|---|---|
| `feba.db.gz` | 2,302 MB | `87b8e150df81f85cfa10650293bb603d` | <https://ndownloader.figshare.com/files/44580595> |
| `code.tar.gz` | 4.6 MB | `74b47706142e59326ecb756b7fa74e76` | <https://ndownloader.figshare.com/files/44580445> |
| `aaseqs.gz` | 43.2 MB | `ec4e093bfade3243d97686a0800d6325` | <https://ndownloader.figshare.com/files/44580544> |

### Schema — READ 2026-08-12, ingestion is de-risked

`code.tar.gz` was fetched (md5 verified) and `feba/lib/db_setup_tables.sql` extracted. 32 tables. The four
that matter, and one fact that removes the main integration risk:

- **`GeneFitness`** = `orgId, locusId, expName, fit, t` — a **continuous fitness value plus a
  t-statistic**. Two-sided *and* significance-weighted; strictly more than the binary essentiality calls
  the current cell uses.
- **`Experiment`** carries **full condition metadata**: `media`, `mediaStrength`, `temperature`, `pH`,
  **`aerobic`**, `liquid`, `shaking`, `expGroup`, `nGenerations`, and **`condition_1..4` with
  `units_1..4` + `concentration_1..4`** — up to four simultaneous defined conditions with concentrations.
  So every experiment is fully characterised, and carbon-source experiments are selectable by `expGroup`.
- **`Gene.sysName`** is documented as *"a locus tag like SO_1446 or **b2338**"*. **Those are E. coli
  b-numbers — the exact join key `dna_decode/fba/conditional_essentiality.py` and the Orth 2011 gold
  standard already use.** The join to iML1515 is direct; no crosswalk needed.
- **`Ortholog`** — cross-organism gene mapping, which is the path to testing catalog transfer across the
  48 organisms rather than assuming it.

Also present and useful later: `SpecificPhenotype` (genes with condition-specific phenotypes — a
ready-made conditional set), `Cofit` / `ConservedCofit` (cofitness), `Organism.taxonomyId` (maps to BiGG).

---

## 2. E. coli genome-wide promoter atlas (Urtecho et al.) ⭐ for design-Q2

**117,556 unique 150 bp sequences, measured promoter activity in vivo.**

| | |
|---|---|
| GEO | **GSE144621** — verified resolving, 60 samples |
| URL | <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144621> |
| Paper | <https://elifesciences.org/reviewed-preprints/92558> · site `ecolipromoterdb.com` |
| Organism | *E. coli* K-12 MG1655 |
| Readout | RNA-seq barcode expression normalised to DNA-seq abundance |

**Processed files are directly downloadable** (no need for raw reads) — verified present:
`GSE144621_U00096.2_frag-rLP5_LB_expression.txt.gz` and `..._M9_expression.txt.gz`, plus tiling/mutagenesis
tables. **Note the LB *and* M9 pair** — that makes it condition-paired, not just a single-condition atlas.

**Why:** ~9× larger than Kosuri (12,563 constructs), same regime — **constructed variation**, which is the
one regime where this project's learned models reliably win. Directly extends Track B, whose promoter arm
is the weak half (0.417 headline / 0.352 within-library, and ranking does **not** transfer across an
unfamiliar design library).

**Gate screen: clean.** G1/G3 clean (direct measurement, constructed library — no ancestry, no sampling
confound).

---

## 3. Goodman/Church/Kosuri 2013 — N-terminal codon bias

**~14,000 synthetic variants of 137 E. coli genes; sfGFP:mCherry expression over a ~200-fold range.**

| | |
|---|---|
| SRA project | **SRP029609** (sample SRS477429; RNA `SRX346948`, DNA `SRX346944`, FlowSeq `SRX346268`) |
| Paper | Science 342:475–479, DOI `10.1126/science.1241934` · PDF at <https://goodman-lab.org/pub_pdfs/10-1126-science-1241934.pdf> |

**Why it complements Kosuri rather than duplicating it:** Kosuri varies **promoter × RBS**. This varies the
**first 33 nt of the coding sequence** across 137 genes at fixed promoter/RBS. Track B currently has no
coding-sequence term at all — this is the missing third axis of "will the host express it?".

The per-variant processed expression table is in the paper's supplementary materials (what re-analyses call
the "Kosuri-All" dataset); the raw reads are the SRA accessions above.

**Gate screen: clean** (constructed, measured).

---

## 4. Strain → product titer — WEAKEST of the four, flagged honestly

| source | what | link |
|---|---|---|
| PLOS ONE 2019 | **~1,200 experimentally realised E. coli cell factories** curated from ~100 papers, with yield/titer/rate + background modifications, referenced to MG1655 | <https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0210558> (open access, PMC6333410) |
| LASER | literature-curated strain designs, E. coli + S. cerevisiae | successor auto-extracted from >15,000 PubMed articles (bioRxiv Dec 2025) |
| Dodecanol DBTL | 60 engineered MG1655 strains, measured dodecanol + full pathway proteomics | ACS Synth Biol 2019, PMID 31072100 |

**Why I rank this last despite it matching Q3 exactly:** it is **literature-curated across ~100 papers**,
so every entry carries a different medium, scale, and fermentation setup. That is a **G2-adjacent
study==class confound** — the strongest predictor of titer may be *which lab measured it*. The project has
been burned by exactly this shape before.

The **60-strain dodecanol set is the cleaner option**: one campaign, one protocol, measured pathway protein
levels *and* titer. Small, but internally consistent — a real validation set for Track A design outputs
rather than a training set.

---

## Checked and NOT recommended (with the reason)

| candidate | verdict |
|---|---|
| **Biolog PM / EcoCyc growth data** (Sarkar 2014, PMC3957686) — 5 independent PM datasets, PM1–4, ~380 conditions | **Licence-walled.** Europe PMC returned `Article with id PMC3957686 is not open access`. The underlying data lives in BioCyc/EcoCyc, which is subscription-gated for most users. Also a methodological caveat: **PM measures respiration, not growth** — EcoCyc does not distinguish them. Would have been ideal for the two-sided carbon benchmark; the Fitness Browser supersedes it and is free. |
| **MtbTnDB** — 64 standardised M. tuberculosis TnSeq screens, conditional essentiality | **Worth a look, not yet verified.** Directly relevant to the TB cell. I did not confirm its download surface or licence this session, so I am not putting it in the recommended list on an unverified basis. |

---

## Suggested download order

1. **`feba.db.gz`** from figshare 25236931 → `D:/` (2.3 GB, CC BY 4.0) — highest value, clears every gate
2. **GSE144621** processed `*_expression.txt.gz` files → `D:/` (small) — design-Q2
3. **Goodman 2013 supplementary table** (per-variant expression) — design-Q2 coding-sequence axis
4. *(optional)* dodecanol 60-strain set — Q3 validation, not training

Nothing here costs money. Nothing here requires a DUA.
