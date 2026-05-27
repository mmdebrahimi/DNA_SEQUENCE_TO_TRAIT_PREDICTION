# Horesh 2021 File F1 Supplementary Audit — per-record label-provenance + H1 + H2 + H3 verdict (2026-05-27)

> Audit memo. NOT a V1 13-column supported memo (no `<!-- memo-schema: ... -->` marker — this is a project-internal analysis sidecar). Sidecar to `research_outputs/ecoli-pathotype-labeled-genome-substrate-survey-2026-05-27.md` Row 17 (which left the Horesh 2021 per-record label-provenance distribution as an honest gap).
>
> **Source file:** `F1_genome_metadata.csv` (4.8 MB, 10,146 rows) from Figshare DOI 10.6084/m9.figshare.13270073 (file ID 25552514, MD5 1e76642e74747a044e9e64bec5becb2d).
> **Paper:** Horesh G, Blackwell GA, Tonkin-Hill G, Corander J, Heinz E, Thomson NR (2021). A comprehensive and high-quality collection of Escherichia coli genomes and their genes. Microb Genom 7(2):000499. PMC8208696.

## Schema discovered

The File F1 CSV has 22 columns: `ID, Assembly_name, PopPUNK, Source, Length_Mbp, Num_CDSs, Num_contigs, ST, MDR, Ab_classes, Pathotype, Phylogroup, Isolation, Year, Country, Continent, Total_AMR_genes, AMR_genes, Total_virulence_genes, Pathotype_Vir_genes, Other_Vir_genes, name_in_presence_absence`.

**Key finding — label-provenance is BUILT INTO the Pathotype column.** Horesh 2021 distinguishes:
- **Independent labels** (NO `(predicted)` suffix): label extracted from source publication
- **Predicted labels** (with `(predicted)` suffix): label derived via ariba + VirulenceFinder DB gene-rule call
- **Not determined**: no label assigned

This is a built-in label-provenance distinction. The substrate-survey memo Row 17 had flagged this as an honest gap — the field DOES exist; it's encoded in the Pathotype-cell suffix rather than a separate column.

## Label-provenance distribution

| Class | Count | % of total | H1 interpretation |
|---|---:|---:|---|
| **Independent** (publication-extracted) | 2,077 | 20.47% | H1-passing per record |
| **Predicted** (gene-rule-derived) | 5,295 | 52.19% | H1-failing per record |
| **Not determined** | 2,774 | 27.34% | No label |

**H1 verdict on full Horesh collection (N=10,146):** **FAIL.** Only 20.47% have independent labels — far below the H1 ≥70% floor.

**H1 verdict on the independent-label subset (N=2,077):** **PASS by construction.** This subset IS the H1-passing pool. Publication attribution: 97.6% (2,028/2,077) cite a specific publication (Kallonen 2017, Salipante 2014, von Mentzer 2014, Ingle 2016, Hazen 2013/2016, Subashchandrabos 2013, etc.).

## Pathotype counts per provenance class

### Full collection (independent + predicted + unknown)

| Canonical pathotype | Total N | Independent | Predicted |
|---|---:|---:|---:|
| EHEC | 3,849 | 6 | 3,843 |
| ExPEC | 2,134 | 1,574 | 560 |
| STEC | 509 | 31 | 478 |
| EPEC | 269 | 269 | 0 |
| aEPEC/EPEC | 226 | 0 | 226 |
| ETEC | 183 | 183 | 0 |
| EAEC | 113 | 2 | 111 |
| EIEC | 57 | 0 | 57 |
| EAEC+STEC hybrid | 19 | 0 | 19 |
| NON-O157EHEC/EPEC | 10 | 10 | 0 |
| COMMENSAL | 2 | 2 | 0 |
| EPEC/ETEC hybrid | 1 | 0 | 1 |

### Independent-label-only subset (the H1-passing pool, N=2,077)

| Pathotype | Independent N | H2 floor (N≥50) | Verdict |
|---|---:|---:|---|
| ExPEC | 1,574 | 50 | **PASS (31× over floor)** |
| EPEC | 269 | 50 | **PASS (5.4×)** |
| ETEC | 183 | 50 | **PASS (3.7×)** |
| STEC | 31 | 50 | FAIL (62%) |
| NON-O157EHEC/EPEC | 10 | 50 | FAIL (20%) |
| EHEC | 6 | 50 | FAIL (12%) |
| EAEC | 2 | 50 | FAIL (4%) |
| COMMENSAL | 2 | 75 | FAIL (3%) |
| EIEC | 0 | 50 | FAIL (0%) — out-of-scope v0 anyway |
| aEPEC/EPEC | 0 | 50 | FAIL — independent label was always lumped under "EPEC" |

**Net H2 verdict from Horesh-only independent subset:** PASSES floor for **ExPEC + EPEC + ETEC** (3 of 6 v0 target classes). FAILS for **STEC, EHEC, EAEC, COMMENSAL** (require Tier-1 supplement from DECA + von Mentzer ETEC reference + curated commensal panel).

## H1-passing subset richness audit

For the 2,077 independent-label records:
- **Year populated:** 2,077 (100%)
- **Country populated:** 2,077 (100%)
- **Isolation field non-Unknown:** 2,074 (99.9%) — Blood / Feces / Urine dominate
- **Publication attribution:** 2,028 (97.6%) — Kallonen 2017 (n=1,319), Salipante 2014 (n=241), von Mentzer 2014 (n=181), Ingle 2016 (n=141), Hazen 2013 (n=70), Hazen 2016 (n=69), PHE (n=31), RefSeq (n=17), etc.
- **Phylogroup populated:** 2,077 (100%; B2 dominant n=1,137, B1 n=289, D n=202, Not Determined n=171)

## H3 fold-construction risk (per-pathotype × top-3 STs on H1-passing subset)

| Pathotype | N | Top-3 STs | Top-3 share | H3 verdict |
|---|---:|---|---:|---|
| ExPEC | 1,574 | ST73 (270) + ST131 (240) + ST95 (170) | 43.2% | **PASS** — ≥5 clade folds feasible without ≤2 dominant STs eating any single fold |
| EPEC | 269 | ST328 (21) + ST517 (21) + ST3 (19) | 22.7% | **PASS** — excellent ST diversity |
| ETEC | 183 | ST443 (17) + ST2368 (11) + ST398 (9) | 20.3% | **PASS** — excellent ST diversity |
| STEC | 30 | ST21 (27) + others (3) | 100% | **FAIL** — 90% of STEC in a single ST; clade-balanced folds infeasible from Horesh alone |
| EHEC | 6 | ST11 dominant | 100% | N too small for H3 anyway |
| EAEC | 2 | ST5655 + ST678 | 100% | N too small for H3 anyway |

**Net H3 verdict:** ExPEC/EPEC/ETEC pass clade-balanced fold-construction feasibility from Horesh alone. STEC/EHEC/EAEC must source from Tier-1 substrate (DECA + von Mentzer ETEC + STEC reference) to achieve fold diversity.

## Source-publication distribution (independent labels only)

| Source | N | Likely pathotype enrichment |
|---|---:|---|
| Kallonen et al. 2017 | 1,319 | ExPEC dominant (UK clinical bacteremia/UTI collection) |
| Salipante et al. 2014 | 241 | ExPEC (clinical isolate collection) |
| von Mentzer et al. 2014 | 181 | ETEC (this is the major ETEC source) |
| Ingle et al. 2016 | 141 | EPEC (Australian infant gut isolates) |
| Hazen et al. 2013/2016 | 139 | Mixed DEC pathotypes (related to DECA collection authors) |
| Public Health England | 31 | Mixed clinical |
| RefSeq | 17 | Type strains |
| Subashchandrabos et al. 2013 | 4 | ExPEC |
| Others | ~5 | (negligible) |

The publication attribution is rich and traceable. Each independent-label record can be cross-referenced to the originating publication's pathotype-assignment methodology (clinical isolation, lab assay, serotype-confirmed, etc.).

## Consolidated verdict

### Re-tiered substrate strategy (post-Horesh-F1 audit)

| Tier | Source | N (H1-passing) | Coverage | Verdict |
|---|---|---:|---|---|
| **1a** | Horesh 2021 independent-label subset | 2,077 | strong on ExPEC + EPEC + ETEC; thin on STEC + EHEC + EAEC + COMMENSAL | **PROMOTED to Tier-1a** — was Tier-2 in substrate-survey; per-record audit upgrades it |
| **1b** | Whittam DECA + von Mentzer 2021 ETEC + N=5 prototype strains | ~100-200 (pending Whittam contact) | strong on STEC + EHEC + EAEC + EIEC (DECA's traditional strengths) | Tier-1b — complements 1a's gaps |
| **2** | Horesh 2021 predicted subset (5,295) | 0 (gene-rule-derived) | large but H1-failing per record | Tier-2 — usable ONLY for assembly QC + clonal-lineage context, NOT pathotype-label evidence |
| **3** | NCBI Pathogen Detection host_disease facet | unknown (pending facet query) | likely strong on ExPEC/UPEC + commensal expansion | Tier-3 — pending Action 4 |
| **4** | EnteroBase BlastFrost-labeled records | (decline) | gene-rule-derived (same markers as v0) | Tier-4 — clonal-lineage + assembly QC only |

### H1 + H2 + H3 verdicts (consolidated)

- **H1:** PASS on the Horesh 1a subset (2,077 records, publication-extracted labels with isolation context). FAIL on the full Horesh collection.
- **H2:** PASS for ExPEC + EPEC + ETEC from Horesh 1a alone. FAILS for STEC + EHEC + EAEC + COMMENSAL — requires Tier-1b supplement (DECA + von Mentzer ETEC + curated commensal panel).
- **H3:** PASS for ExPEC + EPEC + ETEC clade-balanced folds from Horesh 1a (top-3 ST share ≤43%). FAILS for STEC (single-ST 90%) — requires Tier-1b STEC diversification.

### Concrete cohort assembly path (updated)

For a v0 cohort meeting H1 + H2 + H3 floors:

1. **Use the 2,077-record Horesh 1a independent-label subset as the substrate backbone** (publication-extracted labels; rich provenance; 100% Year + Country + Phylogroup populated; 99.9% non-Unknown Isolation).
2. **Sub-select for class balance:** keep all STEC + EPEC + ETEC + EAEC + COMMENSAL + EHEC + NON-O157EHEC/EPEC (~543 strains) + downsample ExPEC to ~150 stratified across ST73 + ST131 + ST95 + others (top-3 share 43%).
3. **Supplement from Tier-1b (DECA + von Mentzer ETEC + prototype strains)** for STEC (need ≥20 more), EHEC (need ≥44 more), EAEC (need ≥48 more), COMMENSAL (need ≥73 more).
4. **Total cohort target:** ~350-500 strains balanced across 6-7 pathotype classes with ≥5 Mash-clade folds.

This is achievable. The substrate-survey concern that "Horesh is hybrid-labeled and partially circular" was correct as written, but per-record auditing reveals the independent subset is large enough to be the substrate backbone — IF supplemented from DECA/von Mentzer for the under-represented pathotypes.

## Open follow-ups

- **Whittam STEC Center direct-contact** (Short-term Action 2-new): still load-bearing for STEC/EHEC/EAEC supplement. Even though Horesh 1a gives 31 STEC + 6 EHEC + 2 EAEC, Tier-1b is needed to reach H2 floors.
- **NCBI Pathogen Detection host_disease facet** (Short-term Action 4-new): COMMENSAL is the weakest class across all substrates (Horesh 1a has N=2; DECA's commensal coverage unknown). Tier-3 may be the primary commensal source.
- **Per-publication methodology audit:** the 2,077 independent-label records cite 8+ distinct publications, each with its own pathotype-assignment methodology. Final H1 confidence could be elevated further by spot-checking 2-3 of the top publications (Kallonen 2017, von Mentzer 2014, Ingle 2016) to confirm their labels are clinically/epidemiologically derived rather than gene-rule-derived in the originating paper.

## Reproducibility note

Computed via Python `csv.DictReader` over `F1_genome_metadata.csv` (file at `C:\Users\Farshad\AppData\Local\Temp\horesh_F1.csv` at audit time). All counts reproducible from the public Figshare deposit (DOI 10.6084/m9.figshare.13270073, file ID 25552514, MD5 1e76642e74747a044e9e64bec5becb2d). The classifier function `classify(pathotype)` returns `independent | predicted | unknown` based on the literal presence of `"(predicted)"` substring in the Pathotype cell. The canonicalization `canon(pathotype)` strips the `"(predicted)"` suffix to group rows by pathotype class regardless of provenance.
