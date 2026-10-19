# F1 — CABBAGE / EMBL-EBI AMR Portal: the label wall is lower than the record says

**Verdict: `INCONCLUSIVE`, awaiting ONE named measurement (G8).** 9 of the 10 rejection gates pass.
This is **not** a cleared candidate and must not be quoted as one.

**What it bears on, stated plainly:** the repo's own strict-MIC 4-drug census (2026-05-18) concluded
**"NO drug clears N=150 per-class at either label-quality bar"**, with `assembly_accession` named as the
structural bottleneck (96% drop on ciprofloxacin). That verdict is load-bearing — it is the measured
basis for *"the binding constraint is LABELS not models"* and for banking the AMR track. **Run through
the repo's own `tier_for_isolate`, this substrate clears N=150 per-class at the STRICT bar on four
cells.** The premise is narrower than the record implies.

## Verified on the real surface, not from the paper

Per the real-surface-first discipline, nothing here is a literature claim:

- FTP listing fetched: `https://ftp.ebi.ac.uk/pub/databases/amr_portal/releases/2026-07/`
- `phenotype.parquet` (8.3 MB) + `phenotype_genotype_merged.parquet` downloaded to
  `D:/dna_decode_cache/amr_portal/`
- **MD5 matched the published `md5sums` EXACTLY for both files** — integrity, not just presence
- Every number below computed from the downloaded table (1,714,486 rows × 34 columns)

It is **8.3 MB and CC-BY 4.0**. This is not an expensive or request-gated acquisition; it is a
seconds-long download.

## The measurement that matters

Strict-MIC tiers from `external_mic_labels.tier_for_isolate` — the same classifier whose BV-BRC run
returned **cipro 17R/4S**:

| cell | isolates (asm+MIC) | strict R | strict S | excluded |
|---|---|---|---|---|
| E. coli × ceftriaxone | 10,568 | **2,802** | **4,212** | 34% |
| Campylobacter jejuni × ciprofloxacin | 5,058 | **1,180** | **3,532** | **7%** |
| E. coli × ciprofloxacin | 11,851 | **509** | **5,218** | 51% |
| Salmonella enterica × ciprofloxacin | 24,131 | **196** | **21,084** | 12% |
| Klebsiella pneumoniae × ciprofloxacin | 4,165 | 402 | **32** | 87% |
| Klebsiella pneumoniae × gentamicin | 4,054 | 132 | 1,307 | 64% |
| E. coli × gentamicin | 13,032 | **67** | 8,192 | 37% |

Four cells clear 150 per class. **Salmonella × cipro is `UNDERPOWERED` on the live report card at 4–5 R
against a bar of 20; here it is 196 R.**

`assembly_ID` is populated on **66.8%** of all rows (ENA GCA preferred, ERZ fallback, keyed on
BioSample) against BV-BRC's ~4% — and that single column is what the 2026-05-18 census identified as the
bottleneck.

## The gate screen (run, not asserted)

`scripts/screen_candidate_gates.py --packet`, `intended_layer=L1_AMR_RS`, on Campylobacter × cipro:

| gate | verdict | basis |
|---|---|---|
| G1 circular label | **pass** | per-row `laboratory_typing_method`, broth dilution 89.8%; AMRFinderPlus annotates GENOTYPE only |
| G2 study==class | **pass (weak)** | 47% top share — **but see the caveat below; the screen flags this itself** |
| G3 sampling-defined | **pass** | an MIC dilution reading with an explicit operator, not a collection context |
| G4 surveillance domination | **pass** | smallest non-surveillance class n=519 (bar 20) |
| G5 assembly attrition | **pass** | 5,058 records with a fetchable assembly (bar 20) |
| G6 MIC censoring | **pass** | **0%** interval-censored-undecidable (`CENSORED_EXCLUDED` = 0) |
| G7 provenance separable | **pass** | country 99.8% / collection_year 99.8% |
| **G8 dedup collapses balance** | **insufficient_data** | **the open gate** — needs `min_effective_lineages` |
| G9 / G10 rule scoreable | **pass** | the deployed QRDR rule is scoreable; 0% off-panel |

## Honest limits — read these before acting

1. **G8 is open and it is the expensive one.** `min_effective_lineages` requires genome clustering over
   the actual assemblies, which must be fetched. **9 of 10 is not a pass.** The repo's own history says
   why this gate matters: every SCORED R class is clonally dominated (≤17 effective lineages), and a
   large isolate count can collapse to a handful of lineages.
2. **Leakage is real and measured: 1,011 of the repo's 1,509 manifest accessions (67.0%) appear in
   CABBAGE.** `PATRIC` *is* BV-BRC, and `NCBI_antibiogram`/`NDARO` are the NCBI-PD family the 10 SCORED
   cells came from. Manageable rather than fatal — **102,699 CABBAGE assemblies (99.0%) are in no repo
   cohort**, so the fail-closed `cohort_manifest` gate costs almost nothing.
3. **G2 is the cohort-wide fallback and the screen says so itself ("weak").** The `database` column is an
   **aggregator label**, semicolon-joined and overlapping (`NCBI_antibiogram;NDARO`), **not a
   BioProject**. The repo's source-diversity bar is per-BioProject, so "7 labels at 47%" is **not that
   measurement**. BioProject-level diversity is **UNMEASURED** and needs BioSample→BioProject resolution.
4. **Useless for the prospective lock.** `collection_year` maxes at **2024**, with **zero** rows ≥2026,
   so nothing here is eligible against the active 2026-08-31 cutoff. A clean negative on that path.
5. **Censoring bites unevenly, and the four clearing cells are the good ones.** E. coli × cipro excludes
   51% (dominant tier `CENSORED_EXCLUDED` = 5,199: a `<=` bound sitting above the S breakpoint);
   Klebsiella × cipro excludes 87% and is S-starved at 32. Campylobacter is the clean one.
6. **E. coli × gentamicin yields only 67 strict R**, so this does **not** serve the open `rmt`/`aac(3)`
   question.
7. **`assembly_ID` is an identifier.** Fetchability is inferred from the portal's own `genomes/` tree
   (indexed by GCA/ERZ with pre-computed GFF annotations) but was **not spot-fetched**.
8. A cleared gate bounds **the label question only** — not regime fit, not reachability, not
   worth-doing. 66.8% assembly coverage is a property of the whole table; the per-cell counts above are
   the ones that matter.

## Why this is not a speculative new source

**The EBI AMR Portal already produced this project's first genuinely-independent TB number** — the
39,193-isolate provenance-disjoint cohort behind RIF sens 0.920 / INH 0.879 raw (lineage-collapsed
0.444 / 0.321). CABBAGE is that portal's expansion to the WHO priority pathogens. The source has cleared
these gates in production once, which is a materially different risk profile from an unexamined deposit.

**And the BV-BRC relationship is the interesting part.** BV-BRC supplies ~45% of CABBAGE (764k entries),
and BV-BRC is the source this project already measured as infeasible. But that failure was **assembly
linkage**, not label availability: 96% of its AST rows lacked a downloadable `assembly_accession`.
CABBAGE re-keys to ENA GCA/ERZ via BioSample — the exact field that failed. **So it plausibly rescues
labels this project already gave up on**, which is a hypothesis the 66.8% figure supports and G8 does
not yet test.

## What this does NOT claim

- It does **not** reopen the embedding tracks. Those are closed on **de-confounded signal-vs-structure**
  grounds, not on label volume, and more labels do not address them. The independent Yu/Barquist 2025
  result already in this harvest reinforces that: more training data does **not** rescue lineage keying.
- It does **not** establish that any new cell would score well. It establishes that **a label substrate
  exists at a scale the record says does not exist**.
- It says nothing about the GENERATIVE north star, whose bottleneck is measured to be the **generator**,
  not the scorer or the label supply.

## The second finding, which is larger than the label counts

**The portal also ships PRE-COMPUTED AMRFinderPlus output, and consuming it naively produces confident
wrong calls.**

Per-assembly at
`genomes/{GCA|ERZ}/<d1>/<d2>/<ACC>/<ACC>_amrfinderplus.tsv.gz` — **~1 KB per genome**. This project runs
AMRFinder through Docker at **~95 s/strain**, so the 5,058-genome Campylobacter cell alone is ~133 hours
of compute the portal has already done. `POINT` rows are present (`Method=POINTP`), so the
`qrdr_point` ciprofloxacin rule is scoreable from it.

**But the header dialect differs and the failure is silent.** The portal writes `Element_symbol` /
`Protein_id` / `Contig_id` (underscores); `eval/amr_rules` reads `r.get("Element symbol")` (spaces). A
missing key yields `None` → `""` → the row contributes no determinant → **a confident SUSCEPTIBLE
instead of an error.** Measured through the deployed `call_resistance` on real portal files:

| genome | drug | as delivered | normalized |
|---|---|---|---|
| GCA_001416835.1 | ciprofloxacin (`gyrA_T86I`) | **S (n=0)** | **R (n=1)** |
| GCA_001416835.1 | tetracycline (`tet(O)`) | **S (n=0)** | **R (n=1)** |
| GCA_001416835.1 | gentamicin (`aph(2'')-If`) | R (n=1) | R (n=1) |

**Two of three drugs were silently wrong**, including `gyrA_T86I` — the determinant this project's own
Campylobacter ciprofloxacin cell is overwhelmingly built on. **Gentamicin was right BY LUCK**: its rule
keys on `Subclass`, a column name containing no underscore. That is what makes the dialect dangerous —
a spot-check on the wrong drug makes the bug invisible.

Same failure CLASS as the TB regression already in the record: a masked-VCF `FILTER` convention mismatch
made a file with 1,382 genuine `GT=1/1` rows parse to zero calls and print a confident S. The remedy
there was to refuse rather than report.

**Shipped: `dna_decode/data/amr_portal.py`** — `header_dialect` / `normalize_header` /
`assert_repo_compatible` (fail-closed) / `prepare` / `genome_url`, with 18 tests. Two design points worth
keeping:

1. The guard is at the **SCHEMA** level, not the result level, and that is strictly better here —
   **zero determinants is a legitimate result** for a susceptible genome, so a "parsed zero rows"
   check would conflate a real negative with a parse failure.
2. `normalize_header` rewrites underscores on **line 1 only**. Determinant symbols in the body
   legitimately contain underscores (`gyrA_T86I`, `rplV_A103V`); rewriting those would corrupt exactly
   the QRDR point mutations the ciprofloxacin rule counts.

**Not claimed:** the `genomes/` tree's coverage was not censused. One of four sampled Campylobacter
accessions (`GCA_900039955.1`) was **absent from the tree**, so coverage is partial and unmeasured.

## Reproduce

```
curl -sL https://ftp.ebi.ac.uk/pub/databases/amr_portal/releases/2026-07/phenotype.parquet \
  -o D:/dna_decode_cache/amr_portal/phenotype.parquet
uv run python scripts/deep_research_screen.py
uv run pytest tests/test_deep_research_candidates.py -q
```

Harvest row: `cabbage_ebi_amr_portal` in `wiki/deep_research_candidates_2026-09-28.json`.
