# HerpesDRG × our HCMV cell: the G1 circularity screen (2026-10-03)

**Verdict: `CONDITIONAL_PASS` — overlap is REAL, SMALL, and mechanically EXCLUDABLE by DOI.**
HerpesDRG is **not** an independent isolate-level label set (it is per-mutation, the CoV-RDB shape), but
it **is** a large held-out-STUDY set, which is precisely what our own HCMV module says an independent
number requires.

## Why this screen had to run first

`dna_decode/data/hcmv_amr.py` sources its catalog **verbatim** from three Chou-lab compilations:

- `PMC3262590` — UL54 fold-change Table 1 + the recombinant-phenotyped BENIGN list
- `PMC5483911` — UL97 codon 591-603 point + deletion fold-changes
- `AAC 10.1128/aac.00922-18` — the UL56 letermovir recombinant EC50 table

and its own docstring already concedes: *"Validation is IN-DISTRIBUTION against the SAME measured
fold-changes the catalog is curated from (a knowledge baseline; an independent number needs held-out
phenotyping studies)."*

HerpesDRG is a systematic review of that same genotype-to-phenotype literature. So the whole question
was whether comparing them is a catalog-vs-catalog restatement (G1 trip) or a genuine held-out read.
**That is decidable by measurement because HerpesDRG records a per-entry source DOI.**

## Measured (`data/raw/herpesdrg/herpesdrg-db.tsv`, 705,930 bytes, MIT)

Filtered to `virus == HCMV` and `status == "A"` — the README's own activity rule.

| quantity | value |
|---|---|
| HCMV active entries | **644** |
| distinct source DOIs | **99** |
| entries from our one exactly-known source DOI (`10.1128/aac.00922-18`) | **12** |
| entries NOT from that DOI | **632 (98.1%)** |
| largest single DOI's share | 35 entries (**5.4%**) |

Genes: UL54 290 · UL97 200 · UL56 106 · UL27 25 · UL89 14 · UL51 9.

Numeric fold-changes are abundant, so the label is quantitative, not a word:

| drug | non-empty | numeric fold-change |
|---|---|---|
| Ganciclovir | 466 | **457** |
| Foscarnet | 292 | **280** |
| Cidofovir | 285 | **277** |
| Letermovir | 155 | **155** |
| Maribavir | 104 | **104** |

`strain_generation`: Recombinant BAC 447 · Marker Transfer 95 · Isolated strain 50 · lab-derived /
drug-induced mutants 51.

## A proxy error I made and am not reporting as a result

My first pass counted studies whose `ref_title` contained "Chou" and got **0 of 644**, which would have
read as "100% non-Chou — no overlap at all". **`ref_title` is the paper title, not the author list**, so
searching it for an author name returns zero by construction — the vacuous-filter shape this repo has on
record. The titles themselves ("Recombinant Phenotyping of Cytomegalovirus UL54 Mutations…") are plainly
Chou-lab papers. The DOI table above replaces it.

## The verdict, with its boundary

**Three distinct readings, and only one is available:**

1. **Independent isolate-level validation — NO.** 592 of 644 entries are engineered single mutants
   (Recombinant BAC / marker transfer). There are no clinical isolates with co-occurring mutations, so
   this cannot test isolate-level performance or epistasis. Same ceiling as the SARS-CoV-2 Mpro cell.
2. **In-distribution restatement — NO, and this is the useful part.** 99 distinct studies against our 3
   sources means the overwhelming majority of HerpesDRG's HCMV content is literature our catalog never
   read. Excluding the overlap is mechanical: filter on `ref_doi`.
3. **HELD-OUT KNOWLEDGE BASELINE — YES, available now.** Score the deployed `hcmv_amr` catalog against
   HerpesDRG entries whose DOI is NOT one of our three sources. That is a stronger tier than the cell's
   current self-declared in-distribution status, and weaker than isolate-level independence. It would
   move the HCMV cell off `NO_FREE_PHENOTYPE` — which is the field this screen was aimed at.

## Honest limits

- **The overlap figure is a LOWER BOUND.** I verified one of our three sources by exact DOI
  (`10.1128/aac.00922-18` → 12 entries). The other two are PMC identifiers (`PMC3262590`,
  `PMC5483911`) whose DOIs I did not resolve, so true overlap is **≥12 entries**, not exactly 12.
  Resolving those two DOIs is the one remaining step before any number is quoted.
- HerpesDRG ships its **own** `call_resistance` interpreter. Using its *interpretation* would be
  circular against our rule; only its **fold-change table** is the label. Do not confuse the two.
- Multiple entries exist per variant across publications and assay methods, so a per-mutation
  aggregation rule (max? median? per-assay?) is a real design choice, not a detail — and the clinical
  binning HerpesDRG itself uses (>15 high / 5–15 moderate / 2–5 low / <2 polymorphism) is **its**
  threshold set, not ours.
- Nothing has been scored. This screen establishes only that scoring is **not** circular by
  construction.

## Recommended next step

Resolve the two PMC DOIs → exclude all overlapping `ref_doi` rows → score the deployed HCMV catalog
against the held-out remainder, reported as a **held-out knowledge baseline** with the engineered-mutant
ceiling stated. That is a new, free, quantitative number for a cell that currently claims none.
