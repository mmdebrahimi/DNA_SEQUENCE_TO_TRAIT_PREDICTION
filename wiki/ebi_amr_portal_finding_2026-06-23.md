# EBI AMR Portal (CABBAGE) — the free measured-AST source that reframes the label wall (2026-06-23)

User surfaced the EBI `/pub/databases` FTP index. Screened against the project's one binding need — a FREE,
INDEPENDENT, isolate-level, MEASURED-phenotype genotype↔phenotype source. **One folder is a potential unlock
for the entire AMR label wall: `amr_portal/`.**

## The find: EBI AMR Portal (MRC CABBAGE project)
`https://ftp.ebi.ac.uk/pub/databases/amr_portal/` · web `https://www.ebi.ac.uk/amr` · "the largest public AMR
dataset in a reconciled, uniform format." Releases monthly (latest 2025-12); Parquet + DuckDB + CSV;
`phenotypes`, `genotypes`, and a `phenotype_genotype_merged` table; genomes indexed by assembly accession.

**Schema (per-isolate, ideal):** `BioSample_ID` (SAMEA…) + `SRA_accession` (ERS…) + `assembly_ID` (GCA…) —
full accession linkage; `organism`/`species`, `antibiotic_name`, `ast_standard` (CLSI/EUCAST),
`laboratory_typing_method`, `measurement`+`measurement_sign`+`measurement_units` (MIC), `resistance_phenotype`
(measured S/I/R), `Updated_phenotype_{CLSI,EUCAST}`, `database` (provenance), `collection_year`, country.

**Coverage (scan of 900k rows):** M. tuberculosis **183k**, Salmonella 73k, Klebsiella 52k, E. coli 28k,
S. aureus 19k, S. pneumoniae 15k, Acinetobacter, N. gonorrhoeae. Drugs: rifampin/isoniazid/ethambutol/PZA
(TB), ciprofloxacin 18k, gentamicin, TMP-SMX, tetracycline 13k, ampicillin, moxifloxacin, amikacin…

**Non-circular ✓ (the gate that killed the TB literature sources):** 674k/900k rows have a measured MIC;
methods = broth dilution 421k / agar dilution 76k / disk diffusion 46k / E-test 27k (all LAB AST). The label
is a wet-lab measurement, not a genotype-derived prediction.

## The make-or-break (honest — don't overclaim)
`database` provenance = `CABBAGE_PubMed_data` (literature) + **NARMS** 190k + **PATRIC** (=BV-BRC) + **NCBI_antibiogram;NDARO** (=NCBI Pathogen Detection) + pathogenwatch + pubMLST + COMPARE_ML. So it AGGREGATES the
exact sources our cohorts come from (NCBI-PD, BV-BRC) → **overlap with our validated cells is real**, and TB
likely overlaps CRyPTIC. It does NOT auto-deliver independence. BUT:
- It carries **BioSample/assembly accessions** → our EXISTING leakage machinery (`dna_decode/eval/cohort_manifest.py`, `biosample_resolver.py`, the external-cohort arm) can carve out the **provenance-disjoint subset
  precisely** — the same gate we already ship.
- It has sources we did NOT use (NARMS, pathogenwatch, pubMLST, COMPARE_ML, literature) → a genuinely
  independent subset almost certainly exists.
- **No DUA, no author contact, no money** — exactly what we could NOT get for the TB gold set (5 sources, all
  author-request/DUA/circular). This makes the independent number MEASURABLE + FREE.

## Why this matters (the two highest-VOI moves it unlocks)
1. **The independent TB number — potentially free, no author contact.** 183k M. tuberculosis rows with measured
   DST. Run our CRyPTIC leakage check (`tb_goldset.assert_independent_aliased`) over the AMR-Portal TB
   BioSamples → the CRyPTIC-disjoint subset with measured RIF/INH = the independent TB gold set we emailed
   authors for. Even a fraction of 183k clears the powering bar. **This could make the Thorpe/India/Ethiopia
   author emails moot.**
2. **External validation for the FROZEN bacterial cells.** 28k E. coli (cipro/tet/…) + Klebsiella + Salmonella
   + S. aureus + N. gonorrhoeae with measured AST. Run the NCBI-PD-disjoint subset through the external-cohort
   revalidation arm (already built) → a 2nd independent number for the deployed cells, and candidate NEW cells
   (Salmonella, N. gonorrhoeae) each with a free measured label.

## Recommended next move
A focused AMR-Portal feasibility + leakage pass (mirrors the CoV-RDB feasibility we just did): download the
2025-12 `phenotype_genotype_merged.parquet` to D:, run the BioSample-level leakage check vs CRyPTIC (TB) and
vs our parquet cohorts (E. coli), and REPORT the provenance-disjoint per-(organism,drug) powering. That tells
us, per cell, whether a free independent number is in reach — before any scoring. Highest-VOI move in the
project right now (it attacks the binding constraint head-on).

## Other EBI folders screened (lower priority)
- `cryptic/` (release_june2022) — our EXISTING TB CRyPTIC source (the in-distribution baseline). Confirms it.
- `ont_tb_eval2022/` — TB ONT reads + MGIT phenotype (Birmingham/Madagascar); raw reads, niche.
- `impc/` — mouse knockout→phenotype (measured; a future NON-AMR organism axis, not now).
- `AllTheBacteria/` + `ENA2018-bacteria-661k/` — massive bacterial GENOME substrates (genomes only, no label).
- `mutfunc/` / `ProtVar/` — variant-effect (curated/predicted → circular-ish; the DMS shape we deprioritized).
- `gwas/` (association, weak), `eva/`, `ega/` (DUA-gated), `opentargets/`, `chembl/` (compound bioactivity, not
  genotype→phenotype) — not the project's shape.

## Provenance
Screened 2026-06-23 from the EBI `/pub/databases` index (user-supplied). AMR Portal release 2025-12
(`releases.yml` latest). Non-circular + provenance scan from streaming `phenotype.csv.gz` (900k rows).
Leakage machinery already in-repo: `dna_decode/eval/{cohort_manifest,biosample_resolver}.py`,
`scripts/external_cohort_*`. Companion: `wiki/tb_goldset_public_source_exhaustion_2026-06-22.md` (the wall this
may break), `wiki/next_independent_label_cell_feasibility_2026-06-23.md` (the same feasibility-first method).
