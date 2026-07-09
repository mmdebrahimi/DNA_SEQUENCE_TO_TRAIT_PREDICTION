# Related data-sources registry (genotype → phenotype substrate)

**Date:** 2026-07-08 · **Branch:** `soraya-j2-fp16fix` · captured by a `/soraya` data-gathering sweep.

Map of the FREE, network-reachable data sources related to the two dna_decode substrates:
- **World-model substrate** (protein/DNA variant-effect — the ESM/foundation-model track): continuous
  fitness (DMS) + binary pathogenicity labels used to test whether a sequence model predicts variant effect.
- **Decoder substrate** (organism genotype → AMR/typing phenotype — the deployed deterministic decoder).

Each row: reachability + size **verified today** (curl), format, capture status, and the headline signal
where a source has been scored. "Captured" = fetched + flows through a working adapter; "Cataloged" =
recipe recorded, not fetched (too big / already decoder-owned).

## World-model substrate (variant-effect)

| Source | Access | Scale (verified) | Format | Status | Headline signal |
|---|---|---|---|---|---|
| **ProteinGym v1.0** | HF `ICML2022/ProteinGym` | 87 substitution assays | per-assay CSV (`mutant`,`DMS_score`) + ref seq | **Captured + scored** (`esm_zeroshot_dms.py`) | ESM2-35M masked-marginal median rho **+0.327** (partial 53/87; full rerun in flight) |
| **ProteinGym v1.1** | Zenodo `13936340` (11 GB zip; 43 MB nested substitutions zip via HTTP range-extract) | **217 assays** (140 free-CPU-tractable ≤1022 aa) | same as v1.0 | **Captured** (`proteingym_v11_fetch.py`); staged at `D:/dna_decode_cache/proteingym/` | tiny-domain (≤120 aa) ESM2-35M median rho **+0.483** |
| **MaveDB** | `api.mavedb.org/api/v1` (search + scores) | **2796 score sets** (~all `protein_coding`; ~38 % ≤400 aa tractable) | JSON score sets; HGVS `hgvs_pro`/`hgvs_nt` variants | **Adapter built + validated on real scores** (`mavedb_cpu_smoke.py` — DNA→protein translate + HGVS parse; mism=0) | **ESM2-35M signal on clean single-missense scans:** PTEN VAMP-seq rho **+0.281**, NUDT15 stability **+0.341** (vs ~0 shuffled). Multi-mutant libraries (UBE2I BarSeq/TileSeq/joint) give chance on their small/noisy single subset — ProteinGym's curation is the value-add. See `wiki/mavedb_real_score_results.json`. |
| **UniProt humsavar** | `ftp.uniprot.org/.../variants/humsavar.txt` | **8.7 MB · 84 845 missense / 13 156 proteins** (LB/B 39 962 · LP/P 33 233 · US 11 650) | fixed-column text (gene, AC, VAR_id, `p.Xxx###Yyy`, category) | **Captured + scored** (`humsavar_fetch.py`) | ESM2-35M discriminates pathogenic vs benign on **MLH1/P40692** (n=128): \|rho\| **0.314**, correct sign (pathogenic → lower ESM score) |
| **ClinVar** | `ftp.ncbi.nlm.nih.gov/.../variant_summary.txt.gz` | **440 MB** clinical variant summary | TSV (gzip) | **Cataloged** (too big for a full local capture on the 20 GB host; subset via NCBI E-utilities or stream-filter to missense) | — |

## Decoder substrate (organism → AMR/typing phenotype)

Mostly already decoder-owned (in-repo catalogs or `D:` caches); listed for completeness of the map.

| Source | Access | Role | Status |
|---|---|---|---|
| **Stanford HIVDB** (PhenoSense fold-change) | `hivdb.stanford.edu` datasets | The first FREE independent isolate-level g→p label the project validated against (NNRTI/NRTI/PI/INSTI/CAI) | Decoder-owned (`data/raw/hiv/`, gitignored) — re-fetchable |
| **Stanford CoV-RDB** | `github.com/hivdb/covid-drdb-payload` | SARS-CoV-2 Mpro measured fold-change | Decoder-owned; today's probed raw path 404'd (repo layout) — catalog only |
| **CRyPTIC** (M. tuberculosis) | Zenodo compendium (`VARIANTS.parquet` + BMD-MIC) | 12 287-isolate TB g→p; WHO-catalogue baseline substrate | Decoder-owned (`D:`) |
| **EBI AMR-Portal** (TB disjoint) | EBI | Provenance-disjoint independent TB DST cohort (39 193 isolates) | Decoder-owned (`data/raw/tb_goldset/`) |
| **BV-BRC / NCBI-PD** | `ftp.bvbrc.org`, NCBI Pathogen Detection | E. coli AMR cohorts (the deployed 6-drug surface) | Decoder-owned |
| **WHO TB catalogue v2 / Napier barcode** | WHO / tbdb | Deterministic TB determinant rule + lineage barcode | Pinned in-repo (`data/raw/`) |

## Baseline predictors & alternative phenotype axes (cataloged)

Reachability verified today; all catalog-only (too big to commit / not variant-level-ESM-scorable).

| Source | Access | Scale (verified) | Role |
|---|---|---|---|
| **AlphaMissense** | `storage.googleapis.com/dm_alphamissense/AlphaMissense_aa_substitutions.tsv.gz` | **1.2 GB** (stream-filter per UniProt) | Per-variant pathogenicity **predictor** — baseline to beat. **Head-to-head (7 proteins, median AUROC):** AM **0.823**. ESM2-**35M** 0.706 (gap +0.117) → ESM2-**650M** **0.811** (gap **+0.012**, near-parity; ESM beats AM on 2/7). **Scale closes the gap — "supervised wins" was a 35M model-size artifact, not a regime property** (`scripts/{humsavar_am_vs_esm,humsavar_am_vs_esm_sweep,kaggle_esm_am_sweep}.py`, `wiki/humsavar_am_vs_esm_sweep{,_650M}.json`). Caveat: AM retains home-field advantage (its training likely overlaps these labels), so ESM reaching parity despite that is the notable result; PROC stays hard for both. |
| **DepMap** (CRISPR gene effect) | figshare (`ndownloader.figshare.com/files/43346616`, 302→S3) | redirect | Gene-**essentiality** fitness across cell lines — a *different* phenotype axis (gene-level, not missense) |
| **gnomAD constraint** | `storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/…tsv` | **95 MB** | Population variant-**tolerance** proxy (o/e, pLI) — gene-level; weak-label for missense tolerance |

## Notes / recipes

- **ProteinGym v1.1 range-extract** — the 11 GB Zenodo zip is not fully downloaded; the substitutions
  sub-archive is pulled via an HTTP `Range` request against the central-directory offset. See
  `scripts/proteingym_v11_fetch.py`.
- **humsavar → ESM scoring** — fetch a protein's sequence (`rest.uniprot.org/uniprotkb/<AC>.fasta`),
  parse its `p.Xxx###Yyy` rows to `(wt,pos,mut)` + label, score masked-marginals, and correlate the ESM
  score against the binary label (worked example: `wiki/humsavar_mlh1_35M_pilot.json`). Proteins >1022 aa
  (most disease genes: ABCA4, SCN5A, RYR1, …) need `--long-mode window`.
- **MaveDB** — the search API is reliable; the per-scoreset `/scores` endpoint 504s intermittently. Prefer
  small (≤400 aa) `protein_coding` score sets for CPU scoring; the adapter translates DNA targets.
- **ClinVar** — for a tractable slice, filter `variant_summary` to `Type == single nucleotide variant` +
  a missense consequence, or query E-utilities per gene, rather than committing the 440 MB dump.
