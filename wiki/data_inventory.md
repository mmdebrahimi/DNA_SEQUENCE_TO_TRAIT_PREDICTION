# Data inventory — what we have, and what already used it (2026-11-22)

**Generated. Do not hand-edit** — run `uv run python scripts/data_inventory.py`.

Read this before claiming a dataset is missing, or that a piece of work would be expensive because the data is not on disk. `consumed_by` / `reported_in` are DERIVED from the real tree, so they answer *has anything already used this?* — the question whose wrong answer causes rediscovery.

- datasets present: **62** (8 with a curated note, 54 undocumented)
- referenced nowhere in code or artifacts: **1**
- **roots not reachable right now: 1** — `D:\dna_decode_cache` (an unplugged drive is not an absent dataset)
- ⚠ curated notes naming datasets that no longer exist: `caur_sra`, `cryptic`, `data files donwload`, `ecoli_sero_asm`, `embeddings`, `essentiality`, `ktype_kaptive`, `models`, `nt-gate-b`, `prospective`, `proteingym`, `refseq`, `revalidation_2026-07-11`, `ss2_checkpoint`, `tb_goldset`, `tb_indep`

## Cohorts / raw inputs (`data/raw`)

| dataset | size | used by | reported in | what it is |
|---|---|---|---|---|
| `campylobacter_provdisjoint_ciprofloxacin` | 0.0MB | 6 | 6 | _(undocumented)_ |
| `ecoli_conditional_essentiality` | 0.0MB | 1 | 2 | _(undocumented)_ |
| `escherichia_coli_shigella_provdisjoint_ceftriaxone` | 0.0MB | 6 | 7 | _(undocumented)_ |
| `escherichia_coli_shigella_provdisjoint_ciprofloxacin` | 0.0MB | 5 | 5 | _(undocumented)_ |
| `escherichia_coli_shigella_provdisjoint_gentamicin` | 0.0MB | 5 | 7 | _(undocumented)_ |
| `escherichia_coli_shigella_provdisjoint_tetracycline` | 0.0MB | 1 | 5 | _(undocumented)_ |
| `klebsiella_provdisjoint_ceftriaxone` | 0.0MB | 2 | 4 | _(undocumented)_ |
| `klebsiella_provdisjoint_ciprofloxacin` | 0.0MB | 10 | 7 | _(undocumented)_ |
| `klebsiella_provdisjoint_gentamicin` | 0.0MB | 2 | 5 | _(undocumented)_ |
| `klebsiella_provdisjoint_meropenem` | 0.0MB | 7 | 6 | _(undocumented)_ |
| `klebsiella_provdisjoint_tetracycline` | 0.0MB | 1 | 4 | _(undocumented)_ |
| `tb_lineage_barcode` | 0.1MB | 7 | 2 | _(undocumented)_ |
| `who_tb_catalogue` | 0.0MB | 10 | 5 | _(undocumented)_ |

## Reference databases (`data/`)

| dataset | size | used by | reported in | what it is |
|---|---|---|---|---|
| `amrfinder_db` | 0.0MB | 19 | 17 | AMRFinderPlus reference database. |
| `amrfinder_runs` | 0.0MB | 45 | 20 | Cached AMRFinderPlus output (main.tsv + mutations.tsv) per assembly. |
| `antimalarial_ref` | 0.0MB | 5 | 2 | _(undocumented)_ |
| `antiviral_ref` | 0.0MB | 2 | 1 | _(undocumented)_ |
| `arabidopsis` | 1.2MB | 17 | 37 | Arabidopsis flowering-time (FT10) phenotype + accession data. |
| `clinvar` | 2.3MB | 14 | 23 | ClinVar variant-classification release. |
| `cyp2d6_psv` | 0.0MB | 2 | 2 | _(undocumented)_ |
| `disinfinder_db` | 0.0MB | 2 | 2 | _(undocumented)_ |
| `ena_wgs` | 0.0MB | 8 | 4 | _(undocumented)_ |
| `external` | 0.0MB | 103 | 248 | External-cohort validation inputs + crosswalks (Oxford, PROBAC, AR Bank, Sci234). |
| `fba` | 0.0MB | **0** | **0** | _(undocumented)_ |
| `forward_ref` | 0.0MB | 11 | 6 | _(undocumented)_ |
| `fungal_ref` | 0.0MB | 14 | 8 | _(undocumented)_ |
| `glm_diffusion_ckpt` | 9.4MB | 2 | 2 | _(undocumented)_ |
| `gono_catalogue` | 0.0MB | 1 | 1 | _(undocumented)_ |
| `hcmv_ref` | 0.0MB | 2 | 1 | _(undocumented)_ |
| `hiv_ref` | 0.0MB | 24 | 5 | Committed HXB2 CDS references (RT / PR / IN / CA). |
| `horesh_2021` | 0.0MB | 1 | **0** | Horesh H4 E. coli pathotype cohort (24 genomes). |
| `horse` | 0.0MB | 20 | 29 | _(undocumented)_ |
| `imputation` | 0.0MB | 21 | 28 | _(undocumented)_ |
| `isaba1_ref` | 0.0MB | 8 | 5 | _(undocumented)_ |
| `j3_abo` | 0.8MB | 3 | 3 | _(undocumented)_ |
| `ktype_db` | 0.0MB | 4 | 4 | _(undocumented)_ |
| `mlst_db` | 0.0MB | 2 | 3 | _(undocumented)_ |
| `mpra` | 95.9MB | 43 | 60 | _(undocumented)_ |
| `mpra_gse108535` | 8.7MB | 3 | 2 | _(undocumented)_ |
| `mpra_human` | 70.8MB | 9 | 8 | _(undocumented)_ |
| `mycoplasma_ref` | 0.0MB | 2 | 2 | _(undocumented)_ |
| `nt_windows` | 0.0MB | 1 | 1 | _(undocumented)_ |
| `pa_clone_durable` | 83.7MB | 2 | 3 | _(undocumented)_ |
| `pathotype_cov_cache` | 0.0MB | 9 | 4 | _(undocumented)_ |
| `pathotype_pergene_cache` | 0.0MB | 9 | 8 | _(undocumented)_ |
| `pgx_1000g` | 0.0MB | 19 | 2 | _(undocumented)_ |
| `pgx_getrm` | 0.0MB | 11 | 12 | _(undocumented)_ |
| `phage_ref` | 0.1MB | 7 | 6 | _(undocumented)_ |
| `plasmidfinder_db` | 0.0MB | 5 | 1 | _(undocumented)_ |
| `pneumo_betalactam_db` | 0.0MB | 2 | 1 | _(undocumented)_ |
| `pneumoserotype_db` | 0.0MB | 5 | 2 | _(undocumented)_ |
| `pointfinder_db` | 0.0MB | 4 | 2 | _(undocumented)_ |
| `promcalc` | 162.0MB | 23 | 25 | _(undocumented)_ |
| `promoteratlas` | 229.6MB | 27 | 33 | _(undocumented)_ |
| `refseq_cache` | 4309.2MB | 24 | 14 | _(undocumented)_ |
| `resfinder_db` | 2.0MB | 6 | 1 | _(undocumented)_ |
| `salmserovar_db` | 0.6MB | 14 | 5 | _(undocumented)_ |
| `sarscov2_ref` | 0.0MB | 12 | 6 | _(undocumented)_ |
| `serotypefinder_db` | 1.4MB | 5 | 1 | _(undocumented)_ |
| `summary_stat_sources` | 1.6MB | 4 | 2 | _(undocumented)_ |
| `unattended` | 0.0MB | 20 | 20 | _(undocumented)_ |
| `virulencefinder_db` | 7.4MB | 14 | 9 | VirulenceFinder allele DB (E. coli-scoped). |

## Datasets nothing references

Either dead weight or a forgotten asset. This is the bucket rediscovery comes from.

- `fba` (ref, 0.0MB) — `C:\Users\Laptop Portable\pythonprojects\dna_decode\data\fba`
