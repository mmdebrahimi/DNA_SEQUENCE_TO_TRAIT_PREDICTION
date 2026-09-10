# Data inventory — what we have, and what already used it (2026-09-10)

**Generated. Do not hand-edit** — run `uv run python scripts/data_inventory.py`.

Read this before claiming a dataset is missing, or that a piece of work would be expensive because the data is not on disk. `consumed_by` / `reported_in` are DERIVED from the real tree, so they answer *has anything already used this?* — the question whose wrong answer causes rediscovery.

- datasets present: **248** (28 with a curated note, 220 undocumented)
- referenced nowhere in code or artifacts: **59**

## Cohorts / raw inputs (`data/raw`)

| dataset | size | used by | reported in | what it is |
|---|---|---|---|---|
| `_biosample_cache.json` | 0.2MB | 1 | **0** | _(undocumented)_ |
| `acinetobacter_meropenem` | 439.6MB | 2 | 1 | _(undocumented)_ |
| `acinetobacter_meropenem_indep` | 506.1MB | 2 | **0** | _(undocumented)_ |
| `acinetobacter_meropenem_run.log` | 0.0MB | **0** | **0** | _(undocumented)_ |
| `ar_bank_caur` | 316.3MB | 3 | 1 | _(undocumented)_ |
| `ar_bank_caur_extval_fluconazole` | 0.0MB | 1 | **0** | _(undocumented)_ |
| `ar_bank_caur_extval_micafungin` | 0.0MB | 1 | **0** | _(undocumented)_ |
| `ar_bank_caur_extval_voriconazole` | 0.0MB | **0** | **0** | _(undocumented)_ |
| `ar_bank_ecoli_extval_ceftriaxone` | 250.8MB | **0** | 1 | _(undocumented)_ |
| `ar_bank_ecoli_extval_ciprofloxacin` | 22.2MB | **0** | **0** | _(undocumented)_ |
| `ar_bank_ecoli_extval_gentamicin` | 41.5MB | **0** | **0** | _(undocumented)_ |
| `ar_bank_efm_kaggle` | 1003.1MB> | 1 | **0** | _(undocumented)_ |
| `ar_bank_enterococcus_faecium_extval_doxycycline` | 0.0MB | 1 | **0** | _(undocumented)_ |
| `ar_bank_enterococcus_faecium_extval_levofloxacin` | 0.0MB | 2 | **0** | _(undocumented)_ |
| `ar_bank_enterococcus_faecium_extval_teicoplanin` | 0.0MB | 1 | **0** | _(undocumented)_ |
| `ar_bank_enterococcus_faecium_extval_vancomycin` | 0.0MB | 1 | **0** | _(undocumented)_ |
| `ar_bank_gono_extval_azithromycin` | 0.0MB | **0** | **0** | _(undocumented)_ |
| `ar_bank_gono_extval_cefixime` | 0.0MB | **0** | **0** | _(undocumented)_ |
| `ar_bank_gono_extval_ceftriaxone` | 0.0MB | **0** | **0** | _(undocumented)_ |
| `ar_bank_gono_extval_ciprofloxacin` | 205.8MB | **0** | **0** | _(undocumented)_ |
| `ar_bank_gono_extval_penicillin` | 0.0MB | **0** | **0** | _(undocumented)_ |
| `ar_bank_gono_extval_tetracycline` | 0.0MB | **0** | **0** | _(undocumented)_ |
| `ar_bank_gono_kaggle` | 78.7MB | 1 | **0** | _(undocumented)_ |
| `ar_bank_gono_kaggle_subset.json` | 0.0MB | **0** | **0** | _(undocumented)_ |
| `ar_bank_gono_smoke` | 8.6MB | 1 | **0** | _(undocumented)_ |
| `ar_bank_kleb_extval_ceftriaxone` | 723.8MB | **0** | **0** | _(undocumented)_ |
| `ar_bank_kleb_extval_ciprofloxacin` | 607.2MB | **0** | **0** | _(undocumented)_ |
| `ar_bank_kleb_extval_gentamicin` | 204.6MB | **0** | **0** | _(undocumented)_ |
| `ar_bank_neisseria_gonorrhoeae_extval_cefixime` | 0.0MB | 1 | **0** | _(undocumented)_ |
| `ar_bank_neisseria_gonorrhoeae_extval_ciprofloxacin` | 0.0MB | 1 | **0** | _(undocumented)_ |
| `ar_bank_staphylococcus_aureus_extval_levofloxacin` | 117.6MB | **0** | **0** | _(undocumented)_ |
| `ar_isolate_bank` | 49.3MB | 5 | **0** | _(undocumented)_ |
| `campy_ncbipd_extval` | 0.0MB | 1 | 1 | _(undocumented)_ |
| `campylobacter_cipro_run.log` | 0.0MB | **0** | **0** | _(undocumented)_ |
| `campylobacter_ciprofloxacin` | 200.8MB | 2 | 3 | _(undocumented)_ |
| `campylobacter_indep_ciprofloxacin` | 202.5MB | 1 | **0** | _(undocumented)_ |
| `campylobacter_provdisjoint_ciprofloxacin` | 272.3MB | 6 | 3 | _(undocumented)_ |
| `carbon_util_bacdive` | 13.6MB | 1 | **0** | _(undocumented)_ |
| `cryptic` | 102.8MB | 22 | 12 | CRyPTIC M. tuberculosis compendium (12,287 isolates, measured BMD MIC). |
| `ecoli_conditional_essentiality` | 0.0MB | 1 | 2 | _(undocumented)_ |
| `ecoli_trn` | 0.1MB | **0** | 1 | _(undocumented)_ |
| `enterobacter_cloacae_cef_run.log` | 0.0MB | **0** | **0** | _(undocumented)_ |
| `enterobacter_cloacae_ceftriaxone` | 210.1MB | **0** | **0** | _(undocumented)_ |
| `escherichia_coli_shigella_provdisjoint_ceftriaxone` | 1030.8MB | 6 | 4 | _(undocumented)_ |
| `escherichia_coli_shigella_provdisjoint_ciprofloxacin` | 1154.8MB | 5 | 2 | _(undocumented)_ |
| `escherichia_coli_shigella_provdisjoint_gentamicin` | 1194.3MB | 5 | 4 | _(undocumented)_ |
| `escherichia_coli_shigella_provdisjoint_tetracycline` | 893.8MB | 1 | 2 | _(undocumented)_ |
| `fungal_ref` | 0.0MB | 13 | 7 | _(undocumented)_ |
| `gono_ncbipd_extval` | 0.0MB | 4 | 1 | _(undocumented)_ |
| `hiv` | 19.6MB | **0** | **0** | _(undocumented)_ |
| `independent_cohort_validate_run.log` | 0.0MB | **0** | **0** | _(undocumented)_ |
| `klebsiella_ceftriaxone` | 1.5MB | **0** | **0** | _(undocumented)_ |
| `klebsiella_cipro` | 1510.8MB | 8 | 6 | _(undocumented)_ |
| `klebsiella_ciprofloxacin` | 0.0MB | **0** | **0** | _(undocumented)_ |
| `klebsiella_gentamicin` | 0.0MB | 1 | 2 | _(undocumented)_ |
| `klebsiella_indep_ciprofloxacin` | 410.7MB | 1 | **0** | _(undocumented)_ |
| `klebsiella_meropenem` | 0.1MB | 1 | **0** | _(undocumented)_ |
| `klebsiella_provdisjoint_ceftriaxone` | 1257.7MB | 2 | 1 | _(undocumented)_ |
| `klebsiella_provdisjoint_ciprofloxacin` | 1195.6MB | 10 | 4 | _(undocumented)_ |
| `klebsiella_provdisjoint_gentamicin` | 1206.3MB | 2 | 2 | _(undocumented)_ |
| `klebsiella_provdisjoint_meropenem` | 1207.8MB | 1 | 2 | _(undocumented)_ |
| `klebsiella_provdisjoint_tetracycline` | 1258.6MB | 1 | 1 | _(undocumented)_ |
| `klebsiella_tetracycline` | 0.6MB | **0** | **0** | _(undocumented)_ |
| `mavedb_full_esm2` | 1.0MB | 1 | **0** | _(undocumented)_ |
| `ncbi_pathogen_xsource` | 430.7MB | 1 | **0** | _(undocumented)_ |
| `oxford` | 36.4MB | 9 | 8 | _(undocumented)_ |
| `pear` | 3.2MB | 5 | 7 | _(undocumented)_ |
| `pneumo_ncbipd_extval` | 0.0MB | 1 | 1 | _(undocumented)_ |
| `pseudomonas_aeruginosa_ciprofloxacin` | 711.4MB | **0** | 2 | _(undocumented)_ |
| `pseudomonas_aeruginosa_meropenem` | 167.6MB | 1 | 1 | _(undocumented)_ |
| `pseudomonas_meropenem_run.log` | 0.0MB | **0** | **0** | _(undocumented)_ |
| `revalidation_2026-07-11` | 5260.7MB | **0** | 1 | A 5.3GB revalidation cohort staged 2026-07-11. |
| `salmonella_cipro_run.log` | 0.0MB | **0** | **0** | _(undocumented)_ |
| `salmonella_ciprofloxacin` | 545.3MB | 2 | 1 | _(undocumented)_ |
| `salmonella_indep_ciprofloxacin` | 546.9MB | **0** | **0** | _(undocumented)_ |
| `sarscov2` | 0.0MB | 8 | 5 | _(undocumented)_ |
| `sci234` | 1.4MB | 4 | 7 | _(undocumented)_ |
| `staph_ncbipd_extval` | 0.0MB | 1 | 1 | _(undocumented)_ |
| `staphylococcus_aureus_oxacillin` | 292.6MB | **0** | **0** | _(undocumented)_ |
| `tb_goldset` | 1.7MB | 15 | 14 | EBI AMR-Portal provenance-disjoint TB cohort (39,193 isolates w/ DST). |
| `tb_indep` | 282.1MB | 10 | 12 | Per-isolate TB independent-run VCFs + results checkpoint (2,845 isolates). |
| `tb_lineage_barcode` | 0.1MB | 7 | 2 | _(undocumented)_ |
| `who_tb_catalogue` | 43.9MB | 8 | 5 | _(undocumented)_ |

## Processed artifacts (`data/processed`)

| dataset | size | used by | reported in | what it is |
|---|---|---|---|---|
| `ecoli_mash_clades.json` | 0.0MB | 1 | **0** | _(undocumented)_ |
| `embeddings` | 2350.8MB | 40 | 103 | HDF5 foundation-model embedding caches (NT, mostly cipro cohorts). |
| `forward_blosum_proteingym.jsonl` | 0.0MB | 1 | **0** | _(undocumented)_ |
| `forward_esm_proteingym.jsonl` | 0.0MB | 1 | **0** | _(undocumented)_ |
| `gate_a_cache.h5` | 3.7MB | **0** | 1 | _(undocumented)_ |
| `gate_a_cohort.parquet` | 0.0MB | **0** | 1 | _(undocumented)_ |
| `gate_b_cohort.parquet` | 0.0MB | 11 | 24 | _(undocumented)_ |
| `gate_b_mini_cef_cohort.parquet` | 0.0MB | 3 | 6 | _(undocumented)_ |
| `gate_b_mini_cohort.parquet` | 0.0MB | 4 | 5 | _(undocumented)_ |
| `gate_b_mini_tet_cohort.parquet` | 0.0MB | 3 | 4 | _(undocumented)_ |
| `gate_b_n40_cipro_cohort.parquet` | 0.0MB | 8 | 15 | _(undocumented)_ |
| `gene_lookup_cache` | 0.0MB | 1 | **0** | _(undocumented)_ |
| `hiv_rt_esm650m_masked_marginals.json` | 0.2MB | 3 | **0** | _(undocumented)_ |
| `mini_cipro_dnabert2_cache.h5` | 0.0MB | **0** | **0** | _(undocumented)_ |
| `mini_cipro_mock_cache.h5` | 52.0MB | **0** | **0** | _(undocumented)_ |
| `mini_cipro_nt_cache.h5` | 140.5MB | 2 | 1 | _(undocumented)_ |
| `models` | 0.8MB | 59 | 138 | Trained classifier pickles + their provenance blocks. |
| `msa_eval_scratch` | 1.8MB | 1 | **0** | _(undocumented)_ |
| `msat_lift_checkpoint.jsonl` | 0.0MB | 1 | 1 | _(undocumented)_ |
| `msat_lift_run20.log` | 0.0MB | **0** | **0** | _(undocumented)_ |
| `pilot_report_phase2_entry_2026-05-13.md` | 0.0MB | **0** | **0** | _(undocumented)_ |
| `plasmid_axis_cache.json` | 0.1MB | 4 | **0** | _(undocumented)_ |
| `plasmid_axis_checkpoint.jsonl` | 0.1MB | 1 | **0** | _(undocumented)_ |
| `prosst_lift_checkpoint.jsonl` | 0.0MB | 7 | 3 | _(undocumented)_ |
| `protein_effect_cache` | 0.3MB | 2 | **0** | _(undocumented)_ |
| `proteingym_inverse_sweep.jsonl` | 0.1MB | 1 | **0** | _(undocumented)_ |
| `proteingym_inverse_sweep_esm.jsonl` | 0.1MB | 1 | **0** | _(undocumented)_ |
| `shared_lineage_ceftriaxone_cohort.parquet` | 0.0MB | **0** | 2 | _(undocumented)_ |
| `shared_lineage_ceftriaxone_dl.parquet` | 0.0MB | **0** | **0** | _(undocumented)_ |
| `shared_lineage_gentamicin_cohort.parquet` | 0.0MB | **0** | 2 | _(undocumented)_ |
| `shared_lineage_gentamicin_dl.parquet` | 0.0MB | **0** | **0** | _(undocumented)_ |
| `shared_lineage_meropenem_cohort.parquet` | 0.0MB | **0** | 2 | _(undocumented)_ |
| `shared_lineage_meropenem_dl.parquet` | 0.0MB | **0** | **0** | _(undocumented)_ |
| `shared_lineage_tetracycline_cohort.parquet` | 0.0MB | **0** | 4 | _(undocumented)_ |
| `shared_lineage_tetracycline_pilot.parquet` | 0.0MB | **0** | **0** | _(undocumented)_ |
| `stage2_n150_cipro_cohort.parquet` | 0.0MB | 15 | 20 | _(undocumented)_ |
| `tb_mic_features_cache.json` | 3.3MB | 2 | **0** | _(undocumented)_ |
| `tb_mic_panel_features_cache.json` | 11.4MB | 1 | **0** | _(undocumented)_ |
| `tb_rif_rpob_cache.json` | 1.2MB | 1 | **0** | _(undocumented)_ |
| `three_way_lift_checkpoint.jsonl` | 0.0MB | 1 | 1 | _(undocumented)_ |
| `virulence_axis_cache.json` | 0.0MB | 3 | **0** | _(undocumented)_ |
| `virulence_axis_checkpoint.jsonl` | 0.0MB | 1 | **0** | _(undocumented)_ |

## Reference databases (`data/`)

| dataset | size | used by | reported in | what it is |
|---|---|---|---|---|
| `amrfinder_db` | 478.2MB | 19 | 17 | AMRFinderPlus reference database. |
| `amrfinder_runs` | 53.2MB | 42 | 16 | Cached AMRFinderPlus output (main.tsv + mutations.tsv) per assembly. |
| `antimalarial_ref` | 0.0MB | 5 | 2 | _(undocumented)_ |
| `antiviral_ref` | 0.0MB | 2 | 1 | _(undocumented)_ |
| `arabidopsis` | 1.7MB | 13 | 11 | Arabidopsis flowering-time (FT10) phenotype + accession data. |
| `clinvar` | 50.1MB | 14 | 17 | ClinVar variant-classification release. |
| `cyp2d6_psv` | 0.0MB | 2 | 2 | _(undocumented)_ |
| `disinfinder_db` | 0.0MB | 2 | 2 | _(undocumented)_ |
| `ena_wgs` | 121.6MB | 8 | 4 | _(undocumented)_ |
| `external` | 4.6MB | 87 | 191 | External-cohort validation inputs + crosswalks (Oxford, PROBAC, AR Bank, Sci234). |
| `fba` | 4.8MB | **0** | **0** | _(undocumented)_ |
| `forward_ref` | 0.0MB | 11 | 5 | _(undocumented)_ |
| `fungal_ref` | 0.0MB | 13 | 7 | _(undocumented)_ |
| `gono_catalogue` | 0.1MB | 1 | 1 | _(undocumented)_ |
| `hcmv_ref` | 0.0MB | 2 | 1 | _(undocumented)_ |
| `hiv_ref` | 0.0MB | 19 | 4 | Committed HXB2 CDS references (RT / PR / IN / CA). |
| `horesh_2021` | 4.6MB | 1 | **0** | Horesh H4 E. coli pathotype cohort (24 genomes). |
| `horse` | 0.0MB | 19 | 13 | _(undocumented)_ |
| `imputation` | 0.0MB | 21 | 26 | _(undocumented)_ |
| `isaba1_ref` | 0.0MB | 8 | 5 | _(undocumented)_ |
| `j3_abo` | 0.8MB | 3 | 3 | _(undocumented)_ |
| `ktype_db` | 0.3MB | 4 | 4 | _(undocumented)_ |
| `mlst_db` | 6.4MB | 2 | 3 | _(undocumented)_ |
| `nt_windows` | 32.9MB | 1 | 1 | _(undocumented)_ |
| `pathotype_cov_cache` | 0.0MB | 9 | 4 | _(undocumented)_ |
| `pathotype_pergene_cache` | 0.0MB | 9 | 8 | _(undocumented)_ |
| `pgx_1000g` | 76.6MB | 14 | 1 | _(undocumented)_ |
| `pgx_getrm` | 0.3MB | 7 | 6 | _(undocumented)_ |
| `phage_ref` | 26.5MB | 7 | 6 | _(undocumented)_ |
| `plasmidfinder_db` | 0.1MB | 5 | 1 | _(undocumented)_ |
| `pneumo_betalactam_db` | 0.0MB | 2 | 1 | _(undocumented)_ |
| `pneumoserotype_db` | 1.6MB | 5 | 2 | _(undocumented)_ |
| `pointfinder_db` | 0.0MB | 3 | 2 | _(undocumented)_ |
| `refseq_cache` | 5476.7MB | 12 | 6 | _(undocumented)_ |
| `resfinder_db` | 2.0MB | 6 | 1 | _(undocumented)_ |
| `salmserovar_db` | 0.6MB | 14 | 5 | _(undocumented)_ |
| `sarscov2_ref` | 0.0MB | 11 | 6 | _(undocumented)_ |
| `serotypefinder_db` | 1.4MB | 5 | 1 | _(undocumented)_ |
| `summary_stat_sources` | 1.6MB | 4 | 2 | _(undocumented)_ |
| `virulencefinder_db` | 7.4MB | 14 | 9 | VirulenceFinder allele DB (E. coli-scoped). |

## External caches (`D:/dna_decode_cache`) — heavy, gitignored, regenerable

| dataset | size | used by | reported in | what it is |
|---|---|---|---|---|
| `_dockertest` | 0.0MB | **0** | **0** | _(undocumented)_ |
| `alphafold` | 50.6MB | 9 | 2 | _(undocumented)_ |
| `alphamissense` | 1182.9MB | 14 | 12 | _(undocumented)_ |
| `arabmagic` | 1.7MB | 1 | **0** | _(undocumented)_ |
| `bakta_out` | 485.0MB | 2 | **0** | _(undocumented)_ |
| `bloom` | 45.7MB | 4 | 1 | _(undocumented)_ |
| `bxd` | 8.1MB | **0** | **0** | _(undocumented)_ |
| `caur_fks1_sra` | 20972.3MB | 1 | **0** | _(undocumented)_ |
| `caur_sra` | 13230.5MB | **0** | **0** | 13.2GB of C. auris SRA reads (reads/ + asm/ + assemblies/). |
| `celegans_ben1` | 0.2MB | 1 | 2 | _(undocumented)_ |
| `clinvar` | 367.4MB | 14 | 17 | ClinVar variant-classification release. |
| `cryptic` | 951.2MB> | 22 | 12 | CRyPTIC M. tuberculosis compendium (12,287 isolates, measured BMD MIC). |
| `darwins_ark` | 23751.5MB | 2 | 1 | _(undocumented)_ |
| `data files donwload` | 3371.8MB | **0** | **0** | 3.4GB directory, name misspelled (`donwload`), contents unclassified. |
| `depmap_pilot` | 1380.3MB | 1 | 1 | _(undocumented)_ |
| `dgrp` | 289.2MB | 1 | 1 | _(undocumented)_ |
| `dog_ref` | 1.3MB | 4 | 2 | _(undocumented)_ |
| `ecoli_sero_asm` | 2915.1MB | 5 | **0** | E. coli serotype cohort assemblies + per-isolate results.jsonl. |
| `ecoli_sero_heldout` | 4749.2MB | 2 | **0** | _(undocumented)_ |
| `embeddings` | 525.0MB | 40 | 103 | HDF5 foundation-model embedding caches (NT, mostly cipro cohorts). |
| `ena_wgs` | 90.9MB | 8 | 4 | _(undocumented)_ |
| `epistasis` | 62.9MB | 38 | 24 | _(undocumented)_ |
| `esm` | 17.4MB | **0** | **0** | _(undocumented)_ |
| `essentiality` | 15.3MB | 54 | 71 | Keio / Fitness Browser conditional-essentiality data. |
| `finngen` | 773.2MB | 1 | 2 | _(undocumented)_ |
| `fitness_browser` | 9210.4MB | 3 | 2 | _(undocumented)_ |
| `fungal_g1` | 15215.7MB | 4 | 4 | _(undocumented)_ |
| `gdsc` | 70.6MB | 1 | 1 | _(undocumented)_ |
| `gemme_work` | 190.3MB | 1 | **0** | _(undocumented)_ |
| `geuvadis` | 18587.5MB | 4 | 1 | _(undocumented)_ |
| `gmap_cli_smoke` | 10.7MB | **0** | **0** | _(undocumented)_ |
| `gnomad` | 1.1MB | 1 | **0** | _(undocumented)_ |
| `gono_asm` | 320.3MB | 1 | **0** | _(undocumented)_ |
| `gtex` | 0.0MB | **0** | **0** | _(undocumented)_ |
| `gwas_catalog` | 750.2MB | 1 | **0** | _(undocumented)_ |
| `hirisplex_src` | 8.0MB | **0** | 1 | _(undocumented)_ |
| `kleb` | 627.0MB | 18 | 7 | _(undocumented)_ |
| `ktype_build` | 0.5MB | 2 | 3 | _(undocumented)_ |
| `ktype_kaptive` | 0.0MB | 1 | **0** | Per-isolate ktype-vs-Kaptive concordance rows (307 genomes). |
| `liftover` | 12.8MB | 10 | 7 | _(undocumented)_ |
| `mavedb` | 0.1MB | 6 | 32 | _(undocumented)_ |
| `msa` | 58.7MB | **0** | **0** | _(undocumented)_ |
| `name_table_variantA_drop_ambiguous.json` | 0.0MB | **0** | **0** | _(undocumented)_ |
| `name_table_variantB.json` | 0.0MB | **0** | **0** | _(undocumented)_ |
| `nt-gate-b` | 3781.4MB | **0** | **0** | 3.8GB NT embedding cache for the gate-B cohort. |
| `opensnp` | 20385.9MB | 6 | 5 | _(undocumented)_ |
| `pan_ukbb` | 1.4MB | 1 | **0** | _(undocumented)_ |
| `pear` | 344.9MB | 5 | 7 | _(undocumented)_ |
| `pgp` | 216.9MB | **0** | **0** | _(undocumented)_ |
| `pgp_uk` | 432.4MB | 1 | 6 | _(undocumented)_ |
| `pneumo_amrfinder_swap` | 0.0MB | 1 | **0** | _(undocumented)_ |
| `pneumo_asm` | 19.0MB | 1 | **0** | _(undocumented)_ |
| `pneumo_gps` | 54.4MB | 5 | **0** | _(undocumented)_ |
| `pneumocat_src` | 647.4MB | **0** | **0** | _(undocumented)_ |
| `pointfinder_concordance` | 0.2MB | 1 | **0** | _(undocumented)_ |
| `precise1k` | 58.7MB | 1 | 8 | _(undocumented)_ |
| `prospective` | 0.0MB | 27 | 50 | Prospective-lock cohort sweeps + accrual state. |
| `prospective_smoke` | 0.0MB | **0** | **0** | _(undocumented)_ |
| `prosst_structures` | 1059.6MB | 1 | 1 | _(undocumented)_ |
| `proteingym` | 8506.3MB | 27 | 12 | ProteinGym DMS substitution benchmark + per-assay reference scores. |
| `refs` | 1111.5MB | 19 | 12 | _(undocumented)_ |
| `refseq` | 18493.8MB> | 56 | 33 | NCBI RefSeq genome cache (FASTA + GFF3 + GBK per accession). |
| `resfinder_collapse` | 1.9MB | 1 | **0** | _(undocumented)_ |
| `salm_asm` | 3599.8MB | 9 | 1 | _(undocumented)_ |
| `salm_test` | 4.8MB | **0** | **0** | _(undocumented)_ |
| `salmonella_antigens.pre_role.fasta` | 0.5MB | **0** | **0** | _(undocumented)_ |
| `salmserovar_cov80_probe.json` | 0.0MB | **0** | 1 | _(undocumented)_ |
| `salmserovar_o_fix_result_run1.json` | 0.0MB | **0** | **0** | _(undocumented)_ |
| `seqsero2_src` | 52.9MB | **0** | **0** | _(undocumented)_ |
| `serovar_table.pre_rename.tsv` | 0.1MB | **0** | **0** | _(undocumented)_ |
| `serovar_table.variantA.tsv` | 0.1MB | **0** | **0** | _(undocumented)_ |
| `spn_pbp_amr` | 443.4MB | **0** | **0** | _(undocumented)_ |
| `ss2` | 4.8MB | **0** | **0** | _(undocumented)_ |
| `ss2_antigens.fasta` | 0.5MB | **0** | 1 | _(undocumented)_ |
| `ss2_checkpoint` | 0.2MB | 4 | 3 | Per-isolate SeqSero2 + salmserovar policy-variant checkpoints (200 isolates). |
| `ss2_ic_dir` | 0.2MB | **0** | 1 | _(undocumented)_ |
| `ss2_initial_conditions.py` | 0.2MB | **0** | **0** | _(undocumented)_ |
| `ss2_run` | 5.0MB | 1 | **0** | _(undocumented)_ |
| `ss2_source.py` | 0.1MB | **0** | 1 | _(undocumented)_ |
| `stage2_kaggle` | 0.9MB | 2 | **0** | _(undocumented)_ |
| `tb_callability` | 0.1MB | 7 | 2 | _(undocumented)_ |
| `tb_indep` | 12496.3MB> | 10 | 12 | Per-isolate TB independent-run VCFs + results checkpoint (2,845 isolates). |
| `torch` | 2483.9MB | 58 | 32 | _(undocumented)_ |

## Datasets nothing references

Either dead weight or a forgotten asset. This is the bucket rediscovery comes from.

- `_dockertest` (cache, 0.0MB) — `D:\dna_decode_cache\_dockertest`
- `bxd` (cache, 8.1MB) — `D:\dna_decode_cache\bxd`
- `caur_sra` (cache, 13230.5MB) — `D:\dna_decode_cache\caur_sra`
- `data files donwload` (cache, 3371.8MB) — `D:\dna_decode_cache\data files donwload`
- `esm` (cache, 17.4MB) — `D:\dna_decode_cache\esm`
- `gmap_cli_smoke` (cache, 10.7MB) — `D:\dna_decode_cache\gmap_cli_smoke`
- `gtex` (cache, 0.0MB) — `D:\dna_decode_cache\gtex`
- `msa` (cache, 58.7MB) — `D:\dna_decode_cache\msa`
- `name_table_variantA_drop_ambiguous.json` (cache, 0.0MB) — `D:\dna_decode_cache\name_table_variantA_drop_ambiguous.json`
- `name_table_variantB.json` (cache, 0.0MB) — `D:\dna_decode_cache\name_table_variantB.json`
- `nt-gate-b` (cache, 3781.4MB) — `D:\dna_decode_cache\nt-gate-b`
- `pgp` (cache, 216.9MB) — `D:\dna_decode_cache\pgp`
- `pneumocat_src` (cache, 647.4MB) — `D:\dna_decode_cache\pneumocat_src`
- `prospective_smoke` (cache, 0.0MB) — `D:\dna_decode_cache\prospective_smoke`
- `salm_test` (cache, 4.8MB) — `D:\dna_decode_cache\salm_test`
- `salmonella_antigens.pre_role.fasta` (cache, 0.5MB) — `D:\dna_decode_cache\salmonella_antigens.pre_role.fasta`
- `salmserovar_o_fix_result_run1.json` (cache, 0.0MB) — `D:\dna_decode_cache\salmserovar_o_fix_result_run1.json`
- `seqsero2_src` (cache, 52.9MB) — `D:\dna_decode_cache\seqsero2_src`
- `serovar_table.pre_rename.tsv` (cache, 0.1MB) — `D:\dna_decode_cache\serovar_table.pre_rename.tsv`
- `serovar_table.variantA.tsv` (cache, 0.1MB) — `D:\dna_decode_cache\serovar_table.variantA.tsv`
- `spn_pbp_amr` (cache, 443.4MB) — `D:\dna_decode_cache\spn_pbp_amr`
- `ss2` (cache, 4.8MB) — `D:\dna_decode_cache\ss2`
- `ss2_initial_conditions.py` (cache, 0.2MB) — `D:\dna_decode_cache\ss2_initial_conditions.py`
- `mini_cipro_dnabert2_cache.h5` (processed, 0.0MB) — `C:\Users\Farshad\PythonProjects\dna_decode\data\processed\mini_cipro_dnabert2_cache.h5`
- `mini_cipro_mock_cache.h5` (processed, 52.0MB) — `C:\Users\Farshad\PythonProjects\dna_decode\data\processed\mini_cipro_mock_cache.h5`
- `msat_lift_run20.log` (processed, 0.0MB) — `C:\Users\Farshad\PythonProjects\dna_decode\data\processed\msat_lift_run20.log`
- `pilot_report_phase2_entry_2026-05-13.md` (processed, 0.0MB) — `C:\Users\Farshad\PythonProjects\dna_decode\data\processed\pilot_report_phase2_entry_2026-05-13.md`
- `shared_lineage_ceftriaxone_dl.parquet` (processed, 0.0MB) — `C:\Users\Farshad\PythonProjects\dna_decode\data\processed\shared_lineage_ceftriaxone_dl.parquet`
- `shared_lineage_gentamicin_dl.parquet` (processed, 0.0MB) — `C:\Users\Farshad\PythonProjects\dna_decode\data\processed\shared_lineage_gentamicin_dl.parquet`
- `shared_lineage_meropenem_dl.parquet` (processed, 0.0MB) — `C:\Users\Farshad\PythonProjects\dna_decode\data\processed\shared_lineage_meropenem_dl.parquet`
- `shared_lineage_tetracycline_pilot.parquet` (processed, 0.0MB) — `C:\Users\Farshad\PythonProjects\dna_decode\data\processed\shared_lineage_tetracycline_pilot.parquet`
- `acinetobacter_meropenem_run.log` (raw, 0.0MB) — `C:\Users\Farshad\PythonProjects\dna_decode\data\raw\acinetobacter_meropenem_run.log`
- `ar_bank_caur_extval_voriconazole` (raw, 0.0MB) — `C:\Users\Farshad\PythonProjects\dna_decode\data\raw\ar_bank_caur_extval_voriconazole`
- `ar_bank_ecoli_extval_ciprofloxacin` (raw, 22.2MB) — `C:\Users\Farshad\PythonProjects\dna_decode\data\raw\ar_bank_ecoli_extval_ciprofloxacin`
- `ar_bank_ecoli_extval_gentamicin` (raw, 41.5MB) — `C:\Users\Farshad\PythonProjects\dna_decode\data\raw\ar_bank_ecoli_extval_gentamicin`
- `ar_bank_gono_extval_azithromycin` (raw, 0.0MB) — `C:\Users\Farshad\PythonProjects\dna_decode\data\raw\ar_bank_gono_extval_azithromycin`
- `ar_bank_gono_extval_cefixime` (raw, 0.0MB) — `C:\Users\Farshad\PythonProjects\dna_decode\data\raw\ar_bank_gono_extval_cefixime`
- `ar_bank_gono_extval_ceftriaxone` (raw, 0.0MB) — `C:\Users\Farshad\PythonProjects\dna_decode\data\raw\ar_bank_gono_extval_ceftriaxone`
- `ar_bank_gono_extval_ciprofloxacin` (raw, 205.8MB) — `C:\Users\Farshad\PythonProjects\dna_decode\data\raw\ar_bank_gono_extval_ciprofloxacin`
- `ar_bank_gono_extval_penicillin` (raw, 0.0MB) — `C:\Users\Farshad\PythonProjects\dna_decode\data\raw\ar_bank_gono_extval_penicillin`
- `ar_bank_gono_extval_tetracycline` (raw, 0.0MB) — `C:\Users\Farshad\PythonProjects\dna_decode\data\raw\ar_bank_gono_extval_tetracycline`
- `ar_bank_gono_kaggle_subset.json` (raw, 0.0MB) — `C:\Users\Farshad\PythonProjects\dna_decode\data\raw\ar_bank_gono_kaggle_subset.json`
- `ar_bank_kleb_extval_ceftriaxone` (raw, 723.8MB) — `C:\Users\Farshad\PythonProjects\dna_decode\data\raw\ar_bank_kleb_extval_ceftriaxone`
- `ar_bank_kleb_extval_ciprofloxacin` (raw, 607.2MB) — `C:\Users\Farshad\PythonProjects\dna_decode\data\raw\ar_bank_kleb_extval_ciprofloxacin`
- `ar_bank_kleb_extval_gentamicin` (raw, 204.6MB) — `C:\Users\Farshad\PythonProjects\dna_decode\data\raw\ar_bank_kleb_extval_gentamicin`
- `ar_bank_staphylococcus_aureus_extval_levofloxacin` (raw, 117.6MB) — `C:\Users\Farshad\PythonProjects\dna_decode\data\raw\ar_bank_staphylococcus_aureus_extval_levofloxacin`
- `campylobacter_cipro_run.log` (raw, 0.0MB) — `C:\Users\Farshad\PythonProjects\dna_decode\data\raw\campylobacter_cipro_run.log`
- `enterobacter_cloacae_cef_run.log` (raw, 0.0MB) — `C:\Users\Farshad\PythonProjects\dna_decode\data\raw\enterobacter_cloacae_cef_run.log`
- `enterobacter_cloacae_ceftriaxone` (raw, 210.1MB) — `C:\Users\Farshad\PythonProjects\dna_decode\data\raw\enterobacter_cloacae_ceftriaxone`
- `hiv` (raw, 19.6MB) — `C:\Users\Farshad\PythonProjects\dna_decode\data\raw\hiv`
- `independent_cohort_validate_run.log` (raw, 0.0MB) — `C:\Users\Farshad\PythonProjects\dna_decode\data\raw\independent_cohort_validate_run.log`
- `klebsiella_ceftriaxone` (raw, 1.5MB) — `C:\Users\Farshad\PythonProjects\dna_decode\data\raw\klebsiella_ceftriaxone`
- `klebsiella_ciprofloxacin` (raw, 0.0MB) — `C:\Users\Farshad\PythonProjects\dna_decode\data\raw\klebsiella_ciprofloxacin`
- `klebsiella_tetracycline` (raw, 0.6MB) — `C:\Users\Farshad\PythonProjects\dna_decode\data\raw\klebsiella_tetracycline`
- `pseudomonas_meropenem_run.log` (raw, 0.0MB) — `C:\Users\Farshad\PythonProjects\dna_decode\data\raw\pseudomonas_meropenem_run.log`
- `salmonella_cipro_run.log` (raw, 0.0MB) — `C:\Users\Farshad\PythonProjects\dna_decode\data\raw\salmonella_cipro_run.log`
- `salmonella_indep_ciprofloxacin` (raw, 546.9MB) — `C:\Users\Farshad\PythonProjects\dna_decode\data\raw\salmonella_indep_ciprofloxacin`
- `staphylococcus_aureus_oxacillin` (raw, 292.6MB) — `C:\Users\Farshad\PythonProjects\dna_decode\data\raw\staphylococcus_aureus_oxacillin`
- `fba` (ref, 4.8MB) — `C:\Users\Farshad\PythonProjects\dna_decode\data\fba`
