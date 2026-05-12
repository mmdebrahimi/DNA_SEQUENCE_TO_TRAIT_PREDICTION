# Architecture

Phase 1 E. coli G2P platform. One-page module map + data-flow overview.

## Layered package structure

```
                              ┌─────────────────────────────────────┐
                              │  config/datasources.yaml            │
                              │  (declarative data-source registry) │
                              └─────────────┬───────────────────────┘
                                            │
                  ┌─────────────────────────┼─────────────────────────┐
                  │                         │                         │
            INGESTION                  FOUNDATION                  EVAL
        dna_decode/data/         dna_decode/models/         dna_decode/eval/
                  │                         │                         │
   ┌──────────────┴──────────────┐ ┌────────┴─────────┐  ┌───────────┴──────────┐
   │ pilot.py    (Step 0.5)      │ │ foundation.py    │  │ cv.py (Step 10)      │
   │ refseq.py   (Step 2)        │ │ (Step 7)         │  │  LOSO + LOMO +       │
   │ annotations.py (Step 3)     │ │  Evo + DNABERT-2 │  │  leave-one-Mash-     │
   │ resistance_db.py (Step 4)   │ │  + NT + GENA-LM  │  │  clade-out           │
   │ ast_data.py (Step 5)        │ │  + MockModel     │  │ metrics.py           │
   │ cohort.py (Step 6)          │ │ cache.py (Step 8)│  │  AUROC/AUPRC/Brier/  │
   │  + candidates_from_bvbrc_ast│ │  HDF5 cache;     │  │  ECE + attribution   │
   │  + save_cohort/load_cohort  │ │  populate() opens│  │  precision           │
   │  + assembly_accession field │ │  HDF5 once       │  │ phylogeny.py         │
   │                             │ │ classifiers.py   │  │  Mash + cluster_by_  │
   │                             │ │  (Step 9)        │  │  ani                 │
   │                             │ │  XGBoost +       │  │ clade_baseline.py    │
   │                             │ │  sigmoid calib   │  │  + validation_gate   │
   │                             │ │ classical_       │  └──────────────────────┘
   │                             │ │  baselines.py    │
   │                             │ │  (Step 18)       │
   │                             │ │  AMRFinder/k-mer │
   │                             │ │  /gene-presence  │
   └─────────────────────────────┘ └──────────────────┘
                                            │
                                            │   gene-level
                                            │   embeddings
                                            ▼
                              ┌─────────────────────────────────────┐
                              │  INTERPRETABILITY                   │
                              │  dna_decode/interp/                 │
                              │  ┌───────────────────────────────┐  │
                              │  │ mutagenesis.py (Step 11)      │  │
                              │  │  - gene_level_mutagenesis     │  │
                              │  │  - saturation_mutagenesis     │  │
                              │  │  - tier_classify (Tier 1-5)   │  │
                              │  │  - _best_tier_across_         │  │
                              │  │    candidates                 │  │
                              │  │  - motif_recovery (Phase 2    │  │
                              │  │    placeholder)               │  │
                              │  └───────────────────────────────┘  │
                              └─────────────────────────────────────┘
                                            │
                                            ▼
                              ┌─────────────────────────────────────┐
                              │  VIZ                                │
                              │  dna_decode/viz/                    │
                              │   browser.py (Step 13)              │
                              │   - render_attribution_plot         │
                              │     (matplotlib PNG)                │
                              │   - export_attribution_tsv          │
                              │  (pygenometracks → Phase 2)         │
                              └─────────────────────────────────────┘
```

## CLI entry points

```
scripts/
├── pilot_gate.py            Step 0.5 — HARD gate before ingestion fires
├── smoke_pipeline.py        Step 15 — synthetic-fixture end-to-end (<60s)
├── pipeline.py              Step 14 — single CLI with 4 subcommands
│                              {ingest, train, predict, attribute}
├── leaderboard.py           Step 17 — shell-loop over pipeline.py train
│                              for (model × drug); writes leaderboard.md
└── quantize_fidelity_check.py  Step 11.5 — 4-bit vs full-precision ISM
                              concordance (one-time validation)
```

## End-to-end data flow

```
1. INGEST       pilot_gate → load_bvbrc_ast → candidates_from_bvbrc_ast
                → build_cohort → save_cohort (parquet) → download_cohort_genomes
                → refseq.download_genome × N (uses assembly_accession)
                → unpacks NCBI ZIP into {genome.fna, annotations.gff3, annotations.gbk}

2. EMBED        EmbeddingCache.populate(model, strain_genomes, annotations)
                → extract_cds_sequences (Step 3) per strain
                → model.embed_batch(...) per gene
                → HDF5 write at /strains/<strain_id>/<gene_id>  (opens file once)

3. TRAIN        load_cohort → cache.bulk_get([(sid, g), ...])
                → aggregate_strain_features (mean-pool)
                → train_xgboost_classifier (sigmoid calibration; minority-class CV)
                → leave_one_strain_out_cv → compute_metrics
                → (optional) train_clade_only_classifier → validation_gate
                → pickle bundle to data/processed/models/<drug>_<model>.pkl

4. ATTRIBUTE    load classifier + cache + annotations + resistance_catalog
                → gene_level_mutagenesis (ISM by embedding-row drop)
                → build_attribution_report (Tier 1-5 via _best_tier_across_candidates)
                → JSON output with per-locus tiers + fractions

5. VIZ          render_attribution_plot → PNG (bar of top-K genes + position panel)
                export_attribution_tsv → TSVs of GeneEffectTable + PositionEffectTable

6. LEADERBOARD  for model in [evo, dnabert2, ...]:
                  for drug in [cipro, ceftriaxone, tet]:
                    subprocess.run(["python", "-m", "scripts.pipeline", "train",
                                    "--drug", drug, "--model", model, ...])
                → aggregate per-model AUROC + clade-gap + tier fractions
                → write data/processed/leaderboard.md

7. (OPTIONAL)   quantize_fidelity_check: compare 4-bit vs full-precision ISM
                attribution tables; top-K intersection + Spearman; GO/NO-GO
```

## Key contracts that must not break

| Contract | Where | Why |
|---|---|---|
| `CandidateStrain.assembly_accession` (NCBI) ≠ `strain_id` (BV-BRC) | `cohort.py` | NCBI Datasets API expects accession; BV-BRC strain IDs would 404. |
| `EmbeddingCache.populate` opens HDF5 once | `cache.py` | Was 2.25M file-opens before Wave 2.5 fix; performance footgun. |
| `build_attribution_report` walks BOTH `gene_id` and `locus_tag` | `mutagenesis.py` | Bakta locus tags never match catalog gene symbols; Wave 3.5 C8 fix. |
| Calibration CV folds bounded by MINORITY class | `classifiers.py` + `classical_baselines.py` | Rare-resistance drugs (low minority class count) crashed before Wave 3.5 C7 fix. |
| `motif_recovery` is a Phase 2 placeholder + warns on call | `mutagenesis.py` | Phase 1 attribution uses Tier 1-5 rubric, not motif recovery. |
| 6 HDF5 end markers in cohort parquet | `cohort.py` save/load round-trip | Wave 3 reads cohort.parquet directly; recomputing build_cohort on every load would be wasteful. |

## Phase 1 success criteria

See `README.md` and `plans/Ecoli_G2P_Platform_Technical_Plan.md` Verification section.

## Future-state pointers

- `TODOS.md` — deferred work + post-Phase-1 polish
- `FUTURE_FEATURES.md` — Phase 2-4+ ideas
- `LESSONS_LEARNED.md` — practical lessons from Phase 1 build
- `plans/Ecoli_G2P_Platform_Technical_Plan.md` Phase 2 Backlog — Attribution Refinement Engine, Captum + diff MLP head, MIC regression, pan-genome graph layer, AlphaFold-inspired arch (Phase 3+ with compute caveat)
