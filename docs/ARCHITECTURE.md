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
   │ mic_tiers.py (2026-05-18)   │ │                  │  │                      │
   │  shared per-drug catalogs:  │ │                  │  │                      │
   │  CLSI/EUCAST breakpoints,   │ │                  │  │                      │
   │  AMRFinder Class filter,    │ │                  │  │                      │
   │  loci_by_mechanism,         │ │                  │  │                      │
   │  classify_tier (HIGH_R/...) │ │                  │  │                      │
   │  classify_gene_symbol +     │ │                  │  │                      │
   │  tolerant prefix match.     │ │                  │  │                      │
   │  Drugs: cipro/cef/tet/gent. │ │                  │  │                      │
   │                             │ │ classifiers.py   │  │  Mash + cluster_by_  │
   │                             │ │  (Step 9)        │  │  ani                 │
   │                             │ │  XGBoost +       │  │ clade_baseline.py    │
   │                             │ │  sigmoid calib   │  │  + validation_gate   │
   │                             │ │ classical_       │  │ loso_kmer.py         │
   │                             │ │  baselines.py    │  │  order-explicit      │
   │                             │ │  (Step 18)       │  │  k-mer + fusion      │
   │                             │ │  AMRFinder/k-mer │  │  LOSO; shared by     │
   │                             │ │  /gene-presence  │  │  smoke gate +        │
   │                             │ │  CONTIG_SEPARATOR│  │  Stage 1 runner      │
   │                             │ │  module constant │  └──────────────────────┘
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
├── quantize_fidelity_check.py  Step 11.5 — 4-bit vs full-precision ISM
│                              concordance (one-time validation)
├── smoke_gate_12strain_cipro.py  EP-framework smoke gate (12-strain mini cohort;
│                              3 variants: NT-XGBoost + k-mer + gene-presence
│                              with INDETERMINATE_IDENTIFIER_OOV guardrail).
│                              EP2 generalization (rename + drug-templated output)
│                              tracked in plans/EP2_Cef_Tet_Smoke_Design_Plan.md.
├── stage1_n40_cipro.py      Phase 2 Stage 1 engineering screen (N=40 cipro; 4
│                              variants gated by max(NT-XGBoost, NT-logreg) − k-mer
│                              ≥ 3 pp; emits verdict + stage2_action separately)
├── probe_nt_cache.py        Cache integrity gate (added 2026-05-15). Thin CLI
│                              wrapper over EmbeddingCache.verify_complete. Run
│                              BEFORE any Stage 1 / EP1 invocation on a populate
│                              that may have been interrupted. Exit 0 = ALL_COMPLETE.
├── build_stage2_n150_cohort.py  Phase 2 Stage 2 cohort builder. Label-stratified
│                              MLST-balanced selection (per-class R/S). Hard-fails
│                              on imbalance + missing MLST. Builds 147-strain cohort
│                              from BV-BRC AST + genome metadata CSV.
├── diagnose_bvbrc_mlst_gaps.py  Layer-by-layer cohort-shortfall diagnostic. Surfaces
│                              whether the bottleneck is MLST, assembly_accession,
│                              AST coverage, or assembly-quality filters.
├── diagnose_gene_presence_auroc.py  AUROC=0.000 root-cause diagnostic (one-time).
├── diagnose_gene_presence_synthetic.py  Synthetic falsifier for the gene-presence
│                              AUROC=0.000 LOSO-base-rate-inversion hypothesis.
├── plot_nt_embeddings_pca_umap.py  Diagnostic 2D projection of NT embeddings;
│                              MLST overlay for lineage-confounding inspection.
├── cipro_attribution_preflight.py  Cipro gene-level ISM attribution audit on
│                              the N=38 NT cache (v1 mean-pool 2026-05-15; v2
│                              expanded loci + signed-positive-delta + frequency
│                              aggregation 2026-05-16). Verdict: STRONG_POSITIVE /
│                              WEAK_POSITIVE / INCONCLUSIVE_MISS. Both v1 and v2
│                              returned INCONCLUSIVE_MISS on cipro mean-pool.
│                              Mean+max preflight v3 refactor deferred (per
│                              plans/Cipro_Decision_Bundle_Technical_Plan.md
│                              Step 2 + Step 4 — closeout falsifier, not
│                              Phase 2 entry-grade).
├── cipro_mic_audit.py        Cipro raw BV-BRC AST/MIC rejoin audit (shipped
│                              2026-05-17). Tiers each strain under CLSI 2024
│                              + EUCAST 14.0 breakpoints; classifies HIGH_R /
│                              HIGH_S / DECISIVE / BORDERLINE / AMBIGUOUS /
│                              CONFLICT / NO_MIC. Output: per-strain mechanism-
│                              class table + cohort signal-quality verdict.
│                              N=38 cipro cohort verdict 2026-05-17: NOISY
│                              (7 HIGH_R / 0 HIGH_S of 40; 9 R have no MIC;
│                              12 S borderline).
├── cipro_mechanism_audit.py  Cipro AMRFinderPlus mechanism audit (shipped
│                              2026-05-17). Per-strain QRDR / plasmid /
│                              regulatory / efflux / porin hit detection
│                              across all R + S strains. Class filter:
│                              CIPRO_RELEVANT_AMR_CLASSES = QUINOLONE_CLASSES |
│                              {MULTIDRUG} — keeps regulatory frameshifts
│                              (marR_V84WfsTer, acrR_S30HfsTer) that come
│                              through with AMRFinder Class=MULTIDRUG. Pinned
│                              synonymous-SNP filter + main.tsv/mutations.tsv
│                              dedupe. N=38 cipro verdict 2026-05-17:
│                              QRDR_DOMINANT (18/20 R have QRDR; 7/20 plasmid;
│                              7/20 S strains carry silent primary mechanism).
├── cipro_mechanism_phenotype_merge.py  Cipro mechanism × MIC merge (shipped
│                              2026-05-17). Joins mechanism audit + MIC audit
│                              into per-strain noise_class table. Strict primary
│                              cipro mechanism = QRDR_target_alteration OR
│                              plasmid_protect_modify only; efflux/regulatory/
│                              porin reported as separate co_resistance_modifiers
│                              column (don't drive noise classification).
│                              mechanism_opacity_flag = True iff HIGH_R + no
│                              primary mechanism (distinguishes tool-incomplete
│                              from labels-noisy). Pre-curated-baseline gate:
│                              RUN_FULL_AND_CLEAN (clean_count >= 20) /
│                              RUN_FULL_ONLY (>= 10) / MECHANISM_DEBUG_BRANCH
│                              (low clean + high opacity) / SUSPEND_CONDITION_4
│                              (low clean + low opacity). N=38 cipro 2026-05-17:
│                              SUSPEND_CONDITION_4 (signal quality 0.17,
│                              opacity_count = 0).
├── cipro_curated_baseline.py  Cipro curated AMR baseline (shipped 2026-05-17;
│                              PIVOT TRIGGER condition 4 test). LR + XGB over
│                              named multi-block feature sets: all / no_POINT /
│                              mechanism_only / POINT_only / kmer_only /
│                              MLST_only / kmer_MLST_only. LOSO. 2-layer verdict:
│                              original_condition_4 (frozen pre-Experiment-2
│                              rule, all-feature AUROC) + amended_condition_4
│                              (load-bearing: no_POINT >= 0.773 OR
│                              mechanism_only >= 0.80). The amended rule isolates
│                              non-textbook genomic signal vs labels-in-genome-
│                              form tautology. Refuses to fire when merge gate's
│                              clean_count < 10. given_suspended_gate field added
│                              for downstream readers.
├── cipro_error_audit.py      Cipro per-strain error audit (placeholder per
│                              Cipro_Decision_Bundle_Technical_Plan Step 10.5;
│                              deferred). Loads Stage 1b NT-LR per-strain
│                              predictions JSON sidecar (when Step 8 ships it);
│                              joins to manifest's noise_class; Fisher exact
│                              label-stratified enrichment test. Pre-conditions
│                              artifact PC1+PC2 required before runtime.
├── bvbrc_strict_mic_4drug_census.py  BV-BRC strict-MIC 4-drug feasibility
│                              census (shipped 2026-05-18; Phase 2 entry).
│                              6-stage pipeline (AST rows -> distinct genomes
│                              -> with MIC -> classification pass -> with
│                              assembly_accession -> passing QC) at TWO label-
│                              quality bars: strict-MIC (HIGH only, 4x safety
│                              margin) + relaxed-MIC (HIGH + DECISIVE). Drugs:
│                              cipro/cef/tet/gent. Imports from mic_tiers.
│                              Headline 2026-05-18: NO drug clears N=150
│                              per-class at either bar; assembly_accession is
│                              the structural bottleneck. Output: wiki/
│                              bvbrc_strict_mic_4drug_census_<date>.{md,json}.
└── run_mechanism_audit_detached.bat  Windows detached batch wrapper for the
                              ~1-hour mechanism audit. Avoids Bash 10-min hard
                              timeout cap. Pattern reusable for any long-running
                              detached job (cf. run_stage1b_detached.bat).

tools/
└── docker_runner.py         Stage 2 bioinformatics-tool runner via Docker Desktop
                              (added 2026-05-15). Single generic run(image, args,
                              mounts, env, capture_output, check, timeout) →
                              CompletedProcess; wraps FileNotFoundError +
                              TimeoutExpired as DockerRunnerError. Replaces
                              .sh wrappers that don't work with Python subprocess
                              on Windows. Routes Mash + AMRFinderPlus + Bakta.
                              9 unit tests at tests/test_docker_runner.py.
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
| `CVResult.strain_ids` is the alignment contract for paired comparisons | `eval/cv.py` | Stage 1 runner's `compute_gate_outcome` validates NT-vs-k-mer strain_ids alignment before computing the gap; mismatch raises (gate-bearing). Fusion strain_ids mismatch suppresses the fusion-outperforms note (diagnostic-only, gate proceeds). |
| `loso_kmer.run_kmer_xgboost_loso` / `run_fusion_loso` respect caller-supplied `strain_ids` verbatim — no internal sort/filter | `eval/loso_kmer.py` | Smoke gate + Stage 1 runner both pass through; preventing alignment drift was the /brainstorm Round 1 finding. Re-raises `ClassifierTrainingError` rather than silent mean-fallback. |
| `CONTIG_SEPARATOR = "N" * 100` module constant in `classical_baselines.py` | `models/classical_baselines.py` | Replaces three magic-string copies; both `loso_kmer` runners import + concatenate via this constant. |
| `_train_baseline_logreg(..., calibrate=False)` returns raw `LogisticRegression` without `CalibratedClassifierCV` | `models/classical_baselines.py` | Calibration at LOSO N≤20 over-corrects to anti-predictive output. Stage 1 + smoke gate pass `calibrate=False` at every gate-bearing call site (pinned by test). |
| `parse_gff3` two-pass parent→CDS gene_symbol propagation | `data/annotations.py` | Bakta-style GFF3 puts `gene=` on parent gene rows; CDS rows link via `Parent=`. Without propagation, gene-presence comparator re-enters INDETERMINATE_IDENTIFIER_OOV after Bakta re-annotation. Same-row `gene=` always wins to preserve RefSeq behavior. |
| `leave_one_*_out_cv` raises on unassigned strain_ids by default | `eval/cv.py` | The legacy silent `__unassigned__` bucketing would warp Stage 2's Mash-clade-out CV (one giant fold for all unassigned). Pass `allow_unassigned=True` only after explicit audit. |
| `build_stage2_n150_cohort` uses per-class label-stratified MLST-balanced selection | `scripts/build_stage2_n150_cohort.py` | Default `build_cohort` selects diversity-only across the merged R+S pool, leaving available R strains on the table at small ceilings (49R/101S → 72R/75S after per-class selection on the same input). |
| BV-BRC cohort R-ceiling bottleneck is `assembly_accession`, NOT MLST | `data/bvbrc_genome.py` + `scripts/diagnose_bvbrc_mlst_gaps.py` | 35,790 of 85,114 raw BV-BRC genome rows lack downloadable NCBI accession. When a cohort comes up short on a class, profile the join layers BEFORE relaxing filters — relaxing MLST would not have helped. |
| `EmbeddingCache.verify_complete` is the consumer-side integrity gate; mean-pool consumers MUST bail unless `report.all_complete` is True | `models/cache.py` + `scripts/probe_nt_cache.py` | `populate(skip_existing=True)` skips at gene-dataset level, not strain level. Stage 1's `cache.list_genes()` + mean-pool admits a strain on ≥1 cached gene → half-flushed strain at crash time becomes a silent landmine. 8 regression tests pin the 4 status buckets (complete / partial / absent / corrupt) + the `all_complete=False` rule. |
| `compute_mash_distances` issues exactly 2 mash invocations (1 sketch + 1 dist), NOT N*(N-1)/2 | `eval/phylogeny.py` | Batched-call discipline (2026-05-15 refactor). At N=147 the prior nested-loop pattern would issue 10,731 invocations; the batched pattern issues 2. Pinned by a call-count regression test in `tests/test_eval_phylogeny.py`. New `use_docker=True` kwarg routes through `tools/docker_runner.run`. |
| `tools/docker_runner.run` wraps `subprocess.run([docker, run, ...])` in argv form; never shell strings | `tools/docker_runner.py` | Avoids Git Bash MSYS path-conversion bug that silently breaks `-v <host>:/<container>` mounts. Python subprocess bypasses shell path-munging by construction. Wraps `FileNotFoundError` (docker not on PATH) + `TimeoutExpired` as `DockerRunnerError`. |
| `cipro_mechanism_audit.py` AMRFinder Class filter keeps MULTIDRUG, not just QUINOLONE | `scripts/cipro_mechanism_audit.py` | Regulatory mutations like `marR_V84WfsTer` + `acrR_S30HfsTer` come through with AMRFinder Class=MULTIDRUG (NOT QUINOLONE). Filtering on QUINOLONE-only would silently drop real cipro-affecting regulatory frameshifts. Same discipline applies to any future cef / tet mechanism audit (cef-relevant: BETA-LACTAM + CARBAPENEM + CEPHALOSPORIN + MULTIDRUG; tet-relevant: TETRACYCLINE + MULTIDRUG). |
| `cipro_mechanism_phenotype_merge.py` strict primary cipro mechanism = QRDR_target_alteration OR plasmid_protect_modify ONLY | `scripts/cipro_mechanism_phenotype_merge.py` | Efflux + regulatory + porin_loss are co-resistance modifiers, not standalone cipro-conferring mechanisms. The merge's `noise_class` is driven by strict primary; co-resistance reported separately. mechanism_opacity_flag separates "AMRFinder is incomplete" from "label is noisy" — both can be true at once, and conflating them loses the remediation signal. |
| `cipro_curated_baseline.py` 2-layer verdict: original_condition_4 (frozen all-feature) + amended_condition_4 (no_POINT >= 0.773 OR mechanism_only >= 0.80) | `scripts/cipro_curated_baseline.py` | The all-feature curated baseline is structurally circular when POINT mutations (gyrA_S83L etc.) dominate — they're essentially labels-in-genome-form. The amended verdict isolates non-textbook genomic signal vs textbook-tautology by gating on the no-POINT and mechanism-only ablation runs. Original verdict preserved for audit-trail discipline only. AMENDED_NO_POINT_GATE_AUROC = max(0.75, STAGE1B_NT_LR_AUROC + 0.10) = 0.773 (consistent with "beat NT by 10pp" framing). |
| Smoke runner output strings templated on `--drug` (2026-05-17) | `scripts/smoke_gate_12strain_cipro.py` | `# Smoke Gate — 12-strain <drug> cohort` heading + `wiki/smoke_gate_12strain_<drug_slug>_<date>.md` output path. Cef + tet smokes 2026-05-17 used this. Script filename rename (smoke_gate_12strain_cipro.py → smoke_gate_12strain.py) deferred as cosmetic tech debt. NT-XGBoost runner falls back to `ast_labels` iteration when `cohort.per_drug_strain_ids[drug]` missing — lets mini cohorts built outside `build_cohort()` (e.g., post-hoc per-drug filters from the cipro N=38 cohort) drive the smoke runner without re-populate. |
| `dna_decode/data/mic_tiers.py` is the single source of truth for per-drug breakpoints + mechanism catalogs; new drug audits MUST use it (not re-hardcode) | `dna_decode/data/mic_tiers.py` + callers | Added 2026-05-18 to prevent drift between `cipro_mic_audit.py` and the 4-drug census. Per-drug data: `DRUG_BREAKPOINTS` (CLSI 2024 + EUCAST 14.0 cipro/cef/tet/gent), `DRUG_AMRFINDER_CLASSES`, `DRUG_LOCI_BY_MECHANISM`, `DRUG_PRIMARY_MECHANISMS`. Helpers: `breakpoints_for(drug)`, `classify_tier(mics, distinct_calls, breakpoints)`, `amrfinder_classes_for(drug)`, `classify_gene_symbol(drug, symbol)` with tolerant prefix-match. `bvbrc_strict_mic_4drug_census.py` uses it; cipro_*.py scripts left as cipro-specific by design. Future drug audits (cef_mechanism_audit, tet_mechanism_audit, gent_mechanism_audit) MUST import from mic_tiers. 76 unit tests at `tests/test_mic_tiers.py`. |
| `pipeline.py predict` emits the v0 JSON + markdown sidecar schema; honest-output discipline (audit_verdict propagates SUSPEND framing, no overclaiming) is a HARD criterion | `scripts/pipeline.py` + `wiki/decoder_v0_ux_and_success_criterion.md` | v0 schema (LOCKED 2026-05-18): `prediction` + `calibrated_probability` + `confidence_tier` (HIGH/MEDIUM/LOW via direct probability compare, not `abs(p-0.5)` — float-precision bug at 0.7/0.3 boundaries) + `top_k_attribution` (gene-level ISM + Tier 1-5 catalog) + `audit_verdict` (cohort_gate_verdict + per-strain noise_class + suspend_gate_fired + verdict_explanation) + `provenance` (model, training_cohort, loso_auroc, trained_on). Train pickle enriched with `training_cohort`, `trained_on`, `n_strains`, `auroc_lomo_clade_out` (2026-05-18). 16 unit tests + 6 E2E integration tests with synthetic fixtures. |

## Phase 1 success criteria

See `README.md` and `plans/Ecoli_G2P_Platform_Technical_Plan.md` Verification section.

## Future-state pointers

- `TODOS.md` — deferred work + post-Phase-1 polish
- `FUTURE_FEATURES.md` — Phase 2-4+ ideas
- `LESSONS_LEARNED.md` — practical lessons from Phase 1 build
- `plans/Ecoli_G2P_Platform_Technical_Plan.md` Phase 2 Backlog — Attribution Refinement Engine, Captum + diff MLP head, MIC regression, pan-genome graph layer, AlphaFold-inspired arch (Phase 3+ with compute caveat)
