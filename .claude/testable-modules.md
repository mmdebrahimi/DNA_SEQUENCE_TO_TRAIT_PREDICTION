# Testable Modules
<!-- Auto-maintained by /test-epilogue — do not edit manually -->

- scripts/audit_cohort.py — Phase 2.5 cohort audit report generator; cohort overview + clade composition + metadata completeness + assembly QC quantiles + AST method breakdown + GO/WARN/NO-GO verdict; argparse + cohort-loader + AST-loader integration tests
- scripts/build_mini_cohort.py — Gate B mini-cohort selector; picks N highest-assembly-quality strains per R/S class from a source cohort parquet (sort: contig_count asc, n50 desc); argparse + ValueError surface for missing class
- scripts/cipro_curated_baseline.py — EP1 audit-infrastructure module-level constants + ABLATION_FEATURE_SETS contract + _load_amrfinder_features JSON parsing. 2-layer verdict (original_condition_4 + amended_condition_4 gated on no_POINT >= 0.773 OR mechanism_only >= 0.80). LR + XGB orchestration over multi-block ablation sets skipped.
- scripts/cipro_mechanism_audit.py — EP1 audit-infrastructure pure-logic: _classify_symbol (gene → mechanism class with tolerant prefix match incl. acrR + bla* + qnr*), _is_synonymous_point (first AA == last AA detection), _parse_amrfinder_outputs (main.tsv POINTX → kind=mutation routing; mutations.tsv class+synonymous filter; cross-TSV dedupe). Docker AMRFinder invocation skipped.
- scripts/cipro_mechanism_phenotype_merge.py — EP1 audit-infrastructure pure-logic: _classify_noise per-row mapping with strict PRIMARY_CIPRO_MECHANISMS = QRDR | plasmid_protect_modify; co_resistance_modifiers separated; mechanism_opacity_flag separates AMRFinder-incomplete from labels-noisy. Pre-curated-baseline gate orchestration skipped.
- scripts/cipro_mic_audit.py — EP1 audit-infrastructure pure-logic: _parse_mic (strips <>= operators, handles NA variants) + _confidence_tier (HIGH_R/HIGH_S/DECISIVE/BORDERLINE/AMBIGUOUS/CONFLICT/NO_MIC classification under CLSI 2024 + EUCAST 14.0 breakpoints). Raw AST CSV load skipped.
- scripts/leaderboard.py — formatting, markdown grouping, bundle-reading helpers (subprocess paths mocked)
- scripts/pipeline.py — CLI dispatcher with ingest/train/predict/attribute subcommands; argparse validation + exit codes
- scripts/populate_cache.py — Phase 2 standalone embedding-cache driver; strain_id ↔ assembly_accession mapping; --allow-mock gate; end-to-end mock populate
- scripts/quantize_fidelity_check.py — top-K intersection + Spearman comparison, GO/NO-GO aggregation, markdown report writer
- scripts/smoke_gate_12strain_cipro.py — 2026-05-17 patches: --drug-templated output strings (write_packet heading + cohort description + default filename slug) + NT-XGBoost runner fallback to ast_labels iteration when cohort.per_drug_strain_ids[drug] missing. Full LOSO + NT cache run is orchestration (skipped).
- scripts/smoke_pipeline.py — end-to-end synthetic-fixture smoke run; AUROC + top-1 attribution assertions
- dna_decode/data/bvbrc_genome.py — Phase 2 BV-BRC Genomes-tab CSV adapter; header normalization (Title Case → lowercase_underscore); organism filter (Species → Genome Name fallback); duplicate-key warning; tolerant int parsing
- dna_decode/viz/browser.py — matplotlib attribution plot rendering + TSV export
