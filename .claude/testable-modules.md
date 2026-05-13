# Testable Modules
<!-- Auto-maintained by /test-epilogue — do not edit manually -->

- scripts/leaderboard.py — formatting, markdown grouping, bundle-reading helpers (subprocess paths mocked)
- scripts/pipeline.py — CLI dispatcher with ingest/train/predict/attribute subcommands; argparse validation + exit codes
- scripts/populate_cache.py — Phase 2 standalone embedding-cache driver; strain_id ↔ assembly_accession mapping; --allow-mock gate; end-to-end mock populate
- dna_decode/data/bvbrc_genome.py — Phase 2 BV-BRC Genomes-tab CSV adapter; header normalization (Title Case → lowercase_underscore); organism filter (Species → Genome Name fallback); duplicate-key warning; tolerant int parsing
- scripts/quantize_fidelity_check.py — top-K intersection + Spearman comparison, GO/NO-GO aggregation, markdown report writer
- scripts/smoke_pipeline.py — end-to-end synthetic-fixture smoke run; AUROC + top-1 attribution assertions
- dna_decode/viz/browser.py — matplotlib attribution plot rendering + TSV export
