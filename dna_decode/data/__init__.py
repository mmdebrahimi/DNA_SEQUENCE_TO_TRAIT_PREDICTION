"""Ingestion + preprocessing + annotation layer.

Modules (populated by later plan steps):
- pilot: Step 0.5 — real-data pilot gate
- refseq: Step 2 — NCBI RefSeq downloader
- annotations: Step 3 — GFF3 / GenBank parser
- resistance_db: Step 4 — CARD + AMRFinder loaders
- ast_data: Step 5 — BV-BRC AST phenotype loader
- cohort: Step 6 — Strain/AST cohort catalog (drug-first)
"""
