"""Evaluation harness: CV strategies + metrics + phylogeny controls.

Modules (populated by later plan steps):
- cv: Step 10 — LOSO + LOMO + leave-one-Mash-clade-out
- phylogeny: Step 10 — Mash/ANI distance + clade clustering
- clade_baseline: Step 10 — clade-only baseline classifier
- metrics: Step 10 — AUROC + AUPRC + attribution-precision + per-clade
- leaderboard: Step 17 — comparative model benchmarking
"""
