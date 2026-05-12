"""Interpretability layer.

Phase 1: ISM-only (gene-level + nucleotide-level saturation mutagenesis).
Captum IG deferred to Phase 2 (requires differentiable MLP head).

Modules (populated by later plan steps):
- mutagenesis: Step 11 — in-silico mutagenesis
"""
