"""Canonical (organism, drug) join key — the SINGLE source of truth.

Three sidecars must agree on how an (organism, drug) cell is identified:
  - the report card (`scripts/build_validation_report_card.py`)
  - the lineage-metrics sidecar (`scripts/compute_lineage_metrics.py`)
  - the SCORED provenance-disjoint JSONs

Before this module each consumer had its own lowercasing `_key`; the disclosure
layer joins lineage metrics onto report-card cells, so a drift in normalization
would silently drop a lineage row. This is the one canonical key (M2).
"""
from __future__ import annotations


def canonical_cell_key(organism: str, drug: str) -> tuple[str, str]:
    """Canonical (organism, drug) cell key: stripped + lowercased."""
    return (organism.strip().lower(), drug.strip().lower())


def cell_key_str(organism: str, drug: str) -> str:
    """String form of the canonical key, for JSON keys / logging: 'org|drug'."""
    org, drug_l = canonical_cell_key(organism, drug)
    return f"{org}|{drug_l}"
