"""Flagellar motility decoder — the AMR/metabolic determinant->phenotype paradigm applied to a
NON-metabolic cell-level trait (the first non-metabolic trait catalog).

Gene presence -> can the cell SWIM? (motile / non-motile), with chemotaxis reported separately.
"""
from .flagellar_catalog import (
    CHEMOTAXIS,
    MOTILITY_MODULES,
    MotilityCall,
    MotilityInputError,
    call_motility,
    catalog_genes,
)

__all__ = [
    "CHEMOTAXIS",
    "MOTILITY_MODULES",
    "MotilityCall",
    "MotilityInputError",
    "call_motility",
    "catalog_genes",
]
