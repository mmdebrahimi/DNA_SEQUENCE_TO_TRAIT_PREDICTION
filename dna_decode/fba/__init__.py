"""FBA metabolic-model cell: a gene edit -> a quantitative cell-level trait.

Mechanistic, constraint-based (flux-balance analysis) decoder over the **iML1515**
E. coli K-12 genome-scale metabolic model (Monk et al. 2017, Nat Biotechnol 35:904) via
`cobrapy`. Unlike the LEARNED regime (a closed negative on organism-level traits — it learns
population structure, not causation), FBA computes phenotype from stoichiometry + known
biochemistry, so it sidesteps population-structure confounding by construction.

This is the first rung where "edit -> quantitative cell-level trait" becomes GENERAL — any of
the model's 1515 genes, not a curated list:

    knock out gene X  ->  predicted growth rate (/h) on a defined medium
                      ->  essential / non-essential call

**Scope (honest):** METABOLIC traits only — growth rate, single/double gene-KO essentiality,
growth on a carbon source. NOT virulence / biofilm / motility / regulation (those are not in a
stoichiometric model). Validated against the free Keio-collection mutant-fitness gold standard
(Baba 2006; Bernstein 2023 method) — see `scripts/fba_keio_validate.py`.
"""
from .model import (
    MODEL_NAME,
    call_essential,
    gene_essentiality,
    knockout_growth,
    load_model,
    resolve_model_path,
    wildtype_growth,
)

__all__ = [
    "MODEL_NAME",
    "call_essential",
    "gene_essentiality",
    "knockout_growth",
    "load_model",
    "resolve_model_path",
    "wildtype_growth",
]
