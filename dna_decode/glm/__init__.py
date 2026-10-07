"""The GLM — the GENERATIVE genome-edit engine. *Name a trait, get the genome edits.*

This package is the product the project is for. It is NOT benchmarked against the deterministic
decoder surface: a lookup table cannot generate, so a head-to-head between them measures the wrong
quantity (recorded as orientation error 0 in `CLAUDE.md`, after four sessions of re-deriving it). The
catalog was scaffolding built so there would be any number at all.

**Why this is buildable without new labels.** The repo's negative-results corpus is real and well
measured, and all of it concerns DISCRIMINATIVE zero-shot scoring of variants in natural populations.
Generation scored by a validated oracle is a different problem: the generator needs an *oracle*, not a
*label set*, and two oracles with real validation already ship here —

* `forward/` — protein variant effect, DMS-validated (Spearman 0.35-0.76 across two beta-lactamases),
* `fba/`     — organism-level growth/flux, Keio accuracy 0.954, per-condition MCC 0.70-0.74,

plus `constraints/`, which is correct BY CONSTRUCTION rather than statistically (reading frame and mass
balance are not claims that can be wrong about a new organism). So: the label wall blocks supervised
genotype->phenotype on natural populations; it does NOT block design-and-score.

**Layering.** `glm` is an INTEGRATION layer and sits at the top of the dependency graph: it may import
`constraints`, `forward`, `fba`, `data`, `eval`; none of those may import `glm`. Keeping that arrow
one-way is what stops a generative experiment from reaching into the frozen decoder surface.

Design + decomposition, written before the prior-art scan so the search checks them rather than
anchoring them: `wiki/glm_architecture_design_2026-10-06.md`, `wiki/glm_decomposition_2026-10-06.md`.
"""
from __future__ import annotations

from dna_decode.glm.edit import (
    DEL,
    EDIT_INTENTS,
    EDIT_KINDS,
    INS,
    REPLACE,
    SUB,
    EditApplicationError,
    EditOverlapError,
    EditSet,
    GenomeEdit,
    RefMismatchError,
    derive_kind,
)

__all__ = [
    "DEL",
    "EDIT_INTENTS",
    "EDIT_KINDS",
    "INS",
    "REPLACE",
    "SUB",
    "EditApplicationError",
    "EditOverlapError",
    "EditSet",
    "GenomeEdit",
    "RefMismatchError",
    "derive_kind",
]
