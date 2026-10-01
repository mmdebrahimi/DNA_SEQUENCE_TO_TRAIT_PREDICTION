"""NON-FROZEN constraint layer — universal biological laws plus clade/regime specialization.

Additive by construction, following the `dna_decode/data/experimental_drug_rules.py` contract: adding a
constraint here touches NO frozen file (`amr_rules.py`, `calibrated_amr_rules.json`, `mic_tiers.py`,
`shipped_decoder_surface.py`, `cohort_manifest.py` stay byte-identical).

Two kinds of rule, and the distinction decides what transfers to an uncharacterised organism:

* **Type A — laws (`axiomatic`).** Universal, organism-independent, true for life we have never sequenced:
  reading frame, mass balance, codon reachability, reference-base agreement. These transfer BY
  CONSTRUCTION and are what a global model needs.
* **Type B — curated associations (`sourced`).** Statements about a specific organism derived from wet-lab
  work on it. They do NOT transfer to an uncharacterised organism; the deployed decoder surface already
  holds these and they are not re-implemented here.

Two ORTHOGONAL specialization axes, because this project measured that taxonomy is not the most predictive
one: the regime axis (natural vs constructed population design, concentrated vs distributed determinant)
discriminated outcomes where taxonomy did not. The regime keys are IMPORTED from `dna_decode.eval.regime`
and never re-declared, so there is one source of truth.

This package is a leaf layer. It imports only `dna_decode.eval.regime`; `forward/`, `typing/` and
`scripts/` import FROM it, never the reverse.
"""
from __future__ import annotations

from dna_decode.constraints.codon_tables import (
    STANDARD,
    STANDARD_STOPS,
    TABLE_IDS,
    TABLES,
    UnknownCladeError,
    UnverifiedTableError,
    stops_for,
    table_for,
    table_meta,
)

__all__ = [
    "STANDARD",
    "STANDARD_STOPS",
    "TABLE_IDS",
    "TABLES",
    "UnknownCladeError",
    "UnverifiedTableError",
    "stops_for",
    "table_for",
    "table_meta",
]
