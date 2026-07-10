"""Canonical parser for NCBI Pathogen Detection's `AST_phenotypes` metadata field.

ONE implementation, shared by `scripts/ncbi_pd_provenance_census.py` and
`scripts/fetch_prospective_cohort.py`, because the field has two traps that bit both.

REAL format, verified against the live PD metadata (2026-07-10):

    "ampicillin=ND,cefazolin=ND,ceftriaxone=R,ciprofloxacin=S,gentamicin=ND"

1. The separator is a **comma**, not a semicolon. (The `;` form appears only in the DERIVED
   `candidates.tsv`, never in the source metadata.)
2. The field is wrapped in **literal double quotes**, and those quotes ride on the FIRST and LAST
   tokens. So the naive membership test `f"{drug}=R" in field.split(",")` can never match a drug that
   happens to sit at either end of the list. Measured over 120k E. coli PD rows, that silently missed
   ciprofloxacin +91, ceftriaxone +62, tetracycline +1 isolates (gentamicin never lands at an end, so 0).
   It can only UNDER-count, never over-count -- the safe direction, but a cell recorded as UNDERPOWERED
   may in fact be powered.

Values other than R/S (`ND` not-determined, `I`, `NS`, `NULL`) are DROPPED: the deterministic decoder
emits a binary R/S call, so those isolates carry no unambiguous ground truth to score against.
"""
from __future__ import annotations

VALID_LABELS = ("R", "S")


def parse_ast_phenotypes(field: str | None, wanted_drugs) -> dict[str, str]:
    """`AST_phenotypes` -> {drug: 'R'|'S'} restricted to `wanted_drugs`. See module docstring."""
    raw = (field or "").strip()
    if not raw or raw == "NULL":
        return {}
    wanted = {d.lower() for d in wanted_drugs}
    out: dict[str, str] = {}
    for part in raw.replace(";", ",").split(","):
        token = part.strip().strip('"').strip("'").strip()
        if "=" not in token:
            continue
        drug, _, val = token.partition("=")
        drug, val = drug.strip().lower(), val.strip().upper()
        if drug in wanted and val in VALID_LABELS:
            out[drug] = val
    return out


def ast_label_for(field: str | None, drug: str) -> str | None:
    """Convenience: the R/S label for ONE drug, or None."""
    return parse_ast_phenotypes(field, {drug}).get(drug.lower())
