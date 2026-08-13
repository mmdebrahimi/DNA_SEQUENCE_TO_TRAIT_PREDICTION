"""Which FBA deletion solves were non-optimal — recorded per CELL, not just counted.

A cobrapy `single_gene_deletion` result carries a `status` column beside `growth`. Every deletion script
in this repo originally read only `growth`, so a non-optimal solve was **indistinguishable from a real
growth value**. That matters more than it sounds, because of how the scripts code a missing growth:

    scripts/fba_regulatory_conditional_test.py:106   d[gid] = 0.0 if g != g else g / wt
    scripts/fba_gapfill_carbon_recheck.py:71         d[gid] = (g != g) or (g < FRAC * wt)

`g != g` is the NaN test, and cobrapy returns NaN exactly when a solve is non-optimal or infeasible. Both
lines therefore code **"the solver failed" as "the gene is essential"** — silently. In an arm that forces
off ~69% of gene-associated reactions before deleting anything, that is not a rounding concern.

**Why cells and not counts.** The first audit added to this codebase counted non-optimal solves per
condition. That found 39 across 15 of 25 carbon conditions — but the count alone cannot answer the
question that matters: are those 39 concentrated in the handful of genes where the model actually
commits to a varying pattern? A single non-optimal cell can create or destroy a 1-of-25 exact-set match,
so a 0.7% global rate does not clear a concentrated subset. This module records the cell ids so the
cross-check is possible.

**Verified tool surface (cobrapy 0.31.1, read from the installed source, not from docs):**
`cobra.flux_analysis.deletion._multi_deletion` builds its result with
`columns=["ids", "growth", "status"]` and documents `status : str — The solution's status.` The `ids`
entry is a **frozenset**, not a bare id, so callers use `next(iter(row["ids"]))`.

Pure and solver-free: everything here operates on an already-computed DataFrame, so it is unit-testable
offline with a synthetic frame — no model, no solver, no 7.4 GB `feba.db`.
"""
from __future__ import annotations

from dataclasses import dataclass, field

OPTIMAL = "optimal"


@dataclass
class DeletionAudit:
    """Solver health for ONE condition's deletion sweep.

    `nonoptimal_cells` and `nan_cells` are deliberately SEPARATE sets. They overlap in practice (a
    non-optimal solve usually yields NaN growth) but they are different failure modes, and collapsing
    them would hide the case that actually worries us: a solve reported `optimal` whose growth is still
    NaN, which no status check would catch.
    """

    condition: str
    n_rows: int = 0
    statuses: dict[str, int] = field(default_factory=dict)
    nonoptimal_cells: set[str] = field(default_factory=set)
    nan_cells: set[str] = field(default_factory=set)
    status_available: bool = True

    @property
    def n_nonoptimal(self) -> int:
        return len(self.nonoptimal_cells)

    @property
    def n_nan_growth(self) -> int:
        return len(self.nan_cells)

    @property
    def n_suspect(self) -> int:
        """Cells that are non-optimal OR NaN — the set an abstention arm would exclude."""
        return len(self.nonoptimal_cells | self.nan_cells)

    def suspect_cells(self) -> set[str]:
        return self.nonoptimal_cells | self.nan_cells

    def to_dict(self) -> dict:
        return {
            "condition": self.condition,
            "n_rows": self.n_rows,
            "status_available": self.status_available,
            "statuses": dict(sorted(self.statuses.items())),
            "n_nonoptimal": self.n_nonoptimal,
            "n_nan_growth": self.n_nan_growth,
            "n_suspect": self.n_suspect,
            "nonoptimal_cells": sorted(self.nonoptimal_cells),
            "nan_cells": sorted(self.nan_cells),
        }


def audit_deletion_frame(res, condition: str) -> DeletionAudit:
    """Audit one cobrapy deletion result frame. Never raises on a shape it does not recognise.

    An older cobrapy without a `status` column must not crash a run that would otherwise succeed, so a
    missing column degrades to `status_available=False` and NaN detection continues. That is an honest
    partial audit, not a silent pass — the flag rides into the artifact.
    """
    audit = DeletionAudit(condition=condition)
    cols = getattr(res, "columns", None)
    audit.status_available = cols is not None and "status" in list(cols)

    for _, row in res.iterrows():
        audit.n_rows += 1
        ids = row["ids"]
        # cobrapy stores a frozenset per row (a multi-deletion could hold several ids); every consumer
        # in this repo is a SINGLE gene deletion, so the sole member is the gene id.
        gid = next(iter(ids)) if not isinstance(ids, str) else ids

        if audit.status_available:
            st = row["status"]
            audit.statuses[str(st)] = audit.statuses.get(str(st), 0) + 1
            if st != OPTIMAL:
                audit.nonoptimal_cells.add(gid)

        g = row["growth"]
        if g != g:                                  # NaN
            audit.nan_cells.add(gid)

    return audit


def merge_audits(audits: dict[str, DeletionAudit]) -> dict:
    """JSON-serialisable rollup across conditions, for an artifact sidecar.

    `suspect_cells` is a flat list of `[gene, condition]` pairs rather than a nested dict so it survives
    a JSON round-trip unchanged and can be loaded straight into a set of tuples by a consumer.
    """
    per_condition = {c: a.to_dict() for c, a in sorted(audits.items())}
    suspect = sorted(
        [g, c] for c, a in audits.items() for g in a.suspect_cells()
    )
    total_rows = sum(a.n_rows for a in audits.values())
    return {
        "n_conditions": len(audits),
        "n_rows_total": total_rows,
        "n_nonoptimal_total": sum(a.n_nonoptimal for a in audits.values()),
        "n_nan_growth_total": sum(a.n_nan_growth for a in audits.values()),
        "n_suspect_total": len(suspect),
        "suspect_fraction": round(len(suspect) / total_rows, 6) if total_rows else 0.0,
        "status_available": all(a.status_available for a in audits.values()),
        "n_conditions_with_nonoptimal": sum(1 for a in audits.values() if a.n_nonoptimal),
        "suspect_cells": suspect,
        "per_condition": per_condition,
    }


def suspect_cell_set(audits: dict[str, DeletionAudit]) -> set[tuple[str, str]]:
    """`{(gene_id, condition)}` — the cells an abstention arm excludes from scoring."""
    return {(g, c) for c, a in audits.items() for g in a.suspect_cells()}
