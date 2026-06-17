"""Step 7 — score the frozen TB cell on an INDEPENDENT post-2023 gold set (deliverable b).

Reuses the Step-3 cell + Step-4 lineage caller + Step-5 collapse, but labels results
`INDEPENDENT_VALIDATION` and scores them SEPARATELY — never merged with the CRyPTIC knowledge-baseline.
BLOCKED-gates to `INDEPENDENT_VALIDATION_BLOCKED_NO_GOLDSET` when no gold set is present.

The gold set is hand-curated per `wiki/tb_independent_goldset_acquisition_2026-06-17.md`; until it lands,
this scorer emits the BLOCKED status (an honest "no independent number yet"), never a fabricated metric.
"""
from __future__ import annotations

from scripts.score_tb_cryptic import BASELINE_LABEL, BLOCKED_LABEL, score_cohort

INDEPENDENT_LABEL = "INDEPENDENT_VALIDATION"
BLOCKED_NO_GOLDSET = "INDEPENDENT_VALIDATION_BLOCKED_NO_GOLDSET"

_INDEPENDENCE_NOTE = (
    "Independent of the WHO-catalogue BUILD: post-2023 isolates (temporal hold-out), since WHO v2 swept "
    "most public pre-2023 TB WGS+pDST. Scored SEPARATELY from the CRyPTIC knowledge-baseline."
)


def score_independent(preds, labels, clusters, *, drug: str) -> dict:
    """Independent-arm scoring. No isolates -> BLOCKED_NO_GOLDSET. Else reuse score_cohort + relabel."""
    if not labels:
        return {"drug": drug, "status": BLOCKED_NO_GOLDSET,
                "reason": "no independent gold set present (hand-curate per the runbook)"}

    out = score_cohort(preds, labels, clusters, drug=drug, cohort_complete=True)
    if out["status"] == BASELINE_LABEL:
        out["status"] = INDEPENDENT_LABEL
        out["independence_note"] = _INDEPENDENCE_NOTE
        out["honesty"] = ("REAL independent validation (post-2023, outside the WHO-catalogue build). "
                          "Small N by design -> read the Wilson CI, not the point estimate.")
    # a BLOCKED_NO_LINEAGE result passes through unchanged (honest)
    elif out["status"] == BLOCKED_LABEL:
        out["independence_note"] = _INDEPENDENCE_NOTE
    return out
