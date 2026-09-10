"""A cell's measured evidence, as one line, at the point of the call.

THE GAP. `dna-amr` has carried an inline trust badge since `trust_surface` shipped: every AMR call
prints its cell's honest tier and headline metric. No other route does. A user running `dna-ktype` gets
a hand-written caveat and nothing else, while that cell's MEASURED evidence -- 0.8788 agreement against
full-locus Kaptive on 307 genomes, verdict NEAR_THE_WZI_METHOD_CEILING -- sits in `cell_registry` and a
report card the user never opens. That is the same failure the doubt layer had and the report card had:
evidence carried only in a registry is not a disclosure.

WHY A SEPARATE MODULE FROM `trust_surface`. `trust_surface.trust_block` is keyed on (drug, organism) and
resolves through the antimicrobial report cards. Most routes here have no drug axis at all -- a capsule
type, a serovar, a variant effect -- so they need the REGISTRY as their source, not a card. Two sources,
two accessors, one honesty standard.

WHAT IT WILL NOT DO:
  * It never invents a number. Every figure comes from the cell's own committed contract.
  * It never upgrades a tier. `KNOWLEDGE_BASELINE` renders as never-measured, which is the honest line
    for a cell nobody has scored -- an absence of measurement, not an absence of doubt.
  * It never replaces a caller's caveat. It is appended beside it, augment-only.

Pure. No I/O beyond importing the committed registry.
"""
from __future__ import annotations

from .cell_registry import cells

# How each tier should READ to someone deciding whether to trust this call. The wording carries the
# limit, not just the label: a tier name alone ("FAITHFUL_TO_TOOL") means nothing at a terminal.
_TIER_GLOSS = {
    "independent_measured": "measured against an INDEPENDENT label (wet-lab or equivalent)",
    "near_independent": "measured against a near-independent label",
    "faithful_to_tool": "agreement with a reference TOOL -- bounds agreement, never correctness",
    "knowledge_baseline": "NEVER MEASURED against a phenotype -- curated-knowledge baseline only",
    "no_free_source": "no free isolate-level phenotype source exists for this cell",
    "not_censused": "shipped but never censused -- no evidence recorded",
}


def _tier_value(cell) -> str:
    t = getattr(cell, "evidence_tier", None)
    return getattr(t, "value", None) or str(t or "").lower()


def cells_for_route(route: str) -> list:
    """Every registered cell served by a CLI route. A route may serve several (drug- or target-keyed)."""
    r = (route or "").strip().lower()
    return [c for c in cells() if str(getattr(c, "route", "")).strip().lower() == r]


def evidence_one_line(route: str, *, max_len: int = 300) -> str | None:
    """One line of measured evidence for a route, or None when the route has no registered cell.

    None means "this route is not in the registry" -- a real answer, and different from "measured and
    clean". A caller that wants the distinction should check `cells_for_route` directly rather than
    reading silence as reassurance.
    """
    found = cells_for_route(route)
    if not found:
        return None
    tiers = {_tier_value(c) for c in found}
    # A route serving several cells can span tiers; report the WEAKEST, because the badge must not
    # imply the strongest evidence applies to whichever cell the caller actually invoked.
    order = ("not_censused", "no_free_source", "knowledge_baseline", "faithful_to_tool",
             "near_independent", "independent_measured")
    weakest = next((t for t in order if t in tiers), sorted(tiers)[0] if tiers else "")
    gloss = _TIER_GLOSS.get(weakest, weakest)

    lead = found[0] if len(found) == 1 else next((c for c in found if _tier_value(c) == weakest), found[0])
    claim = (getattr(lead, "claim_status", "") or "").replace("_", " ").strip()
    span = f" ({len(found)} cells, weakest shown)" if len(tiers) > 1 else ""
    line = f"evidence: {weakest.upper()}{span} -- {gloss}"
    if claim:
        line += f"; {claim}"
    return line if len(line) <= max_len else line[: max_len - 1] + "…"


def routes_missing_evidence_line(rendered_routes) -> list[str]:
    """Registered routes that do NOT yet print an evidence line. The coverage gap, made countable.

    Shipping a partial wiring is fine; shipping it silently is not. This is what lets the guard report
    'N of M routes surface their evidence' instead of a green test over whichever subset was done.
    """
    have = {str(r).strip().lower() for r in rendered_routes}
    registered = {str(getattr(c, "route", "")).strip().lower() for c in cells()}
    registered.discard("")
    return sorted(registered - have)
