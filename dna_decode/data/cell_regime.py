"""Which g->p REGIME does each shipped cell sit in — and a tripwire so a new one must declare it.

WHY. `eval/regime.py` encodes the regime boundary as a function because the boundary has been
mis-stated three times, always by compressing a SCOPED negative into a general one. But that function
only ran where someone remembered to call it. Nothing connected it to the place cells actually enter the
tool, so a decoder could ship in a regime with a recorded negative and nothing would say so.

WHAT WAS MEASURED FIRST, and it changed the design. The intended fix was a `regime` field on
`CellContract`. Counting the construction sites (52) and the regimes actually present made that look like
ceremony -- until the census showed the column is NOT constant: `dna-decode-forward` / `dna-decode-inverse`
are constructed-MOLECULAR (the learned regime that works), `dna-fba` / `dna-essentiality` are
constructed-ORGANISM, and the other 40 routes are curated catalogs. Three regimes, not one. So the mapping
is real information -- it just belongs at the ROUTE level, where it is 44 entries instead of 52 edits of a
mostly-constant column.

THE TRIPWIRE. `regime_for_route` RAISES on a route it does not know. That is the point: a new decoder
cannot register without someone deciding which regime it is in, and a learned decoder aimed at
`natural_organism_zeroshot` is landing in a CLOSED_NEGATIVE. The route set is hand-listed, which is
normally the drift smell -- here the drift is DETECTED by a test rather than silent, and the friction is
the feature (cf. the colour-cell freeze). It is an attention guard, not enforcement.

HONEST SCOPE. This records which regime a cell's METHOD sits in. It does NOT re-validate the cell, and a
`WORKS` regime is not a claim about that particular cell's accuracy -- the evidence tier carries that.
"""
from __future__ import annotations

from dna_decode.eval.regime import REGIMES

CURATED_CATALOG = "curated_catalog_exists"
CONSTRUCTED_MOLECULAR = "constructed_molecular"
CONSTRUCTED_ORGANISM = "constructed_organism_per_condition"


class UndeclaredRegime(KeyError):
    """A route registered a cell without declaring which g->p regime its method sits in."""


# Routes whose method is NOT a curated deterministic catalog. Everything else is one.
ROUTE_REGIME: dict[str, str] = {
    # learned molecular-effect scorers (ESM2 / ProSST / GEMME), DMS-validated
    "dna-decode-forward": CONSTRUCTED_MOLECULAR,
    "dna-decode-inverse": CONSTRUCTED_MOLECULAR,
    # constructed gene deletions -> organism-level growth (iML1515 FBA / conservation-rule essentiality)
    "dna-fba": CONSTRUCTED_ORGANISM,
    "dna-essentiality": CONSTRUCTED_ORGANISM,
}

# The curated-catalog routes as of 2026-09-01. A route absent from BOTH maps raises.
CATALOG_ROUTES: frozenset[str] = frozenset({
    "dna-amr", "dna-clinvar", "dna-hla", "dna-pgx", "dna-pathotype", "dna-phage",
    "dna-plasmid", "dna-serotype", "dna-mlst", "dna-ktype", "dna-salmserovar",
    # CORRECTED 2026-09-04: this read "dna-pneumoserotype", which is not a console script. The old
    # registry shorthand derived the route as dna-<trait> and produced the same wrong name, so two
    # artifacts agreed on a route that does not exist. The real entry point is hyphenated.
    "dna-pneumo-serotype", "dna-resfinder", "dna-pointfinder", "dna-disinfinder",
    "dna-kleb", "dna-motility", "dna-metabolic", "dna-morphology", "dna-pigment",
    "dna-flowering",
    # the frozen 19-cell colour/plumage family
    "dna-alpacacolor", "dna-buffalocolor", "dna-camelcolor", "dna-catcolor", "dna-cattlecolor",
    "dna-coatcolor", "dna-donkeycolor", "dna-foxcolor", "dna-goatcolor", "dna-guineapigcolor",
    "dna-horsecolor", "dna-minkcolor", "dna-mousecolor", "dna-pigcolor", "dna-pigeoncolor",
    "dna-plumage", "dna-rabbitcolor", "dna-roedeercolor", "dna-sheepcolor",
})

# Routes where the LEARNED attempt at this same trait is a recorded negative. The shipped cell is a
# catalog and is unaffected; this exists so nobody re-proposes the learned version as novel.
LEARNED_ATTEMPT_CLOSED: dict[str, str] = {
    "dna-flowering": "Arabidopsis FT10 embeddings: within-group r2 -0.13 vs structure-only spearman "
                     "0.48 -- the embedding learned population structure. "
                     "wiki/phase2_arabidopsis_result_2026-06-12.md",
    "dna-pathotype": "pathotype labels are sampling-defined (G3) and partly tool-derived (G1). "
                     "wiki/negative_results_map_2026-06-13.md",
}

_BY_KEY = {r.key: r for r in REGIMES}


def regime_for_route(route: str) -> str:
    """The regime key for a route. RAISES on an unknown route -- that is the tripwire."""
    if route in ROUTE_REGIME:
        return ROUTE_REGIME[route]
    if route in CATALOG_ROUTES:
        return CURATED_CATALOG
    raise UndeclaredRegime(
        f"route {route!r} registers cells but declares no g->p regime. Add it to ROUTE_REGIME (if its "
        f"method is learned or its variation constructed) or to CATALOG_ROUTES (if it is a curated "
        f"deterministic catalog). Check dna_decode/eval/regime.py first: a learned decoder aimed at a "
        f"natural population with an organism-level endpoint is a CLOSED_NEGATIVE regime.")


def regime_record(route: str) -> dict:
    """Route -> its regime, that regime's measured verdict, and any closed learned attempt."""
    key = regime_for_route(route)
    r = _BY_KEY[key]
    out = {"route": route, "regime": key, "regime_verdict": r.verdict, "evidence": r.evidence,
           "artifact": r.artifact}
    if route in LEARNED_ATTEMPT_CLOSED:
        out["learned_attempt_closed"] = LEARNED_ATTEMPT_CLOSED[route]
    return out


def regime_census(routes) -> dict:
    """Census a set of routes. Raises on the first undeclared one rather than defaulting it."""
    recs = [regime_record(r) for r in sorted(set(routes))]
    counts: dict[str, int] = {}
    for rec in recs:
        counts[rec["regime"]] = counts.get(rec["regime"], 0) + 1
    return {"n_routes": len(recs), "by_regime": counts, "routes": recs,
            "note": "Records which regime each route's METHOD sits in. NOT a re-validation: a WORKS "
                    "regime says the regime has a measured positive, not that this cell is accurate. "
                    "The cell's evidence tier carries that."}
