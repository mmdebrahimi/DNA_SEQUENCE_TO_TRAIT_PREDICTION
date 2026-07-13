"""Pure adapter: committed co-occurrence artifacts -> a {nodes, edges} graph model (2026-07-13).

Reshapes the ALREADY-COMPUTED co-occurrence + cross-axis-lineage JSONs (no new science, no network,
no compute) into the graph the force-directed browser renders. The de-confound is preserved as a
first-class NODE + EDGE attribute so the visual can show it honestly:

  - a NODE's `lineage_status` comes from the cross-axis leave-one-clade-out result
    (`per_gene[gene].generalizes_beyond_lineage` / `clade_concentrated`):
      "generalizes"      -> the association SURVIVES clade-grouped CV (de-confounded signal)
      "lineage_mediated" -> clade-concentrated: the association may be CLONAL population structure
      "untested"         -> no cross-axis entry for this determinant
  - an EDGE is `lineage_mediated=True` iff EITHER endpoint is clade-concentrated — so a co-occurrence
    that could be an artifact of clonal lineage sampling renders dashed, never as a solid/causal link.

Input artifacts (schemas pinned by the census that produced them):
  determinant_cooccurrence_result_*.json   -> per_organism[org].lift_table (edges) + .per_determinant (nodes)
  crossaxis_lineage_deconfound_determinant*.json (optional, per organism) -> per_gene lineage status
"""
from __future__ import annotations

import json
from pathlib import Path

# --- lineage-status vocabulary (the de-confound signal) -------------------------------------------
GENERALIZES = "generalizes"
LINEAGE_MEDIATED = "lineage_mediated"
UNTESTED = "untested"


def _lineage_status(entry: dict | None) -> str:
    """Map one cross-axis per_gene record -> a node lineage_status."""
    if not entry:
        return UNTESTED
    if entry.get("clade_concentrated"):
        return LINEAGE_MEDIATED
    if entry.get("generalizes_beyond_lineage"):
        return GENERALIZES
    return LINEAGE_MEDIATED  # tested, drops below the bar under clade-grouping -> clonally driven


def build_graph(cooc_path: str | Path, organism: str,
                crossaxis_path: str | Path | None = None,
                *, min_cooc: int = 8) -> dict:
    """Return {meta, nodes, edges} for ONE organism from the committed artifacts.

    `min_cooc` drops the long tail of weak co-occurrences (matches the census MIN_COOC default) so the
    graph is legible; it is recorded in meta so the pruning is never silent.
    """
    cooc = json.loads(Path(cooc_path).read_text(encoding="utf-8"))
    per_org = cooc.get("per_organism", {})
    if organism not in per_org:
        raise KeyError(f"organism {organism!r} not in {cooc_path} (have: {sorted(per_org)})")
    org = per_org[organism]

    lineage: dict[str, dict] = {}
    if crossaxis_path is not None:
        cx = json.loads(Path(crossaxis_path).read_text(encoding="utf-8"))
        lineage = cx.get("per_gene", {}) or {}

    per_det = org.get("per_determinant", {}) or {}
    lift_table = org.get("lift_table", {}) or {}

    # --- nodes: every determinant that appears as a co-occurrence endpoint ------------------------
    endpoints: set[str] = set(lift_table)
    for tgt, partners in lift_table.items():
        for p in partners:
            if int(p.get("cooc", 0)) >= min_cooc:
                endpoints.add(p["det"])

    nodes = []
    for det in sorted(endpoints):
        pd = per_det.get(det, {}) or {}
        nodes.append({
            "id": det,
            "organism": organism,
            "prevalence": int(pd.get("n_present", 0)),      # node size
            "auc": pd.get("auc"),                            # imputation strength (may be None)
            "ci_lo": pd.get("ci_lo"), "ci_hi": pd.get("ci_hi"),
            "linked": bool(pd.get("linked", False)),
            "lineage_status": _lineage_status(lineage.get(det)),
        })
    lm_nodes = {n["id"] for n in nodes if n["lineage_status"] == LINEAGE_MEDIATED}

    # --- edges: undirected, deduped on the unordered pair, keep the max cooc/lift -----------------
    seen: dict[frozenset, dict] = {}
    for tgt, partners in lift_table.items():
        for p in partners:
            det = p["det"]
            n_cooc = int(p.get("cooc", 0))
            if n_cooc < min_cooc or det == tgt:
                continue
            key = frozenset((tgt, det))
            if len(key) != 2:
                continue
            prev = seen.get(key)
            if prev is None or n_cooc > prev["cooc"]:
                seen[key] = {
                    "source": tgt, "target": det, "cooc": n_cooc,
                    "lift": p.get("lift"), "p_given_target": p.get("p_given_target"),
                }
    edges = []
    for e in seen.values():
        e = dict(e)
        e["lineage_mediated"] = (e["source"] in lm_nodes) or (e["target"] in lm_nodes)
        edges.append(e)
    edges.sort(key=lambda e: (-e["cooc"], e["source"], e["target"]))

    return {
        "meta": {
            "organism": organism,
            "cooc_artifact": Path(cooc_path).name,
            "crossaxis_artifact": Path(crossaxis_path).name if crossaxis_path else None,
            "verdict": cooc.get("verdict"),
            "n_genomes": org.get("n_genomes"),
            "n_determinants": org.get("n_determinants"),
            "min_cooc": min_cooc,
            "n_nodes": len(nodes), "n_edges": len(edges),
            "n_lineage_mediated_nodes": len(lm_nodes),
            "honest_caveats": cooc.get("honest_caveats"),
        },
        "nodes": nodes,
        "edges": edges,
    }
