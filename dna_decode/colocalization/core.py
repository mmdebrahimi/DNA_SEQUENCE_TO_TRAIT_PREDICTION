"""Pure co-localization core — join resistance genes and plasmid replicons by assembly contig. No I/O.

Input: genes = [{"gene", "contig"}], replicons = [{"replicon", "contig"}] (contigs from the engine's
positions mode). Output: per-gene plasmid-borne candidacy (gene shares a contig with >=1 replicon) +
per-contig grouping. A shared contig is SUGGESTIVE (plasmids often assemble to one contig), not proof of a
single replicon — the caller adds that caveat.
"""
from __future__ import annotations

from collections.abc import Iterable


def colocalize(genes: Iterable[dict], replicons: Iterable[dict]) -> dict:
    """Link acquired resistance genes to plasmid replicons sharing the same contig.

    Returns:
      gene_calls: [{gene, contig, plasmid_borne (bool), replicons_on_contig: [...]}]
      contigs: {contig: {replicons: [...], genes: [...], has_plasmid_marker: bool}}
      summary: counts.
    """
    rep_by_contig: dict[str, list[str]] = {}
    for r in replicons:
        c = r.get("contig")
        if c is None:
            continue
        rep_by_contig.setdefault(c, []).append(r["replicon"])

    gene_calls = []
    contigs: dict[str, dict] = {}
    for g in genes:
        c = g.get("contig")
        reps = sorted(set(rep_by_contig.get(c, []))) if c is not None else []
        plasmid_borne = bool(reps)
        gene_calls.append({"gene": g["gene"], "contig": c, "plasmid_borne": plasmid_borne,
                           "replicons_on_contig": reps})
        if c is not None:
            contigs.setdefault(c, {"replicons": [], "genes": [], "has_plasmid_marker": False})
            contigs[c]["genes"].append(g["gene"])
    for c, reps in rep_by_contig.items():
        contigs.setdefault(c, {"replicons": [], "genes": [], "has_plasmid_marker": False})
        contigs[c]["replicons"] = sorted(set(reps))
        contigs[c]["has_plasmid_marker"] = True
    for c in contigs:
        contigs[c]["genes"] = sorted(set(contigs[c]["genes"]))

    n_borne = sum(1 for g in gene_calls if g["plasmid_borne"])
    return {
        "gene_calls": sorted(gene_calls, key=lambda g: (not g["plasmid_borne"], g["gene"])),
        "contigs": contigs,
        "summary": {"n_genes": len(gene_calls), "n_plasmid_borne": n_borne,
                    "n_chromosomal_or_unplaced": len(gene_calls) - n_borne,
                    "n_replicons": sum(len(v) for v in rep_by_contig.values())},
    }
