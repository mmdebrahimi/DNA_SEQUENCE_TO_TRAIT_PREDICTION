"""`dna-fba` — gene edit -> quantitative cell-level trait via iML1515 FBA.

    dna-fba --gene b0720                 # KO one gene -> growth rate + essential call
    dna-fba --knockout b0720,b0721       # double KO
    dna-fba --wildtype                   # wild-type growth on the current medium
    dna-fba --gene gltA                  # gene NAME also accepted (resolved to b-number)

Mechanistic (constraint-based), deterministic, no learned model. Scope: METABOLIC traits only.
"""
from __future__ import annotations

import argparse
import json
import sys


def _resolve_gene(model, token: str):
    """Accept a b-number (b0720) or a gene NAME (gltA); return the model gene id or None."""
    ids = {g.id for g in model.genes}
    if token in ids:
        return token
    for g in model.genes:
        if getattr(g, "name", None) and g.name.lower() == token.lower():
            return g.id
    return None


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="dna-fba", description=__doc__)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--gene", help="single gene to knock out (b-number or gene name)")
    src.add_argument("--knockout", help="comma-separated genes to knock out together")
    src.add_argument("--wildtype", action="store_true", help="report wild-type growth only")
    ap.add_argument("--frac", type=float, default=0.01,
                    help="essential threshold as fraction of WT growth (default 0.01)")
    ap.add_argument("--model", default=None, help="path to an iML1515 SBML (override)")
    ap.add_argument("--json", action="store_true", help="emit JSON")
    a = ap.parse_args(argv)

    from .model import call_essential, knockout_growth, load_model, wildtype_growth

    try:
        model = load_model(a.model)
    except Exception as e:  # ImportError / FileNotFoundError / parse
        print(f"ERROR: {e}", file=sys.stderr)
        return 2

    wt = wildtype_growth(model)
    rec = {
        "record": "fba-metabolic-trait-v1",
        "model": "iML1515",
        "organism": "Escherichia coli K-12",
        "medium": "model default (glucose M9 aerobic)",
        "wildtype_growth_per_h": round(wt, 4),
        "trait_axis": "growth rate (/h) + gene-KO essentiality",
        "scope": "METABOLIC traits only (growth/essentiality/secretion); NOT virulence/regulation. NOT clinical.",
        "method": "flux-balance analysis (cobrapy); mechanistic, deterministic",
    }

    if a.wildtype:
        rec["query"] = "wildtype"
    else:
        raw = a.gene if a.gene else a.knockout
        tokens = [t.strip() for t in raw.split(",") if t.strip()]
        resolved = []
        for t in tokens:
            gid = _resolve_gene(model, t)
            if gid is None:
                print(f"ERROR: gene '{t}' not in iML1515 (metabolic genes only; "
                      f"non-metabolic genes are out of model scope)", file=sys.stderr)
                return 2
            resolved.append(gid)
        growth = knockout_growth(model, resolved)
        essential = call_essential(growth, wt, a.frac)
        rec["query"] = {"knockout": tokens, "resolved_ids": resolved}
        rec["ko_growth_per_h"] = round(growth, 4)
        rec["growth_fraction_of_wt"] = round(growth / wt, 4) if wt else 0.0
        rec["essential"] = essential
        rec["cell_trait"] = "NON-VIABLE (essential gene)" if essential else "viable"

    if a.json:
        print(json.dumps(rec, indent=2))
        return 0

    # human-readable
    print(f"iML1515 FBA  |  WT growth {rec['wildtype_growth_per_h']} /h (glucose M9 aerobic)")
    if a.wildtype:
        print("  wild-type growth reported.")
    else:
        q = ",".join(rec["query"]["knockout"])
        print(f"  KO {q}  ->  growth {rec['ko_growth_per_h']} /h "
              f"({rec['growth_fraction_of_wt']:.1%} of WT)")
        print(f"  cell-level trait: {rec['cell_trait']}"
              f"  [{'ESSENTIAL' if rec['essential'] else 'NON-ESSENTIAL'}]")
    print(f"  scope: {rec['scope']}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
