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
    ap.add_argument("--synthetic-lethality", action="store_true",
                    help="with --knockout A,B: is the PAIR synthetic-lethal? (both singles viable, double lethal)")
    # compose forward+fba: a POINT MUTATION -> LOF? -> cell trait (requires the protein sequence)
    ap.add_argument("--mutation", default=None,
                    help="a missense in --gene (e.g. D362A); composes forward (LOF?) -> fba (cell trait)")
    seqsrc = ap.add_mutually_exclusive_group()
    seqsrc.add_argument("--protein-seq", default=None, help="protein sequence of --gene (for --mutation)")
    seqsrc.add_argument("--protein-fasta", default=None, help="FASTA file with --gene's protein (for --mutation)")
    ap.add_argument("--forward-method", default="blosum62",
                    choices=["blosum62", "esm2", "prosst", "gemme", "hybrid", "auto"],
                    help="forward scorer for the missense (default blosum62, offline)")
    ap.add_argument("--organism", default=None,
                    help="cross-organism GEM: ecoli(default) | saureus | salmonella | pputida | yeast "
                         "(engine generalizes; only E. coli is Keio-validated). "
                         "P. aeruginosa has NO BiGG model and is refused rather than substituted.")
    ap.add_argument("--model-id", default=None, help="a BiGG model id directly (e.g. iYS854)")
    ap.add_argument("--model", default=None, help="path to an SBML model (override organism/model-id)")
    ap.add_argument("--json", action="store_true", help="emit JSON")
    a = ap.parse_args(argv)

    from .model import call_essential, knockout_growth, load_model, organism_for, wildtype_growth

    try:
        model = load_model(a.model, organism=a.organism, model_id=a.model_id)
    except Exception as e:  # ImportError / FileNotFoundError / parse / unknown organism
        print(f"ERROR: {e}", file=sys.stderr)
        return 2

    # Provenance is read from the LOADED model, never hardcoded -- a record must never claim an
    # organism the run did not actually use (the v0.11.0-v0.12.0 defect).
    loaded_id = a.model_id or getattr(model, "id", None) or "unknown"
    loaded_organism = organism_for(loaded_id)

    wt = wildtype_growth(model)

    # --- compose mode: a POINT MUTATION -> forward (LOF?) -> fba (cell trait) ---
    if a.mutation:
        if not a.gene:
            print("ERROR: --mutation needs --gene (the metabolic gene the missense is in)", file=sys.stderr)
            return 2
        seq = a.protein_seq
        if a.protein_fasta:
            from pathlib import Path
            lines = Path(a.protein_fasta).read_text(encoding="utf-8").splitlines()
            seq = "".join(x.strip() for x in lines if x and not x.startswith(">"))
        if not seq:
            print("ERROR: --mutation needs --protein-seq or --protein-fasta (the gene's protein, for "
                  "forward to score the missense)", file=sys.stderr)
            return 2
        gid = _resolve_gene(model, a.gene)
        if gid is None:
            print(f"ERROR: gene '{a.gene}' not in {loaded_id} (metabolic genes only)", file=sys.stderr)
            return 2
        from .compose import variant_to_cell_trait
        try:
            crec = variant_to_cell_trait(model, gid, seq, a.mutation, method=a.forward_method, frac=a.frac)
        except Exception as e:  # WT-mismatch / forward import / parse
            print(f"ERROR: {e}", file=sys.stderr)
            return 2
        if a.json:
            print(json.dumps(crec, indent=2))
            return 0
        f = crec["forward"]
        print(f"{loaded_id} FBA + forward ({loaded_organism})  |  WT growth {crec['wildtype_growth_per_h']} /h")
        print(f"  missense {a.gene} {crec['mutation']}  ->  forward({f['method']}) "
              f"{f['predicted_effect']} (raw {f['raw_score']})  ->  {crec['lof_call']}")
        if crec["fba_action"] == "conditional":
            print(f"  UNCERTAIN -> reported both ways (forward is a validated RANK, not a calibrated LOF prob):")
            print(f"    if LOF:       {crec['cell_trait_if_LOF']}  (KO growth {crec['ko_growth_per_h_if_LOF']} /h)")
            print(f"    if tolerated: {crec['cell_trait_if_tolerated']}")
        else:
            print(f"  cell-level trait: {crec['cell_trait']}")
        print(f"  scope: {crec['scope']}")
        return 0

    rec = {
        "record": "fba-metabolic-trait-v1",
        "model": loaded_id,
        "organism": loaded_organism,
        # The medium is whatever the loaded GEM's default exchange bounds encode -- for iML1515 that
        # is glucose M9 aerobic, but it is NOT the same for other organisms' models, so do not claim it.
        "medium": ("model default exchange bounds (iML1515: glucose M9 aerobic)"
                   if loaded_id == "iML1515" else f"model default exchange bounds ({loaded_id})"),
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
                print(f"ERROR: gene '{t}' not in {loaded_id} (metabolic genes only; "
                      f"non-metabolic genes are out of model scope)", file=sys.stderr)
                return 2
            resolved.append(gid)
        # synthetic-lethality mode: --knockout A,B --synthetic-lethality
        if a.synthetic_lethality:
            if len(resolved) != 2:
                print("ERROR: --synthetic-lethality needs exactly two genes (--knockout A,B)", file=sys.stderr)
                return 2
            from .model import synthetic_lethality
            sl = synthetic_lethality(model, resolved[0], resolved[1], a.frac)
            sl["genes"] = tokens
            sl["scope"] = f"METABOLIC traits only; {loaded_id} ({loaded_organism}). NOT clinical."
            if a.json:
                print(json.dumps(sl, indent=2))
                return 0
            print(f"{loaded_id} FBA ({loaded_organism})  |  WT growth {sl['wildtype_growth_per_h']} /h")
            print(f"  KO {tokens[0]} alone: {sl['ko_a_growth_per_h']} /h "
                  f"({'essential' if sl['single_a_essential'] else 'viable'})")
            print(f"  KO {tokens[1]} alone: {sl['ko_b_growth_per_h']} /h "
                  f"({'essential' if sl['single_b_essential'] else 'viable'})")
            print(f"  KO {tokens[0]}+{tokens[1]} together: {sl['double_ko_growth_per_h']} /h "
                  f"({'LETHAL' if sl['double_essential'] else 'viable'})")
            print(f"  verdict: {sl['verdict']}")
            print(f"  scope: {sl['scope']}")
            return 0
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
    print(f"{loaded_id} FBA ({loaded_organism})  |  WT growth {rec['wildtype_growth_per_h']} /h "
          f"[{rec['medium']}]")
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
