"""AMR cross-tool concordance — in-package CLI (console `dna-concordance`; also `dna-decode concordance`).

Compares the two INDEPENDENT acquired-gene callers on one genome: AMRFinder (via an existing dna-amr run's
main.tsv — pure, no Docker) vs ResFinder (blastn over its DB). Reports family-level both / amr-only /
resfinder-only + Jaccard agreement. Offline-safe (resfinder side degrades; amr side needs a cached run).

    dna-concordance assembly.fna --amrfinder-run data/amrfinder_runs/GCA_xxx --resfinder-db-dir data/resfinder_db
    dna-concordance --amr-genes blaCTX-M-15,sul1 --resfinder-genes blaCTX-M-15_1_X,tet(A)_6_Y   # gene-set mode
"""
from __future__ import annotations

import argparse
import datetime
import json
import sys
from pathlib import Path

from dna_decode.concordance.core import amr_acquired_genes_from_main, compare
from dna_decode.resfinder.runner import call_resistance_genes, gene_of


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="dna-concordance",
                                 description="AMR cross-tool concordance: AMRFinder (dna-amr) vs ResFinder")
    ap.add_argument("fasta", type=Path, nargs="?", help="genome FASTA (for the ResFinder blastn side)")
    ap.add_argument("--amrfinder-run", type=Path, help="existing AMRFinder run dir (main.tsv) for the amr side")
    ap.add_argument("--amr-genes", help="comma-separated amr gene names (gene-set mode; skips main.tsv)")
    ap.add_argument("--resfinder-db-dir", default="data/resfinder_db", help="dir of ResFinder class .fsa")
    ap.add_argument("--resfinder-genes", help="comma-separated resfinder gene names (gene-set mode; skips blastn)")
    ap.add_argument("--sample-id", default=None)
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--json-only", action="store_true")
    args = ap.parse_args(argv)

    # amr side
    if args.amr_genes is not None:
        amr_genes = {g.strip() for g in args.amr_genes.split(",") if g.strip()}
    elif args.amrfinder_run and (args.amrfinder_run / "main.tsv").exists():
        amr_genes = amr_acquired_genes_from_main(str(args.amrfinder_run / "main.tsv"))
    else:
        print("ERROR: amr side needs --amr-genes OR --amrfinder-run <dir with main.tsv>", file=sys.stderr)
        return 2

    # resfinder side
    res_status = "ok"
    if args.resfinder_genes is not None:
        res_genes = {gene_of(g) if "_" in g else g for g in
                     (x.strip() for x in args.resfinder_genes.split(",")) if g}
    else:
        if not args.fasta or not args.fasta.exists():
            print("ERROR: resfinder side needs a genome FASTA (or --resfinder-genes)", file=sys.stderr)
            return 2
        res_genes = set()
        d = Path(args.resfinder_db_dir)
        dbs = sorted(d.glob("*.fsa")) if d.exists() else []
        if not dbs:
            res_status = "unavailable"
        for db in dbs:
            r = call_resistance_genes(args.fasta, db, drug_class=db.stem)
            if r["status"] != "ok":
                res_status = "unavailable"
                break
            res_genes |= {g["gene"] for g in r["genes"]}

    sample_id = args.sample_id or (args.fasta.stem if args.fasta else "concordance")
    if res_status != "ok":
        rec = {"sample_id": sample_id, "analysis": "amr_concordance", "status": "unavailable",
               "schema": "amr-concordance-v0",
               "reason": "ResFinder side unavailable (no blastn / no DB in --resfinder-db-dir)"}
        if args.out:
            Path(args.out).write_text(json.dumps(rec, indent=2), encoding="utf-8")
        print(json.dumps(rec, indent=2) if args.json_only else f"STATUS: unavailable - {rec['reason']}")
        return 3

    cmp = compare(amr_genes, res_genes)
    rec = {
        "sample_id": sample_id, "analysis": "amr_concordance", "status": "ok",
        "analysis_date": datetime.date.today().isoformat(), "schema": "amr-concordance-v0",
        **{k: cmp[k] for k in ("both", "amr_only", "resfinder_only", "n_both", "n_amr_only",
                               "n_resfinder_only", "n_amr_total", "n_resfinder_total", "agreement")},
        "caveat": ("Family-level comparison (allele-variant suffix stripped) — gene naming differs between "
                   "AMRFinder and ResFinder, so exact-name agreement understates concordance. amr_only / "
                   "resfinder_only are genes one tool's DB calls and the other's does not; neither is ground "
                   "truth. POINT mutations + efflux are excluded (acquired genes only)."),
    }
    if args.out:
        Path(args.out).write_text(json.dumps(rec, indent=2), encoding="utf-8")
    if args.json_only:
        print(json.dumps(rec, indent=2))
    else:
        print(f"sample: {sample_id}  AMR concordance (AMRFinder vs ResFinder, gene-family level)")
        print(f"  agreement (Jaccard): {cmp['agreement']}  | both={cmp['n_both']} "
              f"amr-only={cmp['n_amr_only']} resfinder-only={cmp['n_resfinder_only']}")
        print(f"  both:           {', '.join(cmp['both']) or '(none)'}")
        print(f"  amr-only:       {', '.join(cmp['amr_only']) or '(none)'}")
        print(f"  resfinder-only: {', '.join(cmp['resfinder_only']) or '(none)'}")
        print(f"  {rec['caveat']}")
        if args.out:
            print(f"\n[provenance JSON -> {args.out}]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
