"""Resistance-gene × plasmid co-localization — CLI (console `dna-coloc`; also `dna-decode coloc`).

Runs ResFinder (acquired genes) + PlasmidFinder (replicons) on a genome with the engine's positions mode,
then reports which acquired resistance genes sit on the SAME contig as a plasmid replicon — i.e. likely
plasmid-borne. Offline-safe (status 'unavailable', exit 3). Composes plasmid + resfinder; no new DB.

    dna-coloc assembly.fna --sample-id MY_STRAIN
    dna-coloc X.fna --plasmid-db ... --resfinder-db-dir ...
"""
from __future__ import annotations

import argparse
import datetime
import json
import sys
from pathlib import Path

from dna_decode.colocalization.core import colocalize
from dna_decode.plasmid.runner import (
    PLASMID_COVERAGE_THRESHOLD,
    PLASMID_IDENTITY_THRESHOLD,
    replicon_family,
)
from dna_decode.resfinder.runner import (
    RES_COVERAGE_THRESHOLD,
    RES_IDENTITY_THRESHOLD,
    gene_of,
)
from dna_decode.typing.blast_caller import call_alleles


def _called_with_contig(db, *, identity, coverage):
    """call_alleles with positions -> [(allele_id, contig)] for called hits. Returns (list, status)."""
    res = call_alleles(db[0], db[1], identity_threshold=identity, coverage_threshold=coverage,
                       with_positions=True)
    if res["status"] != "ok":
        return [], res
    out = [(aid, h.get("contig")) for aid, h in res["per_allele"].items() if h["called"]]
    return out, res


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="dna-coloc",
                                 description="Resistance-gene x plasmid-replicon co-localization (same-contig)")
    ap.add_argument("fasta", type=Path, help="genome assembly FASTA")
    ap.add_argument("--plasmid-db", default="data/plasmidfinder_db/enterobacteriales.fsa")
    ap.add_argument("--resfinder-db-dir", default="data/resfinder_db")
    ap.add_argument("--sample-id", default=None)
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--json-only", action="store_true")
    args = ap.parse_args(argv)

    if not args.fasta.exists():
        print(f"ERROR: assembly FASTA not found: {args.fasta}", file=sys.stderr)
        return 2
    sample_id = args.sample_id or args.fasta.stem

    # plasmid replicons (best hit per allele, with contig)
    rep_hits, prep = _called_with_contig((args.fasta, args.plasmid_db),
                                         identity=PLASMID_IDENTITY_THRESHOLD,
                                         coverage=PLASMID_COVERAGE_THRESHOLD)
    if prep["status"] != "ok":
        return _unavail(args, sample_id, f"plasmid side: {prep.get('reason')}")
    replicons = [{"replicon": replicon_family(aid), "contig": c} for aid, c in rep_hits]

    # resfinder acquired genes across every class DB (with contig)
    d = Path(args.resfinder_db_dir)
    dbs = sorted(d.glob("*.fsa")) if d.exists() else []
    if not dbs:
        return _unavail(args, sample_id, f"resfinder side: no class .fsa in {args.resfinder_db_dir}")
    genes = []
    for db in dbs:
        hits, rr = _called_with_contig((args.fasta, db), identity=RES_IDENTITY_THRESHOLD,
                                       coverage=RES_COVERAGE_THRESHOLD)
        if rr["status"] != "ok":
            return _unavail(args, sample_id, f"resfinder side: {rr.get('reason')}")
        genes.extend({"gene": gene_of(aid), "contig": c} for aid, c in hits)

    co = colocalize(genes, replicons)
    rec = {
        "sample_id": sample_id, "analysis": "amr_plasmid_colocalization", "status": "ok",
        "analysis_date": datetime.date.today().isoformat(), "schema": "amr-plasmid-coloc-v0",
        **co,
        "caveat": ("Same-contig co-location of an acquired resistance gene + a plasmid replicon SUGGESTS the "
                   "gene is plasmid-borne (a circularized plasmid often assembles to one contig) but does NOT "
                   "prove a single replicon carries it — fragmented/chimeric assemblies + multi-replicon "
                   "contigs break the inference. 'chromosomal_or_unplaced' = no replicon on that contig "
                   "(could be chromosomal OR a replicon-less plasmid contig). NOT a clinical tool."),
    }
    if args.out:
        Path(args.out).write_text(json.dumps(rec, indent=2), encoding="utf-8")
    if args.json_only:
        print(json.dumps(rec, indent=2))
    else:
        s = co["summary"]
        print(f"sample: {sample_id}  AMR x plasmid co-localization")
        print(f"  {s['n_plasmid_borne']}/{s['n_genes']} acquired genes co-located with a plasmid replicon "
              f"| {s['n_replicons']} replicon hit(s)")
        for g in co["gene_calls"]:
            tag = f"PLASMID-BORNE ({', '.join(g['replicons_on_contig'])})" if g["plasmid_borne"] \
                else "chromosomal/unplaced"
            print(f"  {g['gene']:22} [{g['contig']}]  {tag}")
        if not co["gene_calls"]:
            print("  (no acquired resistance genes detected)")
        print(f"  {rec['caveat']}")
        if args.out:
            print(f"\n[provenance JSON -> {args.out}]")
    return 0


def _unavail(args, sample_id, reason):
    rec = {"sample_id": sample_id, "analysis": "amr_plasmid_colocalization", "status": "unavailable",
           "schema": "amr-plasmid-coloc-v0", "reason": reason}
    if args.out:
        Path(args.out).write_text(json.dumps(rec, indent=2), encoding="utf-8")
    print(json.dumps(rec, indent=2) if args.json_only else f"STATUS: unavailable - {reason}",
          file=sys.stderr if not args.json_only else sys.stdout)
    return 3


if __name__ == "__main__":
    raise SystemExit(main())
