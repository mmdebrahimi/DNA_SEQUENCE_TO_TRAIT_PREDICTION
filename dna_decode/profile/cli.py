"""Unified genome profile — in-package CLI (console `dna-profile`; also `dna-decode profile`).

Runs the assembly-FASTA decoders (pathotype + plasmid + serotype + resfinder) on one genome and emits a
single unified report. Each section is independent + offline-safe (a missing DB / absent blastn marks just
that section 'unavailable'; the others still run).

    dna-profile assembly.fna --sample-id MY_STRAIN
    dna-profile X.fna --plasmid-db ... --serotype-db ... --resfinder-db-dir ... --pathotype-db ...
"""
from __future__ import annotations

import argparse
import datetime
import json
import sys
from pathlib import Path


def _pathotype(fasta, db):
    if not Path(db).exists():
        return {"status": "unavailable", "reason": f"VirulenceFinder DB not found at {db}"}
    try:
        from dna_decode.pathotype.detect import assembly_qc, build_vf_index, detect, parse_fasta
        from dna_decode.pathotype.resolve import resolve_call
        contigs = parse_fasta(fasta)
        det = detect(contigs, build_vf_index(db))
        qc = assembly_qc(contigs)
        call = resolve_call(det["cluster_profile"],
                            partial_clusters=frozenset(det.get("partial_clusters", [])),
                            qc_pass=bool(qc.get("pass", True)))
        return {"status": "ok", "pathotype": call.get("primary"),
                "confidence_tier": call.get("confidence_tier"), "rule_id": call.get("rule_id")}
    except Exception as e:  # offline-safe: never let one section crash the profile
        return {"status": "error", "reason": f"{type(e).__name__}: {str(e)[:120]}"}


def _plasmid(fasta, db):
    from dna_decode.plasmid.runner import call_replicons
    r = call_replicons(fasta, db)
    if r["status"] != "ok":
        return {"status": r["status"], "reason": r.get("reason")}
    return {"status": "ok", "replicons": [x["replicon"] for x in r["replicons"]]}


def _serotype(fasta, db):
    from dna_decode.serotype.runner import call_serotype
    r = call_serotype(fasta, db)
    if r["status"] != "ok":
        return {"status": r["status"], "reason": r.get("reason")}
    return {"status": "ok", "serotype": r.get("serotype"), "o_antigen": r.get("o_antigen"),
            "h_antigen": r.get("h_antigen")}


def _pointfinder(fasta, db_dir):
    from dna_decode.pointfinder.runner import call_point_mutations, parse_overview
    d = Path(db_dir)
    ov = d / "resistens-overview.txt"
    refs = {g: d / f"{g}.fsa" for g in ("gyrA", "parC", "gyrB", "parE")}
    if not ov.exists() or not any(p.exists() for p in refs.values()):
        return {"status": "unavailable", "reason": f"PointFinder DB not found at {db_dir}"}
    r = call_point_mutations(fasta, {g: p for g, p in refs.items() if p.exists()}, parse_overview(ov))
    if r["status"] != "ok":
        return {"status": r["status"], "reason": r.get("reason")}
    return {"status": "ok", "mutations": [m["mutation"] for m in r["mutations"]],
            "resistances": r["resistances"]}


def _resfinder(fasta, db_dir):
    from dna_decode.resfinder.runner import call_resistance_genes
    d = Path(db_dir)
    dbs = sorted(d.glob("*.fsa")) if d.exists() else []
    if not dbs:
        return {"status": "unavailable", "reason": f"no ResFinder class .fsa in {db_dir}"}
    genes = []
    for db in dbs:
        r = call_resistance_genes(fasta, db, drug_class=db.stem)
        if r["status"] != "ok":
            return {"status": "unavailable", "reason": r.get("reason")}
        genes.extend(g["gene"] for g in r["genes"])
    return {"status": "ok", "genes": sorted(set(genes))}


def _amr(main_tsv, drugs, organism, threshold):
    """Per-drug R/S calls from a (cached or freshly-run) AMRFinder main.tsv, each carrying its INLINE trust
    badge (the move-1 trust-surface). Offline-safe: a missing main.tsv marks the section 'unavailable'."""
    if main_tsv is None or not Path(main_tsv).exists():
        return {"status": "unavailable",
                "reason": "no AMRFinder main.tsv (pass --amrfinder-run DIR for a cached run, or "
                          "--run-amrfinder for a Docker run); amr is the only Docker/cached-run section"}
    from dna_decode.data.trust_surface import trust_block
    from dna_decode.eval.amr_rules import call_resistance
    calls = {}
    for drug in drugs:
        try:
            c = call_resistance(Path(main_tsv), drug, threshold, organism=organism)
        except Exception as e:  # one drug failing never sinks the section
            calls[drug] = {"status": "error", "reason": f"{type(e).__name__}: {str(e)[:80]}"}
            continue
        calls[drug] = {"prediction": c["prediction"], "confidence": c.get("confidence"),
                       "n_determinants": c.get("n_determinants"),
                       "validation": trust_block(drug, organism)}
    return {"status": "ok", "organism": organism, "calls": calls}


def _resolve_amr_organism(fasta, explicit, do_identify: bool) -> dict:
    """Decide which organism the AMR section routes on, and record HOW it was decided.

    THREE sources, never collapsed -- a reader must be able to tell a MEASUREMENT from an ASSUMPTION
    from an infrastructure FAULT:
      * `explicit`           -- the user passed --amr-organism. ALWAYS WINS over identification:
                                user authority outranks inference, and an explicit value is also the
                                escape hatch when the router abstains.
      * `identified`         -- `dna-identify` named a supported organism. `organism_assumed` becomes
                                FALSE, because it no longer is one.
      * `assumed_*`          -- the E. coli assumption STANDS, with the reason it stands. An
                                abstention and an unavailable reference are recorded SEPARATELY: the
                                first is a correct biological answer ("not a supported organism"), the
                                second is a fault, and conflating them would let a wedged Docker mount
                                read as a finding.

    Identification is OPT-IN because it needs Docker, and `profile`'s auto-run path (used by
    `dna-decode decode`) is offline-safe. Making it automatic would break that guarantee.
    """
    out = {"organism_requested": explicit, "organism_identify": None}
    if explicit:
        out.update(organism_used=explicit, organism_assumed=False, organism_source="explicit")
        return out
    if not do_identify:
        out.update(organism_used="Escherichia", organism_assumed=True, organism_source="assumed_default")
        return out
    try:
        from dna_decode.data.organism_vocab import routing_token
        from dna_decode.identify.core import decide
        from dna_decode.identify.runner import ReferenceUnavailable, query
        from dna_decode.identify.thresholds import frozen
        import json as _json
        repo = Path(__file__).resolve().parent.parent.parent
        mans = sorted((repo / "wiki").glob("identify_reference_manifest_*.json"))
        if not mans:
            raise ReferenceUnavailable(
                "no identify reference manifest in wiki/; build one with "
                "`uv run python scripts/build_identify_reference.py`")
        labels = {a: o for o, accs in
                  _json.loads(mans[-1].read_text(encoding="utf-8"))["accessions_by_organism"].items()
                  for a in accs}
        res = query(Path(fasta))
        call = decide(res.hits, frozen(), lambda a: labels.get(a), routing_token)
    except ReferenceUnavailable as exc:
        out.update(organism_used="Escherichia", organism_assumed=True,
                   organism_source="assumed_identify_unavailable",
                   organism_identify={"status": "reference_unavailable", "detail": str(exc)[:300]})
        return out
    except Exception as exc:  # any other fault: still a FAULT, never a biological abstention
        out.update(organism_used="Escherichia", organism_assumed=True,
                   organism_source="assumed_identify_failed",
                   organism_identify={"status": "failed",
                                      "detail": f"{type(exc).__name__}: {exc}"[:300]})
        return out

    if call.abstained:
        out.update(organism_used="Escherichia", organism_assumed=True,
                   organism_source="assumed_identify_abstained",
                   organism_identify={"status": "abstained", "reason": call.reason.value,
                                      "nearest_distance": call.nearest_distance})
        return out
    if not call.routing_token:
        # a supported organism with NO AMRFinder -O value (fungal / TB): it routes to its own engine,
        # so the AMR section must NOT pretend an -O value exists
        out.update(organism_used="Escherichia", organism_assumed=True,
                   organism_source="assumed_identified_but_not_amrfinder_routable",
                   organism_identify={"status": "identified", "organism": call.organism,
                                      "routing_token": None,
                                      "note": "AMRFinder has no -O value for this organism"})
        return out
    out.update(organism_used=call.routing_token, organism_assumed=False,
               organism_source="identified",
               organism_identify={"status": "identified", "organism": call.organism,
                                  "routing_token": call.routing_token,
                                  "nearest_distance": call.nearest_distance})
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="dna-profile",
                                 description="Unified genome profile — run all assembly-FASTA decoders")
    ap.add_argument("fasta", type=Path, help="genome assembly FASTA")
    ap.add_argument("--sample-id", default=None)
    ap.add_argument("--pathotype-db", default="data/virulencefinder_db/virulence_ecoli.fsa")
    ap.add_argument("--plasmid-db", default="data/plasmidfinder_db/enterobacteriales.fsa")
    ap.add_argument("--serotype-db", default="data/serotypefinder_db/serotypefinder.fsa")
    ap.add_argument("--resfinder-db-dir", default="data/resfinder_db")
    ap.add_argument("--pointfinder-db-dir", default="data/pointfinder_db/escherichia_coli")
    # AMR R/S section (the only Docker/cached-run section; offline-safe -> 'unavailable' without a source)
    ap.add_argument("--amrfinder-run", type=Path, default=None,
                    help="dir of a CACHED AMRFinder run (uses <dir>/main.tsv; no Docker)")
    ap.add_argument("--run-amrfinder", action="store_true",
                    help="run AMRFinder on the genome via Docker (needs Docker + an AMRFinder DB)")
    ap.add_argument("--amr-organism", default=None,
                    help="organism for AMR routing + the trust badge. If omitted, E. coli is ASSUMED and the "
                         "output is stamped organism_assumed=true (routing + badge are E. coli-oriented -- "
                         "pass this explicitly for a non-E. coli genome)")
    ap.add_argument("--identify-organism", action="store_true",
                    help="resolve the organism from the genome with `dna-identify` instead of ASSUMING "
                         "E. coli (needs Docker + the reference sketch). --amr-organism still wins; an "
                         "abstention keeps the disclosed assumption and records why")
    ap.add_argument("--amr-drugs", default="ciprofloxacin,ceftriaxone,tetracycline,gentamicin",
                    help="comma-separated drugs for the AMR section")
    ap.add_argument("--amr-threshold", type=int, default=None)
    ap.add_argument("--amrfinder-db", type=Path, default=Path("data/amrfinder_db"))
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--json-only", action="store_true")
    args = ap.parse_args(argv)

    if not args.fasta.exists():
        print(f"ERROR: assembly FASTA not found: {args.fasta}", file=sys.stderr)
        return 2
    sample_id = args.sample_id or args.fasta.stem

    # Organism SELF-DISCLOSURE: an omitted --amr-organism means E. coli is ASSUMED. Both the routing AND the
    # trust badge are E. coli-oriented, so the output must NOT look identical to an explicit user choice.
    # With --identify-organism the assumption can be REPLACED by a measurement (`dna-identify`), which is
    # step 1 of the pipeline feeding step 2 instead of a human supplying it.
    org_res = _resolve_amr_organism(args.fasta, args.amr_organism, args.identify_organism)
    amr_org = org_res["organism_used"]
    amr_assumed = org_res["organism_assumed"]

    # Resolve an AMRFinder main.tsv source for the AMR section (cached run > Docker run > none/offline-safe).
    amr_main_tsv = None
    if args.amrfinder_run is not None:
        amr_main_tsv = args.amrfinder_run / "main.tsv"
    elif args.run_amrfinder:
        try:
            from dna_decode.amr.cli import _run_amrfinder_for_genome
            run_dir = _run_amrfinder_for_genome(args.fasta, sample_id, Path("data/amrfinder_runs"),
                                                args.amrfinder_db, organism=amr_org)
            amr_main_tsv = run_dir / "main.tsv"
        except Exception as e:  # offline-safe: a failed Docker run leaves the section 'unavailable'
            print(f"[warn] AMRFinder Docker run failed ({type(e).__name__}: {e}); amr section unavailable",
                  file=sys.stderr)
    amr_drugs = [d.strip() for d in args.amr_drugs.split(",") if d.strip()]

    amr_section = _amr(amr_main_tsv, amr_drugs, amr_org, args.amr_threshold)
    amr_section["organism_requested"] = args.amr_organism   # null when omitted
    amr_section["organism_used"] = amr_org
    amr_section["organism_assumed"] = amr_assumed
    amr_section["organism_source"] = org_res["organism_source"]
    amr_section["organism_identify"] = org_res["organism_identify"]
    decoders = {
        "pathotype": _pathotype(args.fasta, args.pathotype_db),
        "serotype": _serotype(args.fasta, args.serotype_db),
        "plasmid": _plasmid(args.fasta, args.plasmid_db),
        "resfinder": _resfinder(args.fasta, args.resfinder_db_dir),
        "pointfinder": _pointfinder(args.fasta, args.pointfinder_db_dir),
        "amr": amr_section,
    }
    n_ok = sum(1 for d in decoders.values() if d.get("status") == "ok")
    rec = {
        "sample_id": sample_id, "analysis": "genome_profile",
        "analysis_date": datetime.date.today().isoformat(), "schema": "genome-profile-v0",
        "decoders_ok": n_ok, "decoders_total": len(decoders),
        "decoders": decoders,
        "caveat": ("Each section is an independent deterministic curated-DB caller; 'unavailable' = that DB / "
                   "blastn / AMRFinder source is absent, not a negative result. E. coli-oriented (pathotype/"
                   "serotype). The amr section carries each R/S call's INLINE validation tier (trust badge); "
                   "it needs a cached (--amrfinder-run) or Docker (--run-amrfinder) AMRFinder source. "
                   "NOT a clinical tool."),
    }
    if args.out:
        Path(args.out).write_text(json.dumps(rec, indent=2), encoding="utf-8")
    if args.json_only:
        print(json.dumps(rec, indent=2))
    else:
        print(f"sample: {sample_id}  unified genome profile  ({n_ok}/{len(decoders)} decoders ran)")
        p = decoders["pathotype"]
        print(f"  pathotype: {p.get('pathotype') if p['status']=='ok' else '['+p['status']+']'}")
        s = decoders["serotype"]
        print(f"  serotype:  {s.get('serotype') if s['status']=='ok' else '['+s['status']+']'}")
        pl = decoders["plasmid"]
        print(f"  plasmid:   {', '.join(pl.get('replicons', [])) or '(none)' if pl['status']=='ok' else '['+pl['status']+']'}")
        rf = decoders["resfinder"]
        print(f"  resfinder: {', '.join(rf.get('genes', [])) or '(none)' if rf['status']=='ok' else '['+rf['status']+']'}")
        pf = decoders["pointfinder"]
        print(f"  pointfinder: {', '.join(pf.get('mutations', [])) or '(none)' if pf['status']=='ok' else '['+pf['status']+']'}")
        am = decoders["amr"]
        if am["status"] == "ok":
            assumed_note = "  (organism ASSUMED E. coli -- pass --amr-organism for a non-E. coli genome)" \
                if am.get("organism_assumed") else ""
            print(f"  amr ({am['organism']}):{assumed_note}")
            for drug, c in am["calls"].items():
                if "prediction" in c:
                    v = c["validation"]
                    badge = f"{v['tier']}{(' ' + v['headline']) if v.get('headline') else ''}"
                    print(f"    {drug:<14} {c['prediction']:<14} [{c.get('n_determinants')} det]  | {badge}")
                else:
                    print(f"    {drug:<14} [{c.get('status')}]")
        else:
            print(f"  amr: [{am['status']}]")
        print(f"  {rec['caveat']}")
        if args.out:
            print(f"\n[provenance JSON -> {args.out}]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
