"""Full ktype validation: run the wzi caller on the KlebNET-GSP cohort genomes and score vs the MEASURED
serological K-type. This is the genuine wzi-caller-vs-serology number the report card's ceiling pointed at.

PIPELINE per isolate (reuses the project's targeted-map machinery -- single-locus, no whole-genome assembly):
  ERR run accession -> fetch_reads (prefetch+fasterq-dump) -> map_erg11(wzi_ref) (minimap2 -ax sr -> samtools
  consensus over the wzi locus) -> call_ktype(consensus) -> predicted KL -> compare to Serological_Ktype.

The reads are deleted after each isolate (the 447-isolate set would otherwise be ~100+ GB). CHECKPOINTED to
a jsonl (restartable; skips done accessions). `--limit N` runs a pilot. Stage work on D: (large). Docker req.

COHORT: only the 447 isolates with an ERR accession are fetchable (283 'TBD' + 3 unavailable are skipped --
reported, not hidden). CONCORDANCE = predicted KL# vs serological K# (naive numeric match, same honest caveat
as the report card's ceiling: the curated KL<->K equivalence would lift it; this is a lower bound).

Run: `uv run --with openpyxl python scripts/ktype_cohort_validate.py --xlsx <suppl.xlsx> --limit 8`
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from dna_decode.ktype.runner import call_ktype                                   # noqa: E402
from scripts.assemble_sra_cohort import (  # noqa: E402
    MINIMAP2_IMAGE, SAMTOOLS_IMAGE, SRA_IMAGE, fetch_reads, map_erg11,
)

WORK = Path("D:/dna_decode_cache/ktype_build")
WZI_REF = WORK / "wziref" / "wzi_ref.fasta"
DB_DIR = REPO / "data" / "ktype_db"
RESULTS = WORK / "cohort_results.jsonl"


def cohort(xlsx: Path) -> list[dict]:
    import openpyxl
    ws = openpyxl.load_workbook(xlsx, read_only=True)["Supplementary_Table1"]
    rows = list(ws.iter_rows(values_only=True))
    ix = {h: i for i, h in enumerate(rows[0])}
    def c(r, k):
        j = ix.get(k); return r[j] if (j is not None and j < len(r)) else None
    out = []
    for r in rows[1:]:
        if str(c(r, "Included_in_report")).lower() != "include":
            continue
        acc, sero = c(r, "WGS_accession"), c(r, "Serological_Ktype")
        if acc and sero and str(acc).startswith("ERR"):
            out.append({"accession": str(acc).strip(), "serological": str(sero).strip(),
                        "kaptive_kl": str(c(r, "K_locus_Best_match_locus") or "")})
    return out


def _concordant(predicted_kl: str | None, serological: str) -> bool | None:
    if not predicted_kl:
        return None
    pk = predicted_kl.replace("KL", "").replace("K", "").split()[0]
    sero_nums = {x.strip().replace("K", "") for x in serological.split(",")}
    return pk in sero_nums


def _done() -> set[str]:
    if not RESULTS.exists():
        return set()
    return {json.loads(l)["accession"] for l in RESULTS.read_text().splitlines() if l.strip()}


def run_one(iso: dict, threads: int, timeout: float) -> dict:
    acc = iso["accession"]
    reads_dir = WORK / "reads" / acc
    work_dir = WORK / "map" / acc
    cons = work_dir / "wzi_cons.fna"
    try:
        r1, r2 = fetch_reads(acc, reads_dir, SRA_IMAGE, timeout)
        map_erg11(acc, r1, r2, work_dir, cons, WZI_REF, MINIMAP2_IMAGE, SAMTOOLS_IMAGE,
                  threads, min_depth=3, timeout=timeout)
        call = call_ktype(cons, DB_DIR)
        rec = {"accession": acc, "serological": iso["serological"], "kaptive_kl": iso["kaptive_kl"],
               "status": call["status"], "wzi_allele": call.get("wzi_allele"),
               "predicted_kl": call.get("kl_type"),
               "concordant": _concordant(call.get("kl_type"), iso["serological"])}
    except Exception as e:
        rec = {"accession": acc, "serological": iso["serological"], "status": "error",
               "predicted_kl": None, "concordant": None, "reason": f"{type(e).__name__}: {str(e)[:120]}"}
    finally:
        shutil.rmtree(reads_dir, ignore_errors=True)   # bound disk: reads are large + transient
        shutil.rmtree(work_dir, ignore_errors=True)
    return rec


def aggregate() -> dict:
    recs = [json.loads(l) for l in RESULTS.read_text().splitlines() if l.strip()] if RESULTS.exists() else []
    called = [r for r in recs if r.get("predicted_kl")]
    scored = [r for r in called if r.get("concordant") is not None]
    conc = sum(1 for r in scored if r["concordant"])
    return {"n_attempted": len(recs), "n_called": len(called), "n_scored": len(scored),
            "n_concordant": conc, "concordance": round(conc / len(scored), 3) if scored else None,
            "n_error": sum(1 for r in recs if r.get("status") == "error")}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--xlsx", type=Path, default=WORK / "zenodo_suppl.xlsx")
    ap.add_argument("--limit", type=int, default=0, help="run at most N new isolates (0 = all 447)")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--timeout", type=float, default=1800)
    a = ap.parse_args(argv)
    if not WZI_REF.exists() or not (DB_DIR / "wzi.fasta").exists():
        print(f"ERROR: missing wzi ref ({WZI_REF}) or DB ({DB_DIR})", file=sys.stderr); return 2
    isos = cohort(a.xlsx)
    done = _done()
    todo = [i for i in isos if i["accession"] not in done]
    if a.limit:
        todo = todo[:a.limit]
    print(f"[ktype-cohort] {len(isos)} ERR isolates | {len(done)} done | running {len(todo)} this batch")
    RESULTS.parent.mkdir(parents=True, exist_ok=True)
    for n, iso in enumerate(todo, 1):
        rec = run_one(iso, a.threads, a.timeout)
        with RESULTS.open("a", encoding="utf-8") as f:
            f.write(json.dumps(rec) + "\n")
        print(f"  [{n}/{len(todo)}] {rec['accession']}: {rec.get('status')} "
              f"wzi={rec.get('wzi_allele')} KL={rec.get('predicted_kl')} "
              f"sero={rec['serological']} concordant={rec.get('concordant')}")
    agg = aggregate()
    print(f"[ktype-cohort] AGG: {agg}")
    (REPO / "wiki" / "ktype_cohort_validation.json").write_text(
        json.dumps({"schema": "ktype-cohort-validation-v1",
                    "label": "wzi-caller vs MEASURED serological K-type (KlebNET-GSP ERR cohort)",
                    "n_total_ERR_cohort": len(isos), "aggregate": agg,
                    "caveat": "naive KL#==K# concordance (lower bound; curated KL<->K equivalence lifts it). "
                              "wzi single-gene v0 (~94% of full Kaptive). Reads-based targeted wzi mapping."},
                   indent=2), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
