"""Fully-independent arm for the pneumo gene-presence AMR cell: run OUR AMRFinder (not GPS's determinant
calls) on the GPS Poland assemblies -> our own erm/mef/tet determinant calls -> re-score vs measured AST.

This UPGRADES the honesty tier the rule-validation artifact flagged as unmeasured: the gene-presence cell's
0.961/0.932 used GPS pipeline determinant calls (external genotype). This runner produces (a) the
FULLY-INDEPENDENT number (our caller end-to-end) and (b) the OUR-vs-GPS determinant discordance -- the thing
the "AMRFinder ~= GPS BLAST" claim rested on, now MEASURED.

Per isolate: ENA ERS -> assembly (contig.fa.gz) -> AMRFinder (acquired genes, no -O needed for erm/mef/tet)
-> parse main.tsv gene symbols -> our determinant tokens -> pneumo_amr.call_drug -> vs measured disc-zone AST
+ vs GPS determinant call. Checkpointed jsonl; assembly deleted per-isolate (disk-bound). `--limit N` = pilot.

    MSYS_NO_PATHCONV=1 uv run --no-sync python scripts/pneumo_amrfinder_swap.py --limit 12
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
from dna_decode.organism_rules.pneumo_amr import call_drug                  # noqa: E402
from scripts.drug_mechanism_audit import AMRFINDER_DB, AMRFINDER_IMAGE      # noqa: E402
from scripts.pneumo_gps_quellung_validate import _fetch_assembly            # noqa: E402
from tools.docker_runner import run as docker_run                          # noqa: E402

GPS = Path("D:/dna_decode_cache/pneumo_gps")
WORK = Path("D:/dna_decode_cache/pneumo_amrfinder_swap")
DISC = {"erythromycin": {"col": "Erythromycin", "S": 21, "R": 15},
        "tetracycline": {"col": "Tetracycline", "S": 23, "R": 18}}
DET_COL = {"erythromycin": "ERY_Determinant", "tetracycline": "TET_Determinant"}
# AMRFinder gene symbols that map to each drug's gene-presence determinant tokens.
AMR_GENES = {"erythromycin": ("erm", "mef", "msr"), "tetracycline": ("tet",)}


def _zone(v):
    try:
        return float((v or "").strip())
    except ValueError:
        return None


def _measured_rs(drug, z):
    if z is None:
        return None
    bp = DISC[drug]
    return "S" if z >= bp["S"] else "R" if z <= bp["R"] else "I"


def _amrfinder_no_O(fasta: Path, out_dir: Path, timeout: float = 600) -> Path:
    """AMRFinder acquired-gene scan WITHOUT -O (erm/mef/tet are acquired; -O only adds point mutations).
    Mirrors drug_mechanism_audit._run_amrfinder's DB-mount logic (resolve the `latest` symlink)."""
    out_dir.mkdir(parents=True, exist_ok=True)
    main_out = out_dir / "main.tsv"
    db_latest = Path(AMRFINDER_DB) / "latest"
    real_db = db_latest.resolve() if db_latest.exists() else db_latest
    docker_run(
        AMRFINDER_IMAGE,
        ["amrfinder", "-n", f"/in/{fasta.name}", "--database", "/db/latest", "-o", "/out/main.tsv"],
        mounts={str(fasta.parent): "/in:ro", str(real_db): "/db/latest:ro", str(out_dir): "/out"},
        timeout=timeout,
    )
    return main_out


def _our_determinant_tokens(main_tsv: Path) -> list[str]:
    """Parse AMRFinder main.tsv gene symbols -> normalized tokens (erm(B)->ermB) for call_drug matching."""
    if not main_tsv.exists() or main_tsv.stat().st_size == 0:
        return []
    toks = []
    with main_tsv.open(encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            sym = (row.get("Gene symbol") or row.get("Element symbol") or "").strip()
            if sym:
                toks.append(re.sub(r"[^A-Za-z0-9]", "", sym))   # erm(B) -> ermB ; tet(M) -> tetM
    return toks


def cohort():
    m3 = {r["ERR"].strip(): r for r in csv.DictReader(open(GPS / "sd_3.csv", encoding="utf-8")) if r.get("ERR")}
    m4 = {r["Sample_ID"].strip(): r for r in csv.DictReader(open(GPS / "sd_4.csv", encoding="utf-8")) if r.get("Sample_ID")}
    out = []
    for k in set(m3) & set(m4):
        ers = (m3[k].get("ERS") or "").strip()
        # keep isolates with >=1 measured gene-presence drug
        meas = {d: _measured_rs(d, _zone(m3[k].get(DISC[d]["col"]))) for d in DISC}
        if not ers or all(v in (None, "I") for v in meas.values()):
            continue
        out.append({"err": k, "ers": ers, "measured": meas,
                    "gps_det": {d: (m4[k].get(DET_COL[d]) or "") for d in DISC}})
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--limit", type=int, default=0, help="pilot: at most N new isolates (0=all)")
    ap.add_argument("--out", type=Path, default=REPO / "wiki" / "pneumo_amr_amrfinder_swap_pilot.json")
    a = ap.parse_args(argv)

    results = WORK / "swap_results.jsonl"
    results.parent.mkdir(parents=True, exist_ok=True)
    done = ({json.loads(x)["err"] for x in results.read_text().splitlines() if x.strip()}
            if results.exists() else set())
    isos = cohort()
    todo = [c for c in isos if c["err"] not in done]
    if a.limit:
        todo = todo[:a.limit]
    print(f"[amrfinder-swap] {len(isos)} gene-presence isolates | {len(done)} done | running {len(todo)}")

    asm_dir = WORK / "asm"
    asm_dir.mkdir(parents=True, exist_ok=True)   # curl -o needs the dir to exist (else silent fetch fail)
    for n, iso in enumerate(todo, 1):
        fa = asm_dir / f"{iso['ers']}.fna"
        rec = {"err": iso["err"], "ers": iso["ers"]}
        try:
            if not _fetch_assembly(iso["ers"], fa):
                rec["status"] = "assembly_unavailable"
            else:
                main_tsv = _amrfinder_no_O(fa, WORK / "amr" / iso["ers"])
                our_toks = _our_determinant_tokens(main_tsv)
                rec["status"] = "ok"
                rec["our_tokens"] = [t for t in our_toks if any(g in t.lower() for g in ("erm", "mef", "msr", "tet"))]
                for d in DISC:
                    our = call_drug(d, our_toks)
                    gps = call_drug(d, [iso["gps_det"][d]])
                    rec[d] = {"measured": iso["measured"][d], "our_pred": our.prediction if our else None,
                              "gps_pred": gps.prediction if gps else None}
        except Exception as e:
            rec["status"] = "error"; rec["reason"] = f"{type(e).__name__}: {str(e)[:120]}"
        finally:
            import shutil
            fa.unlink(missing_ok=True)                                  # delete the assembly (disk-bound)
            fa.with_suffix(".fa.gz").unlink(missing_ok=True)
            shutil.rmtree(WORK / "amr" / iso["ers"], ignore_errors=True)  # delete AMRFinder scratch
        with results.open("a", encoding="utf-8") as f:
            f.write(json.dumps(rec) + "\n")
        print(f"  [{n}/{len(todo)}] {iso['err']}: {rec['status']} our={rec.get('our_tokens')}")

    # aggregate: fully-independent re-score (our vs measured) + our-vs-GPS determinant agreement
    recs = [json.loads(x) for x in results.read_text().splitlines() if x.strip()]
    agg = {}
    for d in DISC:
        ind_tp = ind_tn = ind_fp = ind_fn = 0
        det_agree = det_total = 0
        for r in recs:
            cell = r.get(d)
            if not cell:
                continue
            m, op, gp = cell["measured"], cell["our_pred"], cell["gps_pred"]
            if op is not None and gp is not None:
                det_total += 1; det_agree += (op == gp)
            if m in ("R", "S") and op in ("R", "S"):
                if m == "R" and op == "R": ind_tp += 1
                elif m == "S" and op == "S": ind_tn += 1
                elif m == "S" and op == "R": ind_fp += 1
                else: ind_fn += 1
        n = ind_tp + ind_tn + ind_fp + ind_fn
        agg[d] = {"fully_independent": {"n": n, "acc": round((ind_tp + ind_tn) / n, 3) if n else None,
                                        "tp": ind_tp, "fp": ind_fp, "tn": ind_tn, "fn": ind_fn},
                  "our_vs_gps_determinant_agreement": {"n": det_total,
                                                       "agree": round(det_agree / det_total, 3) if det_total else None}}
    out = {"schema": "pneumo-amr-amrfinder-swap-v1",
           "label": "FULLY-INDEPENDENT: OUR AMRFinder (not GPS calls) -> our rule -> vs wet-lab disc AST",
           "n_attempted": len(recs), "drugs": agg,
           "note": "our_vs_gps_determinant_agreement MEASURES the previously-unmeasured 'AMRFinder ~= GPS BLAST' "
                   "claim. High agreement => the gene-presence rule-validation number transfers to fully-independent."}
    a.out.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"[amrfinder-swap] AGG: {json.dumps(agg)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
