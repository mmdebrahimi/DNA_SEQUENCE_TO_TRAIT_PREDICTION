"""Independent cross-check on the meropenem FN species claim, via dna-identify (Docker Mash).

THE KILL-TEST. The 2026-10-05 memo publishes, as established fact, that all 16 false negatives in the
`Klebsiella x meropenem` cell are *Klebsiella aerogenes* -- evidence being each assembly's own GenBank
ORGANISM line plus an `ampC-Kaer` determinant fingerprint. Both of those come from the SAME annotation
provenance. `dna-identify` is a different route: Mash sketch distance against a 14-organism closed-set
reference that does NOT contain K. aerogenes and whose thresholds.py records K. aerogenes landing at
0.10849+, above MAX_DISTANCE 0.1036.

PREDICTION, registered before reading any output:
  * the 16 FN should ABSTAIN (out of the supported set) -- or, if called, NOT as klebsiella_pneumoniae
  * the K. pneumoniae true positives should be CALLED klebsiella_pneumoniae

A FN that comes back `klebsiella_pneumoniae` KILLS the published claim for that isolate.
Controls are run in the same pass so a uniformly-abstaining router cannot pass by accident.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[0]))

REPO = Path("C:/Users/Farshad/PythonProjects/dna_decode")
REFSEQ = REPO / "data" / "raw" / "klebsiella_provdisjoint_meropenem" / "refseq"


def call_identify(acc: str) -> dict:
    from dna_decode.identify.cli import main as icli_main
    import contextlib, io
    buf = io.StringIO()
    g = REFSEQ / acc / "genome.fna"
    if not g.exists():
        return {"acc": acc, "error": "genome.fna absent"}
    try:
        with contextlib.redirect_stdout(buf):
            rc = icli_main(["--genome-fasta", str(g), "--json"])
    except SystemExit as e:
        rc = e.code
    except Exception as e:  # noqa: BLE001 - report, never crash the sweep
        return {"acc": acc, "error": f"{type(e).__name__}: {e}"}
    txt = buf.getvalue()
    try:
        rec = json.loads(txt[txt.index("{"): txt.rindex("}") + 1])
    except Exception:
        return {"acc": acc, "error": f"unparsed rc={rc}: {txt.strip()[:160]}"}
    return {"acc": acc, "rc": rc,
            "organism": rec.get("organism"),
            "abstained": rec.get("abstained"),
            "reason": rec.get("abstain_reason"),
            "dist": rec.get("nearest_distance"),
            "nearest": rec.get("nearest_organism") or rec.get("nearest_reference")}


def main() -> int:
    fn = [a for a in (REPO / "fn_accs.tmp").read_text().split() if a]
    tp = [a for a in (REPO / "tp_accs.tmp").read_text().split() if a]
    limit = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    if limit:
        fn, tp = fn[:limit], tp[:limit]

    out = {"false_negatives": [], "true_positive_controls": []}
    for group, accs in (("false_negatives", fn), ("true_positive_controls", tp)):
        for a in accs:
            r = call_identify(a)
            out[group].append(r)
            print(f"  [{group[:2].upper()}] {a}  organism={r.get('organism')} "
                  f"abstained={r.get('abstained')} reason={r.get('reason')} "
                  f"dist={r.get('dist')}{' ERR=' + r['error'] if r.get('error') else ''}", flush=True)

    # ---- the verdict, computed not narrated ----
    fnr = out["false_negatives"]
    killed = [r for r in fnr if r.get("organism") == "klebsiella_pneumoniae"]
    abst = [r for r in fnr if r.get("abstained")]
    errs = [r for r in fnr if r.get("error")]
    ctrl_ok = [r for r in out["true_positive_controls"]
               if r.get("organism") == "klebsiella_pneumoniae"]
    ctrl_n = len(out["true_positive_controls"])

    print("\n--- verdict ---")
    print(f"FN scored={len(fnr) - len(errs)} errors={len(errs)}")
    print(f"FN abstained={len(abst)}  FN called klebsiella_pneumoniae={len(killed)}")
    print(f"controls called klebsiella_pneumoniae={len(ctrl_ok)}/{ctrl_n}")

    if ctrl_n and not ctrl_ok:
        verdict = "INDETERMINATE_CONTROLS_FAILED"
        note = ("the K. pneumoniae controls did NOT resolve to klebsiella_pneumoniae, so the router is "
                "not discriminating here and the FN result carries no information either way")
    elif killed:
        verdict = "CLAIM_CONTRADICTED_FOR_SOME_ISOLATES"
        note = f"{len(killed)} FN were called klebsiella_pneumoniae by an independent route"
    elif len(abst) == len(fnr) - len(errs) and ctrl_ok:
        verdict = "CORROBORATED_FN_ALL_ABSTAIN_CONTROLS_CALLED"
        note = ("every FN is outside the 14-organism supported set while the controls are called -- a "
                "second, method-independent signal that the FN are not the scored species")
    else:
        verdict = "MIXED"
        note = "some FN neither abstained nor were called klebsiella_pneumoniae -- read the rows"
    print(f"verdict: {verdict}")
    print(f"note: {note}")
    out["verdict"] = verdict
    out["note"] = note
    out["prediction_registered_before_run"] = (
        "FN abstain (K. aerogenes is outside the supported 14 and sits above MAX_DISTANCE); "
        "K. pneumoniae controls are called. A FN called klebsiella_pneumoniae kills the claim.")
    (REPO / "identify_crosscheck.tmp.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    print("[-> identify_crosscheck.tmp.json]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
