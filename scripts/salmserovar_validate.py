"""Score the shipped Salmonella serovar caller against a WET-LAB label — and against the tool it mimics.

THE POINT OF THE THREE NUMBERS. `typing:Salmonella:salmserovar` ships at `FAITHFUL_TO_TOOL`, meaning it
has only ever been checked against the reference METHOD. This project's own recorded lesson is that a
policy layer over an external tool must be validated against NAIVE use of that tool on INDEPENDENT
data, because in-cohort agreement only proves the tool works. So three things are measured on the SAME
isolates:

  1. ours vs the wet-lab label     -- the number that would move the cell off FAITHFUL_TO_TOOL
  2. the in-silico tool vs the label -- the incumbent, without which (1) is uninterpretable
  3. ours vs the in-silico tool    -- pure faithfulness, the old number, kept for continuity

(2) is the one that makes this honest. If our caller scores 0.85, that is excellent when the incumbent
scores 0.86 on the same isolates and poor when the incumbent scores 0.97 — and only running both on one
cohort can tell those apart.

RESIDUAL CIRCULARITY IS BOUNDED BY THE DESIGN, not by an assurance. Some public serovar strings really
are tool output pasted into a metadata field, and no filter can prove otherwise per isolate. But since
(1) and (2) are scored against the SAME labels, contamination inflates BOTH — so the DELTA survives it
even where the absolute numbers are optimistic. Quote the delta with more confidence than the level.

Equivalence is decided by `dna_decode.salmserovar.equivalence`, applied identically to both callers so
no leniency can favour one of them.

Needs blastn + network (assembly fetch). Checkpointed to JSONL; restartable; `--limit` runs a pilot.
"""
from __future__ import annotations

import argparse
import collections
import json
import sys
import traceback
from datetime import date as _date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from dna_decode.data.refseq import download_genome, fasta_path  # noqa: E402
from dna_decode.salmserovar.equivalence import equivalent, load_formula_index  # noqa: E402
from dna_decode.salmserovar.runner import call_serovar  # noqa: E402

BLASTN = "C:/Users/Farshad/ncbi-blast/bin/blastn.exe"


def score_pair(pred: str | None, label: str, idx: dict) -> tuple[str, str]:
    """-> (outcome, reason). `no_call` is kept SEPARATE from `miss`: a caller that abstains is not a
    caller that is wrong, and merging them would hide an abstention behind an error rate."""
    if not pred:
        return "no_call", "caller returned no serovar"
    ok, why = equivalent(pred, label, idx)
    return ("hit" if ok else "miss"), why


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cohort-json", type=Path,
                    default=ROOT / "wiki" / "salmserovar_cohort_2026-09-04.json")
    ap.add_argument("--db-dir", type=Path, default=ROOT / "data" / "salmserovar_db")
    ap.add_argument("--asm-dir", type=Path, default=Path("D:/dna_decode_cache/salm_asm"))
    ap.add_argument("--checkpoint", type=Path,
                    default=Path("D:/dna_decode_cache/salm_asm/results.jsonl"))
    ap.add_argument("--limit", type=int, default=0, help="pilot on the first N isolates")
    ap.add_argument("--blastn", default=BLASTN)
    ap.add_argument("--out", type=Path,
                    default=ROOT / "wiki" / f"salmserovar_validation_{_date.today().isoformat()}.json")
    a = ap.parse_args()

    coh = json.loads(a.cohort_json.read_text(encoding="utf-8"))
    isolates = coh["isolates"][:a.limit] if a.limit else coh["isolates"]
    idx = load_formula_index(a.db_dir / "serovar_table.tsv")
    a.asm_dir.mkdir(parents=True, exist_ok=True)

    done: dict[str, dict] = {}
    if a.checkpoint.exists():
        for line in a.checkpoint.read_text(encoding="utf-8").splitlines():
            if line.strip():
                r = json.loads(line)
                done[r["asm_acc"]] = r
    print(f"{len(isolates)} isolates | {len(done)} already checkpointed")

    fh = open(a.checkpoint, "a", encoding="utf-8")
    for n, iso in enumerate(isolates, 1):
        acc = iso["asm_acc"]
        if acc in done:
            continue
        rec = {"asm_acc": acc, "label": iso["serovar_label"],
               "computed": iso.get("computed_serotype")}
        try:
            download_genome(acc, a.asm_dir)
            fa = fasta_path(acc, a.asm_dir)
            if not Path(fa).exists():
                rec["status"] = "assembly_missing"
            else:
                call = call_serovar(fa, a.db_dir, blastn_bin=a.blastn)
                rec["status"] = "ok"
                rec["pred_serovar"] = call.get("serovar")
                rec["pred_formula"] = call.get("antigenic_formula") or call.get("formula")
        except Exception as e:                      # noqa: BLE001 - recorded, never silently dropped
            rec["status"] = f"error:{type(e).__name__}"
            rec["error"] = str(e)[:200]
            rec["trace"] = traceback.format_exc()[-400:]
        fh.write(json.dumps(rec) + "\n")
        fh.flush()
        done[acc] = rec
        if n % 10 == 0 or rec["status"] != "ok":
            print(f"  [{n}/{len(isolates)}] {acc} {rec['status']} "
                  f"pred={rec.get('pred_serovar')!r} label={rec['label']!r}", flush=True)
    fh.close()

    # ---- score the three comparisons on the SAME isolates ---------------------------------------
    rows = [done[i["asm_acc"]] for i in isolates if i["asm_acc"] in done]
    ok_rows = [r for r in rows if r.get("status") == "ok"]
    ours = collections.Counter()
    tool = collections.Counter()
    faith = collections.Counter()
    reasons = collections.Counter()
    disagreements = []

    for r in ok_rows:
        o, why = score_pair(r.get("pred_serovar"), r["label"], idx)
        ours[o] += 1
        reasons[why] += 1
        t, _ = score_pair(r.get("computed"), r["label"], idx)
        tool[t] += 1
        f, _ = score_pair(r.get("pred_serovar"), r.get("computed") or "", idx)
        faith[f] += 1
        if o == "miss" and len(disagreements) < 25:
            disagreements.append({"asm": r["asm_acc"], "label": r["label"],
                                  "ours": r.get("pred_serovar"), "tool": r.get("computed")})

    def rate(c: collections.Counter) -> float | None:
        scored = c["hit"] + c["miss"]
        return (c["hit"] / scored) if scored else None

    n_ok = len(ok_rows)
    print(f"\n=== scored on {n_ok} isolates with a completed call ===")
    for name, c in (("ours vs WET-LAB label", ours), ("in-silico tool vs label", tool),
                    ("ours vs in-silico tool", faith)):
        r_ = rate(c)
        print(f"  {name:<26} hit={c['hit']:<4} miss={c['miss']:<4} no_call={c['no_call']:<4} "
              f"accuracy={r_ if r_ is None else round(r_, 3)}")

    o_r, t_r = rate(ours), rate(tool)
    delta = None if (o_r is None or t_r is None) else o_r - t_r
    if delta is not None:
        print(f"\n  DELTA (ours - incumbent) = {delta:+.3f}   <- the interpretable number")

    status_counts = collections.Counter(r.get("status", "?") for r in rows)
    out = {
        "schema": "salmserovar-validation-v1",
        "date": _date.today().isoformat(),
        "cell": "typing:Salmonella:salmserovar",
        "question": "does the shipped serovar caller match a WET-LAB serovar label, and how does it "
                    "compare to naive use of the in-silico tool it mimics?",
        "cohort": {"n_requested": len(isolates), "n_completed": n_ok,
                   "statuses": dict(status_counts),
                   "n_bioprojects": coh["cohort"]["n_bioprojects"],
                   "largest_source_share": coh["cohort"]["largest_source_share"],
                   "passes_source_diversity_bar": coh["cohort"]["passes_source_diversity_bar"],
                   "n_distinct_serovars": coh["cohort"]["n_distinct_serovars"]},
        "ours_vs_wetlab_label": dict(ours), "ours_accuracy": o_r,
        "insilico_tool_vs_wetlab_label": dict(tool), "tool_accuracy": t_r,
        "ours_vs_insilico_tool": dict(faith), "faithfulness": rate(faith),
        "delta_ours_minus_incumbent": delta,
        "equivalence_reasons": dict(reasons),
        "sample_misses": disagreements,
        "honest_limits": [
            "The label is NCBI-PD's submitter `serovar` restricted to reference public-health labs. It "
            "cannot be proven per isolate that any given string came from slide agglutination; the "
            "provenance filter is a judgment.",
            "Residual circularity is bounded, not eliminated -- but both callers are scored against "
            "the SAME labels, so contamination inflates both and the DELTA survives it.",
            "Equivalence is notation-normalisation plus the committed White-Kauffmann table, applied "
            "identically to both callers. It refuses fuzzy near-misses (Newport/Newbrunswick).",
            "A per-serovar cap flattens the natural prevalence, so this is per-isolate accuracy on a "
            "deliberately diverse mix, NOT population-weighted accuracy.",
            "no_call is reported separately from miss: abstention is not error.",
        ],
    }
    a.out.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
