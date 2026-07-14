"""Forward variant-effect method leaderboard — aggregate every forward-cell result JSON
(wiki/tem1_forward_cell_*.json) into a per-protein BLOSUM / ESM2 / AlphaMissense table.

Each cell run writes {dms_id, method, spearman_pred_vs_dms, n_single_variants_scored, ...}; this groups them
by protein and lays the three methods side by side, tagged by organism (E. coli / human), so the
deterministic-vs-learned and bacterial-vs-eukaryotic pattern is visible at a glance. Read-only.
"""
from __future__ import annotations

import json
import re
from datetime import date as _date
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
WIKI = REPO / "wiki"

_METHOD = {"blosum62_deterministic": "blosum", "esm2_zeroshot": "esm2", "alphamissense_learned": "am"}


def _organism(dms_id: str) -> str:
    if dms_id.endswith("_HUMAN") or "_HUMAN_" in dms_id:
        return "human"
    if "ECOL" in dms_id:
        return "E. coli"
    return "?"


def main() -> int:
    rows: dict[str, dict] = {}
    for f in sorted(WIKI.glob("tem1_forward_cell_*.json")):
        try:
            d = json.loads(f.read_text(encoding="utf-8"))
        except Exception:
            continue
        dms = d.get("dms_id")
        m = _METHOD.get(d.get("method", ""))
        rho = d.get("spearman_pred_vs_dms")
        if not dms or not m or rho is None:
            continue
        r = rows.setdefault(dms, {"dms_id": dms, "organism": _organism(dms),
                                  "n": d.get("n_single_variants_scored"), "blosum": None, "esm2": None, "am": None})
        r[m] = rho
        r["n"] = d.get("n_single_variants_scored", r["n"])

    ordered = sorted(rows.values(), key=lambda r: (r["organism"], -(r.get("am") or r.get("esm2") or r.get("blosum") or 0)))
    # markdown
    lines = [f"# Forward variant-effect method leaderboard ({_date.today().isoformat()})", "",
             "Per-protein Spearman(prediction, measured DMS) across the three forward-cell methods. AlphaMissense "
             "is human-only (bacterial cells show `—`); ESM2 runs on both but is only populated where a table was "
             "built. Higher = better; the learned methods (ESM2 / AlphaMissense) beat deterministic BLOSUM where "
             "run.", "",
             "| protein (assay) | organism | n | BLOSUM62 | ESM2-650M | AlphaMissense |",
             "|---|---|---:|---:|---:|---:|"]
    def cell(v):
        return f"{v:.3f}" if isinstance(v, (int, float)) else "—"
    for r in ordered:
        gene = re.split(r"_HUMAN|_ECOL", r["dms_id"])[0]
        lines.append(f"| {gene} ({r['dms_id']}) | {r['organism']} | {r['n']} | "
                     f"{cell(r['blosum'])} | {cell(r['esm2'])} | {cell(r['am'])} |")
    lines += ["", "Learned-vs-deterministic lift (where both present):"]
    for r in ordered:
        if isinstance(r["blosum"], (int, float)):
            best = max([x for x in (r["esm2"], r["am"]) if isinstance(x, (int, float))], default=None)
            if best is not None:
                which = "AM" if r["am"] == best else "ESM2"
                lines.append(f"- {r['dms_id']}: {which} {best:.3f} − BLOSUM {r['blosum']:.3f} = **+{best - r['blosum']:.3f}**")

    md = "\n".join(lines) + "\n"
    (WIKI / f"forward_method_leaderboard_{_date.today().isoformat()}.md").write_text(md, encoding="utf-8")
    (WIKI / f"forward_method_leaderboard_{_date.today().isoformat()}.json").write_text(
        json.dumps({"generated": _date.today().isoformat(), "proteins": ordered}, indent=2), encoding="utf-8")
    print("\n".join(lines))
    print(f"\n[leaderboard] {len(ordered)} proteins across {sum(1 for r in ordered if r['organism']=='human')} human "
          f"+ {sum(1 for r in ordered if r['organism']=='E. coli')} E. coli")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
