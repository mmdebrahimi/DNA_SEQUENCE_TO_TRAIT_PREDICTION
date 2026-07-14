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

_METHOD = {"blosum62_deterministic": "blosum", "esm2_zeroshot": "esm2", "alphamissense_learned": "am",
           "esm_if_structure": "esm_if"}

# ProteinGym DMS_id = GENE_SPECIES_Author_Year; map the species code (2nd token).
_SPECIES = {"HUMAN": "human", "YEAST": "yeast", "ARATH": "Arabidopsis", "MOUSE": "mouse", "CHICK": "chicken",
            "ECOLI": "E. coli", "ECOLX": "E. coli", "BACSU": "B. subtilis", "PSEAI": "P. aeruginosa",
            "SARS2": "SARS-CoV-2", "CHLRE": "C. reinhardtii", "RHOTO": "R. toruloides"}
# rough kingdom grouping for ordering / narrative
_KINGDOM = {"human": "eukaryote", "yeast": "eukaryote", "Arabidopsis": "eukaryote", "mouse": "eukaryote",
            "chicken": "eukaryote", "C. reinhardtii": "eukaryote", "R. toruloides": "eukaryote",
            "E. coli": "bacterium", "B. subtilis": "bacterium", "P. aeruginosa": "bacterium",
            "SARS-CoV-2": "virus"}


def _organism(dms_id: str) -> str:
    parts = dms_id.split("_")
    code = parts[1] if len(parts) > 1 else ""
    return _SPECIES.get(code, code.lower() or "?")


def _kingdom(org: str) -> str:
    return _KINGDOM.get(org, "?")


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
                                  "n": d.get("n_single_variants_scored"), "blosum": None, "esm2": None,
                                  "am": None, "esm_if": None})
        r[m] = rho
        r["n"] = d.get("n_single_variants_scored", r["n"])

    ordered = sorted(rows.values(), key=lambda r: (_kingdom(r["organism"]), r["organism"],
                                                    -(r.get("am") or r.get("esm2") or r.get("blosum") or 0)))
    # markdown
    lines = [f"# Forward variant-effect method leaderboard ({_date.today().isoformat()})", "",
             "Per-protein Spearman(prediction, measured DMS) across the three forward-cell methods. AlphaMissense "
             "is human-only (bacterial cells show `—`); ESM2 runs on both but is only populated where a table was "
             "built. Higher = better; the learned methods (ESM2 / AlphaMissense) beat deterministic BLOSUM where "
             "run.", "",
             "| protein (assay) | organism | n | BLOSUM62 | ESM2-650M | AlphaMissense | ESM-IF |",
             "|---|---|---:|---:|---:|---:|---:|"]
    def cell(v):
        return f"{v:.3f}" if isinstance(v, (int, float)) else "—"
    for r in ordered:
        gene = r["dms_id"].split("_")[0]
        lines.append(f"| {gene} ({r['dms_id']}) | {r['organism']} | {r['n']} | "
                     f"{cell(r['blosum'])} | {cell(r['esm2'])} | {cell(r['am'])} | {cell(r['esm_if'])} |")
    lines += ["", "Learned-vs-deterministic lift (where both present):"]
    for r in ordered:
        if isinstance(r["blosum"], (int, float)):
            best = max([x for x in (r["esm2"], r["am"], r["esm_if"]) if isinstance(x, (int, float))], default=None)
            if best is not None:
                which = "AM" if r["am"] == best else ("ESM-IF" if r["esm_if"] == best else "ESM2")
                lines.append(f"- {r['dms_id']}: {which} {best:.3f} − BLOSUM {r['blosum']:.3f} = **+{best - r['blosum']:.3f}**")

    md = "\n".join(lines) + "\n"
    (WIKI / f"forward_method_leaderboard_{_date.today().isoformat()}.md").write_text(md, encoding="utf-8")
    (WIKI / f"forward_method_leaderboard_{_date.today().isoformat()}.json").write_text(
        json.dumps({"generated": _date.today().isoformat(), "proteins": ordered}, indent=2), encoding="utf-8")
    print("\n".join(lines))
    from collections import Counter
    by_org = Counter(r["organism"] for r in ordered)
    by_king = Counter(_kingdom(r["organism"]) for r in ordered)
    print(f"\n[leaderboard] {len(ordered)} proteins across {len(by_org)} organisms "
          f"({', '.join(f'{n} {o}' for o, n in by_org.most_common())}); "
          f"kingdoms: {', '.join(f'{n} {k}' for k, n in by_king.most_common())}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
