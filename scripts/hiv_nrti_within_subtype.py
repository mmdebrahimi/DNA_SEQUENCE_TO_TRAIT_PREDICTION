"""NRTI within-subtype transfer check (Wave B Rec 2) — does the consensus-B catalog hold on non-B subtypes?

For a DETERMINISTIC mechanism catalog (the mutation IS the mechanism, on consensus-B numbering), the
de-confound question is a TRANSFER one: does "major NRTI determinant present -> R" work as well outside
subtype B as inside it? Uses the UNFILTERED NRTI dataset (the only one carrying a Subtype column),
restricted to Method=PhenoSense + Type=Clinical for a consistent comparison.

HONEST FINDING up front: the free HIVDB genotype-phenotype data is ~96% subtype B (5143 B vs ~150 non-B
across all clades), so the non-B transfer estimate is UNDER-POWERED — a real limitation of the free data for
this discipline, reported (not hidden). B vs pooled-non-B is the only adequately-powered split.

Catalog = the position-based v0 (deterministic, no data-derivation -> no train/test concern). Label =
PhenoSense fold (independent wet-lab; NOT Sierra). Cutoffs from DRMcv.R. Cite Rhee 2003.
DATA: data/raw/hiv/NRTI_DataSet.Full.txt (gitignored; download .../GenoPhenoDatasets/NRTI_DataSet.Full.txt).
"""
from __future__ import annotations

import json
import sys
from datetime import date as _date
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from dna_decode.data.hiv_amr import call_nrti_from_observed
from scripts.hiv_nnrti_baseline import _confusion
from scripts.hiv_nnrti_validate import _parse_fold, load_rows
from scripts.hiv_nrti_validate import NRTI_DRUGS, NRTI_LOWER_CUTOFF, _NRTI_COL, _observed_nrti

DEFAULT_DATA = REPO / "data" / "raw" / "hiv" / "NRTI_DataSet.Full.txt"


def _group(subtype: str) -> str:
    s = (subtype or "").strip()
    return "B" if s == "B" else ("non-B" if s and s.lower() != "unknown" else "unknown")


def run(path: Path = DEFAULT_DATA) -> dict:
    rows = [r for r in load_rows(path)
            if (r.get("Method", "").strip() == "PhenoSense" and r.get("Type", "").strip() == "Clinical")]
    subtype_counts: dict[str, int] = {}
    for r in rows:
        subtype_counts[r.get("Subtype", "").strip() or "(blank)"] = \
            subtype_counts.get(r.get("Subtype", "").strip() or "(blank)", 0) + 1

    per_drug = {}
    for drug in NRTI_DRUGS:
        col, cutoff = _NRTI_COL[drug], NRTI_LOWER_CUTOFF[drug]
        groups: dict[str, dict] = {"B": {"call": [], "fold": []}, "non-B": {"call": [], "fold": []}}
        for r in rows:
            g = _group(r.get("Subtype", ""))
            if g == "unknown":
                continue
            fold = _parse_fold(r.get(col, ""))
            if fold is None or fold <= 0:
                continue
            call_R = call_nrti_from_observed(drug, {"RT": _observed_nrti(r)}).prediction == "R"
            groups[g]["call"].append(call_R)
            groups[g]["fold"].append(fold)
        out = {}
        for g, d in groups.items():
            if len(d["fold"]) < 15:
                out[g] = {"n": len(d["fold"]), "note": "under-powered (<15)"}
                continue
            call = np.array(d["call"], dtype=bool)
            actual = np.array(d["fold"]) >= cutoff
            out[g] = {"n": len(d["fold"]), "n_R": int(actual.sum()), **_confusion(call, actual)}
        per_drug[drug] = out
    return {
        "artifact": "hiv_nrti_within_subtype_transfer", "schema": "hiv-nrti-subtype-v0",
        "catalog": "position-based v0 (consensus-B numbering)",
        "label_source": "Stanford HIVDB PhenoSense fold (independent wet-lab; NOT Sierra)",
        "filter": "Method=PhenoSense AND Type=Clinical",
        "dataset": str(path), "n_clinical_phenosense": len(rows),
        "subtype_counts": dict(sorted(subtype_counts.items(), key=lambda kv: -kv[1])),
        "honest_caveats": [
            "the free HIVDB gp data is ~96% subtype B -> non-B transfer is UNDER-POWERED (a free-data limit, reported)",
            "per-clade non-B N is tiny; only B-vs-pooled-non-B is adequately powered, and even non-B is small",
            "a similar B vs non-B sens/spec is consistent with transfer; it does NOT prove non-B generalisation at scale",
        ],
        "citation": "Rhee 2003 Nucleic Acids Res 31:298-303; cutoffs from DRMcv.R",
        "per_drug": per_drug,
    }


def render_md(result: dict, generated: str) -> str:
    lines = [f"# HIV NRTI within-subtype transfer check ({generated})", ""]
    lines.append(f"Catalog = {result['catalog']}. Label = {result['label_source']}. "
                 f"Filter = {result['filter']}; N = {result['n_clinical_phenosense']}.")
    lines.append("")
    nonb = sum(v for k, v in result["subtype_counts"].items() if k not in ("B", "Unknown", "(blank)"))
    lines.append(f"**Subtype mix:** B = {result['subtype_counts'].get('B', 0)}, "
                 f"non-B (pooled) = {nonb}. The data is B-dominated -> non-B is under-powered.")
    lines.append("")
    lines.append("| Drug | B n (R) sens/spec/**balacc** | non-B n (R) sens/spec/**balacc** |")
    lines.append("|---|---|---|")
    for drug, g in result["per_drug"].items():
        def cell(x):
            if "balanced_accuracy" not in x:
                return f"{x.get('n')} ({x.get('note','')})"
            return f"{x['n']} ({x['n_R']}) {x['sens']}/{x['spec']}/**{x['balanced_accuracy']}**"
        lines.append(f"| {drug} | {cell(g.get('B', {}))} | {cell(g.get('non-B', {}))} |")
    lines.append("")
    lines.append("## Honest caveats")
    for c in result["honest_caveats"]:
        lines.append(f"- {c}")
    lines.append("")
    lines.append(f"Citation: {result['citation']}.")
    return "\n".join(lines)


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data", type=Path, default=DEFAULT_DATA)
    ap.add_argument("--out-md", type=Path, default=None)
    args = ap.parse_args(argv)
    if not args.data.exists():
        print(f"ERROR: unfiltered dataset not found at {args.data}\n"
              f"  download: curl -L -o {args.data} "
              f"https://hivdb.stanford.edu/download/GenoPhenoDatasets/NRTI_DataSet.Full.txt", file=sys.stderr)
        return 2
    today = _date.today().isoformat()
    result = run(args.data)
    out_md = args.out_md or (REPO / "wiki" / f"hiv_nrti_within_subtype_{today}.md")
    out_md.write_text(render_md(result, today), encoding="utf-8")
    out_md.with_suffix(".json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(render_md(result, today))
    print(f"\n[wrote {out_md} + .json]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
