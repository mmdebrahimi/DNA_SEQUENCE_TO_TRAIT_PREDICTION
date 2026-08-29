"""How SOURCE-CONCENTRATED are the 10 provenance-disjoint SCORED cells?

WHY THIS EXISTS. The report card already discloses that the SCORED cells are CLONALLY dominated — an
over-sampled clone carries the raw metric, so a lineage layer collapses same-label lineages to one vote.
This asks the sibling question one level up: how many distinct **sources** (BioProject, sequencing centre,
country) does each cell actually draw on?

The trigger was concrete. E. coli x ciprofloxacin scores spec 0.700 on the report card. A separate,
accession-disjoint set of 131 isolates scored spec 0.988 with the SAME frozen rule. Checking provenance
refuted the obvious explanation (that the second set was in-distribution and therefore easier): it spans
8 BioProjects while the report-card cohort turned out to be **58 of 60 from one BioProject at one
hospital**. So the estimate that looked pessimistic is measured on roughly one lab.

"Provenance-disjoint" and "provenance-diverse" are different properties. The cells are the first by
construction; this measures whether they are also the second. A cell drawing on one BioProject can be
perfectly disjoint from the tuning data and still be a single-site estimate.

WHAT THIS IS NOT. Not a demotion, not a re-scoring, and not a claim that any published number is wrong.
It is a disclosure measurement in the same spirit as the lineage layer: the numbers stand, and the reader
learns how many independent sources stand behind them.

Run: uv run python scripts/provdisjoint_source_concentration.py
"""
from __future__ import annotations

import csv
import gzip
import io
import json
import sys
import urllib.request
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

# report-card organism -> PD group whose metadata table carries these accessions
PD_GROUP = {
    "escherichia_coli_shigella": "Escherichia_coli_Shigella",
    "klebsiella": "Klebsiella",
    "campylobacter": "Campylobacter",
}
FIELDS = ("bioproject_acc", "sra_center", "geo_loc_name", "collected_by")


def scored_cells() -> list[dict]:
    card = json.loads((ROOT / "wiki" / "decoder_validation_report_card.json").read_text(encoding="utf-8"))
    return [c for c in card["cells"] if c.get("state") == "SCORED"]


def cohort_accessions(organism: str, drug: str) -> set[str]:
    """The exact accessions the cell was scored on, from its cohort selected.tsv."""
    base = ROOT / "data" / "raw" / f"{organism}_provdisjoint_{drug}"
    sel = base / "selected.tsv"
    if not sel.exists():
        return set()
    out = set()
    for line in sel.read_text(encoding="utf-8", errors="replace").splitlines():
        parts = line.split("\t")
        if len(parts) >= 2 and parts[1].strip().upper() in ("R", "S"):
            out.add(parts[0].strip())
    return out


def fetch_provenance(group: str, wanted: set[str]) -> dict[str, dict[str, str]]:
    """{accession: {field: value}} from one PD metadata table. Network."""
    from ncbi_pd_provenance_census import latest_metadata_url

    url = latest_metadata_url(group)
    out: dict[str, dict[str, str]] = {}
    with urllib.request.urlopen(url, timeout=900) as resp:
        raw = gzip.GzipFile(fileobj=resp) if url.endswith(".gz") else resp
        rd = csv.reader(io.TextIOWrapper(raw, encoding="utf-8", errors="replace"), delimiter="\t")
        header = next(rd, None)
        if not header:
            return out
        ix = {c: i for i, c in enumerate(header)}
        ai = ix.get("asm_acc")
        if ai is None:
            return out
        for row in rd:
            if len(row) <= ai:
                continue
            acc = row[ai].strip()
            if acc in wanted:
                out[acc] = {f: (row[ix[f]].strip() if ix.get(f) is not None and len(row) > ix[f] else "")
                            for f in FIELDS}
    return out


def concentration(prov: dict[str, dict[str, str]], field: str) -> dict:
    """Distinct sources and the share held by the largest one. PURE.

    `largest_share` is the number that matters: 8 BioProjects sounds diverse until one of them holds 97%.
    NULL/empty values are counted as UNKNOWN rather than merged into one pseudo-source, which would
    manufacture concentration that is really missing metadata.
    """
    vals = [(v.get(field) or "").strip() or "UNKNOWN" for v in prov.values()]
    known = [v for v in vals if v != "UNKNOWN"]
    c = Counter(known)
    n = len(known)
    return {"n_resolved": len(vals), "n_known": n, "n_unknown": len(vals) - n,
            "distinct": len(c), "largest": c.most_common(1)[0] if c else None,
            "largest_share": round(c.most_common(1)[0][1] / n, 3) if n else None}


def main() -> int:
    cells = scored_cells()
    by_group: dict[str, set[str]] = {}
    for c in cells:
        g = PD_GROUP.get(c["organism"])
        if g:
            by_group.setdefault(g, set()).update(cohort_accessions(c["organism"], c["drug"]))

    prov: dict[str, dict[str, str]] = {}
    errors: dict[str, str] = {}
    for g, want in by_group.items():
        if not want:
            continue
        try:
            prov.update(fetch_provenance(g, want))
        except Exception as exc:
            errors[g] = f"{type(exc).__name__}: {exc}"
            print(f"  {g}: FAILED ({errors[g]})")

    if errors:
        print("INCOMPLETE: at least one metadata table failed; per-cell rows below may be partial.\n")

    rows = []
    print(f"{'cell':46}{'n':>4}{'BioProjects':>13}{'largest share':>15}  dominant source")
    for c in sorted(cells, key=lambda x: (x["organism"], x["drug"])):
        accs = cohort_accessions(c["organism"], c["drug"])
        sub = {a: v for a, v in prov.items() if a in accs}
        if not sub:
            print(f"{c['organism']+' x '+c['drug']:46}{len(accs):>4}{'unresolved':>13}")
            continue
        bp = concentration(sub, "bioproject_acc")
        ctr = concentration(sub, "sra_center")
        dom = bp["largest"][0] if bp["largest"] else "?"
        print(f"{c['organism']+' x '+c['drug']:46}{len(accs):>4}{bp['distinct']:>13}"
              f"{(bp['largest_share'] if bp['largest_share'] is not None else 0):>15.0%}  {dom}")
        rows.append({"organism": c["organism"], "drug": c["drug"], "n_cohort": len(accs),
                     "spec": c.get("spec"), "sens": c.get("sens"),
                     "bioproject": bp, "sra_center": ctr})

    single = [r for r in rows if (r["bioproject"]["largest_share"] or 0) >= 0.80]
    print()
    print(f"cells where ONE BioProject holds >=80% of the cohort: {len(single)}/{len(rows)}")
    for r in single:
        print(f"  {r['organism']} x {r['drug']}  "
              f"{r['bioproject']['largest_share']:.0%} from {r['bioproject']['largest'][0]}")
    print()
    print("These numbers are NOT demotions. The cells remain provenance-DISJOINT from the tuning data;")
    print("this measures whether they are also provenance-DIVERSE, which is a different property and")
    print("was never claimed. A single-BioProject cell is a single-site estimate.")

    (ROOT / "wiki" / "provdisjoint_source_concentration.json").write_text(
        json.dumps({"cells": rows, "errors": errors, "complete": not errors,
                    "n_single_source": len(single)}, indent=2), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
