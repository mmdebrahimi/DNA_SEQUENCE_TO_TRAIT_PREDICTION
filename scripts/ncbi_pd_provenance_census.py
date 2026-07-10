"""Stage-1 PROVENANCE CENSUS (no model) — is a provenance-disjoint NCBI-PD validation subset even powered?

Found via /brainstorm (2026-06-10): the NCBI Pathogen Detection metadata TSV the validator already streams
carries submitter/provenance columns (bioproject_center / collected_by / sra_center). Filtering AST to
submitters OUTSIDE the dominant public-health surveillance networks (NARMS / GenomeTrakr / PulseNet / CDC /
FDA / USDA) yields a FREE, genome-linked, **provenance-disjoint (different-submitter-lab)** subset — a
stronger-than-same-ecosystem (but NOT methodology-independent) validation the decoder has never been run on.

This Stage-1 census does NO modelling + NO genome download: it streams the metadata once per organism and
tallies, per drug, isolates that (a) have a downloadable asm_acc AND (b) carry a `<drug>=R`/`<drug>=S` AST
call, split by submitter-class (ecosystem-EXCLUDED vs OTHER). It reports the OTHER subset's N + R/S balance
so we can decide whether a Stage-2 scoring pass is worth running (gate: >= MIN per class, both classes).

Usage: .venv/Scripts/python.exe scripts/ncbi_pd_provenance_census.py --group Campylobacter --drug ciprofloxacin
       (repeatable; default runs a small tractable panel). Exit 0 always (it's a census, not a gate).
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import urllib.request
from datetime import date as _date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

_SIDECAR = Path(__file__).resolve().parent.parent / "wiki" / "provdisjoint_census_results.json"

# Submitters that are part of the integrated public-health AST ecosystem (NOT provenance-disjoint from the
# BV-BRC/NCBI-PD labels the decoder was tuned/cross-validated on). Case-insensitive substring match across
# bioproject_center + collected_by + sra_center.
ECOSYSTEM = ["narms", "genometrakr", "pulsenet", "cdc", "centers for disease", "fda",
             "food and drug", "usda", "national antimicrobial"]
MIN_PER_CLASS = 20   # Stage-2 powering bar (per the /brainstorm: >=~20/class balanced)

DEFAULT_PANEL = [("Campylobacter", "ciprofloxacin"), ("Salmonella", "ciprofloxacin"),
                 ("Klebsiella", "ciprofloxacin")]


def latest_metadata_url(group: str) -> str:
    base = f"https://ftp.ncbi.nlm.nih.gov/pathogen/Results/{group}/latest_snps/Metadata/"
    with urllib.request.urlopen(base, timeout=120) as r:
        html = r.read().decode("utf-8", "replace")
    m = re.findall(r'href="(PDG[0-9.]+\.metadata\.tsv)"', html)
    if not m:
        raise RuntimeError(f"no PDG metadata at {base}")
    return base + sorted(m)[-1]


from dna_decode.data.pd_ast import parse_ast_phenotypes  # noqa: E402


def is_ecosystem(*cells: str) -> bool:
    blob = " ".join(c for c in cells if c).lower()
    return any(tok in blob for tok in ECOSYSTEM)


def census_group(group: str, drugs: list[str], row_cap: int | None = None) -> dict[str, dict]:
    """Census MANY drugs in ONE pass over the group's PD metadata.

    Streaming the (hundreds-of-MB) metadata once per drug was wasteful; the labels for every drug live in
    the SAME `AST_phenotypes` cell. Returns {drug: result}.

    The label is parsed by the canonical `dna_decode.data.pd_ast.parse_ast_phenotypes`, which strips the
    quotes that wrap the field. The previous `f"{drug}=R" in ast.split(",")` test could never match a drug
    at the FIRST or LAST position (the quote rode on that token) -- an under-count of +91 cipro / +62 cef
    isolates over 120k E. coli rows. Under-counts only, so a `powered` verdict was never falsely claimed.
    """
    url = latest_metadata_url(group)
    wanted = [d.lower() for d in drugs]
    counts = {d: {"other_R": 0, "other_S": 0, "eco_R": 0, "eco_S": 0, "rows": 0, "with_ast": 0}
              for d in wanted}
    other_centers: dict[str, dict[str, int]] = {d: {} for d in wanted}
    n_rows = 0
    capped = False
    with urllib.request.urlopen(url, timeout=300) as resp:  # noqa: S310 (trusted NCBI FTP-over-HTTPS)
        header = resp.readline().decode("utf-8", "replace").rstrip("\n").split("\t")
        idx = {name: header.index(name) for name in
               ("asm_acc", "AST_phenotypes", "bioproject_center", "collected_by", "sra_center")
               if name in header}
        ai, pi = idx.get("asm_acc"), idx.get("AST_phenotypes")
        if ai is None or pi is None:
            return {d: {"group": group, "drug": d, "error": "missing asm_acc/AST_phenotypes columns"}
                    for d in wanted}
        for raw in resp:
            n_rows += 1
            if row_cap and n_rows > row_cap:
                capped = True
                break
            cells = raw.decode("utf-8", "replace").rstrip("\n").split("\t")
            if len(cells) <= max(ai, pi):
                continue
            acc, ast = cells[ai], cells[pi]
            if not acc.startswith(("GCA_", "GCF_")):
                continue
            labels = parse_ast_phenotypes(ast, wanted)
            if not labels:
                continue
            prov = [cells[idx[c]] for c in ("bioproject_center", "collected_by", "sra_center")
                    if c in idx and len(cells) > idx[c]]
            eco = is_ecosystem(*prov)
            key = "eco" if eco else "other"
            center = (cells[idx["bioproject_center"]] if "bioproject_center" in idx
                      and len(cells) > idx["bioproject_center"] else "") or "(blank)"
            for drug, lab in labels.items():
                counts[drug]["with_ast"] += 1
                counts[drug][f"{key}_{lab}"] += 1
                if not eco:
                    other_centers[drug][center] = other_centers[drug].get(center, 0) + 1

    out: dict[str, dict] = {}
    for d in wanted:
        c = counts[d]
        c["rows"] = n_rows
        if capped:
            c["capped"] = True
        nR, nS = c["other_R"], c["other_S"]
        out[d] = {"group": group, "drug": d, **c,
                  "other_total": nR + nS,
                  "powered": nR >= MIN_PER_CLASS and nS >= MIN_PER_CLASS,
                  "top_other_centers": sorted(other_centers[d].items(), key=lambda kv: -kv[1])[:8]}
    return out


def census_one(group: str, drug: str, row_cap: int | None = None) -> dict:
    """Single-drug census (thin wrapper over `census_group`; kept for the existing CLI + tests)."""
    return census_group(group, [drug], row_cap)[drug.lower()]


def census_result_to_sidecar_row(result: dict, today: str, min_per_class: int) -> dict | None:
    """Normalize a census_one() result into a sidecar row (schema provdisjoint-census-results-v1).

    Maps `group` -> `organism` (census_one returns `group`; the sidecar + report-card key on `organism`).
    Returns None for rows that MUST NOT be persisted (M1): an error row, or a row-capped run — these would
    overwrite a prior GOOD powering verdict with degraded data. The caller skips persistence on None.
    """
    if not result or "error" in result or result.get("capped"):
        return None
    if "other_R" not in result or "other_S" not in result:
        return None
    return {"organism": result["group"], "drug": result["drug"],
            "other_R": result["other_R"], "other_S": result["other_S"],
            "powered": bool(result.get("powered")), "date": today}


def upsert_census_result(result: dict, today: str, min_per_class: int, path: Path = _SIDECAR) -> bool:
    """Upsert one census result into the powering sidecar, matched on (organism, drug). Idempotent: a
    re-run for the same organism×drug REPLACES its row, never duplicates. Returns False (skipped) when the
    row is degraded (error/capped) per `census_result_to_sidecar_row`."""
    row = census_result_to_sidecar_row(result, today, min_per_class)
    if row is None:
        return False
    if path.exists():
        doc = json.loads(path.read_text(encoding="utf-8"))
    else:
        doc = {"_schema": "provdisjoint-census-results-v1",
               "_doc": "Machine-readable powering facts from scripts/ncbi_pd_provenance_census.py "
                       "(Stage-1, no-model). Self-persisted per census run; consumed by the report card.",
               "min_per_class": min_per_class, "ecosystem_excluded": ECOSYSTEM, "results": []}
    doc.setdefault("results", [])
    key = (row["organism"].lower(), row["drug"].lower())
    doc["results"] = [r for r in doc["results"]
                      if (r.get("organism", "").lower(), r.get("drug", "").lower()) != key]
    doc["results"].append(row)
    doc["results"].sort(key=lambda r: (r.get("organism", ""), r.get("drug", "")))
    path.write_text(json.dumps(doc, indent=2), encoding="utf-8")
    return True


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--group", help="NCBI PD organism group (e.g. Campylobacter, Salmonella, Klebsiella, Escherichia_coli_Shigella)")
    ap.add_argument("--drug", default="ciprofloxacin")
    ap.add_argument("--row-cap", type=int, default=None, help="cap rows streamed (for huge groups)")
    a = ap.parse_args()
    panel = [(a.group, a.drug)] if a.group else DEFAULT_PANEL
    print(f"Stage-1 provenance census (ecosystem-excluded = {ECOSYSTEM}); powering bar >= {MIN_PER_CLASS}/class\n")
    for grp, drug in panel:
        try:
            r = census_one(grp, drug, a.row_cap)
        except Exception as e:
            print(f"{grp} {drug}: ERROR {type(e).__name__}: {e}"); continue
        if "error" in r:
            print(f"{grp} {drug}: {r['error']}"); continue
        print(f"=== {grp} x {drug} ===")
        print(f"  rows={r['rows']} with_AST(downloadable)={r['with_ast']}  "
              f"ecosystem: {r['eco_R']}R/{r['eco_S']}S  |  OTHER (provenance-disjoint): {r['other_R']}R/{r['other_S']}S")
        print(f"  POWERED (>= {MIN_PER_CLASS}/class both): {r['powered']}")
        if r["top_other_centers"]:
            print(f"  top non-ecosystem submitters: {r['top_other_centers']}")
        if upsert_census_result(r, _date.today().isoformat(), MIN_PER_CLASS):
            print(f"  persisted -> {_SIDECAR.relative_to(Path.cwd()) if _SIDECAR.is_relative_to(Path.cwd()) else _SIDECAR}")
        else:
            print("  NOT persisted (error/capped row — prior powering verdict preserved)")
        print()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
