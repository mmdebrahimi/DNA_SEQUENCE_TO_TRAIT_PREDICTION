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
import re
import sys
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

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


def is_ecosystem(*cells: str) -> bool:
    blob = " ".join(c for c in cells if c).lower()
    return any(tok in blob for tok in ECOSYSTEM)


def census_one(group: str, drug: str, row_cap: int | None = None) -> dict:
    url = latest_metadata_url(group)
    tok_r, tok_s = f"{drug}=R", f"{drug}=S"
    counts = {"other_R": 0, "other_S": 0, "eco_R": 0, "eco_S": 0, "rows": 0, "with_ast": 0}
    other_centers: dict[str, int] = {}
    with urllib.request.urlopen(url, timeout=300) as resp:
        header = resp.readline().decode("utf-8", "replace").rstrip("\n").split("\t")
        idx = {name: header.index(name) for name in
               ("asm_acc", "AST_phenotypes", "bioproject_center", "collected_by", "sra_center")
               if name in header}
        ai, pi = idx.get("asm_acc"), idx.get("AST_phenotypes")
        if ai is None or pi is None:
            return {"group": group, "drug": drug, "error": "missing asm_acc/AST_phenotypes columns"}
        for raw in resp:
            counts["rows"] += 1
            if row_cap and counts["rows"] > row_cap:
                counts["capped"] = True
                break
            cells = raw.decode("utf-8", "replace").rstrip("\n").split("\t")
            if len(cells) <= max(ai, pi):
                continue
            acc, ast = cells[ai], cells[pi]
            if not acc.startswith(("GCA_", "GCF_")) or not ast or ast == "NULL":
                continue
            lab = 1 if tok_r in ast.split(",") else (0 if tok_s in ast.split(",") else None)
            if lab is None:
                continue
            counts["with_ast"] += 1
            prov = [cells[idx[c]] for c in ("bioproject_center", "collected_by", "sra_center")
                    if c in idx and len(cells) > idx[c]]
            eco = is_ecosystem(*prov)
            key = "eco" if eco else "other"
            counts[f"{key}_{'R' if lab else 'S'}"] += 1
            if not eco:
                center = (cells[idx["bioproject_center"]] if "bioproject_center" in idx
                          and len(cells) > idx["bioproject_center"] else "") or "(blank)"
                other_centers[center] = other_centers.get(center, 0) + 1
    nR, nS = counts["other_R"], counts["other_S"]
    powered = nR >= MIN_PER_CLASS and nS >= MIN_PER_CLASS
    top_centers = sorted(other_centers.items(), key=lambda kv: -kv[1])[:8]
    return {"group": group, "drug": drug, **counts,
            "other_total": nR + nS, "powered": powered,
            "top_other_centers": top_centers}


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
        print()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
