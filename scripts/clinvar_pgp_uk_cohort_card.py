#!/usr/bin/env python
"""Build a real-people ClinVar/Mendelian COHORT card from per-individual PGP-UK decode JSONs.

Aggregates N individuals' `*_clinvar*.json` (from scripts/clinvar_decode_vcf.py) into a cohort summary:
per-sample reportable-pathogenic count + benign carrier load, the cohort-wide pathogenic findings (detailed
if any), and the per-gene benign coverage across the panel. HONEST FRAMING: a research demonstration on real
open-consent individuals; a ClinVar classification is a curated allele label, NOT a clinical diagnosis of the
person; 0 pathogenic across a small healthy cohort is the expected ACMG-SF base-rate result. NOT a clinical tool.

Usage:
    uv run python scripts/clinvar_pgp_uk_cohort_card.py --results-dir C:/.../pgp_results \
        --glob "*_clinvar86.json" --panel-note "ACMG SF v3.2 + 5 carrier genes (86 genes, 171,519 variants)" \
        --out-md wiki/clinvar_pgp_uk_cohort_card_2026-07-06.md
"""
from __future__ import annotations

import argparse
import glob
import json
from collections import Counter
from pathlib import Path


def load(results_dir: Path, pattern: str) -> list[dict]:
    out = []
    for f in sorted(glob.glob(str(results_dir / pattern))):
        d = json.loads(Path(f).read_text(encoding="utf-8"))
        if "n_benign" in d and "sample_id" in d:
            out.append(d)
    return out


def build_card(results: list[dict], panel_note: str) -> dict:
    n = len(results)
    per_sample = [{"sample_id": r["sample_id"], "n_pathogenic": r["n_pathogenic"],
                   "n_benign": r["n_benign"]} for r in results]
    all_path = [h for r in results for h in r.get("pathogenic_hits", [])]
    # cohort per-gene benign coverage = sum of each sample's per-gene benign counts
    gene_benign = Counter()
    for r in results:
        for g, c in (r.get("benign_by_gene") or {}).items():
            gene_benign[g] += c
    return {"schema": "clinvar-pgp-uk-cohort-card-v0", "n_individuals": n,
            "samples": [r["sample_id"] for r in results], "panel": panel_note,
            "per_sample": per_sample,
            "cohort_n_pathogenic": sum(r["n_pathogenic"] for r in results),
            "cohort_pathogenic_hits": all_path,
            "benign_carrier_load_range": [min(r["n_benign"] for r in results),
                                          max(r["n_benign"] for r in results)] if results else [0, 0],
            "top_genes_by_benign_coverage": dict(gene_benign.most_common(15)),
            "cohort": "PGP-UK (Personal Genome Project UK), ENA PRJEB17529, open-consent GRCh37",
            "honest_tier": ("research demonstration on real open-consent individuals; curated ClinVar allele "
                            "classifications, NOT a clinical diagnosis; 0 reportable pathogenic across a small "
                            "healthy cohort is the expected ACMG-SF base-rate result. NOT a clinical tool.")}


def render_md(card: dict) -> str:
    n = card["n_individuals"]
    lo, hi = card["benign_carrier_load_range"]
    L = [f"# ClinVar/Mendelian real-people COHORT card — {n} PGP-UK individuals (2026-07-06)", "",
         f"_{card['honest_tier']}_", "",
         f"**Cohort:** {card['cohort']} · **N = {n}** · **Panel:** {card['panel']}", "",
         "## Per-individual", "",
         "| sample | reportable pathogenic (P/LP) | benign carrier variants (B/LB) |",
         "|---|---|---|"]
    for s in card["per_sample"]:
        L.append(f"| {s['sample_id']} | {s['n_pathogenic']} | {s['n_benign']} |")
    L += ["",
          f"**Cohort reportable-pathogenic findings: {card['cohort_n_pathogenic']}** "
          f"(0 is the expected result — ACMG-SF pathogenic carriage is ~1–2% per person, so 0 across N={n} "
          f"healthy participants is the most likely outcome). Benign carrier load per person: "
          f"{lo}–{hi} variants across the panel."]
    if card["cohort_pathogenic_hits"]:
        L += ["", "### Pathogenic hits (detail)", ""]
        for h in card["cohort_pathogenic_hits"]:
            L.append(f"- {h['gene']} {h['chrom']}:{h['pos']} {h['ref']}>{h['alt']} — {h['significance']} "
                     f"({h['stars']}★) · {h.get('disease','')}")
    L += ["", "## Top genes by benign-variant coverage (cohort sum)", "",
          "| gene | benign hits (Σ across cohort) |", "|---|---|"]
    for g, c in card["top_genes_by_benign_coverage"].items():
        L.append(f"| {g} | {c} |")
    L += ["", "_Accuracy-vs-truth is not applicable here (PGP-UK ships no ClinVar truth); this is real-world "
          "deployment coverage of the deterministic Mendelian decoder. NOT a clinical tool._", ""]
    return "\n".join(L)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Build a PGP-UK real-people ClinVar cohort card.")
    ap.add_argument("--results-dir", type=Path, required=True)
    ap.add_argument("--glob", default="*_clinvar*.json")
    ap.add_argument("--panel-note", default="committed ClinVar panel")
    ap.add_argument("--out-md", type=Path, required=True)
    ap.add_argument("--out-json", type=Path, default=None)
    args = ap.parse_args(argv)
    results = load(args.results_dir, args.glob)
    if not results:
        print(f"no {args.glob} in {args.results_dir}")
        return 1
    card = build_card(results, args.panel_note)
    args.out_md.write_text(render_md(card), encoding="utf-8")
    (args.out_json or args.out_md.with_suffix(".json")).write_text(json.dumps(card, indent=2), encoding="utf-8")
    print(f"[clinvar cohort card N={card['n_individuals']} -> {args.out_md}]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
