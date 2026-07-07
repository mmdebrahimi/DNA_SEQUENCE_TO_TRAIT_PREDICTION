#!/usr/bin/env python
"""Build a real-people PGx COHORT card from per-individual PGP-UK decode JSONs (scripts/pgx_decode_pgp_uk.py).

Aggregates N individuals' `*_pgx.json` results into a cohort summary: the per-gene diplotype + metabolizer
phenotype distribution, observed star-allele counts (a small-N allele-frequency readout), and the
VKORC1/SLCO1B1 genotype spread. Emits both a markdown card + a JSON.

HONEST FRAMING: a deployment/robustness demonstration on real independent-cohort individuals; the allele
frequencies are a tiny-N observation (not a population estimate), and there is no GeT-RM truth here (accuracy
lives in the concordance cells). NOT a clinical tool.

Usage:
    uv run python scripts/pgx_pgp_uk_cohort_card.py --results-dir C:/.../pgp_results --out-md wiki/pgx_pgp_uk_cohort_card_2026-07-06.md
"""
from __future__ import annotations

import argparse
import glob
import json
import re
from collections import Counter
from pathlib import Path

DIPLO_GENES = ["cyp2c19", "cyp2c9", "cyp2d6", "cyp3a5", "tpmt", "cyp2b6"]
SNP_GENES = ["vkorc1", "slco1b1"]


def load_results(results_dir: Path) -> list[dict]:
    out = []
    for f in sorted(glob.glob(str(results_dir / "*_pgx.json"))):
        d = json.loads(Path(f).read_text(encoding="utf-8"))
        if "results" in d and "sample_id" in d:
            out.append(d)
    return out


def build_card(results: list[dict]) -> dict:
    n = len(results)
    per_gene = {}
    for g in DIPLO_GENES:
        diplos, phenos, alleles = Counter(), Counter(), Counter()
        for r in results:
            cell = r["results"].get(g, {})
            dp = cell.get("diplotype")
            if dp:
                diplos[dp] += 1
                for a in re.findall(r"\*\w+", dp):
                    alleles[a] += 1
            ph = cell.get("phenotype_abbrev") or ("(calling-only)" if g == "cyp2c8" else "-")
            phenos[ph] += 1
        per_gene[g] = {"diplotypes": dict(diplos.most_common()),
                       "phenotypes": dict(phenos.most_common()),
                       "observed_allele_counts": dict(alleles.most_common())}
    for g in SNP_GENES:
        gts = Counter()
        for r in results:
            gts[r["results"].get(g, {}).get("genotype", "?")] += 1
        per_gene[g] = {"genotypes": dict(gts.most_common())}
    return {"schema": "pgx-pgp-uk-cohort-card-v0", "n_individuals": n,
            "samples": [r["sample_id"] for r in results], "per_gene": per_gene,
            "cohort": "PGP-UK (Personal Genome Project UK), ENA PRJEB17529, open-consent GRCh37",
            "honest_tier": ("deployment/robustness demonstration on real independent-cohort individuals; "
                            "tiny-N observed allele counts (NOT a population estimate); no GeT-RM truth here "
                            "(accuracy lives in the concordance cells). NOT a clinical tool.")}


def render_md(card: dict) -> str:
    n = card["n_individuals"]
    L = [f"# PGx real-people COHORT card — {n} PGP-UK individuals ({card.get('date','2026-07-06')})", "",
         f"_{card['honest_tier']}_", "",
         f"**Cohort:** {card['cohort']} · **N = {n}** · samples: {', '.join(card['samples'])}", ""]
    # per-individual table
    L += ["## Per-individual calls", "",
          "| sample | " + " | ".join(g.upper() for g in DIPLO_GENES + SNP_GENES) + " |",
          "|" + "---|" * (len(DIPLO_GENES) + len(SNP_GENES) + 1)]
    for s in card["samples"]:
        L.append("| " + s + " | " + " | ".join(
            _cell(card, g, s) for g in DIPLO_GENES + SNP_GENES) + " |")
    # per-gene distributions
    L += ["", "## Per-gene distribution (across the cohort)", ""]
    for g in DIPLO_GENES:
        pg = card["per_gene"][g]
        L.append(f"- **{g.upper()}** — phenotypes: {pg['phenotypes']}; diplotypes: {pg['diplotypes']}; "
                 f"observed allele counts: {pg['observed_allele_counts']}")
    for g in SNP_GENES:
        L.append(f"- **{g.upper()}** — genotypes: {card['per_gene'][g]['genotypes']}")
    L += ["", "_**CYP2D6 caveat (load-bearing honesty):** the CYP2D6 call here is a SNP-proxy diplotype from "
          "a called VCF — it CANNOT see the structural alleles (*5 deletion / *xN duplication / *13/*36/*68 "
          "hybrids), so every CYP2D6 cell carries `cnv_hybrid_unassessed`. The copy-number half is resolvable "
          "only from a BAM/CRAM (dna_decode.pgx.cyp2d6_structural, 26/26 on 1000G CRAMs); PGP-UK ships VCFs, "
          "not reads._",
          "", "_Accuracy-vs-truth lives in the GeT-RM concordance cells (`wiki/pgx_report_card.md`); this "
          "card is real-world deployment coverage, not a new accuracy number._", ""]
    return "\n".join(L)


# the per-individual cells are re-derived from the JSONs at render time (kept out of the card dict for size)
_RESULTS_CACHE: dict = {}


def _cell(card, gene, sample):
    r = _RESULTS_CACHE.get(sample, {}).get("results", {}).get(gene, {})
    if gene in SNP_GENES:
        return r.get("genotype", "?")
    dp = r.get("diplotype") or "null"
    ab = r.get("phenotype_abbrev") or "-"
    return f"{dp} {ab}"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Build a PGP-UK real-people PGx cohort card.")
    ap.add_argument("--results-dir", type=Path, required=True)
    ap.add_argument("--out-md", type=Path, required=True)
    ap.add_argument("--out-json", type=Path, default=None)
    args = ap.parse_args(argv)
    results = load_results(args.results_dir)
    if not results:
        print(f"no *_pgx.json in {args.results_dir}")
        return 1
    for r in results:
        _RESULTS_CACHE[r["sample_id"]] = r
    card = build_card(results)
    args.out_md.write_text(render_md(card), encoding="utf-8")
    (args.out_json or args.out_md.with_suffix(".json")).write_text(json.dumps(card, indent=2), encoding="utf-8")
    print(f"[cohort card N={card['n_individuals']} -> {args.out_md}]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
