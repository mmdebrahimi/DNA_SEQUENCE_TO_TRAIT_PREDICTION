"""Mine the GWAS Catalog associations for STRONG SINGLE-SNP trait associations — the product-aligned use.

The GWAS Catalog full table is polygenic summary stats (the deterministic decoder does not consume it
wholesale). Its ONE product-aligned use: filter for strong, single-locus, genome-wide-significant
associations that could SEED new deterministic single-locus cells (the lactase/earwax/eye-colour pattern).
Emits a committed candidate shortlist ranked by effect size. Streams the 716 MB TSV on D: (never loads it
all). NOT a decoder integration — a candidate-discovery capture.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
DEFAULT_TSV = Path("D:/dna_decode_cache/gwas_catalog/gwas-catalog-download-associations-alt-full.tsv")
_RSID = re.compile(r"^rs\d+$")
P_MAX = 5e-8          # genome-wide significance (the CLEAN filter — p is reliable)
OR_MIN = 2.0          # strong effect — the oligogenic/near-Mendelian tail
# GWAS Catalog "OR or BETA" is a DIRTY free field (mixes OR / beta / malformed values up to ~1e14). Real
# single-SNP trait ORs are rarely > ~20 (except ABO-class). Cap to exclude data-entry artifacts so the
# shortlist is genuine strong-effect biology, not OR-field noise. (verify-in-batch fix, 2026-07-04.)
OR_MAX = 20.0


def _f(x: str) -> float | None:
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def mine(tsv: Path, out: Path, top: int = 400) -> dict:
    cand = []
    n_rows = n_single = 0
    with open(tsv, encoding="utf-8", errors="replace") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        idx = {h: i for i, h in enumerate(hdr)}
        need = ["SNPS", "DISEASE/TRAIT", "MAPPED_GENE", "P-VALUE", "OR or BETA", "STRONGEST SNP-RISK ALLELE"]
        if any(k not in idx for k in need):
            return {"error": f"missing columns; have {hdr[:8]}"}
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) <= idx["OR or BETA"]:
                continue
            n_rows += 1
            snp = c[idx["SNPS"]].strip()
            if not _RSID.match(snp):                          # single rsID only (drop haplotypes/interactions)
                continue
            n_single += 1
            p = _f(c[idx["P-VALUE"]])
            orb = _f(c[idx["OR or BETA"]])
            if p is None or p > P_MAX or orb is None or orb < OR_MIN or orb > OR_MAX:
                continue
            cand.append({"rsid": snp, "trait": c[idx["DISEASE/TRAIT"]][:70],
                         "gene": c[idx["MAPPED_GENE"]][:30], "p_value": p, "or_or_beta": orb,
                         "risk_allele": c[idx["STRONGEST SNP-RISK ALLELE"]][:20]})
    cand.sort(key=lambda r: -r["or_or_beta"])
    top_cand = cand[:top]
    out.parent.mkdir(parents=True, exist_ok=True)
    cols = ["rsid", "trait", "gene", "or_or_beta", "p_value", "risk_allele"]
    with open(out, "w", encoding="utf-8", newline="") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in top_cand:
            fh.write("\t".join(str(r[k]).replace("\t", " ") for k in cols) + "\n")
    return {"n_rows_scanned": n_rows, "n_single_snp": n_single, "n_strong_candidates": len(cand),
            "n_written": len(top_cand), "out": str(out),
            "filter": f"single rsID + P<={P_MAX} + OR/beta>={OR_MIN}",
            "top5": [{"rsid": r["rsid"], "trait": r["trait"], "gene": r["gene"], "or": r["or_or_beta"]}
                     for r in top_cand[:5]]}


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tsv", type=Path, default=DEFAULT_TSV)
    ap.add_argument("--out", type=Path, default=REPO / "data" / "summary_stat_sources" / "gwas_single_snp_candidates.tsv")
    ap.add_argument("--top", type=int, default=400)
    a = ap.parse_args(argv)
    if not a.tsv.exists():
        print(f"ERROR: GWAS associations TSV not found at {a.tsv}", file=sys.stderr)
        return 2
    res = mine(a.tsv, a.out, a.top)
    import json
    print(json.dumps(res, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
