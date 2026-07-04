"""Capture the bounded INDEX/MANIFEST of the polygenic summary-stat sources (PGS Catalog / Pan-UKBB / FinnGen).

User-directed (2026-07-04): capture GWAS Catalog + FinnGen + PGS + Pan-UKBB + OpenGWAS. These are polygenic
prior-layer sources the deterministic determinant decoder does NOT directly consume (honest caveat recorded
in wiki/free_data_capture_registry_2026-07-04.md) — but the user wants them AVAILABLE. Judgment on HOW: the
full per-phenotype summary-stat bulk is TB-scale + not decoder-consumable, so we capture each source's
bounded, useful INDEX (which traits/scores/SNPs exist) to a committed compact TSV, with the full-bulk fetch
recipe. Raw indexes cache to D: (C: disk-tight). GWAS Catalog (associations table) + OpenGWAS (auth-gated)
are handled separately.
"""
from __future__ import annotations

import csv
import gzip
import io
import json
import sys
import urllib.request
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
OUT = REPO / "data" / "summary_stat_sources"
DCACHE = Path("D:/dna_decode_cache")


def _get(url: str, timeout: int = 60) -> bytes | None:
    try:
        with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "dna_decode/1.0"}), timeout=timeout) as r:
            return r.read()
    except Exception as e:
        print(f"  WARN fetch failed {url}: {e}", file=sys.stderr)
        return None


def capture_pgs(out: Path, page: int = 250) -> dict:
    """PGS Catalog: paginate /rest/score/all -> committed index (id, name, trait, n_variants, ftp)."""
    rows, offset, total = [], 0, None
    while True:
        b = _get(f"https://www.pgscatalog.org/rest/score/all?limit={page}&offset={offset}")
        if not b:
            break
        d = json.loads(b)
        total = d.get("count", total)
        for s in d.get("results", []):
            rows.append({"pgs_id": s.get("id"), "name": (s.get("name") or "")[:60],
                         "trait": (s.get("trait_reported") or "")[:80],
                         "n_variants": (s.get("variants_number") or ""),
                         "ftp_scoring_file": s.get("ftp_scoring_file") or ""})
        offset += page
        if not d.get("next") or offset > (total or 0):
            break
    _write_tsv(out, ["pgs_id", "name", "trait", "n_variants", "ftp_scoring_file"], rows)
    return {"source": "PGS Catalog", "n": len(rows), "total_reported": total, "out": str(out),
            "bulk_recipe": "per-score harmonized scoring files at each row's ftp_scoring_file (PGS Catalog FTP)"}


def capture_panukbb(out: Path) -> dict:
    """Pan-UKBB: fetch the phenotype manifest (bgz; gzip-readable) -> committed compact index."""
    raw = DCACHE / "pan_ukbb" / "phenotype_manifest.tsv.bgz"
    if not raw.exists():
        b = _get("https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_release/phenotype_manifest.tsv.bgz", timeout=300)
        if b:
            raw.parent.mkdir(parents=True, exist_ok=True)
            raw.write_bytes(b)
    if not raw.exists():
        return {"source": "Pan-UKBB", "n": 0, "note": "manifest fetch failed"}
    txt = gzip.open(raw, "rt", encoding="utf-8", errors="replace")
    rdr = csv.DictReader(txt, delimiter="\t")
    keep = ["trait_type", "phenocode", "description", "pheno_sex", "n_cases_full_cohort_both_sexes"]
    rows = [{k: (r.get(k) or "")[:80] for k in keep} for r in rdr]
    _write_tsv(out, keep, rows)
    return {"source": "Pan-UKBB", "n": len(rows), "out": str(out),
            "bulk_recipe": "per-phenotype multi-ancestry sumstats via each row's manifest aws/gs path (TB-scale)"}


def capture_finngen(out: Path) -> dict:
    """FinnGen: the R12 manifest already on D: -> committed compact index."""
    src = DCACHE / "finngen" / "finngen_R12_manifest.tsv"
    if not src.exists():
        return {"source": "FinnGen", "n": 0, "note": "R12 manifest not on D: — fetch from the FinnGen public release"}
    rdr = csv.DictReader(src.open(encoding="utf-8", errors="replace"), delimiter="\t")
    keep = ["phenocode", "phenotype", "category", "num_cases", "num_controls"]
    rows = [{k: (r.get(k) or "")[:80] for k in keep} for r in rdr]
    _write_tsv(out, keep, rows)
    return {"source": "FinnGen R12", "n": len(rows), "out": str(out),
            "bulk_recipe": "per-endpoint GWAS sumstats at gs://finngen-public-data-r12/summary_stats/ (TB-scale)"}


def _write_tsv(out: Path, cols: list[str], rows: list[dict]):
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", encoding="utf-8", newline="") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r.get(c, "")).replace("\t", " ").replace("\n", " ") for c in cols) + "\n")


def main(argv=None) -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    res = [capture_pgs(OUT / "pgs_catalog_index.tsv"),
           capture_panukbb(OUT / "pan_ukbb_phenotype_index.tsv"),
           capture_finngen(OUT / "finngen_r12_endpoint_index.tsv")]
    for r in res:
        print(f"  {r['source']:12} n={r.get('n')} -> {r.get('out', r.get('note'))}")
    (OUT / "capture_summary.json").write_text(json.dumps(res, indent=2), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
