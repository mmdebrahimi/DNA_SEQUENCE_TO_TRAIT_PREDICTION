"""Ingest FinnGen (Finnish biobank) GWAS summary statistics -> a compact HOST-genetics locus-priors table.

FinnGen publishes free, no-auth GWAS summary stats (R12: 2,470 endpoints) at a public GCS bucket. This is a
HOST (human) genetics resource -- it maps human variants -> disease-endpoint association. For dna_decode it is
a locus-PRIORS / feature layer for a possible future human-variant arm; it does NOT feed the pathogen decoder
(that is organism/pathogen genetics -- a different axis). Banked, honestly scoped, NOT overclaimed.

Streams a per-endpoint .gz (schema: #chrom pos ref alt rsids nearest_genes pval mlogp beta sebeta af_* ),
filters to genome-wide-significant variants (default p < 5e-8), clumps to the lead variant per nearest_gene,
and emits a small priors table. The 810 MB per-endpoint .gz is gitignored-class (D: cache); the emitted
priors table is tiny + committable.
"""
from __future__ import annotations

import csv
import gzip
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
GW_SIG_MLOGP = 7.30103  # -log10(5e-8)


def ingest(gz_path: Path, endpoint: str, n_cases=None, n_controls=None, mlogp_min=GW_SIG_MLOGP) -> dict:
    """Stream the summary-stat gz; keep GW-sig variants; clump to lead-variant-per-gene."""
    lead_by_gene: dict[str, dict] = {}
    n_total = n_sig = 0
    with gzip.open(gz_path, "rt") as fh:
        header = fh.readline().lstrip("#").rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in fh:
            n_total += 1
            f = line.rstrip("\n").split("\t")
            try:
                mlogp = float(f[idx["mlogp"]])
            except (ValueError, IndexError):
                continue
            if mlogp < mlogp_min:
                continue
            n_sig += 1
            gene = f[idx["nearest_genes"]] or "NA"
            rec = {"variant": f"{f[idx['chrom']]}:{f[idx['pos']]}:{f[idx['ref']]}:{f[idx['alt']]}",
                   "rsid": f[idx["rsids"]], "gene": gene, "pval": f[idx["pval"]], "mlogp": round(mlogp, 2),
                   "beta": float(f[idx["beta"]]), "direction": "risk" if float(f[idx["beta"]]) > 0 else "protective"}
            if gene not in lead_by_gene or mlogp > lead_by_gene[gene]["mlogp"]:
                lead_by_gene[gene] = rec
    leads = sorted(lead_by_gene.values(), key=lambda r: -r["mlogp"])
    return {"_schema": "finngen-locus-priors-v1", "endpoint": endpoint, "release": "R12",
            "resource": "FinnGen (host/human genetics) -- free, no-auth GWAS summary statistics",
            "scope_note": ("HOST-genetics locus priors for a future human-variant arm; does NOT feed the "
                           "pathogen/organism decoder (different genetics axis). Banked, not load-bearing."),
            "n_cases": n_cases, "n_controls": n_controls, "n_variants": n_total,
            "n_genome_wide_sig": n_sig, "gw_sig_threshold": "p < 5e-8", "n_gene_loci": len(leads),
            "lead_loci_by_gene": leads}


def _manifest_row(manifest: Path, endpoint: str) -> dict | None:
    with open(manifest, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["phenocode"] == endpoint:
                return r
    return None


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gz", type=Path, required=True, help="per-endpoint summary-stat .gz")
    ap.add_argument("--endpoint", required=True)
    ap.add_argument("--manifest", type=Path, default=None, help="R12 manifest tsv (for case/control counts)")
    ap.add_argument("--out", type=Path, default=None)
    a = ap.parse_args(argv)
    nc = ncc = None
    if a.manifest and a.manifest.exists():
        row = _manifest_row(a.manifest, a.endpoint)
        if row:
            nc, ncc = row.get("num_cases"), row.get("num_controls")
    res = ingest(a.gz, a.endpoint, nc, ncc)
    out = a.out or REPO / "wiki" / f"finngen_priors_{a.endpoint}.json"
    out.write_text(json.dumps(res, indent=2), encoding="utf-8")
    print(f"{a.endpoint}: {res['n_variants']} variants, {res['n_genome_wide_sig']} GW-sig, "
          f"{res['n_gene_loci']} gene-loci -> {out}")
    for r in res["lead_loci_by_gene"][:8]:
        print(f"   {r['gene']:12} {r['variant']:22} {r['rsid']:14} mlogp={r['mlogp']} {r['direction']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
