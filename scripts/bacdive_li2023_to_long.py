"""Convert the Li et al. 2023 (PMC10729968) BacDive carbon-utilization release
into the long-format export `dna_decode/data/bacdive.py` expects.

Source data (public, OSF project jwkr7 — https://osf.io/jwkr7/):
  - bacdive_growth_data_final.csv : WIDE matrix, rows = BacDive id, cols = 58
    carbon sources, cells in {0.0, 1.0, blank}. The paper's finalized growth data.
  - bacdive_16S.fna : FASTA whose headers carry `>bacdive_<id> ENA|acc|<Species> 16S...`
    — the ONLY species label in the release (the CSVs have no organism column).

This adapter melts the wide matrix to long (strain_id × carbon_source × utilization),
joins the species from the 16S headers, and writes a long CSV with the columns the
loader maps tolerantly. assembly_accession is left EMPTY — the release links genomes
via proGenomes biosample IDs (in bacdive.pk, 298 MB), not GCF/GCA; resolving those is
a downstream step only worth doing if the count/balance feasibility gate passes.

Run:
  uv run python scripts/bacdive_li2023_to_long.py \
    --growth data/raw/carbon_util_bacdive/bacdive_growth_data_final.csv \
    --fasta16s data/raw/carbon_util_bacdive/bacdive_16S.fna \
    --out data/raw/carbon_util_bacdive/bacdive_carbon_long.csv
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd

_HDR = re.compile(r">bacdive_(\d+)\s+\S+\|.*?\|(.+?)\s+(?:partial\s+)?16S", re.IGNORECASE)
_HDR2 = re.compile(r">bacdive_(\d+)\s+.*?\|(.+?)\s+16S", re.IGNORECASE)


def species_map_from_16s(fasta_path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    with open(fasta_path, encoding="utf-8", errors="replace") as f:
        for ln in f:
            if not ln.startswith(">"):
                continue
            m = _HDR.match(ln) or _HDR2.match(ln)
            if m:
                out[m.group(1)] = m.group(2).strip()
    return out


def convert(growth_path: Path, fasta_path: Path, out_path: Path) -> dict:
    sp = species_map_from_16s(fasta_path)
    g = pd.read_csv(growth_path, dtype=str)
    g.columns = [c.strip() for c in g.columns]
    idcol = g.columns[0]
    g[idcol] = g[idcol].astype(str).str.strip()
    carbons = [c for c in g.columns if c != idcol]

    rows: list[dict] = []
    for _, r in g.iterrows():
        sid = r[idcol]
        organism = sp.get(sid, "")
        for c in carbons:
            v = r[c]
            if v is None or str(v).strip() == "":
                continue
            try:
                iv = int(float(v))
            except ValueError:
                continue
            if iv not in (0, 1):
                continue
            rows.append({
                "strain_id": f"bacdive_{sid}",
                "carbon_source": c.strip().lower(),
                "utilization": iv,
                "organism": organism,
                "assembly_accession": "",
            })
    out = pd.DataFrame(rows, columns=["strain_id", "carbon_source", "utilization",
                                      "organism", "assembly_accession"])
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, index=False)
    n_strains = out["strain_id"].nunique()
    n_ecoli = out[out["organism"].str.contains("Escherichia coli", case=False, na=False)]["strain_id"].nunique()
    return {
        "long_rows": len(out),
        "distinct_strains": n_strains,
        "distinct_carbon_sources": out["carbon_source"].nunique(),
        "strains_with_species": int((out["organism"].str.len() > 0).groupby(out["strain_id"]).any().sum()),
        "ecoli_strains": n_ecoli,
        "out_path": str(out_path),
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--growth", required=True, type=Path)
    ap.add_argument("--fasta16s", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args(argv)
    summary = convert(args.growth, args.fasta16s, args.out)
    for k, v in summary.items():
        print(f"{k}: {v}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
