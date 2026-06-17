"""Stage a bounded, balanced CRyPTIC VCF subset into the local cache (acquisition-only).

The TB-decoder plan needs a bounded VCF subset on disk for the eventual v1b run (the full 12,287
cohort ~2.6 GB is too large for this host's chronically-full C:). This is a THIN downloader — it
ONLY populates `data/raw/cryptic/vcf_cache/` (gitignored), reusing the proven fetch + cache logic in
`scripts/cryptic_feasibility_probe.py`. It deliberately contains NO genotype/parse logic (that is the
plan's Stage-0 `tb_vcf.py` parser module). Restartable + idempotent: cached files are skipped.

Default = a balanced 150 R / 150 S rifampicin HIGH-quality set (~300 VCFs, ~65 MB). RIF-R isolates are
frequently MDR, so INH-R gets incidental coverage too. Bounded subset per the plan's Open Question B
default (exact final size is a ratification call).

Run:
  uv run python scripts/stage_tb_vcf_subset.py                  # 150/class RIF
  uv run python scripts/stage_tb_vcf_subset.py --drug RIF --per-class 150
Exit: 0 = done (or already complete), 2 = reuse table absent.
"""
from __future__ import annotations

import argparse
import sys
from datetime import date as _date
from pathlib import Path

# Reuse the feasibility probe's table loader + fetch/cache (single source of truth for the FTP path).
from scripts.cryptic_feasibility_probe import load_rows, fetch_vcf, REUSE_CSV


def pick_balanced(rows: list[dict], drug_code: str, per_class: int) -> list[dict]:
    """First `per_class` HIGH-quality R and S isolates (with a VCF path) for a drug code."""
    ph, q = f"{drug_code}_BINARY_PHENOTYPE", f"{drug_code}_PHENOTYPE_QUALITY"
    picked: dict[str, list[dict]] = {"R": [], "S": []}
    for row in rows:
        call = (row.get(ph) or "").strip().upper()
        qual = (row.get(q) or "").strip().upper()
        vcf = (row.get("VCF") or "").strip()
        if call in ("R", "S") and qual == "HIGH" and vcf and len(picked[call]) < per_class:
            picked[call].append(row)
        if len(picked["R"]) >= per_class and len(picked["S"]) >= per_class:
            break
    return picked["R"] + picked["S"]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--drug", default="RIF", help="CRyPTIC drug code (RIF, INH, ...)")
    ap.add_argument("--per-class", type=int, default=150)
    a = ap.parse_args()
    if not REUSE_CSV.exists():
        print(f"missing {REUSE_CSV}")
        return 2
    rows = load_rows()
    sample = pick_balanced(rows, a.drug, a.per_class)
    print(f"{_date.today().isoformat()} staging {len(sample)} {a.drug} VCFs "
          f"({a.per_class}/class HIGH-quality) into the cache ...")
    ok = miss = 0
    for i, row in enumerate(sample, 1):
        data = fetch_vcf(row["VCF"].strip())
        if data is None:
            miss += 1
        else:
            ok += 1
        if i % 25 == 0 or i == len(sample):
            print(f"  {i}/{len(sample)}  ok={ok} miss={miss}", flush=True)
    print(f"done: {ok} cached, {miss} fetch-miss.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
