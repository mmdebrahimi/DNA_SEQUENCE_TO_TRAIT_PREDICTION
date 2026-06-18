"""Populate the FULL per-drug CRyPTIC cohort (masked + regeno VCFs) to D: for the v1b baseline.

Ratified B: the baseline label needs the full per-drug PREVALENCE-PRESERVING HIGH-quality eligible
cohort (NOT a balanced subset). Ratified full-regeno (C3): callability needs the per-isolate regeno VCF
(~177 MB each) -> the full cohort is ~1.6 TB, so it MUST land on D: (4.4 TB free), never C:.

Eligible = HIGH phenotype quality + R/S + has both a VCF and a REGENOTYPED_VCF path. ALL eligible are
taken (prevalence-preserving). Restartable / skip-existing (the D: cache dedups by mangled rel-path).

Sizing: masked ~217 KB/isolate (~2 GB/cohort, ~hours); regeno ~177 MB/isolate (~1.6 TB/cohort, ~days).
Run masked first (fast -> enables a fuller v1b run with the callability caveat), regeno as the long job.

  uv run python scripts/populate_tb_cohort.py --drug rifampicin --kind masked
  uv run python scripts/populate_tb_cohort.py --drug rifampicin --kind regeno   # the ~1.6 TB job
"""
from __future__ import annotations

import argparse
from datetime import date as _date
from pathlib import Path

from dna_decode.organism_rules import tb_vcf

DEFAULT_CACHE = tb_vcf.DEFAULT_CACHE  # D:/dna_decode_cache/cryptic


def eligible_cohort(drug_code: str) -> list[dict]:
    """All HIGH-quality R/S isolates with both VCF + REGENOTYPED_VCF paths (prevalence-preserving)."""
    ph, q = f"{drug_code}_BINARY_PHENOTYPE", f"{drug_code}_PHENOTYPE_QUALITY"
    out = []
    for r in tb_vcf.reuse_rows():
        call = (r.get(ph) or "").strip().upper()
        qual = (r.get(q) or "").strip().upper()
        masked, regeno = tb_vcf.vcf_paths_for(r)
        if call in ("R", "S") and qual == "HIGH" and masked and regeno:
            out.append(r)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--drug", default="rifampicin")
    ap.add_argument("--kind", choices=["masked", "regeno", "both"], default="masked")
    ap.add_argument("--cache", default=str(DEFAULT_CACHE), help="D: cohort cache root")
    ap.add_argument("--max", type=int, default=0, help="cap isolates (0 = all eligible)")
    a = ap.parse_args()

    code = {"rifampicin": "RIF", "isoniazid": "INH"}[a.drug]
    cohort = eligible_cohort(code)
    if a.max:
        cohort = cohort[: a.max]
    n_r = sum(1 for r in cohort if (r.get(f"{code}_BINARY_PHENOTYPE") or "").upper() == "R")
    cache = Path(a.cache)
    kinds = ["masked", "regeno"] if a.kind == "both" else [a.kind]
    print(f"{_date.today().isoformat()} {a.drug}: {len(cohort)} eligible HIGH-quality isolates "
          f"({n_r}R/{len(cohort) - n_r}S, prevalence-preserving) -> {cache} kinds={kinds}")
    if "regeno" in kinds:
        print(f"  WARNING: regeno is ~177 MB/isolate -> ~{len(cohort) * 0.177:.0f} GB. Ensure D: space.")

    for kind in kinds:
        ok = miss = 0
        for i, r in enumerate(cohort, 1):
            masked, regeno = tb_vcf.vcf_paths_for(r)
            rel = masked if kind == "masked" else regeno
            data = tb_vcf.fetch_vcf(rel, kind, cache_dir=cache)
            ok += data is not None
            miss += data is None
            if i % 50 == 0 or i == len(cohort):
                print(f"  [{kind}] {i}/{len(cohort)} ok={ok} miss={miss}", flush=True)
        print(f"  [{kind}] done: {ok} cached, {miss} miss.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
