"""`dna-hla` — HLA drug-hypersensitivity carriage from a VCF (GRCh38), via validated tag SNPs.

    dna-hla sample.vcf --allele b5701                 # HLA-B*57:01 / abacavir (the only shipped tag)
    dna-hla cohort.vcf --allele b5701 --sample NA12878 --json-only

Deterministic tag-SNP -> HLA-allele-carriage -> CPIC drug action. Sibling of dna-pgx / dna-clinvar.
NOT full sequence-based HLA typing; NOT a clinical tool.

SHIPPED SURFACE IS ONE ALLELE. B*58:01/allopurinol and A*31:01/carbamazepine were MEASURED against real
1000G HLA truth (2026-07-06) and DEMOTED: rs9263726 misses 39% of B*58:01 carriers (sens 0.61 / PPV 0.18 --
unsafe for an SJS/TEN screen) and rs1061235 is not paneled on 1000G at all (sens 0.0). They are a documented
negative in `catalog._UNVALIDATED_TAGS`, not a routable cell. This docstring used to show `--allele b5801`
and the parser description used to claim "abacavir/allopurinol/carbamazepine" -- an advertisement of exactly
the two screens the project's own validation rejected.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from dna_decode.hla import HLA_ALLELES
from dna_decode.hla.caller import call_hla


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        prog="dna-hla",
        description="Deterministic HLA drug-hypersensitivity carriage via a validated tag SNP, from a VCF "
                    "(GRCh38). Shipped surface is HLA-B*57:01 / abacavir ONLY (rs2395029; validated vs 1000G "
                    "HLA truth: sens 0.979 / spec 0.992 / PPV 0.855, n=1103). The allopurinol (B*58:01) and "
                    "carbamazepine (A*31:01) tags FAILED that same validation and are not shipped. "
                    "NOT sequence-based typing; NOT clinical.")
    ap.add_argument("vcf", type=Path, help="VCF (GRCh38)")
    ap.add_argument("--allele", default="b5701", choices=list(HLA_ALLELES),
                    help="HLA allele to screen (default b5701 = HLA-B*57:01 / abacavir)")
    ap.add_argument("--sample", default=None, help="sample COLUMN name in a multi-sample VCF (default: first)")
    ap.add_argument("--sample-id", default=None, help="label for the report (default: VCF stem)")
    ap.add_argument("--out", type=Path, default=None, help="write provenance JSON here")
    ap.add_argument("--json-only", action="store_true")
    args = ap.parse_args(argv)

    if not args.vcf.exists():
        print(f"ERROR: VCF not found: {args.vcf}", file=sys.stderr)
        return 2

    rec = call_hla(args.vcf, args.allele, sample=args.sample)
    rec["sample_id"] = args.sample_id or args.sample or args.vcf.stem
    if args.out:
        Path(args.out).write_text(json.dumps(rec, indent=2), encoding="utf-8")
    if args.json_only:
        print(json.dumps(rec, indent=2))
    else:
        print(f"sample: {rec['sample_id']}  {rec['allele']} ({rec['assembly']})  drug: {rec['drug']}")
        print(f"  tag {rec['tag_rsid']} {rec['position']} {rec['tag_ref_alt']}  GT={rec['tag_gt']}  "
              f"copies={rec['tag_copies']}  -> {rec['zygosity']}")
        print(f"  RISK: {rec['risk_call']}  [{rec['reaction']}; proxy tier: {rec['proxy_tier']}]")
        if rec["flags"]:
            print(f"  flags: {', '.join(rec['flags'])}")
        print(f"  {rec['caveat']}")
        if args.out:
            print(f"\n[provenance JSON -> {args.out}]")
    return 0 if rec["status"] in ("ok", "assumed_reference") else 3


if __name__ == "__main__":
    raise SystemExit(main())
