"""Human pharmacogenomics (CYP2C19) -- in-package CLI (console `dna-pgx`; also `dna-decode pgx`).

Phased VCF (GRCh38) -> CYP2C19 diplotype + CPIC metabolizer phenotype + provenance. Deterministic
variant->catalog caller (sibling of dna-amr's HIV/TB target-site cells). Pure-stdlib VCF parse, no Docker.

    dna-pgx sample.vcf --sample-id MY_SAMPLE
    dna-pgx cohort.vcf --sample NA12878      # pick a sample column by name
    dna-pgx sample.vcf --json-only
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from dna_decode.pgx.runner import call_cyp2c19


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        prog="dna-pgx",
        description="Deterministic human pharmacogenomics: CYP2C19 diplotype + CPIC metabolizer phenotype "
                    "from a phased VCF (GRCh38). NOT a clinical tool.")
    ap.add_argument("vcf", type=Path, help="phased VCF (GRCh38)")
    ap.add_argument("--sample-id", default=None, help="label for the report (default: VCF stem)")
    ap.add_argument("--sample", default=None, help="sample COLUMN name in a multi-sample VCF (default: first)")
    ap.add_argument("--out", type=Path, default=None, help="write provenance JSON here")
    ap.add_argument("--json-only", action="store_true")
    args = ap.parse_args(argv)

    if not args.vcf.exists():
        print(f"ERROR: VCF not found: {args.vcf}", file=sys.stderr)
        return 2

    rec = call_cyp2c19(args.vcf, sample_id=args.sample_id, sample_column=args.sample)

    if args.out:
        Path(args.out).write_text(json.dumps(rec, indent=2), encoding="utf-8")
    if args.json_only:
        print(json.dumps(rec, indent=2))
    else:
        print(f"sample: {rec['sample_id']}  gene: {rec['gene']} ({rec['assembly']})")
        if rec["status"] != "ok":
            print(f"STATUS: {rec['status']} - {rec.get('reason')}")
        else:
            print(f"DIPLOTYPE: {rec['diplotype']}   PHENOTYPE: {rec['phenotype']} "
                  f"({rec['phenotype_abbrev']})   [phasing: {rec['phasing']}]")
            for vc in rec["variant_calls"]:
                state = ("absent" if not vc["found"]
                         else "no-call" if vc["no_call"]
                         else {0: "ref", 1: "het", 2: "hom-alt"}.get(vc["alt_count"], "?"))
                print(f"  {vc['star']:4} {vc['rsid']:12} chr10:{vc['pos']}  {state}"
                      + (f"  GT={vc['gt']}" if vc["gt"] else ""))
            if rec["flags"]:
                print(f"  flags: {', '.join(rec['flags'])}")
            print(f"  {rec['caveat']}")
        if args.out:
            print(f"\n[provenance JSON -> {args.out}]")
    return 0 if rec["status"] == "ok" else 3


if __name__ == "__main__":
    raise SystemExit(main())
