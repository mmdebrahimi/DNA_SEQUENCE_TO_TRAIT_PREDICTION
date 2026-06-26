"""Human pharmacogenomics (CYP2C19) -- in-package CLI (console `dna-pgx`; also `dna-decode pgx`).

Phased VCF (GRCh38) -> CYP2C19 diplotype + CPIC metabolizer phenotype + provenance. Deterministic
variant->catalog caller (sibling of dna-amr's HIV/TB target-site cells). Pure-stdlib VCF parse, no Docker.

    dna-pgx sample.vcf --sample-id MY_SAMPLE                 # CYP2C19 (default)
    dna-pgx sample.vcf --gene cyp2c9                         # CYP2C9 (activity-score)
    dna-pgx sample.vcf --gene vkorc1                         # VKORC1 warfarin sensitivity
    dna-pgx cohort.vcf --sample NA12878 --json-only
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from dna_decode.pgx import PGX_GENES
from dna_decode.pgx.runner import call_cyp2c19, call_cyp2c9


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        prog="dna-pgx",
        description="Deterministic human pharmacogenomics: CYP2C19/CYP2C9 diplotype + CPIC metabolizer "
                    "phenotype, or VKORC1 warfarin sensitivity, from a phased VCF (GRCh38). NOT a clinical tool.")
    ap.add_argument("vcf", type=Path, help="phased VCF (GRCh38)")
    ap.add_argument("--gene", default="cyp2c19", choices=list(PGX_GENES),
                    help="gene to call (default cyp2c19)")
    ap.add_argument("--sample-id", default=None, help="label for the report (default: VCF stem)")
    ap.add_argument("--sample", default=None, help="sample COLUMN name in a multi-sample VCF (default: first)")
    ap.add_argument("--out", type=Path, default=None, help="write provenance JSON here")
    ap.add_argument("--json-only", action="store_true")
    args = ap.parse_args(argv)

    if not args.vcf.exists():
        print(f"ERROR: VCF not found: {args.vcf}", file=sys.stderr)
        return 2

    # VKORC1 is a single-SNP genotype->sensitivity call (not a diplotype) -> its own shape.
    if args.gene == "vkorc1":
        from dna_decode.pgx.vkorc1 import call_vkorc1
        rec = call_vkorc1(args.vcf, sample=args.sample)
        rec["sample_id"] = args.sample_id or args.sample or args.vcf.stem
        if args.out:
            Path(args.out).write_text(json.dumps(rec, indent=2), encoding="utf-8")
        if args.json_only:
            print(json.dumps(rec, indent=2))
        else:
            print(f"sample: {rec['sample_id']}  gene: VKORC1 ({rec['assembly']})  {rec['position']}")
            print(f"  rs9923231 {rec['genomic_ref_alt']}  GT={rec['genomic_gt']}  "
                  f"cDNA {rec['cdna_genotype']}  -> {rec['warfarin_sensitivity']} ({rec['dose_category']})")
            if rec["flags"]:
                print(f"  flags: {', '.join(rec['flags'])}")
            print(f"  {rec['caveat']}")
        return 0

    rec = (call_cyp2c9 if args.gene == "cyp2c9" else call_cyp2c19)(
        args.vcf, sample_id=args.sample_id, sample_column=args.sample)

    if args.out:
        Path(args.out).write_text(json.dumps(rec, indent=2), encoding="utf-8")
    if args.json_only:
        print(json.dumps(rec, indent=2))
    else:
        print(f"sample: {rec['sample_id']}  gene: {rec['gene']} ({rec['assembly']})")
        if rec["status"] != "ok":
            print(f"STATUS: {rec['status']} - {rec.get('reason')}")
        elif rec["phenotype_status"] == "phenotype_withheld":
            hits = ", ".join(f"{h['rsid']}->{h['implies']}" for h in rec["sentinel_hits"])
            print(f"PHENOTYPE WITHHELD (non-core allele detected: {hits})")
            print(f"  core-proxy diplotype: {rec['core_proxy_diplotype']}  (NOT a final call -- a "
                  f"non-core star allele the core SNP set cannot resolve is present)")
        else:
            print(f"DIPLOTYPE: {rec['diplotype']}   PHENOTYPE: {rec['phenotype']} "
                  f"({rec['phenotype_abbrev']})   [phasing: {rec['phasing']}; "
                  f"confidence: {rec['phenotype_confidence']}]")
            if rec["phenotype_status"] == "phase_ambiguous":
                print(f"  PHASE-AMBIGUOUS: alternate resolution = {rec['alternate_diplotype']} "
                      f"({rec['alternate_phenotype']}); reported call assumes trans (the standard).")
        if rec["status"] == "ok":
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
    # exit nonzero when the phenotype is NOT safe to consume (parse failure OR withheld). phase_ambiguous
    # keeps the (flagged) call -> exit 0.
    return 0 if (rec["status"] == "ok" and rec["phenotype_status"] != "phenotype_withheld") else 3


if __name__ == "__main__":
    raise SystemExit(main())
