"""`dna-clinvar` — first-class CLI for the deterministic human Mendelian (ClinVar) decoder.

Point at ONE individual VCF -> curated ClinVar germline classifications (P/LP + B/LB) for every variant the
person carries, with gene / disease / gold-star review level. Fail-closed: not-in-panel -> INDETERMINATE
(absence != benign). Build-agnostic (pass --panel for a GRCh37 panel on a GRCh37 VCF). NOT a clinical tool.

    dna-clinvar sample.vcf.gz --sample-id S1                 # committed ACMG-SF-v3.2 panel (GRCh38)
    dna-clinvar sample.vcf.gz --panel data/clinvar/clinvar_panel_grch37.tsv   # GRCh37 individual VCF
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from dna_decode.clinvar.decode import decode_vcf
from dna_decode.data.clinvar import ClinVarDecoder

_HONEST = ("research demonstration; curated ClinVar allele classifications for carried variants, NOT a "
           "clinical diagnosis of the individual. NOT a clinical tool.")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        prog="dna-clinvar",
        description="Deterministic human Mendelian decoder: curated ClinVar germline pathogenicity (P/LP + "
                    "B/LB) over the committed panel, from a VCF. NOT a clinical tool.")
    ap.add_argument("vcf", type=Path, help="individual VCF (.vcf or .vcf.gz)")
    ap.add_argument("--sample-id", default=None)
    ap.add_argument("--panel", type=Path, default=None, help="ClinVar panel TSV/TSV.gz (default: committed GRCh38)")
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--json-only", action="store_true")
    args = ap.parse_args(argv)

    if not args.vcf.exists():
        print(f"ERROR: VCF not found: {args.vcf}")
        return 2
    decoder = ClinVarDecoder.from_tsv(args.panel) if args.panel else ClinVarDecoder.from_tsv()
    rep = decode_vcf(args.vcf, decoder, sample=args.sample_id)
    rep["schema"] = "clinvar-vcf-decode-v0"
    rep["panel"] = str(args.panel) if args.panel else "committed GRCh38 clinvar_panel (ACMG SF v3.2 + carriers)"
    rep["honest_tier"] = _HONEST

    if args.out:
        args.out.write_text(json.dumps(rep, indent=2), encoding="utf-8")
    if args.json_only:
        print(json.dumps({k: rep[k] for k in ("sample_id", "n_pathogenic", "n_benign",
                                              "n_indeterminate_not_in_panel", "pathogenic_hits")}, indent=2))
    else:
        print(f"sample: {rep['sample_id']}  (panel: {rep['panel']})")
        print(f"  reportable pathogenic (P/LP): {rep['n_pathogenic']}   benign carrier (B/LB): {rep['n_benign']}"
              f"   not-in-panel (abstained): {rep['n_indeterminate_not_in_panel']}")
        for h in rep["pathogenic_hits"]:
            print(f"  PATHOGENIC  {h['gene']} {h['chrom']}:{h['pos']} {h['ref']}>{h['alt']} — "
                  f"{h['significance']} ({h['stars']}★) · {h.get('disease','')}")
        print(f"  {_HONEST}")
        if args.out:
            print(f"\n[provenance -> {args.out}]")
    return 0 if (rep["n_pathogenic"] + rep["n_benign"]) >= 1 else 1


if __name__ == "__main__":
    raise SystemExit(main())
