"""`dna-decode pigment` — visible-trait pigmentation decoder (v0 = IrisPlex eye colour).

    dna-decode pigment --genotypes rs12913832=GG,rs1800407=TT,rs12896399=GG,rs16891982=GG,rs1393350=GG,rs12203592=CC
    dna-pigment --genotypes rs12913832=AA,... --json

Deterministic multinomial-logistic model over 6 curated SNPs (published coefficients) -> P(blue/intermediate/
brown) + a category call. Pure-python, wheel-only, offline. Scope: benign visible-trait genetics, NOT forensic.
Hair/skin (full HIrisPlex-S 41-SNP) + VCF input + openSNP scoring = documented v0.1 follow-ons.
"""
from __future__ import annotations

import argparse
import json
import sys


def _parse_genotypes(spec: str) -> dict:
    out = {}
    for tok in spec.split(","):
        tok = tok.strip()
        if not tok:
            continue
        if "=" not in tok:
            raise ValueError(f"bad genotype token {tok!r}; expected rsID=GT (e.g. rs12913832=GG)")
        rsid, gt = tok.split("=", 1)
        out[rsid.strip()] = gt.strip()
    return out


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    ap = argparse.ArgumentParser(
        prog="dna-decode pigment",
        description="Visible-trait pigmentation decoder (v0 = IrisPlex eye colour, 6 SNPs, deterministic).",
        epilog="scope: benign visible-trait genetics, NOT a forensic/surveillance tool.",
    )
    ap.add_argument("--trait", choices=["eye", "hair", "skin"], default="eye",
                    help="which pigmentation trait to predict (default eye = IrisPlex 6-SNP). hair/skin use "
                         "the HIrisPlex-S models recovered + held-out-validated from the erasmusmc webtool.")
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--genotypes",
                     help="comma-separated rsID=GT. eye: the 6 IrisPlex SNPs; hair/skin: the HIrisPlex-S "
                          "panel (22/41 SNPs) — e.g. rs12913832=GG,rs1805007=CC,...")
    src.add_argument("--vcf", help="a genome VCF (.vcf/.vcf.gz); the trait's SNPs are extracted by rsID")
    ap.add_argument("--allow-missing", action="store_true",
                    help="impute a missing SNP as x=0 + cap confidence low (biased; use knowingly)")
    ap.add_argument("--json", action="store_true", help="emit the result as JSON")
    args = ap.parse_args(argv)

    from dna_decode.pigment import MissingGenotypeError, predict_eye_color
    from dna_decode.pigment.vcf_input import genotypes_from_vcf

    try:
        if args.trait == "eye":
            genos = genotypes_from_vcf(args.vcf) if args.vcf else _parse_genotypes(args.genotypes)
            d = predict_eye_color(genos, allow_missing=args.allow_missing).as_dict()
        else:
            from dna_decode.pigment.hirisplex_models import load_hirisplex_models
            from dna_decode.pigment.multinomial import predict as mn_predict
            model = load_hirisplex_models()[f"{args.trait}_colour"]
            if args.vcf:
                genos = genotypes_from_vcf(args.vcf, [(s.rsid, s.counted_allele) for s in model.snps])
            else:
                genos = _parse_genotypes(args.genotypes)
            d = mn_predict(model, genos, allow_missing=args.allow_missing).as_dict()
    except FileNotFoundError as e:
        print(f"error: cannot read --vcf: {e}", file=sys.stderr)
        return 2
    except (ValueError, MissingGenotypeError) as e:
        print(f"error: {e}", file=sys.stderr)
        return 2

    if args.json:
        print(json.dumps(d, indent=2))
        return 0

    if args.trait == "eye":
        print("pigmentation decode — eye colour (IrisPlex)")
        print(f"  call: {d['call'].upper()}   confidence: {d['confidence']}")
        print(f"  P(blue)={d['p_blue']:.3f}  P(intermediate)={d['p_intermediate']:.3f}  P(brown)={d['p_brown']:.3f}")
        print(f"  counted alleles: {d['counted_alleles']}")
        for n in d["notes"]:
            print(f"  note: {n}")
        print("  [deterministic Walsh-2011 IrisPlex coefficients (curated); reference-integrity biology-checked]")
    else:
        print(f"pigmentation decode — {args.trait} colour (HIrisPlex-S)")
        print(f"  call: {d['call'].upper()}   confidence: {d['confidence']}")
        print("  probabilities: " + "  ".join(f"P({k})={v:.3f}" for k, v in d["probabilities"].items()))
        for n in d["notes"]:
            print(f"  note: {n}")
        print("  [multinomial model recovered + held-out-validated from the HIrisPlex-S webtool (2026-07-30)]")
    print("  [scope: benign visible-trait genetics, NOT a forensic/surveillance tool]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
