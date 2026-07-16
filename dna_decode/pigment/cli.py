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
    ap.add_argument("--genotypes", required=True,
                    help="comma-separated rsID=GT for the 6 IrisPlex SNPs (rs12913832, rs1800407, rs12896399, "
                         "rs16891982, rs1393350, rs12203592), e.g. rs12913832=GG,rs1800407=TT,...")
    ap.add_argument("--allow-missing", action="store_true",
                    help="impute a missing non-HERC2 SNP as x=0 + cap confidence low (biased; use knowingly)")
    ap.add_argument("--json", action="store_true", help="emit the result as JSON")
    args = ap.parse_args(argv)

    from dna_decode.pigment import MissingGenotypeError, predict_eye_color

    try:
        genos = _parse_genotypes(args.genotypes)
        res = predict_eye_color(genos, allow_missing=args.allow_missing)
    except (ValueError, MissingGenotypeError) as e:
        print(f"error: {e}", file=sys.stderr)
        return 2

    if args.json:
        print(json.dumps(res.as_dict(), indent=2))
        return 0

    d = res.as_dict()
    print("pigmentation decode — eye colour (IrisPlex)")
    print(f"  call: {d['call'].upper()}   confidence: {d['confidence']}")
    print(f"  P(blue)={d['p_blue']:.3f}  P(intermediate)={d['p_intermediate']:.3f}  P(brown)={d['p_brown']:.3f}")
    print(f"  counted alleles: {d['counted_alleles']}")
    for n in d["notes"]:
        print(f"  note: {n}")
    print("  [deterministic Walsh-2011 IrisPlex coefficients (curated); reference-integrity biology-checked]")
    print("  [scope: benign visible-trait genetics, NOT a forensic/surveillance tool]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
