"""`dna-decode catcolor` — cat coat-colour decoder (curated-catalog + X-linked epistasis, v0).

    dna-decode catcolor --loci O=o,A=a/a,B=B/B,D=D/D                 # black (non-orange male)
    dna-decode catcolor --loci O=O/o,B=B/B --sex female             # tortoiseshell (X-linked mosaic)
    dna-decode catcolor --loci O=O/o,W=ws/w,B=B/B                    # calico (tortie + white spotting)
    dna-decode catcolor --loci O=o,C=cs/cs,B=B/B                     # seal-point Siamese
    dna-decode catcolor --loci W=W/w,O=o/o,B=b/b                     # white (dominant white, epistatic)

Deterministic epistatic resolution over the OMIA-curated loci (W/KIT, O/ARHGAP36 X-linked, A/ASIP, B/TYRP1,
D/MLPH, C/TYR). The O locus takes ONE allele for a male (hemizygous, e.g. O=O) or TWO for a female (e.g.
O=O/o); --sex may be given but is inferred from the O zygosity when omitted. Pure-python, wheel-only, offline.
Calls the COLOUR + major pattern (tortie/calico/pointed/bicolor) — NOT tabby sub-pattern, shade, or spotting
extent (those ABSTAIN). KNOWLEDGE_BASELINE. Scope: benign companion-animal visible-trait genetics.
"""
from __future__ import annotations

import argparse
import json
import sys


def _parse_loci(spec: str) -> dict:
    out = {}
    for tok in spec.split(","):
        tok = tok.strip()
        if not tok:
            continue
        if "=" not in tok:
            raise ValueError(f"bad locus token {tok!r}; expected LOCUS=a1/a2 (or O=O for a male)")
        loc, gt = tok.split("=", 1)
        out[loc.strip()] = gt.strip()
    return out


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    ap = argparse.ArgumentParser(
        prog="dna-decode catcolor",
        description="Cat coat-colour decoder (curated-catalog + X-linked epistasis over W/O/A/B/D/C, deterministic).",
        epilog="scope: benign companion-animal visible-trait genetics. Tabby sub-pattern / silver / wideband "
               "ABSTAIN via --present. The O locus is X-linked (1 allele = male, 2 = female). NOT a human tool.",
    )
    ap.add_argument("--loci", required=True,
                    help="comma-separated LOCUS=a1/a2 for W/O/A/B/D/C, e.g. O=O/o,B=B/B,D=D/D,A=a/a. The O "
                         "(X-linked) locus takes ONE allele for a male (O=O) or TWO for a female (O=O/o)")
    ap.add_argument("--sex", choices=["male", "female"], default=None,
                    help="cat sex (optional — inferred from the O-locus zygosity when omitted)")
    ap.add_argument("--present", default="",
                    help="comma-separated v0-UNMODELLED loci present (TA tabby-subpattern, I silver, WB "
                         "wideband, KA karpati/roan) — the affected appearance axis ABSTAINS")
    ap.add_argument("--json", action="store_true", help="emit the result as JSON")
    args = ap.parse_args(argv)

    from dna_decode.pigment.cat_coat import CatInputError, call_cat_coat

    try:
        loci = _parse_loci(args.loci)
        present = [p.strip() for p in args.present.split(",") if p.strip()]
        res = call_cat_coat(loci, sex=args.sex, present_loci=present)
    except (ValueError, CatInputError) as e:
        print(f"error: {e}", file=sys.stderr)
        return 2

    if args.json:
        print(json.dumps(res.as_dict(), indent=2))
        return 0

    d = res.as_dict()
    print("coat-colour decode — cat (W/O/A/B/D/C epistasis; O is X-linked)")
    print(f"  coat colour: {d['coat_color'].upper()}   confidence: {d['confidence']}")
    print(f"  base: {d['base_color']}   sex-basis: {d['sex_basis']}"
          + ("   tortoiseshell: yes" if d['is_tortoiseshell'] else "")
          + (f"   colorpoint: {d['colorpoint']}" if d['colorpoint'] else "")
          + (f"   white: {d['white_pattern']}" if d['white_pattern'] else "")
          + ("   dilute: yes" if d['dilute'] else ""))
    print(f"  per-locus: {d['per_locus']}")
    for ax in d["abstains_on"]:
        print(f"  ABSTAIN: {ax}")
    for n in d["notes"]:
        print(f"  note: {n}")
    print("  [deterministic curated OMIA loci (W/KIT, O/ARHGAP36, A/ASIP, B/TYRP1, D/MLPH, C/TYR); biology-checked]")
    print("  [KNOWLEDGE_BASELINE: curated catalog, no free per-individual validation substrate]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
