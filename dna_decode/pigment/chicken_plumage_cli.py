"""`dna-decode plumage` — chicken plumage-colour decoder (curated-catalog + Z-linked epistasis, v0).

    dna-decode plumage --loci E=E/E                                  # extended black
    dna-decode plumage --loci E=E/E,B=B --sex female                 # barred hen (Z-linked, ZW hemizygous)
    dna-decode plumage --loci E=E/E,B=B/B --sex male                 # Barred Plymouth Rock cock
    dna-decode plumage --loci E=EWh/EWh,B=B/b+                       # wheaten -> barring barely shows (canvas)
    dna-decode plumage --loci E=E/E,I=I/i+                           # dominant white (Leghorn)
    dna-decode plumage --loci E=E/E,BL=Bl/bl+                        # blue (Andalusian)

Deterministic epistatic resolution over the OMIA-curated loci (E/MC1R, B/CDKN2A Z-linked, S/SLC45A2 Z-linked,
I/PMEL17, Bl blue, lav/MLPH, C/TYR). The Z-linked B and S loci take ONE allele for a FEMALE (ZW hemizygous,
e.g. B=B) or TWO for a MALE (ZZ, e.g. B=B/b+) — REVERSED from mammals; --sex may be given but is inferred from
the Z-locus zygosity. Pure-python, wheel-only, offline. Calls the eumelanin CANVAS + major modifiers
(barred/silver/blue/lavender/white) — NOT fine feather pattern/lacing/pencilling (those ABSTAIN).
KNOWLEDGE_BASELINE. Scope: benign livestock visible-trait genetics.
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
            raise ValueError(f"bad locus token {tok!r}; expected LOCUS=a1/a2 (or B=B for a female)")
        loc, gt = tok.split("=", 1)
        out[loc.strip()] = gt.strip()
    return out


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    ap = argparse.ArgumentParser(
        prog="dna-decode plumage",
        description="Chicken plumage-colour decoder (curated-catalog + Z-linked epistasis over E/B/S/I/Bl/lav/c).",
        epilog="scope: benign livestock visible-trait genetics. Fine pattern (Columbian/mottling/pencilling) "
               "ABSTAIN via --present. The B/S loci are Z-linked (1 allele = female ZW, 2 = male ZZ). NOT a human tool.",
    )
    ap.add_argument("--loci", required=True,
                    help="comma-separated LOCUS=a1/a2 for E/B/S/I/BL/LAV/C, e.g. E=E/E,B=B/b+,S=S/s+. The "
                         "Z-linked B and S loci take ONE allele for a FEMALE (ZW, B=B) or TWO for a MALE (ZZ)")
    ap.add_argument("--sex", choices=["male", "female"], default=None,
                    help="chicken sex (optional — inferred from the Z-linked B/S zygosity; birds are ZW so a "
                         "FEMALE is hemizygous, REVERSED from mammals)")
    ap.add_argument("--present", default="",
                    help="comma-separated v0-UNMODELLED loci present (CO Columbian, MO mottling, PG pattern/"
                         "pencilling, DB dark-brown, SP spangling) — the affected appearance axis ABSTAINS")
    ap.add_argument("--json", action="store_true", help="emit the result as JSON")
    args = ap.parse_args(argv)

    from dna_decode.pigment.chicken_plumage import ChickenInputError, call_chicken_plumage

    try:
        loci = _parse_loci(args.loci)
        present = [p.strip() for p in args.present.split(",") if p.strip()]
        res = call_chicken_plumage(loci, sex=args.sex, present_loci=present)
    except (ValueError, ChickenInputError) as e:
        print(f"error: {e}", file=sys.stderr)
        return 2

    if args.json:
        print(json.dumps(res.as_dict(), indent=2))
        return 0

    d = res.as_dict()
    print("plumage decode — chicken (E/B/S/I/Bl/lav/c epistasis; B/S are Z-linked)")
    print(f"  plumage: {d['plumage'].upper()}   confidence: {d['confidence']}")
    print(f"  eumelanin canvas: {d['eumelanin_canvas']}   sex-basis: {d['sex_basis']}"
          + ("   barred: yes" if d['barred'] else "")
          + ("   silver: yes" if d['silver'] else "")
          + (f"   dilution: {d['dilution']}" if d['dilution'] else "")
          + (f"   white: {d['white_type']}" if d['white_type'] else ""))
    print(f"  per-locus: {d['per_locus']}")
    for ax in d["abstains_on"]:
        print(f"  ABSTAIN: {ax}")
    for n in d["notes"]:
        print(f"  note: {n}")
    print("  [deterministic curated OMIA loci (E/MC1R, B/CDKN2A, S/SLC45A2, I/PMEL17, Bl, lav/MLPH, C/TYR); biology-checked]")
    print("  [KNOWLEDGE_BASELINE: curated catalog, no free per-individual validation substrate]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
