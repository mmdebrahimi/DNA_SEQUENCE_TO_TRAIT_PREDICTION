"""`dna-decode pigeoncolor` — pigeon plumage-colour decoder (curated-catalog + Z-linked epistasis, v0).

    dna-decode pigeoncolor --loci B=B+/B+,C=+/+                     # blue bar (wild rock pigeon)
    dna-decode pigeoncolor --loci B=BA/B+,C=C/+                     # ash-red checker
    dna-decode pigeoncolor --loci B=B+/B+,E=e/e                     # recessive red (SOX10, epistatic over B)
    dna-decode pigeoncolor --loci B=B+/B+,D=d/d                     # dun (dilute blue)
    dna-decode pigeoncolor --loci B=BA --sex female                # ash-red hen (Z-linked, ZW hemizygous)

Deterministic epistatic resolution over the Shapiro-lab molecularly-confirmed loci (B/TYRP1 Z-linked,
E/SOX10 recessive-red, D/SLC45A2 Z-linked dilute, C/NDP wing pattern). The Z-linked B and D loci take ONE
allele for a FEMALE (ZW hemizygous) or TWO for a MALE (ZZ) — REVERSED from mammals. Pure-python, wheel-only,
offline. Calls base colour + dilute + wing pattern — NOT modifiers/shade (those ABSTAIN). KNOWLEDGE_BASELINE.
Scope: benign hobby/livestock visible-trait genetics.
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
            raise ValueError(f"bad locus token {tok!r}; expected LOCUS=a1/a2 (or B=B+ for a female)")
        loc, gt = tok.split("=", 1)
        out[loc.strip()] = gt.strip()
    return out


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    ap = argparse.ArgumentParser(
        prog="dna-decode pigeoncolor",
        description="Pigeon plumage-colour decoder (curated-catalog + Z-linked epistasis over B/E/D/C, deterministic).",
        epilog="scope: benign hobby/livestock visible-trait genetics. Modifiers (spread/grizzle/almond) ABSTAIN "
               "via --present. The B/D loci are Z-linked (1 allele = female ZW, 2 = male ZZ). NOT a human tool.",
    )
    ap.add_argument("--loci", required=True,
                    help="comma-separated LOCUS=a1/a2 for B/E/D/C, e.g. B=BA/B+,E=E+/e,D=D/d,C=C/+. The Z-linked "
                         "B and D loci take ONE allele for a FEMALE (ZW, B=B+) or TWO for a MALE (ZZ)")
    ap.add_argument("--sex", choices=["male", "female"], default=None,
                    help="pigeon sex (optional — inferred from the Z-linked B/D zygosity; birds are ZW so a "
                         "FEMALE is hemizygous, REVERSED from mammals)")
    ap.add_argument("--present", default="",
                    help="comma-separated v0-UNMODELLED loci present (S spread, G grizzle, AL almond, IN indigo) "
                         "— the affected appearance axis ABSTAINS")
    ap.add_argument("--json", action="store_true", help="emit the result as JSON")
    args = ap.parse_args(argv)

    from dna_decode.pigment.pigeon_plumage import PigeonInputError, call_pigeon_plumage

    try:
        loci = _parse_loci(args.loci)
        present = [p.strip() for p in args.present.split(",") if p.strip()]
        res = call_pigeon_plumage(loci, sex=args.sex, present_loci=present)
    except (ValueError, PigeonInputError) as e:
        print(f"error: {e}", file=sys.stderr)
        return 2

    if args.json:
        print(json.dumps(res.as_dict(), indent=2))
        return 0

    d = res.as_dict()
    print("plumage decode — pigeon (B/E/D/C epistasis; B/D are Z-linked)")
    print(f"  plumage: {d['plumage'].upper()}   confidence: {d['confidence']}")
    print(f"  base: {d['base_color']}   sex-basis: {d['sex_basis']}"
          + ("   dilute: yes" if d['dilute'] else "")
          + (f"   wing pattern: {d['wing_pattern']}" if d['wing_pattern'] != "n/a" else "")
          + ("   recessive-red: yes" if d['is_recessive_red'] else ""))
    print(f"  per-locus: {d['per_locus']}")
    for ax in d["abstains_on"]:
        print(f"  ABSTAIN: {ax}")
    for n in d["notes"]:
        print(f"  note: {n}")
    print("  [deterministic curated loci (B/TYRP1, E/SOX10, D/SLC45A2, C/NDP; Domyan 2014 / Vickrey 2018); biology-checked]")
    print("  [KNOWLEDGE_BASELINE: curated catalog, no free per-individual validation substrate]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
