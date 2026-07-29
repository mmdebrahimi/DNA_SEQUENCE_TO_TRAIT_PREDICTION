"""`dna-decode coatcolor` — dog coat-colour decoder (curated-catalog + epistasis, v0).

    dna-decode coatcolor --loci E=e/e,B=B/B,D=D/D                       # yellow Labrador
    dna-decode coatcolor --loci E=E/E,K=KB/KB,B=b/b,D=d/d               # Weimaraner (isabella)
    dna-decode coatcolor --loci E=E/E,K=ky/ky,A=at/at,B=B/B,D=D/D       # black-and-tan
    dna-coatcolor --loci E=E/e,B=B/b,D=D/D --present M,S --json         # merle+spotting -> abstains on those axes

Deterministic epistatic resolution over the five classic solid-colour loci (E/K/A/B/D). Pure-python,
wheel-only, offline. Calls the COLOUR (pigment type + eumelanin colour + distribution) — NOT shade,
coat length, or spotting extent (those ABSTAIN). Scope: benign companion-animal visible-trait genetics.
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
            raise ValueError(f"bad locus token {tok!r}; expected LOCUS=a1/a2 (e.g. E=e/e)")
        loc, gt = tok.split("=", 1)
        out[loc.strip()] = gt.strip()
    return out


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    ap = argparse.ArgumentParser(
        prog="dna-decode coatcolor",
        description="Dog coat-colour decoder (curated-catalog epistasis over E/K/A/B/D loci, deterministic).",
        epilog="scope: benign companion-animal visible-trait genetics. Pattern loci (merle/spotting/ticking) "
               "ABSTAIN via --present. NOT a human/forensic tool.",
    )
    ap.add_argument("--loci", required=True,
                    help="comma-separated LOCUS=a1/a2 for E/K/A/B/D, e.g. E=e/e,K=ky/ky,A=at/at,B=B/b,D=D/d "
                         "(E is required — the pigment-type switch)")
    ap.add_argument("--present", default="",
                    help="comma-separated PATTERN loci present but not modelled in v0 (M,S,T,H) — the affected "
                         "appearance axis ABSTAINS instead of a wrong solid-colour call")
    ap.add_argument("--json", action="store_true", help="emit the result as JSON")
    args = ap.parse_args(argv)

    from dna_decode.pigment.dog_coat import CoatInputError, call_coat_color

    try:
        loci = _parse_loci(args.loci)
        present = [p.strip() for p in args.present.split(",") if p.strip()]
        res = call_coat_color(loci, present_loci=present)
    except (ValueError, CoatInputError) as e:
        print(f"error: {e}", file=sys.stderr)
        return 2

    if args.json:
        print(json.dumps(res.as_dict(), indent=2))
        return 0

    d = res.as_dict()
    print("coat-colour decode — dog (E/K/A/B/D epistasis)")
    print(f"  coat colour: {d['coat_color'].upper()}   confidence: {d['confidence']}")
    print(f"  pigment type: {d['pigment_type']}"
          + (f"   eumelanin: {d['eumelanin_color']}" if d['eumelanin_color'] else "")
          + f"   distribution: {d['distribution']}")
    print(f"  per-locus: {d['per_locus']}")
    for ax in d["abstains_on"]:
        print(f"  ABSTAIN: {ax}")
    for n in d["notes"]:
        print(f"  note: {n}")
    print("  [deterministic curated OMIA loci (E/MC1R, K/CBD103, A/ASIP, B/TYRP1, D/MLPH); biology-checked]")
    print("  [scope: coat COLOUR only — not shade/length/spotting; benign companion-animal genetics]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
