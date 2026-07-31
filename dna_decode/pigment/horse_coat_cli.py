"""`dna-decode horsecolor` — horse coat-colour decoder (curated-catalog + epistasis, v0).

    dna-decode horsecolor --loci E=E/e,A=A/a                          # bay
    dna-decode horsecolor --loci E=e/e,A=a/a,CR=Cr/N                  # palomino (chestnut + 1 cream)
    dna-decode horsecolor --loci E=E/E,A=a/a,D=D/nd1                  # grullo (black dun)
    dna-decode horsecolor --loci E=E/e,A=A/a,G=G/n                    # grey (born bay, greys with age)
    dna-decode horsecolor --loci E=E/e,A=A/a --present Z,TO --json    # silver+tobiano -> abstain on those axes

Deterministic epistatic resolution over the five OMIA-curated loci (E/MC1R, A/ASIP, CR/SLC45A2, D/TBX3,
G/STX17). Pure-python, wheel-only, offline. Calls the COLOUR (base + cream/dun dilution + grey) — NOT
sooty/flaxen shade, spotting extent, or champagne/silver/pearl (those ABSTAIN). KNOWLEDGE_BASELINE (curated
catalog; no free per-individual validation substrate). Scope: benign livestock/companion visible-trait genetics.
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
            raise ValueError(f"bad locus token {tok!r}; expected LOCUS=a1/a2 (e.g. E=E/e)")
        loc, gt = tok.split("=", 1)
        out[loc.strip()] = gt.strip()
    return out


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    ap = argparse.ArgumentParser(
        prog="dna-decode horsecolor",
        description="Horse coat-colour decoder (curated-catalog epistasis over E/A/CR/D/G loci, deterministic).",
        epilog="scope: benign livestock/companion visible-trait genetics. Dilution/pattern loci "
               "(champagne/silver/pearl/roan/tobiano/appaloosa) ABSTAIN via --present. NOT a human tool.",
    )
    ap.add_argument("--loci", required=True,
                    help="comma-separated LOCUS=a1/a2 for E/A/CR/D/G, e.g. E=E/e,A=A/a,CR=Cr/N,D=nd1/nd1,G=n/n "
                         "(E is required — the pigment-type switch; CR uses Cr/N, D uses D/nd1/nd2, G uses G/n)")
    ap.add_argument("--present", default="",
                    help="comma-separated v0-UNMODELLED loci present (CH champagne, Z silver, PRL pearl, RN "
                         "roan, TO tobiano/overo, LP appaloosa, STY sooty/flaxen) — the affected axis ABSTAINS")
    ap.add_argument("--json", action="store_true", help="emit the result as JSON")
    args = ap.parse_args(argv)

    from dna_decode.pigment.horse_coat import HorseInputError, call_horse_coat

    try:
        loci = _parse_loci(args.loci)
        present = [p.strip() for p in args.present.split(",") if p.strip()]
        res = call_horse_coat(loci, present_loci=present)
    except (ValueError, HorseInputError) as e:
        print(f"error: {e}", file=sys.stderr)
        return 2

    if args.json:
        print(json.dumps(res.as_dict(), indent=2))
        return 0

    d = res.as_dict()
    print("coat-colour decode — horse (E/A/CR/D/G epistasis)")
    print(f"  coat colour: {d['coat_color'].upper()}   confidence: {d['confidence']}")
    print(f"  base: {d['base_color']}   pigment: {d['pigment_type']}"
          + (f"   dilutions: {', '.join(d['dilutions'])}" if d['dilutions'] else "")
          + ("   greying-with-age: yes" if d['greying_with_age'] else ""))
    print(f"  per-locus: {d['per_locus']}")
    for ax in d["abstains_on"]:
        print(f"  ABSTAIN: {ax}")
    for n in d["notes"]:
        print(f"  note: {n}")
    print("  [deterministic curated OMIA loci (E/MC1R, A/ASIP, CR/SLC45A2, D/TBX3, G/STX17); biology-checked]")
    print("  [KNOWLEDGE_BASELINE: curated catalog, no free per-individual validation substrate]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
