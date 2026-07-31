"""`dna-decode {rabbit,mouse,cattle,pig,sheep}color` — mammalian coat-colour CLIs over the shared engine.

    dna-decode rabbitcolor --loci A=A/a,B=B/b,C=C/C,D=D/d,E=E/e     # rabbit A-E series
    dna-decode mousecolor  --loci A=a/a,B=b/b,C=C/C,P=p/p,E=E/E     # mouse (+ pink-eye)
    dna-decode cattlecolor --loci E=ED/e                           # cattle (MC1R + PMEL dilution)
    dna-decode pigcolor    --loci KIT=I/i+,E=e/e                    # pig (KIT dominant white + MC1R)
    dna-decode sheepcolor  --loci A=AWt/a,E=E+/E+                   # sheep (ASIP + MC1R)

Each is a thin wrapper over dna_decode.pigment.mammal_color (the shared Extension/Agouti/Brown/Dilute/Albino
epistatic engine + per-organism OMIA catalogs). Pure-python, wheel-only, offline, deterministic. Calls colour
+ agouti pattern — NOT fine pattern/spotting/shade. KNOWLEDGE_BASELINE. Benign livestock/lab-animal genetics.
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


def _main(organism: str, argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    from dna_decode.pigment.mammal_color import MAMMAL_CATALOGS, MammalInputError, call_mammal_color

    cat = MAMMAL_CATALOGS[organism]
    loci_help = ", ".join(f"{s} ({cat.loci[s].gene})" for s in cat.loci)
    ap = argparse.ArgumentParser(
        prog=f"dna-decode {organism}color",
        description=f"{organism.capitalize()} coat-colour decoder (curated OMIA epistatic catalog, deterministic).",
        epilog=f"loci: {loci_help}. scope: benign livestock/lab-animal visible-trait genetics; calls colour + "
               "agouti pattern, not fine pattern/spotting/shade. NOT a human tool.",
    )
    ap.add_argument("--loci", required=True,
                    help=f"comma-separated LOCUS=a1/a2 for {organism} loci ({', '.join(cat.loci)})")
    ap.add_argument("--json", action="store_true", help="emit the result as JSON")
    args = ap.parse_args(argv)

    try:
        res = call_mammal_color(cat, _parse_loci(args.loci))
    except (ValueError, MammalInputError) as e:
        print(f"error: {e}", file=sys.stderr)
        return 2

    if args.json:
        print(json.dumps(res.as_dict(), indent=2))
        return 0

    d = res.as_dict()
    print(f"coat-colour decode — {organism} (curated OMIA A/B/C/D/E epistasis)")
    print(f"  coat colour: {d['coat_color'].upper()}   confidence: {d['confidence']}")
    print(f"  pigment: {d['pigment_type']}"
          + (f"   eumelanin: {d['base_eumelanin']}" if d['base_eumelanin'] else "")
          + f"   pattern: {d['pattern']}"
          + (f"   dilutions: {', '.join(d['dilutions'])}" if d['dilutions'] else ""))
    print(f"  per-locus: {d['per_locus']}")
    for n in d["notes"]:
        print(f"  note: {n}")
    print(f"  [deterministic curated OMIA loci ({', '.join(cat.loci[s].gene for s in cat.loci)}); biology-checked]")
    print("  [KNOWLEDGE_BASELINE: curated catalog, no free per-individual validation substrate]")
    return 0


def rabbit_main(argv=None) -> int:
    return _main("rabbit", argv)


def mouse_main(argv=None) -> int:
    return _main("mouse", argv)


def cattle_main(argv=None) -> int:
    return _main("cattle", argv)


def pig_main(argv=None) -> int:
    return _main("pig", argv)


def sheep_main(argv=None) -> int:
    return _main("sheep", argv)


def goat_main(argv=None) -> int:
    return _main("goat", argv)


def alpaca_main(argv=None) -> int:
    return _main("alpaca", argv)


def guineapig_main(argv=None) -> int:
    return _main("guineapig", argv)


def fox_main(argv=None) -> int:
    return _main("fox", argv)


def donkey_main(argv=None) -> int:
    return _main("donkey", argv)


def buffalo_main(argv=None) -> int:
    return _main("buffalo", argv)


def camel_main(argv=None) -> int:
    return _main("camel", argv)


def mink_main(argv=None) -> int:
    return _main("mink", argv)


def roedeer_main(argv=None) -> int:
    return _main("roedeer", argv)


if __name__ == "__main__":
    import sys as _s
    org = _s.argv[1] if len(_s.argv) > 1 else "rabbit"
    raise SystemExit(_main(org, _s.argv[2:]))
