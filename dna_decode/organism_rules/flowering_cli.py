"""`dna-decode flowering` — Arabidopsis flowering-habit decoder (FRI/FLC curated causal loci).

    dna-decode flowering --fri Col   --flc Col         # -> summer_annual_early (the Col-0 rapid cycler)
    dna-decode flowering --fri Sf-2  --flc Col         # -> winter_annual_late  (the Col-FRI line)
    dna-decode flowering --fri Sf-2  --flc Da\\(1\\)-12  # -> early via the FLC route (functional FRI!)
    dna-decode flowering --list-alleles

v0 takes ALLELE CALLS (the wheel-only pattern the HIV/fungal cells use); genome-mode (detect the FRI-Col stop /
FRI-Ler start-codon deletion from sequence) needs verified variant coordinates — a v0.1 follow-on. Deterministic,
pure-python, offline. Scope: flowering HABIT/direction (~40-70% of long-day variation), NOT days-to-flower.
"""
from __future__ import annotations

import argparse
import json
import sys


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    ap = argparse.ArgumentParser(
        prog="dna-decode flowering",
        description="Arabidopsis flowering-habit decoder — deterministic FRI/FLC curated-locus rule "
                    "(late iff functional FRI AND strong FLC).",
        epilog="scope: HABIT/direction (~40-70% of long-day variation), NOT quantitative days-to-flower.",
    )
    ap.add_argument("--fri", help="FRI allele name (Col/Ler/Cvi/Wil-2/Sf-2/H51) or a status "
                                  "(functional/lof/unknown)")
    ap.add_argument("--flc", help="FLC allele name (Col/Ler/Van-0/Bur-0/Da(1)-12/Shakhdara) or a status "
                                  "(functional/weak/lof/unknown)")
    ap.add_argument("--list-alleles", action="store_true", help="print the curated allele catalog + exit")
    ap.add_argument("--json", action="store_true", help="emit the call as JSON")
    args = ap.parse_args(argv)

    from dna_decode.organism_rules.arabidopsis_flowering import (
        FLC_ALLELES, FRI_ALLELES, FloweringInputError, call_flowering_habit,
    )

    if args.list_alleles:
        print("curated causal-allele catalog (Arabidopsis flowering habit)\n")
        for locus, table in (("FRI", FRI_ALLELES), ("FLC", FLC_ALLELES)):
            print(f"  {locus}:")
            for name, meta in table.items():
                lesion = f" — {meta['lesion']}" if meta["lesion"] else ""
                print(f"    {name:12} {meta['status']:11}{lesion}")
                print(f"                 source: {meta['source']}")
            print()
        return 0

    if not args.fri or not args.flc:
        print("error: both --fri and --flc are required (or use --list-alleles)", file=sys.stderr)
        return 2

    try:
        call = call_flowering_habit(args.fri, args.flc)
    except FloweringInputError as e:
        print(f"error: {e}", file=sys.stderr)
        return 2

    if args.json:
        print(json.dumps(call.as_dict(), indent=2))
        return 0

    d = call.as_dict()
    print("Arabidopsis flowering-habit decode (FRI/FLC curated loci)")
    print(f"  habit: {d['habit'].upper()}   confidence: {d['confidence']}")
    print(f"  vernalization required: {d['vernalization_required']}")
    print(f"  FRI: {d['fri_status']}   FLC: {d['flc_status']}   rule: {d['rule']}")
    for n in d["notes"]:
        print(f"  note: {n}")
    print(f"  [{d['scope_limit']}]")
    print("  [cannot see: " + d["undetectable_mechanisms"][0] + "]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
