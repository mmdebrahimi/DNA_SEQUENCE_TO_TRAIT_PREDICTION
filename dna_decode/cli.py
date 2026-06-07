"""Unified `dna-decode` tool entry — one command, the trait decoders underneath.

Turns the two validated deterministic decoders into a single coherent tool (per the project north star:
"AI DNA decoder tool, not papers"). Thin dispatcher — each subcommand delegates argv verbatim to the
existing decoder `main(argv)`, so the per-decoder CLIs (`dna-amr`, `dna-pathotype`) stay independently
usable and their behavior/output is unchanged.

    dna-decode amr --drug ciprofloxacin --amrfinder-run data/amrfinder_runs/GCA_xxx.x
    dna-decode amr --drug ceftriaxone  --genome-fasta X.fna --sample-id X   # needs Docker + AMRFinder DB
    dna-decode pathotype assembly.fna --sample-id MY_STRAIN
    dna-decode list        # what this tool decodes + per-trait validation status
    dna-decode --version

Deterministic mechanism-feature decoders — NOT embeddings (the frozen-genome-embedding thesis was tested
to a decisive FAIL on every reachable de-confounded substrate; see plans/AMR_embedding_niche_decision_2026-06-05.md
+ CHANGELOG 0.3.0). NOT a clinical decision tool.
"""
from __future__ import annotations

import argparse
import sys

# Per-trait capability registry: subcommand -> (delegate dotted-path main, one-line capability + validation).
TRAITS = {
    "amr": {
        "summary": "antibiotic resistance R/S (cipro/cef/tet/gent/meropenem) - E.coli/Klebsiella/Pseudomonas/S.aureus (--organism)",
        "validation": "cipro 0.925 (held-out 0.862, cross-source 1.0) | cef 0.933 | gent 0.945 | tet 0.833 | mero 0.867; cross-organism (see capstone)",
    },
    "pathotype": {
        "summary": "E. coli pathotype (EPEC/EHEC/ETEC/UPEC/EAEC/...) compatibility call + abstention",
        "validation": "VirulenceFinder-marker resolver; ExPEC recall 0.917; rest documented scope-limit",
    },
}


def _version() -> str:
    try:
        from importlib.metadata import version
        return version("dna_decode")
    except Exception:
        return "unknown"


def _delegate(trait: str, rest: list[str]) -> int:
    if trait == "amr":
        from dna_decode.amr.cli import main as amr_main
        return amr_main(rest)
    if trait == "pathotype":
        from dna_decode.pathotype.cli import main as patho_main
        return patho_main(rest)
    raise ValueError(f"unknown trait: {trait}")


def _print_list() -> int:
    print(f"dna-decode {_version()} - deterministic genotype->phenotype decoders\n")
    for name, meta in TRAITS.items():
        print(f"  {name:10} {meta['summary']}")
        print(f"  {'':10} validation: {meta['validation']}")
    print("\nrun `dna-decode <trait> --help` for a decoder's options.")
    return 0


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    ap = argparse.ArgumentParser(
        prog="dna-decode",
        description="Unified DNA trait decoder — deterministic, interpretable, mechanism-feature based.",
        epilog="traits: " + ", ".join(TRAITS) + ".  `dna-decode list` for validation status.",
    )
    ap.add_argument("--version", action="version", version=f"dna-decode {_version()}")
    sub = ap.add_subparsers(dest="trait", metavar="{amr,pathotype,list}")
    # Register thin pass-through subparsers; real arg parsing happens in each decoder's main().
    for name, meta in TRAITS.items():
        sub.add_parser(name, add_help=False, help=meta["summary"])
    sub.add_parser("list", help="show what this tool decodes + per-trait validation status")

    # Split argv at the subcommand so the rest passes through verbatim (incl. --help) to the decoder.
    if not argv:
        ap.print_help()
        return 0
    trait = argv[0]
    if trait in ("-h", "--help"):
        ap.print_help()
        return 0
    if trait == "--version":
        print(f"dna-decode {_version()}")
        return 0
    if trait == "list":
        return _print_list()
    if trait not in TRAITS:
        ap.error(f"unknown trait {trait!r}; choose from {', '.join(TRAITS)} (or `list`)")
    return _delegate(trait, argv[1:])


if __name__ == "__main__":
    raise SystemExit(main())
