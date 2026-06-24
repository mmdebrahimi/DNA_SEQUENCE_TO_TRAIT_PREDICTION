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
        "summary": "antibiotic resistance R/S - bacterial (cipro/cef/tet/gent/meropenem; E.coli/Klebsiella/Pseudomonas/S.aureus) + FUNGAL azole/echinocandin (fluconazole/voriconazole/caspofungin/micafungin; C. auris) via --drug",
        "validation": "bacterial: cipro 0.925 (held-out 0.862, cross-source 1.0) | cef 0.933 | gent 0.945 | tet 0.833 | mero 0.867; cross-organism (capstone). fungal C. auris fluconazole G1: sens 1.0 across clades, label-limited spec (wiki/fungal_ep7_g1_closeout_2026-06-08)",
    },
    "pathotype": {
        "summary": "E. coli pathotype (EPEC/EHEC/ETEC/UPEC/EAEC/...) compatibility call + abstention",
        "validation": "VirulenceFinder-marker resolver; ExPEC recall 0.917; rest documented scope-limit",
    },
    "plasmid": {
        "summary": "plasmid Inc-replicon typing (IncF/IncH/IncI/IncX/IncN/... via PlasmidFinder allele DB) - composes with amr (is the resistance plasmid-borne?)",
        "validation": "deterministic PlasmidFinder-blastn caller (identity 95 / coverage 60); faithful-to-tool, not an independent baseline; offline-safe degrade",
    },
    "serotype": {
        "summary": "E. coli O:H serotype (wzx/wzy/wzm/wzt O-antigen + fliC H-antigen via SerotypeFinder allele DB)",
        "validation": "deterministic SerotypeFinder-blastn caller (identity 85 / coverage 60); faithful-to-tool; O?/H? when a locus is unresolved; offline-safe degrade",
    },
    "resfinder": {
        "summary": "acquired AMR genes (ResFinder allele DB) - an INDEPENDENT cross-tool check vs amr (AMRFinder DB)",
        "validation": "deterministic ResFinder-blastn caller (identity 90 / coverage 60); caller_is_independent_baseline=True (acquired genes only, no point-mutations/efflux); offline-safe degrade",
    },
    "pointfinder": {
        "summary": "chromosomal AMR point mutations (PointFinder; v0 E. coli FQ QRDR gyrA/parC/gyrB/parE) - INDEPENDENT vs amr's AMRFinder POINT",
        "validation": "deterministic blastn + codon-position lookup vs resistens-overview; caller_is_independent_baseline=True; epistasis recorded not enforced; offline-safe degrade",
    },
    "disinfinder": {
        "summary": "biocide/disinfectant resistance genes (DisinFinder; qac/form... quaternary-ammonium + formaldehyde) - often plasmid-borne (pair with coloc)",
        "validation": "deterministic DisinFinder-blastn caller (identity 90 / coverage 60); faithful-to-tool; offline-safe degrade",
    },
    "mlst": {
        "summary": "MLST sequence type (PubMLST; v0 E. coli Achtman adk/fumC/gyrB/icd/mdh/purA/recA) - exact-allele -> profile -> ST",
        "validation": "deterministic blastn 100/100 exact-allele + PubMLST profile lookup; novel/incomplete -> ST not guessed; `dna-mlst --fetch-db` to install; offline-safe degrade",
    },
    "ktype": {
        "summary": "Klebsiella K-antigen (capsule) type via the wzi allele scheme (BIGSdb Pasteur, Kleborate-bundled) - the serotype sibling",
        "validation": "deterministic wzi-blastn caller (identity 90 / coverage 80); faithful-to-tool; wzi->K ~94% NOT one-to-one; full-locus Kaptive more accurate; offline-safe degrade",
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
    if trait == "plasmid":
        from dna_decode.plasmid.cli import main as plasmid_main
        return plasmid_main(rest)
    if trait == "serotype":
        from dna_decode.serotype.cli import main as serotype_main
        return serotype_main(rest)
    if trait == "resfinder":
        from dna_decode.resfinder.cli import main as resfinder_main
        return resfinder_main(rest)
    if trait == "pointfinder":
        from dna_decode.pointfinder.cli import main as pointfinder_main
        return pointfinder_main(rest)
    if trait == "disinfinder":
        from dna_decode.disinfinder.cli import main as disinfinder_main
        return disinfinder_main(rest)
    if trait == "mlst":
        from dna_decode.mlst.cli import main as mlst_main
        return mlst_main(rest)
    if trait == "ktype":
        from dna_decode.ktype.cli import main as ktype_main
        return ktype_main(rest)
    if trait == "concordance":
        from dna_decode.concordance.cli import main as concordance_main
        return concordance_main(rest)
    if trait == "profile":
        from dna_decode.profile.cli import main as profile_main
        return profile_main(rest)
    if trait == "coloc":
        from dna_decode.colocalization.cli import main as coloc_main
        return coloc_main(rest)
    raise ValueError(f"unknown trait: {trait}")


# Cross-decoder ANALYSES (compose the decoders; NOT new traits/DBs — kept out of TRAITS so the
# decoder registry contract stays the 5-decoder set).
ANALYSES = {
    "concordance": "AMR cross-tool concordance (AMRFinder vs ResFinder acquired-gene calls)",
    "profile": "unified genome profile - run all assembly-FASTA decoders in one report",
    "coloc": "resistance-gene x plasmid co-localization (is this acquired AMR gene plasmid-borne?)",
}


def _print_list() -> int:
    print(f"dna-decode {_version()} - deterministic genotype->phenotype decoders\n")
    for name, meta in TRAITS.items():
        print(f"  {name:11} {meta['summary']}")
        print(f"  {'':11} validation: {meta['validation']}")
    print("\nanalyses (compose the decoders):")
    for name, summary in ANALYSES.items():
        print(f"  {name:11} {summary}")
    print("\nrun `dna-decode <trait|analysis> --help` for options.")
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
    for name, summary in ANALYSES.items():
        sub.add_parser(name, add_help=False, help=summary)
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
    if trait not in TRAITS and trait not in ANALYSES:
        ap.error(f"unknown subcommand {trait!r}; traits: {', '.join(TRAITS)}; "
                 f"analyses: {', '.join(ANALYSES)} (or `list`)")
    return _delegate(trait, argv[1:])


if __name__ == "__main__":
    raise SystemExit(main())
