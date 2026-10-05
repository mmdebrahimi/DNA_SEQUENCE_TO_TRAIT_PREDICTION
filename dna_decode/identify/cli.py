"""`dna-identify` -- given a genome, name which SUPPORTED organism it is, or abstain.

    dna-identify --genome-fasta path/to/genome.fna
    dna-identify --genome-fasta X.fna --json
    dna-identify --list-organisms

Step 1 of the decoder pipeline: every other trait REQUIRES `--organism`, so this is the piece that
lets a genome route itself. The printed `routing_token` is the AMRFinder `-O` value, which is ALSO
what `call_resistance(organism=...)` resolves through (verified: the genus-prefix fallback makes one
token serve both), so it can be passed straight to the next decoder.

CLOSED SET. The question is "which of the supported organisms is this, if any?", never "what is
this?". An organism outside the set ABSTAINS. Measured: 44 of 48 out-of-set genomes abstain; the 4
that do not are Klebsiella congeners whose distance is BELOW the in-set ceiling, a named limit no
threshold can fix (see `thresholds.py`).

Thresholds are NOT restated here -- the argparse defaults derive from `thresholds.py`. A CLI that
re-declares a validated constant silently overrides it the moment the constant moves, which has
already happened once in this repo (a `--coverage default=80.0` kept running a threshold replaced
five days earlier, worth -0.225 accuracy).
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from dna_decode.data.cell_evidence_line import evidence_one_line
from dna_decode.data.organism_vocab import BY_CANONICAL, routing_token, supported_canonicals
from dna_decode.identify import core
from dna_decode.identify.runner import DEFAULT_REFERENCE, ReferenceUnavailable, query
from dna_decode.identify.thresholds import AMBIGUITY_MARGIN, MAX_DISTANCE, frozen


def _manifest_labels(manifest: Path | None) -> dict[str, str]:
    """accession -> canonical organism, from the committed reference manifest."""
    repo = Path(__file__).resolve().parent.parent.parent
    p = manifest
    if p is None:
        cands = sorted((repo / "wiki").glob("identify_reference_manifest_*.json"))
        if not cands:
            raise ReferenceUnavailable(
                "no reference manifest in wiki/; build one with "
                "`uv run python scripts/build_identify_reference.py`")
        p = cands[-1]
    m = json.loads(p.read_text(encoding="utf-8"))
    return {a: o for o, accs in m["accessions_by_organism"].items() for a in accs}


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="dna-identify",
        description="Closed-set organism identification from an assembled genome (Mash sketch "
                    "distance). Abstains rather than guess on an organism outside the supported set.",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genome-fasta", type=Path, help="assembled genome (FASTA)")
    ap.add_argument("--reference", type=Path, default=DEFAULT_REFERENCE,
                    help="Mash reference sketch (default: the D: cache)")
    ap.add_argument("--manifest", type=Path, default=None,
                    help="reference manifest (default: newest wiki/identify_reference_manifest_*.json)")
    # DERIVED from thresholds.py -- never restated. See the module docstring.
    ap.add_argument("--max-distance", type=float, default=MAX_DISTANCE,
                    help=f"abstain above this Mash distance (default {MAX_DISTANCE}, derived)")
    ap.add_argument("--ambiguity-margin", type=float, default=AMBIGUITY_MARGIN,
                    help=f"abstain when two organisms are within this distance "
                         f"(default {AMBIGUITY_MARGIN}, derived)")
    ap.add_argument("--top-k", type=int, default=5, help="hits to retain in the record")
    ap.add_argument("--json", action="store_true", help="emit the record as JSON")
    ap.add_argument("--list-organisms", action="store_true",
                    help="print the supported closed set and exit")
    return ap


def main(argv: list[str] | None = None) -> int:
    ap = build_parser()
    a = ap.parse_args(argv)

    if a.list_organisms:
        # COVERAGE, not just vocabulary. The vocabulary declares more organisms than the reference
        # sketch actually contains, and an organism with no reference genomes can only ever ABSTAIN.
        # Printing the vocabulary alone would advertise coverage that does not exist -- the `--help`
        # over-claim pattern (`dna-hla` advertised two SNP tags its own measurement had demoted).
        try:
            covered = set(_manifest_labels(a.manifest).values())
        except (ReferenceUnavailable, OSError, KeyError, json.JSONDecodeError):
            covered = set()
        allc = supported_canonicals()
        print(f"closed set: {len(covered)} of {len(allc)} vocabulary organisms have reference "
              f"coverage and are IDENTIFIABLE; the rest can only ABSTAIN.")
        print("\nIDENTIFIABLE (in the reference sketch):")
        for c in allc:
            if c not in covered:
                continue
            tok = BY_CANONICAL[c].amrfinder_organism
            print(f"   {c:<28} -O {tok}" if tok
                  else f"   {c:<28} (no AMRFinder -O value; routes to its own engine)")
        missing = [c for c in allc if c not in covered]
        if missing:
            print("\nDECLARED but NOT in the reference -- a genome of these ABSTAINS, it is NOT "
                  "identified:")
            for c in missing:
                print(f"   {c}")
            print("   (add reference genomes for them and rebuild with "
                  "scripts/build_identify_reference.py)")
        return 0

    if not a.genome_fasta:
        ap.error("--genome-fasta is required (or use --list-organisms)")

    th = (frozen() if (a.max_distance == MAX_DISTANCE and a.ambiguity_margin == AMBIGUITY_MARGIN)
          else core.Thresholds(max_distance=a.max_distance, ambiguity_margin=a.ambiguity_margin))

    try:
        labels = _manifest_labels(a.manifest)
        res = query(a.genome_fasta, a.reference)
    except ReferenceUnavailable as exc:
        # NOT an abstention record: an infrastructure fault must not read as "organism unsupported".
        print(f"REFERENCE UNAVAILABLE: {exc}", file=sys.stderr)
        return 3

    call = core.decide(res.hits, th, lambda acc: labels.get(acc), routing_token, top_k=a.top_k)
    rec = call.as_dict()
    rec["query"] = str(a.genome_fasta)
    rec["reference"] = str(a.reference)

    if a.json:
        print(json.dumps(rec, indent=2))
        return 0

    if call.abstained:
        print(f"ABSTAIN ({call.reason.value})")
        if call.nearest_distance is not None:
            print(f"  nearest reference distance : {call.nearest_distance:.5f}")
        if call.runner_up_organism:
            print(f"  runner-up                  : {call.runner_up_organism} "
                  f"({call.runner_up_distance:.5f})")
        print("  the supported set does not contain this organism, or it sits between two of them.")
        print("  NOT identified -- this is a closed-set router, not open-world taxonomy.")
        return 0

    print(f"ORGANISM : {call.organism}")
    print(f"-O value : {call.routing_token}" if call.routing_token
          else "-O value : (none -- routes to this organism's own engine, not AMRFinder)")
    print(f"distance : {call.nearest_distance:.5f}")
    if call.runner_up_organism:
        print(f"runner-up: {call.runner_up_organism} ({call.runner_up_distance:.5f})")
    print("scope    : closed-set, species-level, assembled genomes only; a Klebsiella congener "
          "(variicola/michiganensis) is a KNOWN mis-call class")
    # The cell's evidence tier, from its committed registry contract. BESIDE the scope line, not
    # replacing it: a caveat states the method's limits, this states what was measured. Evidence
    # carried only in the registry is not a disclosure. Wired 2026-10-04 -- this route shipped
    # earlier the same day WITHOUT it, and the exhaustive-wiring guard is what surfaced that.
    _ev = evidence_one_line("dna-identify")
    if _ev:
        print(f"  {_ev}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
