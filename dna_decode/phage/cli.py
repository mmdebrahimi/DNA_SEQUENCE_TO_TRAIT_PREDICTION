"""`dna-decode phage` — bacteriophage genome -> host-receptor-class decoder (the first non-AMR cell).

    # wheel-only, offline (catalogue lookup by NCBI lineage):
    dna-decode phage --lineage Vequintavirus
    dna-phage --lineage Tequintavirus,Demerecviridae --json

    # genome-homology transfer (needs a reference set + blastn; BASEL genomes fetched on demand):
    dna-decode phage --genome-fasta my_phage.fna

Two paths, one honest scope (RECEPTOR-CLASS only, NOT the full host-range matrix):
  - `--lineage`  : pure catalogue lookup (the receptor the paper assigns to that genus/family). Wheel-only,
    no BLAST, no downloads. RBP-variable clades (T-even, Drexlerviridae) return INDETERMINATE by design
    (receptor varies within the clade) rather than a fabricated single receptor.
  - `--genome-fasta` : genome-homology receptor TRANSFER — the query inherits the receptor of its nearest
    BLAST neighbour among a labelled reference set (default = the BASEL collection). Needs blastn + the
    reference genomes; degrades to an actionable message (never a fabricated call) when either is absent.

Tier: KNOWLEDGE_BASELINE / in-distribution (catalogue from Maffei 2021 PLOS Biology 3001424, scored on
BASEL). NOT a clinical/biocontrol decision tool.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

_SCOPE = ("scope: receptor-CLASS only (NOT the full phage x strain host-range matrix, which is polygenic); "
          "tier KNOWLEDGE_BASELINE / in-distribution (BASEL, Maffei 2021 PLOS Biology 3001424).")
_DEFAULT_MANIFEST = "data/phage_ref/basel_manifest.tsv"
_DEFAULT_GENOME_DIR = "data/phage_ref/basel"


def _lineage_result(lineage: list[str]) -> dict:
    from dna_decode.data.phage_receptor import label_receptor_for_lineage, receptor_for_lineage
    tr = receptor_for_lineage(lineage)
    label = label_receptor_for_lineage(lineage)
    if tr is None:
        return {"mode": "lineage", "lineage": lineage, "status": "INDETERMINATE",
                "receptor": None, "reason": "no catalogued receptor for this lineage (uncatalogued taxon)"}
    if label is None:  # catalogued but RBP-variable -> abstain, don't fabricate a single receptor
        return {"mode": "lineage", "lineage": lineage, "status": "INDETERMINATE", "receptor": None,
                "taxon": tr.taxon, "rank": tr.rank, "clade_conserved": False,
                "receptors_in_clade": list(tr.receptors),
                "reason": f"receptor varies within {tr.taxon} by receptor-binding protein ({', '.join(tr.receptors)}) "
                          "-- not clade-conserved, so no single receptor is assigned"}
    return {"mode": "lineage", "lineage": lineage, "status": "CALLED", "receptor": label,
            "taxon": tr.taxon, "rank": tr.rank, "clade_conserved": True, "note": tr.note}


def _genome_result(genome_fasta: str, manifest: str, genome_dir: str) -> dict:
    from dna_decode.phage.receptor_caller import _load_manifest, call_receptor
    if not Path(genome_fasta).exists():
        return {"mode": "genome", "status": "ERROR", "receptor": None,
                "reason": f"cannot read --genome-fasta: {genome_fasta}"}
    if not Path(manifest).exists() or not Path(genome_dir).exists():
        return {"mode": "genome", "status": "INDETERMINATE", "receptor": None,
                "reason": f"reference set not found ({manifest} / {genome_dir}); the BASEL genomes are "
                          "regenerable via `uv run python scripts/fetch_basel_genomes.py` (gitignored, ~7MB). "
                          "Or use `--lineage <genus>` for the wheel-only catalogue path."}
    refs, receptors = _load_manifest(manifest, genome_dir)
    if not refs:
        return {"mode": "genome", "status": "INDETERMINATE", "receptor": None,
                "reason": "reference manifest yielded no labelled phages"}
    call = call_receptor(genome_fasta, refs, receptors)
    return {"mode": "genome", "status": call.status, "receptor": call.predicted_receptor,
            "nearest_reference": call.nearest_label, "percent_identity": call.percent_identity,
            "method": call.method, "reason": call.reason,
            "n_reference_phages": len(refs), "reason_if_indeterminate": call.reason or None}


def _read_first_protein(fasta: str) -> str | None:
    seq = []
    started = False
    with open(fasta, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if started:
                    break
                started = True
                continue
            if started:
                seq.append(line.strip())
    return "".join(seq) or None


def _rbp_result(rbp_fasta: str) -> dict:
    from dna_decode.phage.rbp_caller import call_rbp_from_protein
    if not Path(rbp_fasta).exists():
        return {"mode": "rbp", "status": "ERROR", "receptor": None,
                "reason": f"cannot read --rbp-fasta: {rbp_fasta}"}
    prot = _read_first_protein(rbp_fasta)
    if not prot:
        return {"mode": "rbp", "status": "ERROR", "receptor": None,
                "reason": "no protein sequence found in --rbp-fasta"}
    call = call_rbp_from_protein(prot)
    return {"mode": "rbp", "status": call.status, "receptor": call.predicted_receptor,
            "nearest_reference": call.nearest_phage, "similarity": call.similarity,
            "method": "rbp_kmer_transfer_v0", "reason": call.reason}


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    ap = argparse.ArgumentParser(
        prog="dna-decode phage",
        description="Bacteriophage genome/lineage -> host-receptor-class decoder (first non-AMR cell).",
        epilog=_SCOPE)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--lineage",
                     help="comma-separated NCBI taxa, genus-first (e.g. Vequintavirus or "
                          "Tequintavirus,Demerecviridae) -> catalogue receptor lookup (wheel-only, offline)")
    src.add_argument("--genome-fasta",
                     help="a phage genome FASTA -> genome-homology receptor transfer vs a reference set (needs blastn)")
    src.add_argument("--rbp-fasta",
                     help="a phage tail-fiber RECEPTOR-BINDING-PROTEIN FASTA -> RBP k-mer transfer (wheel-only, "
                          "offline; covers the RBP-variable mixed clades Tsx/OmpC/FhuA/OmpA that --lineage abstains on)")
    ap.add_argument("--reference-manifest", default=_DEFAULT_MANIFEST,
                    help=f"labelled reference manifest TSV (default {_DEFAULT_MANIFEST})")
    ap.add_argument("--reference-dir", default=_DEFAULT_GENOME_DIR,
                    help=f"dir of <accession>.fna reference genomes (default {_DEFAULT_GENOME_DIR})")
    ap.add_argument("--json", action="store_true", help="emit the result as JSON")
    args = ap.parse_args(argv)

    if args.lineage:
        res = _lineage_result([t.strip() for t in args.lineage.split(",") if t.strip()])
    elif args.rbp_fasta:
        res = _rbp_result(args.rbp_fasta)
    else:
        res = _genome_result(args.genome_fasta, args.reference_manifest, args.reference_dir)
    res["scope"] = _SCOPE

    if args.json:
        print(json.dumps(res, indent=2))
    else:
        st = res["status"]
        if st == "CALLED":
            print(f"receptor: {res['receptor']}  ({res['mode']} path)")
            if res.get("nearest_reference"):
                score = res.get("percent_identity") if res.get("percent_identity") is not None else res.get("similarity")
                label = "pident" if res.get("percent_identity") is not None else "rbp k-mer similarity"
                print(f"  nearest reference phage: {res['nearest_reference']} ({label} {score})")
            if res.get("taxon"):
                print(f"  taxon: {res['taxon']} ({res.get('rank')})")
        else:
            print(f"{st}: {res.get('reason') or 'no call'}")
        print(_SCOPE)
    return 0 if res["status"] in ("CALLED", "INDETERMINATE") else 2


if __name__ == "__main__":
    raise SystemExit(main())
