"""`dna-decode kleb` / `dna-kleb` — Klebsiella phage depolymerase -> capsule KL-type (fetch-only cell).

    dna-decode kleb --depolymerase-fasta my_dpo.faa            # ranked top-K KL-types
    dna-kleb --depolymerase-fasta my_dpo.faa --top-k 5 --json

Predicts which Klebsiella capsule type(s) a phage depolymerase (enzymatic tail-spike domain) targets, as a
RANKED shortlist (depolymerases are promiscuous). FETCH-ONLY: the package bundles NO data — build the local
reference first with `scripts/fetch_dpotropisearch.py --accept-license` (Zenodo CC-BY / repo non-commercial —
verify your use). Without a reference, the cell degrades to an actionable INDETERMINATE + the fetch command.

Tier: KNOWLEDGE_BASELINE / in-distribution (prophage-LCA labels; clonality-corrected LOO top-1 ~0.45 /
top-5 ~0.60, +0.49 over null). NOT a clinical tool.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

_SCOPE = ("scope: capsule KL-type only; in-distribution (DpoTropiSearch prophage-LCA labels); "
          "clonality-corrected LOO top-1 ~0.45 / top-5 ~0.60. NOT clinical.")
_FETCH_HINT = ("no local reference — build it (fetch-only, no bundled data): "
               "`uv run --with biopython python scripts/fetch_dpotropisearch.py --accept-license` "
               "(Zenodo CC-BY / repo Decapsulate Non-Commercial License — verify your use), "
               "or set $DPO_KLEB_REFERENCE.")


def _read_first_protein(fasta: str) -> str | None:
    seq, started = [], False
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


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    ap = argparse.ArgumentParser(
        prog="dna-decode kleb",
        description="Klebsiella phage depolymerase -> capsule KL-type (fetch-only cross-organism cell).",
        epilog=_SCOPE)
    ap.add_argument("--depolymerase-fasta", required=True,
                    help="a phage depolymerase (tail-spike enzymatic domain) protein FASTA")
    ap.add_argument("--reference", help="path to the local KL-type reference (else auto-resolve / $DPO_KLEB_REFERENCE)")
    ap.add_argument("--top-k", type=int, default=5, help="number of ranked KL-types to return (default 5)")
    ap.add_argument("--json", action="store_true", help="emit the result as JSON")
    args = ap.parse_args(argv)

    from dna_decode.kleb.depolymerase_caller import call_kltype, load_reference, resolve_reference

    if not Path(args.depolymerase_fasta).exists():
        print(f"error: cannot read --depolymerase-fasta: {args.depolymerase_fasta}", file=sys.stderr)
        return 2
    ref_path = resolve_reference(args.reference)
    if ref_path is None:
        res = {"status": "INDETERMINATE", "ranked_kltypes": [], "reason": _FETCH_HINT, "scope": _SCOPE}
        print(json.dumps(res, indent=2) if args.json else f"INDETERMINATE: {_FETCH_HINT}\n{_SCOPE}")
        return 0
    prot = _read_first_protein(args.depolymerase_fasta)
    if not prot:
        print("error: no protein sequence found in --depolymerase-fasta", file=sys.stderr)
        return 2
    kmers, kltype = load_reference(ref_path)
    call = call_kltype(prot, kmers, kltype, top_k=args.top_k)
    res = {"status": call.status, "ranked_kltypes": list(call.ranked_kltypes),
           "top_similarity": call.top_similarity, "method": call.method,
           "reference": str(ref_path), "n_reference": len(kmers),
           "reason": call.reason, "scope": _SCOPE}
    if args.json:
        print(json.dumps(res, indent=2))
    else:
        if call.status == "CALLED":
            print(f"predicted capsule KL-type(s), ranked: {', '.join(call.ranked_kltypes)}")
            print(f"  top domain k-mer similarity: {round(call.top_similarity, 3)}  (ref: {len(kmers)} domains)")
        else:
            print(f"INDETERMINATE: {call.reason}")
        print(_SCOPE)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
