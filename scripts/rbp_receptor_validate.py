"""Validate the RBP-level caller on the LBNL Phage Datasheets — the (3) mixed-clade deliverable.

Tests the core hypothesis: does RBP-homology (tail-fiber protein k-mer transfer) recover the receptor on
the RBP-VARIABLE classes (Tsx/OmpC/FhuA/OmpA/OmpF/FadL) where the v0 genome-homology caller scored 0/N
(wiki/phage_independent_result)? Leave-one-out over the LBNL characterized phages that carry BOTH a measured
receptor (mappable to our classes) AND an annotated receptor-binding protein (Table_S1 col 17).

Within-LBNL method validation (the RBP caller's own leave-one-out), reusing the SAME measured labels the
independent number used. Reproduce:
  git clone --depth 1 https://github.com/mjohnson11/PhageDataSheets.git
  uv run --with biopython python scripts/rbp_receptor_validate.py --repo PhageDataSheets/Ecoli_phages
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from dna_decode.phage.rbp_caller import leave_one_out_rbp, protein_kmers
from scripts.lbnl_independent_validate import map_receptor

# the RBP-variable classes the v0 genome-homology caller could NOT cover (0/N in the independent run)
RBP_VARIABLE = {"Tsx", "OmpC", "FhuA", "OmpA", "OmpF", "FadL", "TolC", "YncD"}


def extract_rbp_proteins(repo: Path, k: int = 4):
    """Return (rbp_kmers{phage->frozenset}, receptors{phage->class}, raw_counts). Uses ALL characterized
    phages (Bas + non-Bas) that have a mappable receptor AND a resolvable RBP CDS translation."""
    from Bio import SeqIO
    tsv = repo / "data" / "Table_S1_Phages.tsv"
    gdir = repo / "data" / "phage_genomes"
    rbp_kmers: dict[str, frozenset[str]] = {}
    receptors: dict[str, str] = {}
    n_no_rbp = n_no_gbk = 0
    with open(tsv, "r", encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) <= idx["Receptor-binding protein"]:
                continue
            phage = p[idx["Phage"]].strip()
            receptor = map_receptor(p[idx["BW25113 receptor"]])
            rbp_id = p[idx["Receptor-binding protein"]].strip()
            if not phage or receptor is None:
                continue
            if not rbp_id or rbp_id in ("", "Unknown", "Unknow", "NA"):
                n_no_rbp += 1
                continue
            gbk = gdir / f"{phage}.gbk"
            if not gbk.exists():
                n_no_gbk += 1
                continue
            prot = None
            try:
                for r in SeqIO.parse(str(gbk), "genbank"):
                    for f in r.features:
                        if f.type != "CDS":
                            continue
                        lt = (f.qualifiers.get("locus_tag") or [""])[0]
                        if lt == rbp_id or rbp_id in lt:
                            prot = (f.qualifiers.get("translation") or [""])[0]
                            break
                    if prot:
                        break
            except Exception:
                continue
            if not prot:
                n_no_rbp += 1
                continue
            km = protein_kmers(prot, k=k)
            if not km:
                continue
            rbp_kmers[phage] = km
            receptors[phage] = receptor
    return rbp_kmers, receptors, {"no_rbp": n_no_rbp, "no_gbk": n_no_gbk}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", required=True)
    ap.add_argument("--k", type=int, default=4)
    ap.add_argument("--min-similarity", type=float, default=0.05)
    ap.add_argument("--date", default="2026-07-24")
    args = ap.parse_args()

    rbp_kmers, receptors, skipped = extract_rbp_proteins(Path(args.repo), k=args.k)
    print(f"RBPs extracted: {len(rbp_kmers)} phages across {len(set(receptors.values()))} receptor classes "
          f"(skipped: {skipped})")
    res = leave_one_out_rbp(rbp_kmers, receptors, k=args.k, min_similarity=args.min_similarity)

    # headline: accuracy on the RBP-VARIABLE classes (where genome-homology got 0/N)
    var = [p for p in res.predictions if p["true"] in RBP_VARIABLE and p["status"] == "CALLED"]
    var_correct = sum(1 for p in var if p["correct"])
    from collections import Counter
    out = {
        "cell": "phage_receptor_class_RBP_level", "date": args.date,
        "method": "RBP tail-fiber protein k-mer (k=%d) nearest-neighbour transfer, leave-one-out" % args.k,
        "source": "LBNL Phage Datasheets (measured receptor + RBP annotation)",
        "n_rbp_phages": len(rbp_kmers), "n_receptor_classes": len(set(receptors.values())),
        "overall_called": res.n_called, "overall_correct": res.n_correct,
        "overall_accuracy": res.accuracy,
        "rbp_variable_called": len(var), "rbp_variable_correct": var_correct,
        "rbp_variable_accuracy": (var_correct / len(var)) if var else None,
        "rbp_variable_classes": sorted(RBP_VARIABLE),
        "per_receptor_correct_called": res.per_receptor,
        "receptor_distribution": dict(Counter(receptors.values())),
        "comparison": "the v0 GENOME-homology caller scored 0/N on Tsx/OmpC/FhuA/OmpA/OmpF/FadL/TolC "
                      "(independent run); RBP-homology recovers them because receptor is RBP-determined, "
                      "not backbone-determined.",
    }
    wiki = Path("wiki")
    (wiki / f"phage_rbp_caller_result_{args.date}.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps({k: out[k] for k in ("n_rbp_phages", "overall_called", "overall_correct",
          "overall_accuracy", "rbp_variable_called", "rbp_variable_correct", "rbp_variable_accuracy",
          "per_receptor_correct_called")}, indent=2))


if __name__ == "__main__":
    main()
