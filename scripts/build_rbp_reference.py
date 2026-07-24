"""Build the committed RBP reference FASTA for the wheel-only `dna-decode phage --rbp-fasta` caller.

Extracts the annotated tail-fiber receptor-binding protein of each characterized phage (Table_S1 col 17)
from the LBNL Phage Datasheets and writes data/phage_ref/rbp_reference.faa with headers
`>{phage}|{receptor}`. This makes the RBP-level caller self-contained (no repo clone at run time).

ATTRIBUTION: RBP protein sequences + receptor labels are from the LBNL/Arkin-Mutalik Phage Datasheets
(Moriniere et al. 2026; github.com/mjohnson11/PhageDataSheets), redistributed under its MIT License.

Usage:  uv run --with biopython python scripts/build_rbp_reference.py --repo <PhageDataSheets/Ecoli_phages>
"""
from __future__ import annotations

import argparse
from pathlib import Path

from scripts.lbnl_independent_validate import map_receptor

_HEADER = (
    "; dna-decode phage RBP reference — tail-fiber receptor-binding proteins + measured receptor class.\n"
    "; Source: LBNL/Arkin-Mutalik Phage Datasheets (Moriniere et al. 2026; github.com/mjohnson11/"
    "PhageDataSheets), MIT License. Receptor measured by genome-wide genetic screens on E. coli K-12 BW25113.\n"
    "; Header format: >{phage}|{receptor_class}\n"
)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", required=True)
    ap.add_argument("--out", default="data/phage_ref/rbp_reference.faa")
    args = ap.parse_args()

    from Bio import SeqIO
    repo = Path(args.repo)
    tsv = repo / "data" / "Table_S1_Phages.tsv"
    gdir = repo / "data" / "phage_genomes"
    entries: list[tuple[str, str, str]] = []
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
            if not phage or receptor is None or not rbp_id or rbp_id in ("Unknown", "Unknow", "NA"):
                continue
            gbk = gdir / f"{phage}.gbk"
            if not gbk.exists():
                continue
            prot = None
            try:
                for r in SeqIO.parse(str(gbk), "genbank"):
                    for f in r.features:
                        if f.type == "CDS" and ((f.qualifiers.get("locus_tag") or [""])[0] == rbp_id
                                                or rbp_id in (f.qualifiers.get("locus_tag") or [""])[0]):
                            prot = (f.qualifiers.get("translation") or [""])[0]
                            break
                    if prot:
                        break
            except Exception:
                continue
            if prot:
                entries.append((phage, receptor, prot))
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", encoding="utf-8") as fh:
        fh.write(_HEADER)
        for phage, receptor, prot in sorted(entries):
            fh.write(f">{phage}|{receptor}\n{prot}\n")
    print(f"wrote {len(entries)} RBP reference proteins -> {out} ({out.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
