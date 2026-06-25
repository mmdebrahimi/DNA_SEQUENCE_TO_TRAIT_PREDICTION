"""Build the pneumococcus cps-reference DB (`cps_references.fasta`) for `dna-pneumo-serotype` from PneumoCaT.

PneumoCaT (`phe-bioinformatics/PneumoCaT`) ships a Stage-1 all-serotype capsular locus reference at
`streptococcus-pneumoniae-ctvdb/reference.fasta` (95 sequences, one per serotype, headers = the serotype
name: `04`, `06A`, `19F`, ...). This script rewrites it under the `serotype__<ST>__01` header convention the
caller expects. The DB is gitignored-class (NOT committed); this script makes the build reproducible.

    git clone --depth 1 https://github.com/phe-bioinformatics/PneumoCaT <src>
    uv run python scripts/build_pneumo_cps_db.py --pneumocat-src <src> --out data/pneumoserotype_db

Provenance: PneumoCaT Stage-1 capsular reference (Kapatai et al. 2016, PeerJ; PHE).
"""
from __future__ import annotations

import argparse
from pathlib import Path


def build(ref_fasta: Path, out_fasta: Path) -> int:
    recs: list[tuple[str, str]] = []
    hdr, seq = None, []
    for line in ref_fasta.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if hdr is not None:
                recs.append((hdr, "".join(seq)))
            hdr, seq = line[1:].strip(), []
        else:
            seq.append(line.strip())
    if hdr is not None:
        recs.append((hdr, "".join(seq)))
    out_fasta.parent.mkdir(parents=True, exist_ok=True)
    with out_fasta.open("w", encoding="utf-8") as f:
        for st, s in recs:
            f.write(f">serotype__{st}__01\n")
            for i in range(0, len(s), 70):
                f.write(s[i:i + 70] + "\n")
    return len(recs)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pneumocat-src", type=Path, required=True, help="a PneumoCaT clone")
    ap.add_argument("--out", type=Path, default=Path("data/pneumoserotype_db"))
    a = ap.parse_args(argv)
    ref = a.pneumocat_src / "streptococcus-pneumoniae-ctvdb" / "reference.fasta"
    if not ref.exists():
        print(f"ERROR: PneumoCaT Stage-1 reference not found at {ref}")
        return 2
    n = build(ref, a.out / "cps_references.fasta")
    print(f"wrote {n} serotype references -> {a.out / 'cps_references.fasta'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
