"""Build the real Salmonella serovar DB (data/salmserovar_db/) for `dna-salmserovar` from SeqSero2.

Two artifacts, both derived from a SeqSero2 clone:
  * `serovar_table.tsv` (O<tab>H1<tab>H2<tab>Serovar) — the White-Kauffmann-Le Minor scheme, from SeqSero2's
    `bin/Initial_Conditions.py` parallel lists (`phaseO`/`phase1`/`phase2`/`sero`; ~2578 formulas).
  * `salmonella_antigens.fasta` (headers `<axis>__<antigen>__<id>`, axis in {O,H1,H2}) — reformatted from
    SeqSero2's `seqsero2_db/H_and_O_and_specific_genes.fasta` (O-group wzx/wzy alleles + fliC=H1 + fljB=H2).

The DB is gitignored-class (NOT committed); this script makes the build reproducible.

    git clone --depth 1 https://github.com/denglab/SeqSero2 <src>
    uv run python scripts/build_salmserovar_db.py --seqsero2-src <src> --out data/salmserovar_db
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path


def build_table(ic_dir: Path, out_tsv: Path) -> int:
    sys.path.insert(0, str(ic_dir))
    import Initial_Conditions as IC  # noqa: E402
    table: dict[tuple[str, str, str], str] = {}
    for o, h1, h2, sv in zip(IC.phaseO, IC.phase1, IC.phase2, IC.sero):
        o, h1, h2, sv = (str(x).strip() for x in (o, h1, h2, sv))
        h2 = h2 or "-"
        if o and sv and (o, h1, h2) not in table:   # first-wins: subspecies-I named serovars come first
            table[(o, h1, h2)] = sv
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with out_tsv.open("w", encoding="utf-8") as f:
        f.write("O\tH1\tH2\tSerovar\n")
        for (o, h1, h2), sv in table.items():
            f.write(f"{o}\t{h1}\t{h2}\t{sv}\n")
    return len(table)


def _axis_antigen(header: str) -> tuple[str, str] | None:
    """SeqSero2 header -> (axis, antigen). O-4_wzx__id -> (O,4); fliC_i_...__id -> (H1,i);
    fljB_1,2_...__id -> (H2,1,2). Specific genes (oafA/sdf-I/tyr/gntR) -> None (v0 uses O/H1/H2 only)."""
    first = header.split("_")[0]
    if first.startswith(("O-", "O:")):
        return "O", first[2:]                       # 'O-4'->'4' ; 'O-9,46'->'9,46' ; 'O:22'->'22'
    parts = header.split("_")
    if first == "fliC" and len(parts) >= 2:
        return "H1", parts[1]
    if first == "fljB" and len(parts) >= 2:
        return "H2", parts[1]
    return None


def build_fasta(src_fasta: Path, out_fasta: Path) -> int:
    recs, hdr, seq = [], None, []
    for line in src_fasta.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if hdr is not None:
                recs.append((hdr, "".join(seq)))
            hdr, seq = line[1:].strip(), []
        else:
            seq.append(line.strip())
    if hdr is not None:
        recs.append((hdr, "".join(seq)))
    n = 0
    out_fasta.parent.mkdir(parents=True, exist_ok=True)
    with out_fasta.open("w", encoding="utf-8") as f:
        for i, (h, s) in enumerate(recs):
            aa = _axis_antigen(h)
            if aa is None or not s:
                continue
            axis, antigen = aa
            f.write(f">{axis}__{antigen}__{i}\n")
            for j in range(0, len(s), 70):
                f.write(s[j:j + 70] + "\n")
            n += 1
    return n


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--seqsero2-src", type=Path, required=True, help="a SeqSero2 clone")
    ap.add_argument("--out", type=Path, default=Path("data/salmserovar_db"))
    a = ap.parse_args(argv)
    ic_dir = a.seqsero2_src / "bin"
    src_fasta = a.seqsero2_src / "seqsero2_db" / "H_and_O_and_specific_genes.fasta"
    if not (ic_dir / "Initial_Conditions.py").exists() or not src_fasta.exists():
        print(f"ERROR: SeqSero2 layout not found under {a.seqsero2_src}")
        return 2
    n_table = build_table(ic_dir, a.out / "serovar_table.tsv")
    n_fasta = build_fasta(src_fasta, a.out / "salmonella_antigens.fasta")
    print(f"wrote {n_table} serovar formulas -> {a.out / 'serovar_table.tsv'}")
    print(f"wrote {n_fasta} antigen alleles  -> {a.out / 'salmonella_antigens.fasta'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
