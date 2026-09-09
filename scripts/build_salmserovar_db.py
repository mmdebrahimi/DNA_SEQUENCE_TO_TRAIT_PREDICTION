"""Build the real Salmonella serovar DB (data/salmserovar_db/) for `dna-salmserovar` from SeqSero2.

Two artifacts, both derived from a SeqSero2 clone:
  * `serovar_table.tsv` (O<tab>H1<tab>H2<tab>Serovar) — the White-Kauffmann-Le Minor scheme, from SeqSero2's
    `bin/Initial_Conditions.py` parallel lists (`phaseO`/`phase1`/`phase2`/`sero`; ~2578 formulas).
  * `salmonella_antigens.fasta` (headers `<axis>__<antigen>__<id>__<seqsero2-header>`, axis in
    {O,H1,H2,SP}) — reformatted from SeqSero2's `seqsero2_db/H_and_O_and_specific_genes.fasta`
    (O-group wzx/wzy alleles + fliC=H1 + fljB=H2 + the two tyr genes SeqSero2's O logic consults).

THE 4TH FIELD IS THE FIX FOR A MEASURED DEFECT. The original schema was `<axis>__<antigen>__<id>`,
which has no field for the SEMANTIC QUALIFIER SeqSero2 attaches to six of its O entries
(`not_in_3,10`, `wzy_partial`, `wbaV`). Dropping it made the caller read a DIFFERENTIAL MARKER as a
standalone positive allele — measured at 60% hit-rate against a 2% true prevalence — and flattened the
O9-vs-O9,46 discriminator so plain O9 was unemittable. See
`wiki/salmserovar_db_semantics_2026-09-09.md` + `wiki/seqsero2_o_algorithm_reading_2026-09-09.md`.

The DB is gitignored-class (NOT committed); this script makes the build reproducible.

    git clone --depth 1 https://github.com/denglab/SeqSero2 <src>
    uv run python scripts/build_salmserovar_db.py --seqsero2-src <src> --out data/salmserovar_db

The antigen FASTA can also be rebuilt alone from a copy of SeqSero2's DB file (e.g. one extracted
from the pinned biocontainer), which needs no clone and leaves `serovar_table.tsv` untouched:

    uv run python scripts/build_salmserovar_db.py --antigens-fasta <H_and_O_...fasta> \
        --out data/salmserovar_db
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


# Specific-gene entries SeqSero2's *O* decision procedure consults (`call_O_and_H_type`, its
# `Special_dict` branch): tyrosinase splits O-9 from O-2 INSIDE the wbaV branch. Every other specific
# gene (oafA / sdf-I / gntR) feeds H-antigen or subspecies logic we do not implement, and is still
# dropped. These are carried on their OWN axis so they can never be selected as an O/H1/H2 antigen --
# `runner.parse_axis_antigen` returns None for any axis outside {O, H1, H2}.
SPECIAL_O_GENES = ("tyr-O-9", "tyr-O-2")


def _axis_antigen(header: str) -> tuple[str, str] | None:
    """SeqSero2 header -> (axis, antigen). O-4_wzx__id -> (O,4); fliC_i_...__id -> (H1,i);
    fljB_1,2_...__id -> (H2,1,2); tyr-O-9__609 -> (SP,tyr-O-9). Other specific genes
    (oafA/sdf-I/gntR) -> None (they feed logic this caller does not implement)."""
    first = header.split("_")[0]
    if first.startswith(("O-", "O:")):
        return "O", first[2:]                       # 'O-4'->'4' ; 'O-9,46'->'9,46' ; 'O:22'->'22'
    parts = header.split("_")
    if first == "fliC" and len(parts) >= 2:
        return "H1", parts[1]
    if first == "fljB" and len(parts) >= 2:
        return "H2", parts[1]
    for gene in SPECIAL_O_GENES:
        if header.startswith(gene + "_"):
            return "SP", gene
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
            if axis in ("O", "SP"):
                # 4th field = the ORIGINAL SeqSero2 header, verbatim. Every branch condition in
                # `call_O_and_H_type` names a literal key (`O-9,46_wbaV__1002`), so the role cannot be
                # reduced to a tidy enum without re-specifying the algorithm. It embeds `__` itself;
                # the reader recovers it with `"__".join(parts[3:])`, which is exact.
                if any(c.isspace() for c in h):
                    # blastn truncates `sseqid` at the first whitespace, so a header with a space
                    # would silently lose its role and read as an ordinary positive allele -- the
                    # exact failure this field exists to end. Refuse rather than emit a lossy entry.
                    raise ValueError(f"O/SP header contains whitespace, cannot round-trip: {h!r}")
                f.write(f">{axis}__{antigen}__{i}__{h}\n")
            else:
                # H1/H2 headers are left EXACTLY as they were. The role is consumed only by the O
                # decision procedure, and some fljB headers carry whitespace (which blastn truncates)
                # -- so adding a field there would be both useless and lossy. Keeping them
                # byte-identical also means the H axes CANNOT change, which is what the frozen
                # acceptance bar requires of them.
                f.write(f">{axis}__{antigen}__{i}\n")
            for j in range(0, len(s), 70):
                f.write(s[j:j + 70] + "\n")
            n += 1
    return n


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--seqsero2-src", type=Path, help="a SeqSero2 clone (builds BOTH artifacts)")
    ap.add_argument("--antigens-fasta", type=Path,
                    help="a copy of SeqSero2's H_and_O_and_specific_genes.fasta; rebuilds ONLY the "
                         "antigen FASTA and leaves serovar_table.tsv untouched")
    ap.add_argument("--out", type=Path, default=Path("data/salmserovar_db"))
    a = ap.parse_args(argv)
    if not a.seqsero2_src and not a.antigens_fasta:
        print("ERROR: one of --seqsero2-src / --antigens-fasta is required")
        return 2
    if a.antigens_fasta:
        if not a.antigens_fasta.exists():
            print(f"ERROR: {a.antigens_fasta} not found")
            return 2
        n_fasta = build_fasta(a.antigens_fasta, a.out / "salmonella_antigens.fasta")
        print(f"wrote {n_fasta} antigen alleles  -> {a.out / 'salmonella_antigens.fasta'}")
        print("serovar_table.tsv NOT rebuilt (pass --seqsero2-src for that)")
        return 0
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
