"""`dna-motility` -- flagellar motility from gene presence (the first non-metabolic trait catalog).

    dna-motility --genes flhD,flhC,fliA,fliC,motA,motB,fliF,fliG,flhA,fliI   # -> MOTILE
    dna-motility --genes fliC,motA,motB                                       # -> NON-MOTILE (no master flhDC)
    dna-motility --feature-table GCF_..._feature_table.txt.gz                 # read present symbols genome-wide

Gene presence -> can the cell SWIM? Deterministic, offline, KNOWLEDGE_BASELINE. NOT clinical.
"""
from __future__ import annotations

import argparse
import gzip
import json
import sys
from pathlib import Path


def _feature_table_symbols(path: str) -> set[str]:
    op = gzip.open if str(path).endswith(".gz") else open
    syms: set[str] = set()
    with op(path, "rt", encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            isym = header.index("symbol")
        except ValueError:
            return syms
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) > isym and p[0] == "CDS" and p[isym]:
                syms.add(p[isym])
    return syms


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="dna-motility", description=__doc__)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--genes", help="comma-separated PRESENT gene symbols (e.g. flhD,flhC,fliC,motA,motB)")
    src.add_argument("--feature-table", help="NCBI feature_table.txt(.gz): read present symbols genome-wide")
    src.add_argument("--list", action="store_true", help="list the catalogued flagellar/chemotaxis modules")
    ap.add_argument("--json", action="store_true", help="emit JSON")
    a = ap.parse_args(argv)

    from dna_decode.motility.flagellar_catalog import (
        CHEMOTAXIS, MOTILITY_MODULES, MotilityInputError, call_motility)

    if a.list:
        for name, (comb, genes, role) in MOTILITY_MODULES.items():
            print(f"  {name:16} [{comb}] {', '.join(genes):22} {role}")
        print(f"  {'chemotaxis':16} [{CHEMOTAXIS[0]}] {', '.join(CHEMOTAXIS[1]):22} {CHEMOTAXIS[2]} (NOT a motility gate)")
        return 0

    if a.feature_table:
        if not Path(a.feature_table).exists():
            print(f"error: feature table not found: {a.feature_table}", file=sys.stderr)
            return 2
        present = _feature_table_symbols(a.feature_table)
    else:
        present = {g.strip() for g in a.genes.split(",") if g.strip()}

    try:
        call = call_motility(present)
    except MotilityInputError as e:
        print(f"error: {e}", file=sys.stderr)
        return 2

    rec = {
        "record": "motility-trait-v1",
        "trait": "flagellar swimming motility",
        "verdict": call.verdict,
        "motile": call.motile,
        "module_status": call.module_status,
        "missing_modules": list(call.missing_modules),
        "chemotaxis_competent": call.chemotaxis_competent,
        "note": call.note,
        "scope": "gene-presence -> swim/no-swim DIRECTION (not speed); KNOWLEDGE_BASELINE; NOT clinical.",
    }
    if a.json:
        print(json.dumps(rec, indent=2))
        return 0
    print(f"flagellar motility: {call.verdict}")
    for name, ok in call.module_status.items():
        print(f"  {'OK ' if ok else 'MISSING '}{name}")
    print(f"  chemotaxis-competent: {call.chemotaxis_competent}")
    print(f"  note: {call.note}")
    print(f"  scope: {rec['scope']}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
