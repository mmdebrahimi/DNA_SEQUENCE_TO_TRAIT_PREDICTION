"""AlphaMissense method for the forward variant-effect predictor — the HUMAN/eukaryotic upgrade.

AlphaMissense (Cheng et al. 2023, Science) ships a precomputed pathogenicity score in [0,1] for every
possible missense variant of the HUMAN proteome (higher = more pathogenic = more DAMAGING — the OPPOSITE
polarity to BLOSUM/ESM). It does NOT cover non-human proteins (that is why the bacterial cell can't use it).

`am_table_for_mutants` builds a {DMS-numbered mutation -> AM pathogenicity} table by looking up the cached
`am_filtered.tsv` for a protein's UniProt accession, applying the DMS->UniProt position OFFSET. predict_effect
then reads it as raw_score = 1 - AM (so higher = benign = preserved, consistent with the BLOSUM/ESM sign),
and tiers via AlphaMissense's own published class thresholds (benign <= 0.34, pathogenic >= 0.564).
"""
from __future__ import annotations

from pathlib import Path

# AlphaMissense published class thresholds (Cheng et al. 2023).
AM_BENIGN_MAX = 0.34
AM_PATHOGENIC_MIN = 0.564

_AM_CACHE: dict[str, dict[str, float]] = {}


def load_am_for_uniprot(am_tsv: Path, uniprot: str) -> dict[str, float]:
    """{variant_1letter (UniProt numbering, e.g. 'V133A') -> AM pathogenicity} for one UniProt accession.
    am_filtered.tsv is headerless: uniprot \\t protein_variant \\t am_pathogenicity \\t am_class."""
    key = f"{am_tsv}::{uniprot}"
    if key in _AM_CACHE:
        return _AM_CACHE[key]
    out: dict[str, float] = {}
    with open(am_tsv, encoding="utf-8") as fh:
        for ln in fh:
            p = ln.rstrip("\n").split("\t")
            if len(p) >= 3 and p[0] == uniprot:
                try:
                    out[p[1]] = float(p[2])
                except ValueError:
                    pass
    _AM_CACHE[key] = out
    return out


def am_table_for_mutants(am_by_variant: dict[str, float], offset: int, mutants) -> dict[str, float]:
    """{DMS-numbered mutation 'wt{pos}alt' -> AM pathogenicity}. UniProt pos = DMS pos + offset; a mutation
    with no AM entry is omitted (its DMS variant is not AlphaMissense-covered)."""
    table: dict[str, float] = {}
    for m in mutants:
        m = m.strip()
        if ":" in m or len(m) < 3 or not m[1:-1].isdigit():
            continue
        wt, pos, alt = m[0], int(m[1:-1]), m[-1]
        am = am_by_variant.get(f"{wt}{pos + offset}{alt}")
        if am is not None:
            table[m] = am
    return table


def am_tier(am: float) -> str:
    """AlphaMissense class -> forward-cell tier (benign->preserved, pathogenic->damaging, else uncertain)."""
    if am <= AM_BENIGN_MAX:
        return "preserved"
    if am >= AM_PATHOGENIC_MIN:
        return "damaging"
    return "uncertain"
