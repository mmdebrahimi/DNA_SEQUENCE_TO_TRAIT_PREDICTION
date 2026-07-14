"""Position-novelty self-awareness flag — the g5 refined survivor, GENERALIZED (2026-07-13).

`/innovate` round-1 refined survivor **g5-selfaware-flag**, validated for HIV NNRTI only in
`scripts/hiv_blindspot_position_novelty.py` (`FLAG_RECOVERS_BLINDSPOT`, median lift 3.98). This module
generalizes it into a CELL-AGNOSTIC detector usable across every target-site catalog (HIV RT/PR/IN/CA,
SARS-CoV-2 Mpro, fungal ERG11/FKS, …).

THE FLAG (deterministic, no model): given an observed genotype + a curated DRM catalog, `position_novel`
fires when the genotype carries a substitution AT a catalogued DRM position whose SPECIFIC substitution
is NOT itself catalogued. It is a "the catalog call may be INCOMPLETE here" self-awareness signal —
NOT a resistance prediction. It never emits R/S; it flags where a susceptible-by-absence catalog call is
least trustworthy (a novel substitution sits on a known resistance residue).

Faithful to the HIV logic: the HIV script filters observed substitutions to catalogued positions
(`_observed_rt_mutations` iterates `NNRTI_POSITIONS`) then fires if any is not itself catalogued. Here
the same two steps are explicit + catalog-driven, so ANY cell's DRM set drives the flag.

Pure + dependency-free. Catalog adapters lazy-import the committed cell catalogs on demand.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field

# a standard single-residue substitution: <WT><pos><MUT> (MUT may be a stop `*`). Complex delins
# (e.g. "VF125AL") deliberately do NOT parse -> excluded from position logic (not a point substitution).
_SUB = re.compile(r"^([A-Z])(\d+)([A-Z*])$")


def parse_substitution(s: str) -> tuple[str, int, str] | None:
    """'K103N' -> ('K', 103, 'N'); non-standard / delins -> None. Pure."""
    m = _SUB.match(s.strip())
    if not m:
        return None
    return m.group(1), int(m.group(2)), m.group(3)


def catalog_positions(catalog_drms) -> set[int]:
    """The set of residue POSITIONS covered by a DRM catalog (parseable single substitutions only)."""
    out: set[int] = set()
    for m in catalog_drms:
        p = parse_substitution(m)
        if p:
            out.add(p[1])
    return out


@dataclass
class PositionNoveltyResult:
    position_novel: bool
    novel_substitutions: list[str] = field(default_factory=list)   # at a catalogued position, NOT catalogued
    catalogued_hits: list[str] = field(default_factory=list)       # observed subs that ARE catalogued
    at_catalog_positions: list[str] = field(default_factory=list)  # observed subs sitting on a catalogued residue
    n_catalog_positions: int = 0

    def as_dict(self) -> dict:
        return {
            "position_novel": self.position_novel,
            "novel_substitutions": self.novel_substitutions,
            "catalogued_hits": self.catalogued_hits,
            "at_catalog_positions": self.at_catalog_positions,
            "n_catalog_positions": self.n_catalog_positions,
        }


def position_novelty(observed_substitutions, catalog_drms) -> PositionNoveltyResult:
    """Compute the position-novelty flag for one genotype against one catalog. Pure, cell-agnostic.

    `observed_substitutions` + `catalog_drms` are iterables of `<WT><pos><MUT>` strings. The flag fires
    when the genotype carries a substitution at a catalogued position whose specific substitution is not
    itself in the catalog.
    """
    catalog = set(catalog_drms)
    positions = catalog_positions(catalog)
    at_pos, novel, hits = [], [], []
    for s in observed_substitutions:
        p = parse_substitution(s)
        if p is None or p[1] not in positions:
            continue
        at_pos.append(s)
        if s in catalog:
            hits.append(s)
        else:
            novel.append(s)
    return PositionNoveltyResult(
        position_novel=bool(novel),
        novel_substitutions=sorted(set(novel)),
        catalogued_hits=sorted(set(hits)),
        at_catalog_positions=sorted(set(at_pos)),
        n_catalog_positions=len(positions),
    )


# --- cell catalog registry (lazy — adapters pull the committed DRM sets on demand) ----------------
def catalog_drms_for(cell: str) -> frozenset[str]:
    """Return the curated DRM set for a named cell. Lazy-imports the committed catalog. Raises on unknown.

    Cells (v0): 'hiv-nnrti-rt' / 'sarscov2-mpro' / 'fungal-fluconazole-erg11' /
    'fungal-voriconazole-erg11'. Extend by adding an adapter here (the flag itself is cell-agnostic).
    """
    cell = cell.strip().lower()
    if cell == "hiv-nnrti-rt":
        from ..data.hiv_amr import NNRTI_RT_MAJOR_DRMS
        return frozenset(NNRTI_RT_MAJOR_DRMS)
    if cell == "sarscov2-mpro":
        from ..data.sarscov2_amr import MPRO_MAJOR_DRMS
        return frozenset(MPRO_MAJOR_DRMS)
    if cell.startswith("fungal-") and cell.endswith("-erg11"):
        drug = cell[len("fungal-"):-len("-erg11")]
        from ..data.fungal_amr import FUNGAL_RESISTANCE_MUTATIONS
        by_gene = FUNGAL_RESISTANCE_MUTATIONS.get(drug, {})
        return frozenset(by_gene.get("ERG11", set()))
    raise KeyError(f"unknown cell {cell!r}; known: hiv-nnrti-rt / sarscov2-mpro / fungal-<drug>-erg11")


KNOWN_CELLS = ("hiv-nnrti-rt", "sarscov2-mpro", "fungal-fluconazole-erg11", "fungal-voriconazole-erg11")


def flag_for_cell(observed_substitutions, cell: str) -> PositionNoveltyResult:
    """Position-novelty flag for a genotype against a NAMED cell's committed catalog."""
    return position_novelty(observed_substitutions, catalog_drms_for(cell))
