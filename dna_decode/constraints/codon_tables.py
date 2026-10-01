"""Clade-aware genetic-code tables — ONE canonical source for the whole repo.

WHY THIS EXISTS. Before this module the genetic code was hardcoded as universal in five separate places
(`forward/genome_edit.py`, `forward/inverse.py`, `typing/codon_map.py`, `scripts/forward_inverse_roundtrip.py`,
`scripts/fungal_erg11_caller.py`). All three in-package copies were verified byte-identical standard tables
with `TGA/TAA/TAG -> '*'`, i.e. NCBI translation table 1. That is **wrong for real clades**: Mycoplasma and
Spiroplasma read UGA as tryptophan, vertebrate mitochondria read AGA/AGG as stop, ciliates read UAA/UAG as
glutamine. Pointed at a Mycoplasma gene, a table-1-only caller reports a normal tryptophan as a nonsense
mutation. Not a live defect for the shipped cells (E. coli, Klebsiella, HIV, SARS-CoV-2, Candida, TB all use
standard assignments) -- a latent assumption that becomes a defect the moment the organism set widens.

PROVENANCE: EVERY table is stored as the four VERBATIM NCBI strings (`AAs` / `Base1` / `Base2` / `Base3`)
exactly as printed by https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi, and the 64-codon mapping is
DERIVED from them by `_decode`. Nothing here is hand-typed from memory -- writing reference biology from
recall is the fabrication hazard this project refuses elsewhere (a volunteered ECOFF value was rejected for
exactly that reason). A table added without `source_url` + `verbatim_quote` is REFUSED by `table_for`.

REFUSAL over guessing: `table_for` raises `UnknownCladeError` for an unrecognised clade rather than quietly
returning the standard table. Applying the wrong table corrupts every call downstream, so a declared unknown
is strictly safer than a plausible default.
"""
from __future__ import annotations

from dataclasses import dataclass

_BASES = "TCAG"
NCBI_SOURCE_URL = "https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi"


class UnknownCladeError(KeyError):
    """Raised when a clade has no declared translation table. Never defaults to standard."""


class UnverifiedTableError(RuntimeError):
    """Raised when a table lacks `source_url` + `verbatim_quote`. Provenance is not optional."""


@dataclass(frozen=True)
class TranslationTable:
    """One NCBI genetic code, stored as its verbatim printed strings plus derived provenance."""

    table_id: int
    name: str
    aas: str            # the verbatim `AAs` line
    base1: str          # the verbatim `Base1` line
    base2: str          # the verbatim `Base2` line
    base3: str          # the verbatim `Base3` line
    source_url: str
    differences_from_standard: str

    @property
    def verbatim_quote(self) -> str:
        """The verbatim NCBI block this table was transcribed from."""
        return (f"AAs  = {self.aas}\nBase1  = {self.base1}\n"
                f"Base2  = {self.base2}\nBase3  = {self.base3}")

    @property
    def verified(self) -> bool:
        return bool(self.source_url) and bool(self.aas)


# The four strings per table are VERBATIM from NCBI (fetched 2026-10-01). Do not edit by hand; re-fetch.
_B1 = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
_B2 = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
_B3 = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"

TABLES: dict[int, TranslationTable] = {
    1: TranslationTable(
        1, "The Standard Code",
        "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        _B1, _B2, _B3, NCBI_SOURCE_URL, "Not applicable (this is the standard)."),
    2: TranslationTable(
        2, "The Vertebrate Mitochondrial Code",
        "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
        _B1, _B2, _B3, NCBI_SOURCE_URL,
        "AGA: Ter (*) vs. Arg (R); AGG: Ter (*) vs. Arg (R); AUA: Met (M) vs. Ile (I); "
        "UGA: Trp (W) vs. Ter (*)"),
    4: TranslationTable(
        4, "The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code",
        "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        _B1, _B2, _B3, NCBI_SOURCE_URL, "UGA: Trp (W) vs. Ter (*)"),
    6: TranslationTable(
        6, "The Ciliate, Dasycladacean and Hexamita Nuclear Code",
        "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        _B1, _B2, _B3, NCBI_SOURCE_URL,
        "UAA: Gln (Q) vs. Ter (*); UAG: Gln (Q) vs. Ter (*)"),
    11: TranslationTable(
        11, "The Bacterial, Archaeal and Plant Plastid Code",
        "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        _B1, _B2, _B3, NCBI_SOURCE_URL,
        "Internal assignments are the same as in the standard code (differs only in permitted start "
        "codons, which this layer does not model)."),
}


def _decode(t: TranslationTable) -> dict[str, str]:
    """Build the 64-codon map from the verbatim NCBI strings. Derivation, not transcription."""
    if not t.verified:
        raise UnverifiedTableError(
            f"translation table {t.table_id} lacks source_url/verbatim_quote; refusing to use it")
    for nm, s in (("AAs", t.aas), ("Base1", t.base1), ("Base2", t.base2), ("Base3", t.base3)):
        if len(s) != 64:
            raise ValueError(f"table {t.table_id}: {nm} line is {len(s)} chars, expected 64")
    out: dict[str, str] = {}
    for i in range(64):
        out[t.base1[i] + t.base2[i] + t.base3[i]] = t.aas[i]
    if len(out) != 64:
        raise ValueError(f"table {t.table_id}: decoded {len(out)} distinct codons, expected 64")
    return out


# Clade key -> NCBI table id. Nested keys use dots; specialization is resolved by the registry.
TABLE_IDS: dict[str, int] = {
    "standard": 1,
    "eukaryote": 1,
    "eukaryote.nuclear": 1,
    "virus": 1,                      # viruses are translated by the host's machinery
    "bacteria": 11,
    "archaea": 11,
    "plant.plastid": 11,
    "bacteria.mycoplasma": 4,
    "bacteria.spiroplasma": 4,
    "mitochondrion.mold": 4,
    "mitochondrion.protozoan": 4,
    "mitochondrion.coelenterate": 4,
    "mitochondrion.vertebrate": 2,
    "eukaryote.ciliate": 6,
    "eukaryote.dasycladacean": 6,
    "eukaryote.hexamita": 6,
}


def table_for(clade: str) -> dict[str, str]:
    """The 64-codon map for `clade`. REFUSES an unrecognised clade — never defaults to standard."""
    key = (clade or "").strip()
    if key not in TABLE_IDS:
        raise UnknownCladeError(
            f"no translation table declared for clade {clade!r}. Declare it in TABLE_IDS with an NCBI "
            f"table id rather than assuming the standard code. Known: {sorted(TABLE_IDS)}")
    return _decode(TABLES[TABLE_IDS[key]])


def table_meta(clade: str) -> TranslationTable:
    """The provenance record behind `clade`'s table."""
    key = (clade or "").strip()
    if key not in TABLE_IDS:
        raise UnknownCladeError(f"no translation table declared for clade {clade!r}")
    return TABLES[TABLE_IDS[key]]


#: The canonical standard code. Every in-repo codon table aliases THIS object.
STANDARD: dict[str, str] = _decode(TABLES[1])

#: Convenience aliases for the clade tables that differ from standard.
VERTEBRATE_MITOCHONDRIAL: dict[str, str] = _decode(TABLES[2])
MYCOPLASMA_SPIROPLASMA: dict[str, str] = _decode(TABLES[4])
CILIATE_NUCLEAR: dict[str, str] = _decode(TABLES[6])
BACTERIAL_PLASTID: dict[str, str] = _decode(TABLES[11])

#: Stop codons UNDER THE STANDARD CODE ONLY. The name carries the scope deliberately: an unqualified
#: `STOPS` was the clade-blind trap this module exists to remove, surviving in its one derived constant --
#: `TGA` is in it but is TRYPTOPHAN under table 4, and `AGA`/`AGG` are absent from it but ARE stops under
#: table 2. Use `stops_for(clade)` whenever a clade is known.
STANDARD_STOPS: frozenset[str] = frozenset(c for c, aa in STANDARD.items() if aa == "*")


def stops_for(clade: str) -> frozenset[str]:
    """The stop codons for `clade`'s genetic code. Refuses an unknown clade, like `table_for`."""
    return frozenset(c for c, aa in table_for(clade).items() if aa == "*")
