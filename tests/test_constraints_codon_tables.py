"""Guards for the canonical clade-aware codon tables (plan Step 1, 2026-10-01)."""

from __future__ import annotations

import pytest

from dna_decode.constraints import codon_tables as T

# The PRE-CHANGE values, pinned as literals. Three in-package tables were verified byte-identical
# standard tables before this plan ran; these two assignments are the ones a clade table changes.
PRE_CHANGE_TGA = "*"
PRE_CHANGE_AGA = "R"


def test_standard_is_64_codons_with_the_three_canonical_stops():
    assert len(T.STANDARD) == 64
    assert sorted(T.STANDARD_STOPS) == ["TAA", "TAG", "TGA"]
    # clade-aware accessor: table 4 loses TGA, table 2 gains AGA/AGG
    assert "TGA" not in T.stops_for("bacteria.mycoplasma")
    assert {"AGA", "AGG"} <= T.stops_for("mitochondrion.vertebrate")
    assert sorted(T.stops_for("standard")) == ["TAA", "TAG", "TGA"]


def test_standard_matches_the_pre_change_literals():
    """Pure-refactor pin: consolidation must not move a single assignment."""
    assert T.STANDARD["TGA"] == PRE_CHANGE_TGA
    assert T.STANDARD["AGA"] == PRE_CHANGE_AGA


def test_every_table_is_DERIVED_from_its_verbatim_ncbi_strings():
    """Nothing here may be hand-typed from memory. The 64-codon map is decoded from the four verbatim
    NCBI lines, so the stored provenance and the usable data cannot drift apart."""
    for tid, meta in T.TABLES.items():
        assert len(meta.aas) == 64, tid
        assert len(meta.base1) == len(meta.base2) == len(meta.base3) == 64, tid
        assert meta.source_url.startswith("https://www.ncbi.nlm.nih.gov/"), tid
        assert "AAs  =" in meta.verbatim_quote and "Base3  =" in meta.verbatim_quote, tid
        decoded = T._decode(meta)
        assert len(decoded) == 64, tid
        # the decode must be reproducible from the quote alone
        assert decoded[meta.base1[0] + meta.base2[0] + meta.base3[0]] == meta.aas[0], tid


def test_the_derived_differences_match_ncbis_own_stated_differences():
    """THE CROSS-CHECK that the transcription is faithful: NCBI prints both an `AAs` line and a prose
    'Differences from Standard Code' list. We decode the former; if the two disagree, the transcription
    is wrong. Table 11 is the control -- NCBI states its internal assignments equal the standard."""
    diff = lambda tab: sorted(c for c in T.STANDARD if tab[c] != T.STANDARD[c])
    assert diff(T.MYCOPLASMA_SPIROPLASMA) == ["TGA"]                      # UGA: Trp vs Ter
    assert diff(T.VERTEBRATE_MITOCHONDRIAL) == ["AGA", "AGG", "ATA", "TGA"]  # four stated differences
    assert diff(T.CILIATE_NUCLEAR) == ["TAA", "TAG"]                      # UAA/UAG: Gln vs Ter
    assert T.BACTERIAL_PLASTID == T.STANDARD                              # control


def test_the_clade_tables_are_load_bearing_not_decorative():
    """If a clade table never changes an assignment it buys nothing. These are the specific reassignments
    that would make a table-1-only caller report a false nonsense mutation."""
    assert T.MYCOPLASMA_SPIROPLASMA["TGA"] == "W" and T.STANDARD["TGA"] == "*"
    assert T.VERTEBRATE_MITOCHONDRIAL["AGA"] == "*" and T.STANDARD["AGA"] == "R"
    assert T.CILIATE_NUCLEAR["TAA"] == "Q" and T.STANDARD["TAA"] == "*"


def test_table_for_REFUSES_an_unknown_clade_rather_than_defaulting():
    """A wrong table corrupts every downstream call, so a declared unknown is safer than a plausible
    default. Non-vacuous: a known clade must succeed in the same test."""
    assert T.table_for("bacteria.mycoplasma")["TGA"] == "W"
    with pytest.raises(T.UnknownCladeError):
        T.table_for("nonsense_clade")
    with pytest.raises(T.UnknownCladeError):
        T.table_for("")


def test_an_unsourced_table_is_REFUSED_and_the_guard_is_non_vacuous():
    """Provenance is not optional. Proven both directions: a sourced table decodes, an unsourced one
    raises -- so the check cannot pass by never firing."""
    good = T.TABLES[1]
    assert len(T._decode(good)) == 64
    unsourced = T.TranslationTable(
        99, "fabricated", good.aas, good.base1, good.base2, good.base3, "", "none")
    with pytest.raises(T.UnverifiedTableError):
        T._decode(unsourced)


def test_a_malformed_table_is_rejected_rather_than_silently_truncated():
    good = T.TABLES[1]
    short = T.TranslationTable(
        98, "short", good.aas[:60], good.base1, good.base2, good.base3, T.NCBI_SOURCE_URL, "none")
    with pytest.raises(ValueError, match="expected 64"):
        T._decode(short)


def test_table_meta_carries_provenance_for_every_declared_clade():
    for clade in T.TABLE_IDS:
        meta = T.table_meta(clade)
        assert meta.verified is True, clade
        assert meta.name and meta.differences_from_standard, clade
