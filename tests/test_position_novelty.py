"""Tests for the generalized position-novelty self-awareness flag (pure, offline, CI-safe)."""
from __future__ import annotations

import pytest

from dna_decode.eval import position_novelty as PN


def test_parse_substitution():
    assert PN.parse_substitution("K103N") == ("K", 103, "N")
    assert PN.parse_substitution("E166*") == ("E", 166, "*")     # stop
    assert PN.parse_substitution(" G190A ") == ("G", 190, "A")   # whitespace-tolerant
    assert PN.parse_substitution("VF125AL") is None              # complex delins -> not a point sub
    assert PN.parse_substitution("junk") is None


def test_catalog_positions():
    assert PN.catalog_positions({"K103N", "K103S", "Y181C", "VF125AL"}) == {103, 181}


def test_flag_fires_on_novel_at_catalogued_position():
    catalog = {"K103N", "Y181C", "G190A"}
    # V179D is NOT catalogued and NOT at a catalogued position (179 not in {103,181,190}) -> no fire
    r = PN.position_novelty({"V179D"}, catalog)
    assert r.position_novel is False and r.at_catalog_positions == []
    # K103R IS at a catalogued position (103) but is NOT catalogued (only K103N is) -> FIRES
    r = PN.position_novelty({"K103R"}, catalog)
    assert r.position_novel is True and r.novel_substitutions == ["K103R"] and r.catalogued_hits == []
    # K103N is catalogued -> hit, not novel; flag does NOT fire on it alone
    r = PN.position_novelty({"K103N"}, catalog)
    assert r.position_novel is False and r.catalogued_hits == ["K103N"] and r.novel_substitutions == []


def test_mixed_genotype():
    catalog = {"K103N", "Y181C"}
    # carries a catalogued DRM (Y181C), a novel-at-position (K103R), and an off-position sub (V179D)
    r = PN.position_novelty({"Y181C", "K103R", "V179D"}, catalog)
    assert r.position_novel is True
    assert r.catalogued_hits == ["Y181C"]
    assert r.novel_substitutions == ["K103R"]
    assert r.at_catalog_positions == ["K103R", "Y181C"]   # V179D excluded (off catalogued positions)
    assert r.n_catalog_positions == 2


def test_faithful_to_hiv_logic():
    # equivalence to the HIV script's rule: filter observed to catalogued positions, fire if any not catalogued.
    catalog = {"K103N", "K103S", "V106A", "V106M"}
    # K103 has two catalogued substitutions; a K103T is novel-at-a-multiply-catalogued-position
    r = PN.position_novelty({"K103T"}, catalog)
    assert r.position_novel is True and r.novel_substitutions == ["K103T"]
    # both catalogued substitutions present -> no novelty
    r = PN.position_novelty({"K103N", "V106A"}, catalog)
    assert r.position_novel is False


def test_result_as_dict_no_phenotype_claim():
    r = PN.position_novelty({"K103R"}, {"K103N"})
    d = r.as_dict()
    assert set(d) == {"position_novel", "novel_substitutions", "catalogued_hits",
                      "at_catalog_positions", "n_catalog_positions"}
    # the flag never emits an R/S / prediction field
    assert "prediction" not in d and "resistance" not in d


def test_catalog_registry_multi_cell():
    # generalization proof: the registry serves >=2 distinct-kingdom cells from committed catalogs
    hiv = PN.catalog_drms_for("hiv-nnrti-rt")
    mpro = PN.catalog_drms_for("sarscov2-mpro")
    assert "K103N" in hiv and PN.catalog_positions(hiv)           # HIV RT NNRTI
    assert PN.catalog_positions(mpro)                             # SARS-CoV-2 Mpro
    with pytest.raises(KeyError):
        PN.catalog_drms_for("nonexistent-cell")


def test_flag_for_cell_real_catalogs():
    # a novel substitution at HIV RT position 103 fires against the real committed HIV catalog
    r = PN.flag_for_cell({"K103R"}, "hiv-nnrti-rt")
    assert r.position_novel is True and r.novel_substitutions == ["K103R"]
    # SARS-CoV-2 Mpro: a catalogued E166 substitution is a hit, not novel (uses the real Mpro catalog)
    mpro = PN.catalog_drms_for("sarscov2-mpro")
    e166 = next((m for m in mpro if m.startswith("E166")), None)
    if e166:
        r = PN.flag_for_cell({e166}, "sarscov2-mpro")
        assert r.position_novel is False and r.catalogued_hits == [e166]
