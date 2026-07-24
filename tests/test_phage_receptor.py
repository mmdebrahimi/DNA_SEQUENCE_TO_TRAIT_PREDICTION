"""Tests for the phage host-receptor catalog (dna_decode/data/phage_receptor.py).

Offline + deterministic: catalog integrity, VERBATIM-grounded taxon->receptor lookups, the
honest INDETERMINATE on an uncatalogued taxon, and reference-phage self-consistency. No network.
"""
from __future__ import annotations

from dna_decode.data import phage_receptor as pr


def test_receptor_classes_are_frozen_and_named():
    # every receptor named in the BASEL Results must be in the label space
    for r in ("FhuA", "BtuB", "YncD", "TolC", "LptD", "LamB", "FepA",
              "OmpC", "OmpF", "OmpA", "Tsx", "FadL", "NfrA", "LPS_core", "ECA"):
        assert pr.is_receptor_class(r), f"{r} missing from RECEPTOR_CLASSES"
    assert not pr.is_receptor_class("NotAReceptor")


def test_every_taxon_primary_is_a_valid_receptor_class():
    # no catalog entry may point at a receptor outside the label space (anti-typo / anti-fabrication)
    for table in (pr.FAMILY_RECEPTOR, pr.GENUS_RECEPTOR):
        for tr in table.values():
            assert pr.is_receptor_class(tr.primary), f"{tr.taxon} primary {tr.primary} invalid"
            assert tr.primary in tr.receptors
            for r in tr.receptors:
                assert pr.is_receptor_class(r), f"{tr.taxon} lists invalid receptor {r}"


def test_verbatim_receptor_assignments():
    # spot-check the load-bearing VERBATIM assignments from Maffei 2021
    assert pr.primary_receptor_for_taxon("Tequatrovirus") == "OmpC"      # T4-like
    assert pr.primary_receptor_for_taxon("Tequintavirus") == "BtuB"      # T5-like
    assert pr.primary_receptor_for_taxon("Teseptimavirus") == "LPS_core" # T7-like
    assert pr.primary_receptor_for_taxon("Enquatrovirus") == "NfrA"      # N4-like
    assert pr.primary_receptor_for_taxon("Lambdavirus") == "LamB"        # lambda
    assert pr.primary_receptor_for_taxon("Drexlerviridae") == "FhuA"     # family rank
    assert pr.primary_receptor_for_taxon("Demerecviridae") == "BtuB"     # T5-like family
    assert pr.primary_receptor_for_taxon("Schitoviridae") == "NfrA"


def test_uncatalogued_taxon_is_indeterminate_not_guessed():
    assert pr.receptor_for_taxon("Caudoviricetes") is None       # too high a rank
    assert pr.receptor_for_taxon("SomeNovelGenus") is None
    assert pr.primary_receptor_for_taxon("Unknownvirus") is None


def test_lineage_is_first_match_in_order_loader_passes_genus_first():
    # receptor_for_lineage returns the FIRST catalogued taxon in the given order; the loader passes
    # [genus, family] so a catalogued genus wins over its family fallback.
    tr = pr.receptor_for_lineage(["Tequintavirus", "Demerecviridae"])   # genus first
    assert tr is not None and tr.primary == "BtuB" and tr.rank == "genus"
    # a phage whose genus is uncatalogued falls back to the family assignment
    tr2 = pr.receptor_for_lineage(["SomeDrexlerGenus", "Drexlerviridae"])
    assert tr2 is not None and tr2.taxon == "Drexlerviridae" and tr2.primary == "FhuA"
    # nothing catalogued -> None
    assert pr.receptor_for_lineage(["Novelgenus", "Novelviridae"]) is None


def test_label_only_for_clade_conserved_taxa():
    # clade-conserved clades yield a training label...
    assert pr.label_receptor_for_lineage(["Vequintavirus", "unclassified"]) == "ECA"
    assert pr.label_receptor_for_lineage(["Tequintavirus", "Demerecviridae"]) == "BtuB"
    assert pr.label_receptor_for_lineage(["Teseptimavirus", "Autotranscriptaviridae"]) == "LPS_core"
    assert pr.label_receptor_for_lineage(["Augustepiccardvirus", "Drexlerviridae"]) == "LptD"
    # ...but RBP-variable clades (T-even, Drexlerviridae) return None, NOT a fabricated single label
    assert pr.label_receptor_for_lineage(["Tequatrovirus", "Straboviridae"]) is None
    assert pr.label_receptor_for_lineage(["SomeDrexlerGenus", "Drexlerviridae"]) is None
    assert pr.label_receptor_for_lineage(["Novelgenus", "unclassified"]) is None


def test_rbp_variable_clades_are_flagged():
    assert pr.GENUS_RECEPTOR["Tequatrovirus"].clade_conserved is False
    assert pr.FAMILY_RECEPTOR["Drexlerviridae"].clade_conserved is False
    assert pr.FAMILY_RECEPTOR["Straboviridae"].clade_conserved is False
    assert pr.GENUS_RECEPTOR["Vequintavirus"].clade_conserved is True


def test_reference_phages_are_self_consistent():
    assert len(pr.REFERENCE_PHAGES) >= 6
    for rp in pr.REFERENCE_PHAGES:
        assert pr.is_receptor_class(rp.receptor), f"{rp.name} receptor {rp.receptor} invalid"
        assert rp.accession.startswith("NC_"), f"{rp.name} accession looks wrong"
        # the reference phage's genus_hint must resolve to a receptor consistent with its label
        tr = pr.receptor_for_taxon(rp.genus_hint)
        assert tr is not None, f"{rp.name} genus_hint {rp.genus_hint} not catalogued"
        assert rp.receptor in tr.receptors, f"{rp.name} {rp.receptor} not in {rp.genus_hint} receptors"
