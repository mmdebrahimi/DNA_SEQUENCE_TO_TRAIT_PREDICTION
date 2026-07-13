"""Tests for the gene->sequence lookup layer (pure record logic offline; live UniProt fetch network-guarded)."""
from __future__ import annotations

import os

import pytest

from dna_decode.protein_effect import gene_lookup as G


def _hit(acc, seq, gene="gyrA", org="Escherichia coli (strain K12)"):
    return {"primaryAccession": acc, "sequence": {"value": seq}, "organism": {"scientificName": org},
            "entryType": "UniProtKB reviewed (Swiss-Prot)",
            "proteinDescription": {"recommendedName": {"fullName": {"value": "DNA gyrase subunit A"}}}}


def test_record_single_hit_assembles_provenance():
    rec = G._record("gyrA", "83333", [_hit("P0AES4", "MSDLA")], reviewed=True)
    assert rec["accession"] == "P0AES4" and rec["sequence"] == "MSDLA" and rec["length"] == 5
    assert rec["reviewed"] is True and rec["protein_name"] == "DNA gyrase subunit A"
    assert "canonical" in rec["provenance"]["numbering"].lower() and rec["provenance"]["organism_id"] == "83333"


def test_record_not_found_and_ambiguous():
    with pytest.raises(G.GeneLookupNotFound):
        G._record("nope", "83333", [], reviewed=True)
    with pytest.raises(G.GeneLookupAmbiguous):
        G._record("dup", "83333", [_hit("P1", "MA"), _hit("P2", "MB")], reviewed=True)


def test_cache_path_shape(tmp_path):
    p = G._cache_path("GyrA", "83333", tmp_path)
    assert p.name == "83333_gyra.json"          # lowercased gene, organism-prefixed


def test_pinned_organism_is_ecoli_k12():
    assert G.ECOLI_K12["organism_id"] == "83333"


@pytest.mark.skipif(os.environ.get("DNA_DECODE_ALLOW_NET") != "1",
                    reason="live UniProt fetch — set DNA_DECODE_ALLOW_NET=1 to run")
def test_live_fetch_gyra(tmp_path):
    rec = G.fetch_protein_sequence("gyrA", cache_dir=tmp_path)
    assert rec["accession"] == "P0AES4" and rec["length"] == 875
    assert rec["sequence"][82] == "S"           # gyrA S83 (the QRDR position) — numbering face-check
