"""Offline tests for the clinical-gene landscape census (scripts/clinical_gene_landscape_census.py).

No MaveDB/ClinVar/AM network — pins the offset-application, curated gene->UniProt fallback, and the state
classifier via a monkeypatched census_gene path. The AUROC math is covered by test_clinical_variant_effect.py.
"""
from __future__ import annotations

import scripts.clinical_gene_landscape_census as census
from scripts.clinical_gene_landscape_census import CURATED_UNIPROT, census_gene


def test_curated_uniprot_covers_every_clinical_gene():
    # every gene in the census list must have a curated UniProt fallback (so none is lost to a missing MaveDB id)
    for g in census.CLINICAL_GENES:
        assert g in CURATED_UNIPROT, f"{g} missing a curated UniProt fallback"
    # spot-check a few canonical accessions
    assert CURATED_UNIPROT["TP53"] == "P04637"
    assert CURATED_UNIPROT["LDLR"] == "P01130"
    assert CURATED_UNIPROT["F9"] == "P00740"


def test_census_gene_uses_curated_when_mavedb_id_absent(monkeypatch):
    # MaveDB assay has NO uniprot -> curated fallback fills it; offset stays from the assay
    monkeypatch.setattr(census, "fetch_dms_offset", lambda urn, off: {("R", 175, "H"): -2.0, ("P", 72, "R"): 0.5})
    monkeypatch.setattr(census, "fetch_clinvar_missense", lambda g, use_cache=True: {
        ("R", 175, "H"): "PATH", ("P", 72, "R"): "BENIGN"})
    monkeypatch.setattr(census, "load_am", lambda up: {("R", 175, "H"): 0.99, ("P", 72, "R"): 0.05})
    rec = census_gene("TP53", {"uniprot": None, "offset": 0, "urn": "urn:x", "n_variants": 100})
    assert rec["uniprot"] == "P04637"
    assert rec["uniprot_source"] == "curated"


def test_census_gene_no_uniprot_when_absent_and_uncurated(monkeypatch):
    rec = census_gene("ZZZ_UNKNOWN", {"uniprot": None, "offset": 0, "urn": "urn:x", "n_variants": 10})
    assert rec["state"] == "NO_UNIPROT_ID"


def test_census_gene_offset_is_applied_to_the_join(monkeypatch):
    # DMS numbered 1-based in a domain; offset 486 must shift positions before the ClinVar/AM join.
    # ClinVar/AM are keyed at UniProt positions (mavedb_pos + 486).
    def fake_dms(urn, off):
        return {(w, p + off, a): s for (w, p, a), s in {("A", 1, "V"): 0.1, ("G", 2, "R"): 0.2}.items()}
    monkeypatch.setattr(census, "fetch_dms_offset", fake_dms)
    # a benign-poor overlap -> SINGLE_CLASS (both mapped to UniProt 487/488)
    monkeypatch.setattr(census, "fetch_clinvar_missense", lambda g, use_cache=True: {
        ("A", 487, "V"): "PATH", ("G", 488, "R"): "PATH"})
    monkeypatch.setattr(census, "load_am", lambda up: {("A", 487, "V"): 0.9, ("G", 488, "R"): 0.8})
    rec = census_gene("MLH1", {"uniprot": "P40692", "offset": 486, "urn": "urn:x", "n_variants": 5000})
    # both joined at the shifted positions -> proves the offset was applied
    assert rec["n_joined"] == 2
    assert rec["state"] == "SINGLE_CLASS"   # 2 path / 0 benign
    assert rec["n_path"] == 2 and rec["n_benign"] == 0
