"""Offline tests for the gnomAD frequency-benign source + gnomAD-benign census (no network).

Pins the hgvsp parser, the AF-threshold filter, the exclude-contradictory-variant join, and the
circularity-flagged AM reporting.
"""
from __future__ import annotations

from scripts.gnomad_benign import parse_hgvsp, fetch_gnomad_benign, DEFAULT_AF_MIN
import scripts.gnomad_benign as gb
import scripts.clinical_gnomad_benign_census as census


def test_parse_hgvsp_single_missense():
    assert parse_hgvsp("p.Gly324Ser") == ("G", 324, "S")
    assert parse_hgvsp("p.Asp342Asn") == ("D", 342, "N")


def test_parse_hgvsp_rejects_non_missense():
    assert parse_hgvsp("p.Trp212Ter") is None      # nonsense
    assert parse_hgvsp("p.=") is None               # synonymous
    assert parse_hgvsp("p.Gly324=") is None
    assert parse_hgvsp("c.100A>G") is None          # not protein
    assert parse_hgvsp("") is None


def test_fetch_gnomad_benign_af_threshold(monkeypatch, tmp_path):
    # feed a cached raw list; only AF>=af_min missense survive
    raw = [{"hgvsp": "p.Gly324Ser", "af": 0.005}, {"hgvsp": "p.Asp342Asn", "af": 0.00005},
           {"hgvsp": "p.Ala391Thr", "af": 0.08}, {"hgvsp": "p.Trp10Ter", "af": 0.9}]
    monkeypatch.setattr(gb, "GNOMAD_CACHE", tmp_path)
    (tmp_path / "LDLR.json").write_text(__import__("json").dumps(raw), encoding="utf-8")
    out = fetch_gnomad_benign("LDLR", af_min=1e-4)
    assert ("G", 324, "S") in out and ("A", 391, "T") in out
    assert ("D", 342, "N") not in out              # below threshold
    assert all(len(k) == 3 for k in out)           # nonsense p.Trp10Ter dropped by the parser


def test_default_af_min_is_acmg_adjacent():
    assert DEFAULT_AF_MIN == 1e-4


def test_census_gnomad_excludes_contradictory_and_flags_am_circular(monkeypatch):
    # positives = ClinVar-path; negatives = gnomAD-benign; a variant that is BOTH is dropped.
    dms = {("R", i, "H"): float(-i) for i in range(1, 21)}          # 20 pathogenic-ish
    dms.update({("I", i, "V"): float(i) for i in range(21, 41)})     # 20 benign-ish
    monkeypatch.setattr(census, "fetch_dms_offset", lambda urn, off: dms)
    monkeypatch.setattr(census, "fetch_clinvar_missense", lambda g, use_cache=True:
                        {("R", i, "H"): "PATH" for i in range(1, 21)})
    monkeypatch.setattr(census, "fetch_gnomad_benign", lambda g, af_min: {("I", i, "V"): 0.01 for i in range(21, 41)})
    monkeypatch.setattr(census, "load_am", lambda up: {**{k: 0.9 for k in dms}})  # AM covers all
    rec = census.census_gene_gnomad("TESTG", {"uniprot": "P00000", "offset": 0, "urn": "urn:x"}, DEFAULT_AF_MIN)
    assert rec["state"] == "SCORED"
    assert rec["n_path_clinvar"] == 20 and rec["n_benign_gnomad"] == 20
    assert "dms_ceiling_auroc" in rec and "blosum_floor_auroc" in rec
    # AM is reported under the CIRCULAR-flagged key, never a plain am_auroc
    assert "am_auroc_CIRCULAR" in rec
    assert "am_auroc" not in rec


def test_census_gnomad_still_underpowered_when_benign_sparse(monkeypatch):
    monkeypatch.setattr(census, "fetch_dms_offset", lambda urn, off: {("R", i, "H"): float(i) for i in range(1, 40)})
    monkeypatch.setattr(census, "fetch_clinvar_missense", lambda g, use_cache=True:
                        {("R", i, "H"): "PATH" for i in range(1, 40)})
    monkeypatch.setattr(census, "fetch_gnomad_benign", lambda g, af_min: {("R", 1, "H"): 0.01})  # 1 benign only
    monkeypatch.setattr(census, "load_am", lambda up: {})
    rec = census.census_gene_gnomad("KRAS", {"uniprot": "P01116", "offset": 0, "urn": "urn:x"}, DEFAULT_AF_MIN)
    assert rec["state"] == "STILL_UNDERPOWERED"
