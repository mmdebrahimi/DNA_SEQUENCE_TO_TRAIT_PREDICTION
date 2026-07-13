"""Tests for the rung-2 + rung-3 fused query (phenotype routing offline; molecular rung on a tiny synthetic logp)."""
from __future__ import annotations

import pytest

from dna_decode.protein_effect import integration as IQ


def test_hiv_phenotype_known_drm():
    # K103N is a catalogued NNRTI DRM -> known resistance, RT determinant, resistant drugs listed
    ph = IQ.known_phenotype("K103N", gene="RT", organism="HIV-1")
    assert ph["is_known_resistance_mutation"] is True
    assert "RT:K103N" in ph["determinant"] and "efavirenz" in ph["resistant_drugs"]
    assert "hiv_amr" in ph["catalog"]


def test_hiv_phenotype_non_drm_is_negative():
    # a benign RT polymorphism not in the major-DRM set -> known_resistance False (catalog present, no hit)
    ph = IQ.known_phenotype("K102Q", gene="RT", organism="HIV-1")
    assert ph["is_known_resistance_mutation"] is False and ph["resistant_drugs"] is None


def test_unsupported_organism_honest_fallback():
    # E. coli has NO in-repo mutation-level catalog -> None + an honest note (never fabricated)
    ph = IQ.known_phenotype("S83L", gene="gyrA", organism="Escherichia coli")
    assert ph["catalog"] is None and ph["is_known_resistance_mutation"] is None
    assert "AMRFinder" in ph["caveat"] and "FULL GENOME" in ph["caveat"]


def test_hiv_gene_drug_routing():
    assert set(IQ._hiv_drugs_for_gene("PR")) and all(
        __import__("dna_decode.data.hiv_amr", fromlist=["x"]).gene_for_hiv_drug(d) == "PR"
        for d in IQ._hiv_drugs_for_gene("PR"))
    assert IQ._hiv_drugs_for_gene("RT")           # RT has NNRTI + NRTI drugs


def test_load_logp_both_formats(tmp_path):
    import json
    raw = tmp_path / "raw.json"; raw.write_text(json.dumps({"1": {"A": -1.0}, "2": {"A": -2.0}}))
    wrapped = tmp_path / "w.json"; wrapped.write_text(json.dumps({"sequence": "AA", "logp": {"1": {"A": -1.0}}}))
    assert IQ.load_logp(raw) == {1: {"A": -1.0}, 2: {"A": -2.0}}
    assert IQ.load_logp(wrapped) == {1: {"A": -1.0}}          # unwraps the predictor format


def test_integrated_query_fuses_both_rungs():
    # tiny synthetic logp for a 3-residue "RT-like" seq; K103N is out of range so use position 2
    seq = "MKA"
    logp = {2: {a: -3.0 for a in IQ.P.AA}}
    logp[2]["K"] = -0.2; logp[2]["N"] = -1.0
    r = IQ.integrated_query("K2N", gene="RT", organism="HIV-1", sequence=seq, logp=logp)
    assert r["schema"] == "integrated-mutation-query-v1"
    assert "damage_llr" in r["molecular_effect_rung2"]
    assert "known_phenotype_rung3" in r and "honest_framing" in r
    # K2N isn't K103N so it's not a catalogued DRM -> phenotype negative, but molecular rung still present
    assert r["known_phenotype_rung3"]["is_known_resistance_mutation"] is False
