"""Offline tests for the deployable forward runner + CLI wiring.

The heavy computes (_esm2_logp_table / _prosst_table_from_structure / _gemme_table) are monkeypatched to
cheap stand-ins, so NO model download / no torch forward pass / no Docker runs here. This pins:
  - method resolution + input-aware `auto`
  - honest graceful degradation (provenance: method_requested/method_used/degraded/degrade_reason)
  - --no-degrade hard error
  - CLI --capabilities + strong-method wiring via precomputed tables
"""
from __future__ import annotations

import json

from dna_decode.forward import deploy
from dna_decode.forward.deploy import predict_effect_deployable
import dna_decode.forward.cli as fcli

SEQ = "MKVLAAGGWY"                      # 10-aa toy protein
ALL = {"torch": True, "transformers": True, "prosst": True, "docker": True}
NONE = {"torch": False, "transformers": False, "prosst": False, "docker": False}


def _cheap_esm_table(seq):
    # pos-table shape {pos:{aa:logp}} covering every position/residue
    aa = "ACDEFGHIKLMNPQRSTVWY"
    return {i + 1: {a: (0.5 if a == seq[i] else -0.5) for a in aa} for i in range(len(seq))}


def _cheap_variant_table(seq, offset=0.0):
    aa = "ACDEFGHIKLMNPQRSTVWY"
    return {f"{seq[i]}{i+1}{a}": (hash((i, a)) % 100) / 100.0 + offset
            for i in range(len(seq)) for a in aa if a != seq[i]}


def test_esm2_with_supplied_table_no_compute(monkeypatch):
    # supplying esm_table must NOT trigger the heavy compute
    monkeypatch.setattr(deploy, "_esm2_logp_table", lambda *a, **k: (_ for _ in ()).throw(AssertionError("computed!")))
    d = predict_effect_deployable(SEQ, "M1L", method="esm2", esm_table=_cheap_esm_table(SEQ), caps=ALL)
    assert d["method_used"] == "esm2"
    assert d["degraded"] is False and d["method_requested"] == "esm2"


def test_auto_picks_esm2_when_only_seq_and_esm_deps(monkeypatch):
    called = {}
    monkeypatch.setattr(deploy, "_esm2_logp_table", lambda seq, **k: called.setdefault("y", _cheap_esm_table(seq)))
    # ESM deps present, no structure, no msa -> auto must resolve to esm2 (NOT hybrid/prosst/gemme)
    d = predict_effect_deployable(SEQ, "M1L", method="auto",
                                  caps={"torch": True, "transformers": True, "prosst": True, "docker": True})
    assert d["method_used"] == "esm2"
    assert "y" in called


def test_auto_degrades_to_blosum_when_no_deps():
    d = predict_effect_deployable(SEQ, "M1L", method="auto", caps=NONE)
    assert d["method_used"] == "blosum62"


def test_requested_esm2_degrades_when_no_torch():
    d = predict_effect_deployable(SEQ, "M1L", method="esm2", caps=NONE, degrade=True)
    assert d["method_used"] == "blosum62"
    assert d["degraded"] is True
    assert "esm2" in d["degrade_reason"] and "blosum62" in d["degrade_reason"]


def test_no_degrade_raises():
    try:
        predict_effect_deployable(SEQ, "M1L", method="esm2", caps=NONE, degrade=False)
    except RuntimeError as e:
        assert "not runnable" in str(e)
        return
    raise AssertionError("expected RuntimeError under --no-degrade")


def test_prosst_needs_structure_input_even_with_deps():
    # deps present but NO structure supplied -> prosst not computable -> degrade
    d = predict_effect_deployable(SEQ, "M1L", method="prosst", caps=ALL, degrade=True)
    assert d["degraded"] is True
    assert "structure" in d["degrade_reason"]


def test_hybrid_with_two_supplied_tables(monkeypatch):
    monkeypatch.setattr(deploy, "_esm2_logp_table", lambda *a, **k: (_ for _ in ()).throw(AssertionError("computed!")))
    d = predict_effect_deployable(SEQ, "M1L", method="hybrid", esm_table=_cheap_esm_table(SEQ),
                                  prosst_table=_cheap_variant_table(SEQ), caps=ALL)
    assert d["method_used"] == "hybrid"
    assert d["degraded"] is False


def test_gemme_computes_from_msa(monkeypatch):
    monkeypatch.setattr(deploy, "_gemme_table", lambda seq, muts, msa: _cheap_variant_table(seq))
    d = predict_effect_deployable(SEQ, "M1L", method="gemme", msa="fake.a3m", caps=ALL)
    assert d["method_used"] == "gemme"


# ---- CLI wiring ----

def test_cli_capabilities_exits_0(capsys):
    rc = fcli.main(["--capabilities"])
    assert rc == 0
    out = capsys.readouterr().out
    assert "capability preflight" in out


def test_cli_blosum_unchanged_path(capsys):
    rc = fcli.main(["--mutation", "M1L", "--protein-seq", SEQ, "--method", "blosum62"])
    assert rc == 0
    out = capsys.readouterr().out
    assert "forward variant-effect decode" in out and "blosum62" in out


def test_cli_esm2_via_supplied_table(tmp_path, capsys):
    tbl = tmp_path / "esm.json"
    tbl.write_text(json.dumps({str(k): v for k, v in _cheap_esm_table(SEQ).items()}))
    rc = fcli.main(["--mutation", "M1L", "--protein-seq", SEQ, "--method", "esm2",
                    "--esm-table", str(tbl), "--json"])
    assert rc == 0
    d = json.loads(capsys.readouterr().out)
    assert d["method_used"] == "esm2" and d["degraded"] is False


def test_cli_mutation_required_without_capabilities():
    try:
        fcli.main(["--protein-seq", SEQ])
    except SystemExit as e:
        assert e.code == 2
        return
    raise AssertionError("expected SystemExit when --mutation missing")


def test_cli_hybrid_via_supplied_tables(tmp_path, capsys):
    e = tmp_path / "e.json"; p = tmp_path / "p.json"
    e.write_text(json.dumps({str(k): v for k, v in _cheap_esm_table(SEQ).items()}))
    p.write_text(json.dumps(_cheap_variant_table(SEQ)))
    rc = fcli.main(["--mutation", "M1L", "--protein-seq", SEQ, "--method", "hybrid",
                    "--esm-table", str(e), "--prosst-table", str(p), "--json"])
    assert rc == 0
    d = json.loads(capsys.readouterr().out)
    assert d["method_used"] == "hybrid"
