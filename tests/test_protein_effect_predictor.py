"""Tests for the rung-2 mutation-effect predictor (pure functions; ESM path needs no coverage here)."""
from __future__ import annotations

import pytest

from dna_decode.protein_effect import predictor as P


def test_parse_mutation_ok_and_errors():
    assert P.parse_mutation("A123V") == ("A", 123, "V")
    assert P.parse_mutation("a1v") == ("A", 1, "V")             # case-insensitive by design
    for bad in ("A123", "123V", "AV", "B10C", "A10Z"):          # malformed / non-standard residue
        with pytest.raises(P.MutationParseError):
            P.parse_mutation(bad)
    with pytest.raises(P.MutationParseError):
        P.parse_mutation("A10A")                                 # no-op


def test_apply_edit_deterministic():
    seq = "MAQVK"
    e = P.apply_edit(seq, "Q3E")
    assert e["wt"] == "Q" and e["pos"] == 3 and e["mut"] == "E"
    assert e["mut_seq"] == "MAEVK" and e["wt_seq"] == "MAQVK"


def test_apply_edit_rejects_mismatch_and_range():
    with pytest.raises(P.MutationParseError):
        P.apply_edit("MAQVK", "A3E")        # position 3 is Q, not A
    with pytest.raises(P.MutationParseError):
        P.apply_edit("MAQVK", "K9A")        # position out of range


def test_damage_llr_frozen_sign_contract():
    # logP(wt) - logP(mut): a mutation to a LOWER-prob residue is MORE damaging (positive, larger)
    logp = {5: {"A": -1.0, "V": -5.0, "G": -0.5}}
    assert P.damage_llr(logp, 5, "A", "V") == pytest.approx(4.0)    # A(-1) -> V(-5): damaging
    assert P.damage_llr(logp, 5, "A", "G") == pytest.approx(-0.5)   # A(-1) -> G(-0.5): tolerated (negative)


def test_position_percentile_and_direction_hint():
    # WT=A; V is the least likely (most damaging) of all substitutions -> percentile ~1.0 -> deleterious
    logp = {1: {a: -1.0 for a in P.AA}}
    logp[1]["A"] = -0.1        # WT most likely
    logp[1]["V"] = -9.0        # V very unlikely -> most damaging
    pct = P.position_percentile(logp, 1, "A", "V")
    assert pct == pytest.approx(1.0) and P.direction_hint(pct) == "likely-deleterious"
    # a mutation to a residue nearly as likely as WT -> low percentile -> tolerated
    logp[1]["G"] = -0.15
    assert P.direction_hint(P.position_percentile(logp, 1, "A", "G")) == "likely-tolerated"


def test_predict_output_schema_and_honesty():
    seq = "MAQVK"
    logp = {3: {a: -3.0 for a in P.AA}}
    logp[3]["Q"] = -0.2; logp[3]["E"] = -1.0
    out = P.predict(seq, "Q3E", logp)
    assert out["schema"] == "protein-mutation-effect-v1"
    assert out["sequence_change"]["certain"] is True and out["sequence_change"]["position"] == 3
    assert out["damage_llr_definition"].startswith("logP(wt) - logP(mut)")
    assert "NOT a per-mutation probability" in out["honest_caveat"]
    assert out["provenance"]["model"] == P.MODEL
    assert out["direction_hint"] in ("likely-deleterious", "uncertain", "likely-tolerated")
