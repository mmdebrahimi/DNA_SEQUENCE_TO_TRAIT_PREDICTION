"""Compose forward+fba: point-mutation -> cell-trait. Pure decision logic + real-model smoke (slow).

The slow tests use a SYNTHETIC protein sequence (forward scores whatever seq is supplied; the KO uses
the real iML1515 model gene) so no network is needed. Damaging/conservative BLOSUM62 calls are pinned.
"""
from __future__ import annotations

import pytest

from dna_decode.fba.compose import decide_fba_action

# controlled residues: M1 D2 D3 D4 I5 I6 I7 K8 K9 K10 W11 W12 W13
_SEQ = "MDDDIIIKKKWWW"


# ---- pure: forward.predicted_effect -> FBA action ----

def test_decide_action_mapping():
    assert decide_fba_action("damaging") == "knockout"
    assert decide_fba_action("preserved") == "wildtype"
    assert decide_fba_action("uncertain") == "conditional"
    assert decide_fba_action("abstain") == "conditional"
    assert decide_fba_action("anything-unknown") == "conditional"  # never force a call


# ---- real-model smoke (needs cobra + iML1515) ----

@pytest.fixture(scope="module")
def model():
    pytest.importorskip("cobra")
    from dna_decode.fba.model import load_model
    return load_model()


@pytest.mark.slow
def test_damaging_missense_at_essential_gene_is_nonviable(model):
    # b0720 (gltA) is essential on glucose M9; D2W is BLOSUM62 -4 -> damaging -> LOF
    from dna_decode.fba.compose import variant_to_cell_trait
    rec = variant_to_cell_trait(model, "b0720", _SEQ, "D2W")
    assert rec["forward"]["predicted_effect"] == "damaging"
    assert rec["lof_call"] == "LOF"
    assert rec["fba_action"] == "knockout"
    assert rec["essential_if_lost"] is True
    assert "NON-VIABLE" in rec["cell_trait"]


@pytest.mark.slow
def test_conservative_missense_is_tolerated_wildtype(model):
    # I5L is BLOSUM62 +2 -> preserved -> tolerated -> no metabolic change
    from dna_decode.fba.compose import variant_to_cell_trait
    rec = variant_to_cell_trait(model, "b0720", _SEQ, "I5L")
    assert rec["forward"]["predicted_effect"] == "preserved"
    assert rec["lof_call"] == "TOLERATED"
    assert rec["fba_action"] == "wildtype"
    assert "viable" in rec["cell_trait"]


@pytest.mark.slow
def test_damaging_missense_at_nonessential_gene_is_viable_altered(model):
    # b4025 (pgi) is NON-essential on glucose M9; a damaging missense -> LOF -> viable, altered flux
    from dna_decode.fba.compose import variant_to_cell_trait
    rec = variant_to_cell_trait(model, "b4025", _SEQ, "D2W")
    assert rec["fba_action"] == "knockout"
    assert rec["essential_if_lost"] is False
    assert "viable" in rec["cell_trait"]


@pytest.mark.slow
def test_uncertain_missense_reports_conditional_both_ways(model):
    # D2T is BLOSUM62 -1 -> uncertain -> DO NOT force; report both ways
    from dna_decode.fba.compose import variant_to_cell_trait
    rec = variant_to_cell_trait(model, "b0720", _SEQ, "D2T")
    assert rec["forward"]["predicted_effect"] == "uncertain"
    assert rec["fba_action"] == "conditional"
    assert "cell_trait_if_LOF" in rec and "cell_trait_if_tolerated" in rec
    assert "cell_trait" not in rec  # no forced single call
