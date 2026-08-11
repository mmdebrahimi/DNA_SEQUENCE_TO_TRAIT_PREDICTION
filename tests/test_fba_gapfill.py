"""Model gap-finding + repair — pure logic (wheel-only) + the real measured sucrose case."""
from __future__ import annotations

import pytest

from dna_decode.fba.gapfill import DeadEnd, find_dead_ends


# ---- pure: the structural diagnostic ----

def test_metabolite_that_is_made_but_never_used_is_a_dead_end():
    # R1 makes X and nothing consumes it -> X cannot carry steady-state flux, so R1 is dead weight.
    des = find_dead_ends({"R1": {"A": -1, "X": 1}})
    assert DeadEnd("X", "no_consumer", ("R1",)) in des


def test_metabolite_that_is_used_but_never_made_is_a_dead_end():
    des = find_dead_ends({"R1": {"Y": -1, "B": 1}})
    kinds = {d.metabolite: d.kind for d in des}
    assert kinds["Y"] == "no_producer"


def test_a_metabolite_with_both_sides_is_not_a_dead_end():
    des = find_dead_ends({"R1": {"A": -1, "X": 1}, "R2": {"X": -1, "B": 1}})
    assert "X" not in {d.metabolite for d in des}


def test_a_REVERSIBLE_reaction_counts_as_both_producer_and_consumer():
    """The correctness crux. A reversible reaction can run either way, so scoring it one-directionally
    invents dead ends that are not dead -- and this diagnostic is only useful if its hits are real."""
    des = find_dead_ends({"R1": {"A": -1, "X": 1}, "R2": {"X": -1, "B": 1}}, reversible={"R2"})
    assert "X" not in {d.metabolite for d in des}
    # ... and with R2 irreversible in the *other* direction, X is still fine (R2 consumes it)
    assert "X" not in {d.metabolite for d in find_dead_ends({"R1": {"A": -1, "X": 1}, "R2": {"X": -1}})}


def test_reversible_producer_alone_still_leaves_the_partner_dead():
    # a lone reversible reaction makes A both produced and consumed, so nothing is a dead end
    assert find_dead_ends({"R1": {"A": -1, "X": 1}}, reversible={"R1"}) == []


def test_zero_coefficients_are_ignored():
    des = find_dead_ends({"R1": {"A": -1, "X": 1, "Z": 0.0}})
    assert "Z" not in {d.metabolite for d in des}


def test_empty_model_has_no_dead_ends():
    assert find_dead_ends({}) == []


def test_dead_end_is_serialisable():
    d = DeadEnd("suc6p_c", "no_consumer", ("SUCptspp",)).as_dict()
    assert d == {"metabolite": "suc6p_c", "kind": "no_consumer", "reactions": ["SUCptspp"]}


# ---- real model: the measured sucrose gap, found and repaired ----

@pytest.mark.slow
def test_finds_and_repairs_the_measured_sucrose_gap():
    """The end-to-end claim, on the real models.

    iML1515 predicts NO growth on sucrose, but BW25113 has a sucrose carbon-source experiment in the
    Wetmore/Keio RB-TnSeq set (that assay only runs sources the organism grows on), so this is a measured
    FALSE NEGATIVE. Three things must hold, and the first is the interesting one: the structural diagnostic
    has to surface `suc6p_c` WITHOUT being told anything about sucrose.
    """
    pytest.importorskip("cobra")
    from dna_decode.fba.gapfill import orphan_uptake_targets, propose_repair, verify_repair
    from dna_decode.fba.model import load_model

    m = load_model()
    donor = load_model(model_id="iYS1720")  # Salmonella -- can catabolise sucrose

    # 1. label-free: the dead end is found by structure alone
    orphans = {d.metabolite for d in orphan_uptake_targets(m)}
    assert "suc6p_c" in orphans, "the structural diagnostic must surface the gap unprompted"

    # 2. the proposal
    p = propose_repair(m, donor, "EX_sucr_e")
    assert p["status"] == "ok" and p["growth_before"] == 0.0
    proposed = [c["id"] for c in p["candidates"][0]]
    assert "FFSD" in proposed

    # 3. the repair is measured, and does NOT leak into unrelated carbon sources
    v = verify_repair(m, donor, proposed, "EX_sucr_e",
                      specificity_exchanges=("EX_glc__D_e", "EX_xyl__D_e", "EX_cellb_e"))
    assert v["repaired"] is True
    assert v["growth_after"] > 1.0 > v["growth_before"]
    assert v["specificity_unchanged"], f"repair leaked onto {v['specificity_changed_exchanges']}"
