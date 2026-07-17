"""Tests for the shipped `dna-decode inverse` cell.

Pins the honesty rails (it RANKS, never doses), the coordinate gate, and the tie/diversity behaviour that
the real blaTEM data forced.
"""
from __future__ import annotations

import pytest

from dna_decode.forward.inverse import (
    UNSUPPORTED_CLAIMS,
    InverseError,
    enumerate_candidates,
    propose_edits,
)

TEM1_HEAD = "MSIQHFRVALIPFFAAFCLPVFA"


# ---- the candidate space ---------------------------------------------------------------------------------

def test_without_a_cds_every_substitution_is_a_candidate():
    c, space = enumerate_candidates("MK")
    assert space.endswith("no_CDS")
    assert len(c) == 2 * 19                       # 2 residues x 19 alternatives


def test_with_a_cds_only_single_nt_reachable_edits_are_candidates():
    c, space = enumerate_candidates("M", "ATG")
    assert space == "single_nt_accessible_from_CDS"
    alts = {m[-1] for m, *_ in c}
    assert "L" in alts and "V" in alts
    assert "W" not in alts                        # 2 nt away from ATG


def test_a_cds_that_does_not_translate_to_the_protein_fails_loudly():
    """Coordinate gate: never propose an edit in the wrong frame."""
    with pytest.raises(InverseError, match="does not translate"):
        enumerate_candidates("M", "GGG")          # GGG is Gly, not Met


def test_a_short_cds_fails_loudly():
    with pytest.raises(InverseError, match="too short"):
        enumerate_candidates("MK", "ATG")


def test_nonsense_is_never_proposed():
    c, _ = enumerate_candidates("Y", "TAT")       # TAT -> TAA/TAG are stops, 1 nt away
    assert all(m[-1] != "*" for m, *_ in c)


def test_a_non_protein_sequence_is_refused():
    with pytest.raises(InverseError, match="standard amino acids"):
        enumerate_candidates("MKZ1")


# ---- the proposal ----------------------------------------------------------------------------------------

def test_percentile_zero_proposes_more_damaging_edits_than_percentile_one():
    lo = propose_edits(TEM1_HEAD, 0.0, top_k=3)
    hi = propose_edits(TEM1_HEAD, 1.0, top_k=3)
    assert max(p.score for p in lo.proposals) <= min(p.score for p in hi.proposals)


def test_an_out_of_range_percentile_is_refused():
    for bad in (-0.1, 1.1):
        with pytest.raises(InverseError, match="must be in"):
            propose_edits(TEM1_HEAD, bad)


def test_esm2_without_a_table_is_refused_rather_than_silently_falling_back():
    """A silent fallback to blosum62 would misreport which oracle produced the proposal."""
    with pytest.raises(InverseError, match="esm_table"):
        propose_edits(TEM1_HEAD, 0.5, method="esm2")


def test_an_unknown_method_is_refused():
    with pytest.raises(InverseError, match="unknown method"):
        propose_edits(TEM1_HEAD, 0.5, method="magic")


# ---- diversity (forced by the real BLOSUM tie structure) -------------------------------------------------

def test_diverse_default_returns_one_edit_per_residue():
    r = propose_edits(TEM1_HEAD, 0.02, top_k=5)
    assert r.proposals
    assert len({p.pos for p in r.proposals}) == len(r.proposals)


def test_non_diverse_may_repeat_a_residue():
    """The plain window is available but is NOT the default: BLOSUM's 7 distinct scores over ~1.9k blaTEM
    candidates make it return k shots at the same residue."""
    r = propose_edits(TEM1_HEAD, 0.02, top_k=5, diverse=False)
    assert len({p.pos for p in r.proposals}) <= len(r.proposals)


def test_a_large_tie_group_is_surfaced_in_the_notes():
    r = propose_edits(TEM1_HEAD, 0.02, top_k=5)
    assert any("TIE" in n for n in r.notes), r.notes


# ---- the honesty rails -----------------------------------------------------------------------------------

def test_the_output_never_claims_a_magnitude():
    d = propose_edits(TEM1_HEAD, 0.5, top_k=3).as_dict()
    assert d["claim"].startswith("predicted-damage RANK")
    assert any("MAGNITUDE" in c for c in d["does_not_support"])
    assert any("CLINICAL" in c for c in d["does_not_support"])
    for p in d["proposals"]:
        assert "effect" not in p and "fold_change" not in p     # no dose field can be misread
        assert "score_percentile" in p


def test_every_proposal_ships_its_evidence_and_scope():
    d = propose_edits(TEM1_HEAD, 0.5).as_dict()
    assert d["regime"] == "B_molecular"
    assert d["research_use_only"] is True
    assert d["evidence"]["rank_inverse_beats_no_oracle_null"] == "4/4 usable proteins"
    assert len(d["does_not_support"]) == len(UNSUPPORTED_CLAIMS)


def test_best_of_k_guidance_ships_with_every_call():
    """Top-1 is ~4x worse than best-of-5; a user who assays only the first proposal must be told."""
    r = propose_edits(TEM1_HEAD, 0.5, top_k=5)
    assert any("assay all k" in n for n in r.notes)


# ---- CLI -------------------------------------------------------------------------------------------------

def test_cli_json_roundtrip(capsys):
    from dna_decode.forward.inverse_cli import main

    rc = main(["--protein-seq", TEM1_HEAD, "--target-percentile", "0.1", "--top-k", "3", "--json"])
    assert rc == 0
    import json
    d = json.loads(capsys.readouterr().out)
    assert len(d["proposals"]) == 3
    assert d["tool"] == "dna_decode.forward.inverse"


def test_cli_rejects_a_bad_percentile(capsys):
    from dna_decode.forward.inverse_cli import main

    assert main(["--protein-seq", TEM1_HEAD, "--target-percentile", "5"]) == 2


def test_cli_is_routable_from_the_unified_entrypoint(capsys):
    import dna_decode.cli as uni

    assert "inverse" in uni.TRAITS
    rc = uni.main(["inverse", "--protein-seq", TEM1_HEAD, "--target-percentile", "0.5", "--json"])
    assert rc == 0
