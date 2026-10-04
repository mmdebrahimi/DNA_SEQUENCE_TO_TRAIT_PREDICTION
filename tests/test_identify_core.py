"""Tests for the identification decision core (Step 3). Fully offline -- no Docker, no D:, no sketch.

The decision rules are where a router becomes dangerous or safe, so the emphasis is on the ABSTAIN
paths and on the two traps this repo has on record: resolving a tie by sort order, and letting a
broken invocation read as a confident biological answer.
"""
from __future__ import annotations

import pytest

from dna_decode.identify.core import (
    AbstainReason,
    Hit,
    IdentifyCall,
    Thresholds,
    accession_from_reference_id,
    decide,
    parse_mash_dist,
)

TH = Thresholds(max_distance=0.05, ambiguity_margin=0.005)

# accession -> organism, standing in for the reference manifest
LABELS = {
    "GCA_ec1": "escherichia_coli",
    "GCA_ec2": "escherichia_coli",
    "GCA_kp1": "klebsiella_pneumoniae",
    "GCA_cj1": "campylobacter",
    "GCA_unlabelled": None,
}
TOKENS = {
    "escherichia_coli": "Escherichia",
    "klebsiella_pneumoniae": "Klebsiella_pneumoniae",
    "campylobacter": "Campylobacter",
}


def _org(acc): return LABELS.get(acc)
def _tok(org): return TOKENS.get(org)
def _ref(acc): return f"/d/dna_decode_cache/refseq/{acc}/genome.fna"
def _hit(acc, dist): return Hit(_ref(acc), "q.fna", dist, 0.0, "900/1000")
def _decide(hits, th=TH): return decide(hits, th, _org, _tok)


# --------------------------------------------------------------------------- parsing

def test_parse_reads_the_verified_five_field_format_and_sorts_nearest_first():
    out = "\n".join([
        f"{_ref('GCA_kp1')}\tq.fna\t0.0400\t0.0\t500/1000",
        f"{_ref('GCA_ec1')}\tq.fna\t0.0010\t0.0\t980/1000",
    ])
    hits = parse_mash_dist(out)
    assert [h.distance for h in hits] == [0.0010, 0.0400]
    assert hits[0].shared_hashes == "980/1000"


@pytest.mark.parametrize("junk", ["", "   ", "not\ttab\tseparated", "a\tb\tnotafloat\t0\t1/2",
                                  "only\tthree\tfields"])
def test_malformed_rows_are_skipped_not_raised(junk):
    assert parse_mash_dist(junk) == []


def test_crlf_rows_parse():
    """The sketch inputs were bitten by CRLF once already; the parser must not inherit that."""
    row = f"{_ref('GCA_ec1')}\tq.fna\t0.001\t0.0\t980/1000\r\n"
    assert len(parse_mash_dist(row)) == 1


@pytest.mark.parametrize("ref,expect", [
    ("/d/dna_decode_cache/refseq/GCA_000023665.1/genome.fna", "GCA_000023665.1"),
    ("D:\\dna_decode_cache\\refseq\\GCA_1.1\\genome.fna", "GCA_1.1"),
    ("GCA_1.1/genome.fna", "GCA_1.1"),
    ("genome.fna", ""),
    ("", ""),
])
def test_accession_comes_from_the_parent_directory(ref, expect):
    assert accession_from_reference_id(ref) == expect


# --------------------------------------------------------------------------- abstain paths

def test_no_hits_abstains_NO_HIT():
    c = _decide([])
    assert c.abstained and c.reason is AbstainReason.NO_HIT and c.organism is None


def test_a_garbage_stream_abstains_rather_than_raising():
    """A broken container must not produce an exception that aborts a cohort run, nor a confident
    call. It lands on NO_HIT -- which is why the RUNNER, not this function, is responsible for
    distinguishing a failed invocation (REFERENCE_UNAVAILABLE) from a genuine no-match."""
    c = _decide(parse_mash_dist("ERROR: could not open reference\n"))
    assert c.abstained and c.reason is AbstainReason.NO_HIT


def test_hits_that_map_to_nothing_abstain_UNLABELLED_REFERENCE():
    c = _decide([_hit("GCA_unlabelled", 0.001)])
    assert c.abstained and c.reason is AbstainReason.UNLABELLED_REFERENCE
    assert c.nearest_distance == 0.001


def test_a_distant_nearest_neighbour_abstains_ABOVE_MAX_DISTANCE():
    """The out-of-set case: a genome whose closest reference is far away is not evidence of anything.
    This is the rule that makes the router closed-set rather than a nearest-label guesser."""
    c = _decide([_hit("GCA_ec1", 0.30)])
    assert c.abstained and c.reason is AbstainReason.ABOVE_MAX_DISTANCE


def test_two_organisms_within_the_margin_abstain_AMBIGUOUS_TOP2():
    c = _decide([_hit("GCA_ec1", 0.0100), _hit("GCA_kp1", 0.0120)])
    assert c.abstained and c.reason is AbstainReason.AMBIGUOUS_TOP2
    assert c.runner_up_organism == "klebsiella_pneumoniae"


def test_an_exact_tie_between_two_organisms_ABSTAINS_and_is_not_resolved_by_order():
    """THE DOCUMENTED TIE TRAP. Equal distances must not be settled by sort position -- the same
    inputs in the other order must give the same answer."""
    a = _decide([_hit("GCA_ec1", 0.01), _hit("GCA_kp1", 0.01)])
    b = _decide([_hit("GCA_kp1", 0.01), _hit("GCA_ec1", 0.01)])
    assert a.abstained and b.abstained
    assert a.reason is b.reason is AbstainReason.AMBIGUOUS_TOP2


# --------------------------------------------------------------------------- call paths

def test_a_clear_nearest_organism_is_called_with_its_routing_token():
    c = _decide([_hit("GCA_ec1", 0.0005), _hit("GCA_kp1", 0.0400)])
    assert not c.abstained
    assert c.organism == "escherichia_coli" and c.routing_token == "Escherichia"
    assert c.nearest_distance == 0.0005


def test_two_reference_genomes_of_the_SAME_organism_are_not_ambiguity():
    """Load-bearing on this corpus. E. coli is 68% of the reference, so the top two ROWS are almost
    always both E. coli; comparing raw rows instead of best-per-organism would abstain constantly."""
    c = _decide([_hit("GCA_ec1", 0.0010), _hit("GCA_ec2", 0.0011), _hit("GCA_kp1", 0.0400)])
    assert not c.abstained and c.organism == "escherichia_coli"
    assert c.runner_up_organism == "klebsiella_pneumoniae", (
        "the runner-up must be the next DIFFERENT organism, not the second row")


def test_a_single_organism_reference_hit_still_calls():
    c = _decide([_hit("GCA_cj1", 0.002)])
    assert not c.abstained and c.organism == "campylobacter"


def test_unlabelled_rows_are_skipped_without_blocking_a_labelled_hit_behind_them():
    c = _decide([_hit("GCA_unlabelled", 0.0001), _hit("GCA_ec1", 0.0009)])
    assert not c.abstained and c.organism == "escherichia_coli"


# --------------------------------------------------------------------------- invariants

def test_a_call_can_never_be_both_an_organism_and_an_abstention():
    with pytest.raises(ValueError):
        IdentifyCall("escherichia_coli", "Escherichia", True, AbstainReason.NO_HIT)
    with pytest.raises(ValueError):
        IdentifyCall(None, None, False)
    with pytest.raises(ValueError):
        IdentifyCall(None, None, True, None)


@pytest.mark.parametrize("bad", [0.0, -0.1, 1.5])
def test_nonsensical_thresholds_are_refused(bad):
    with pytest.raises(ValueError):
        Thresholds(max_distance=bad, ambiguity_margin=0.001)


def test_negative_ambiguity_margin_is_refused():
    with pytest.raises(ValueError):
        Thresholds(max_distance=0.05, ambiguity_margin=-0.001)


def test_as_dict_states_what_the_cell_does_NOT_support():
    """The closed-set scope travels with every record. A consumer must not read an abstention as
    'unidentifiable organism' nor a call as open-world taxonomy."""
    d = _decide([_hit("GCA_ec1", 0.0005)]).as_dict()
    assert d["closed_set"] is True
    joined = " ".join(d["does_not_support"]).lower()
    assert "open-world" in joined and "subspecies" in joined and "mixtures" in joined
    assert d["abstain_reason"] is None and d["routing_token"] == "Escherichia"


def test_thresholds_are_a_PARAMETER_so_behaviour_moves_with_them():
    """Non-vacuity for the seam: the same hits must decide differently under different thresholds,
    proving no value is baked into the module."""
    hits = [_hit("GCA_ec1", 0.02)]
    assert not _decide(hits, Thresholds(0.05, 0.005)).abstained
    assert _decide(hits, Thresholds(0.01, 0.005)).abstained
