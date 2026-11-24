"""Tests for the extracted uniform-admissibility gate.

The load-bearing test is `test_the_committed_refusal_is_reproduced`: the whole justification for
this module is that it is an EXTRACTION of `scripts/glm_beyond_kmer.py:228-245`, so if it does not
reproduce the committed 2026-11-22 refusal the extraction changed behaviour and nothing downstream
can be trusted. Everything else exists to stop the predicate being accidentally constant.
"""
from __future__ import annotations

import math

import pytest

from dna_decode.glm.eval_admissibility import (
    ADMISSIBLE,
    REFUSED_BOTH_BELOW_UNIFORM,
    REFUSED_MODEL_BELOW_UNIFORM,
    REFUSED_NULL_BELOW_UNIFORM,
    STATES,
    UNIFORM_CE_NATS_DNA,
    admissibility,
)

# The committed numbers, from wiki/glm_beyond_kmer_2026-11-22.json, primary cell
# left_only|markov_ecoli plus the left_only|uniform cell's paired delta.
COMMITTED_NULL_CE = 1.4291
COMMITTED_MODEL_CE = 1.4098
COMMITTED_UNIFORM_DELTA = -0.0235


def test_uniform_constant_is_ln4():
    assert UNIFORM_CE_NATS_DNA == math.log(4)
    assert round(UNIFORM_CE_NATS_DNA, 4) == 1.3863


# ---- the extraction pin: the committed refusal ------------------------------------------------

def test_the_committed_refusal_is_reproduced():
    """Both predictors scored WORSE than uniform, so the set cannot carry a verdict."""
    got = admissibility(null_ce_nats=COMMITTED_NULL_CE, model_ce_nats=COMMITTED_MODEL_CE)
    assert not got.admissible
    assert got.state == REFUSED_BOTH_BELOW_UNIFORM
    # Both margins negative is the specific condition behind EVAL_SET_UNINFORMATIVE_BELOW_UNIFORM's
    # "neither beats a random guess" wording.
    assert got.null_margin_nats < 0
    assert got.model_margin_nats < 0


def test_the_committed_refusal_is_identical_via_the_paired_delta_form():
    """The scorer tested the ROUNDED paired delta; that path must reach the same state."""
    got = admissibility(null_ce_nats=COMMITTED_NULL_CE,
                        model_uniform_delta_nats=COMMITTED_UNIFORM_DELTA)
    assert got.state == REFUSED_BOTH_BELOW_UNIFORM
    assert got.model_margin_nats == COMMITTED_UNIFORM_DELTA


def test_the_two_model_forms_agree_on_the_committed_numbers():
    """Guards the equivalence the refactor relies on: uniform CE is exactly ln4 per window, so the
    paired delta against uniform IS ln4 - model_ce up to the scorer's 4-dp rounding."""
    assert COMMITTED_UNIFORM_DELTA == pytest.approx(
        UNIFORM_CE_NATS_DNA - COMMITTED_MODEL_CE, abs=5e-5)


# ---- non-vacuity: the predicate is not constant ------------------------------------------------

def test_a_predictable_set_admits():
    got = admissibility(null_ce_nats=1.30, model_ce_nats=1.25)
    assert got.admissible
    assert got.state == ADMISSIBLE
    assert got.null_margin_nats > 0 and got.model_margin_nats > 0


def test_every_state_is_reachable():
    """A state nothing can emit is the vacuous-guard shape this repo has paid for."""
    reached = {
        admissibility(1.30, 1.25).state,
        admissibility(1.45, 1.25).state,
        admissibility(1.30, 1.45).state,
        admissibility(1.45, 1.45).state,
    }
    assert reached == set(STATES)


def test_only_the_null_failing_is_its_own_state():
    got = admissibility(null_ce_nats=1.45, model_ce_nats=1.25)
    assert got.state == REFUSED_NULL_BELOW_UNIFORM


def test_only_the_model_failing_is_its_own_state():
    got = admissibility(null_ce_nats=1.30, model_ce_nats=1.45)
    assert got.state == REFUSED_MODEL_BELOW_UNIFORM


# ---- the null-only pre-scoring screen ----------------------------------------------------------

def test_null_only_admits_when_the_null_beats_uniform():
    got = admissibility(null_ce_nats=1.30)
    assert got.admissible
    assert got.model_assessed is False
    assert got.model_ce_nats is None and got.model_margin_nats is None


def test_null_only_refuses_when_the_null_loses_to_uniform():
    got = admissibility(null_ce_nats=1.45)
    assert not got.admissible
    assert got.state == REFUSED_NULL_BELOW_UNIFORM
    assert got.model_assessed is False


def test_null_only_can_never_emit_a_model_state():
    """An unassessed model must not disqualify a set -- the screen is weaker on purpose."""
    for null_ce in (1.0, 1.30, UNIFORM_CE_NATS_DNA, 1.45, 2.0):
        assert admissibility(null_ce_nats=null_ce).state in (
            ADMISSIBLE, REFUSED_NULL_BELOW_UNIFORM)


def test_an_admitting_null_only_screen_does_not_survive_a_failing_model():
    """The cheap screen admitting is NOT a promise the full gate will admit."""
    assert admissibility(1.30).admissible
    assert not admissibility(1.30, model_ce_nats=1.45).admissible


# ---- boundary: exactly at uniform ---------------------------------------------------------------

def test_a_null_exactly_at_uniform_refuses():
    """Strict `<` on the margin. Pinned so a future float change cannot silently flip it."""
    got = admissibility(null_ce_nats=UNIFORM_CE_NATS_DNA)
    assert not got.admissible
    assert got.state == REFUSED_NULL_BELOW_UNIFORM
    assert got.null_margin_nats == 0.0


def test_a_model_exactly_at_uniform_refuses():
    got = admissibility(null_ce_nats=1.30, model_ce_nats=UNIFORM_CE_NATS_DNA)
    assert got.state == REFUSED_MODEL_BELOW_UNIFORM


def test_a_zero_paired_delta_refuses():
    got = admissibility(null_ce_nats=1.30, model_uniform_delta_nats=0.0)
    assert got.state == REFUSED_MODEL_BELOW_UNIFORM


# ---- the alphabet-size keyword -------------------------------------------------------------------

def test_a_different_alphabet_size_is_expressible_without_editing_the_module():
    """A 20-letter alphabet: ln(20)=2.9957. A null at 2.5 beats it; the DNA constant would not."""
    got = admissibility(null_ce_nats=2.5, uniform=math.log(20))
    assert got.admissible
    assert not admissibility(null_ce_nats=2.5).admissible


# ---- the module says nothing about any model -----------------------------------------------------

def test_the_result_carries_no_verdict_about_a_model():
    """Admissibility is necessary, not sufficient; a confound verdict must not live here."""
    d = admissibility(1.30, 1.25).to_dict()
    blob = " ".join(str(v) for v in d.values()).upper()
    for forbidden in ("BEYOND_KMER", "CONFOUND", "ORTHOLOG", "LEAK"):
        assert forbidden not in blob


def test_to_dict_is_json_serialisable():
    import json
    json.loads(json.dumps(admissibility(1.30, 1.25).to_dict()))
    json.loads(json.dumps(admissibility(1.45).to_dict()))


# ---- BOTH model forms at once: the shape the scorer actually calls ------------------------------
# `glm_beyond_kmer.py` passes model_ce_nats AND model_uniform_delta_nats together, because it wants
# the cross-entropy recorded in the artifact while the DECISION is taken on the rounded delta. No
# test exercised that call shape, so the documented precedence was unenforced: if precedence
# silently flipped to the cross-entropy, the extraction would stop being bit-identical to the
# inline original exactly in the near-boundary cases the delta form exists to handle.

def test_the_delta_takes_precedence_over_the_cross_entropy_when_both_are_given():
    """Documented precedence: the delta decides, the cross-entropy is only carried through."""
    got = admissibility(null_ce_nats=1.30, model_ce_nats=1.45, model_uniform_delta_nats=0.01)
    # 1.45 is WORSE than uniform, so a cross-entropy-led decision would refuse here.
    assert got.state == ADMISSIBLE
    assert got.model_margin_nats == 0.01
    # ...and the cross-entropy is still reported, unrounded, for the artifact.
    assert got.model_ce_nats == 1.45
    assert got.model_assessed is True


def test_the_reverse_disagreement_refuses_on_the_delta_too():
    """The mirror case, so the test above cannot pass by the predicate being constant."""
    got = admissibility(null_ce_nats=1.30, model_ce_nats=1.25, model_uniform_delta_nats=-0.01)
    assert got.state == REFUSED_MODEL_BELOW_UNIFORM
    assert got.model_ce_nats == 1.25 and got.model_margin_nats == -0.01


def test_a_negative_delta_with_a_failing_null_is_still_the_both_state():
    """The committed call shape at the committed numbers reaches REFUSED_BOTH through the delta."""
    got = admissibility(null_ce_nats=COMMITTED_NULL_CE, model_ce_nats=COMMITTED_MODEL_CE,
                        model_uniform_delta_nats=COMMITTED_UNIFORM_DELTA)
    assert got.state == REFUSED_BOTH_BELOW_UNIFORM
