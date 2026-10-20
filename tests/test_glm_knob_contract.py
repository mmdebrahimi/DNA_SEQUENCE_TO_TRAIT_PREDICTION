"""The pre-run contract on the GLM sampling knobs.

These exist because a sweep ran 21 minutes of compute over a knob that did nothing, and the
resulting invariance was reported as a finding. The guard that was missing is cheap and is
pinned here: a swept variable must be shown to move the output BEFORE the sweep is trusted.
"""
from __future__ import annotations

import pytest

torch = pytest.importorskip("torch")
import torch.nn.functional as F  # noqa: E402
from transformers import TopKLogitsWarper  # noqa: E402

from dna_decode.glm.knob_contract import (  # noqa: E402
    TOP_K_KMER_FLOOR,
    SamplingContractError,
    assert_knob_wired,
    requested_non_default,
    validate_sampling,
    verify_knob_live,
)


# ---------------------------------------------------------------- fail-fast validation
@pytest.mark.parametrize("t", [0.0, -1.0, -0.5])
def test_a_non_positive_temperature_is_refused(t):
    with pytest.raises(SamplingContractError, match="temperature must be > 0"):
        validate_sampling(t, 0)


def test_NaN_temperature_is_refused():
    with pytest.raises(SamplingContractError, match="NaN"):
        validate_sampling(float("nan"), 0)


@pytest.mark.parametrize("t", [0.7, 1.0, 1.5, 2.0])
def test_valid_temperatures_pass(t):
    validate_sampling(t, 0)


def test_top_k_below_the_floor_is_refused_and_the_message_explains_KMER_semantics():
    with pytest.raises(SamplingContractError) as e:
        validate_sampling(1.0, 4)
    msg = str(e.value)
    assert "below the measured floor" in msg
    # The whole point: a user typing 4 means "4 bases". The message must say otherwise.
    assert "not the 4 bases" in msg and "UNREACHABLE" in msg


@pytest.mark.parametrize("k", [0, TOP_K_KMER_FLOOR, TOP_K_KMER_FLOOR + 1, 512])
def test_disabled_or_at_or_above_the_floor_passes(k):
    validate_sampling(1.0, k)


# ---------------------------------------------------------------- the floor is DERIVED
def _zero_cells(scores, top_k, k=6, V=4096):
    """Count (position, base) cells whose marginal is exactly 0 -- i.e. unreachable bases.

    Replicates the remote _BPLogitsProcessor marginalisation: warp, softmax over the k-mer
    logits, then sum k-mer probability into the base each k-mer carries at each position.
    """
    idx = torch.zeros(k, V, dtype=torch.long)
    for km in range(V):
        rem = km
        for pos in range(k - 1, -1, -1):
            idx[pos, km] = rem % 4
            rem //= 4
    warped = TopKLogitsWarper(top_k)(torch.zeros((1, 1), dtype=torch.long), scores.clone())
    probs = F.softmax(warped[:, :V].float(), dim=-1)
    return sum(1 for pos in range(k) for nt in range(4)
               if float(probs[:, idx[pos] == nt].sum()) <= 1e-9)


def test_the_floor_is_the_SMALLEST_value_with_no_reachability_loss():
    """Pins TOP_K_KMER_FLOOR to a re-derived measurement, not to a remembered number.

    Lowering the constant without re-measuring fails here. This is the difference between a
    derived threshold and an asserted one -- the habit this whole module exists to fix.
    """
    conds = [(s, sharp) for s in (0, 1, 2, 3, 4) for sharp in (1.0, 2.0, 4.0)]

    def worst(top_k):
        out = 0
        for seed, sharp in conds:
            torch.manual_seed(seed)
            out = max(out, _zero_cells(torch.randn(1, 4096) * sharp, top_k))
        return out

    assert worst(TOP_K_KMER_FLOOR) == 0, "the floor itself loses reachability; re-derive it"
    assert worst(TOP_K_KMER_FLOOR - 1) > 0, (
        f"top_k={TOP_K_KMER_FLOOR - 1} shows no reachability loss, so the floor is set "
        f"HIGHER than measured -- lower it to the derived boundary rather than padding it")


def test_the_refused_region_really_does_lose_reachability():
    """Non-vacuity for the refusal: the value a user would naively type is genuinely broken."""
    torch.manual_seed(0)
    assert _zero_cells(torch.randn(1, 4096) * 2.0, 4) > 0


# ---------------------------------------------------------------- the inert-knob guard
def test_a_requested_knob_with_NO_warper_is_refused():
    """The exact 2026-10-19 regression: a stamped setting that installed nothing."""
    with pytest.raises(SamplingContractError, match="INERT"):
        assert_knob_wired(0.7, 0, [])


def test_the_default_knob_with_no_warper_is_FINE():
    """T=1.0 installs nothing by design; refusing that would break every default run."""
    assert_knob_wired(1.0, 0, [])


def test_a_requested_knob_WITH_a_warper_passes():
    assert_knob_wired(0.7, 0, [object()])


@pytest.mark.parametrize("t,k,expected", [(1.0, 0, False), (0.7, 0, True),
                                          (1.0, 50, True), (0.7, 50, True)])
def test_non_default_detection(t, k, expected):
    assert requested_non_default(t, k) is expected


# ---------------------------------------------------------------- the end-to-end check
def test_verify_knob_live_CATCHES_a_knob_that_ignores_its_argument():
    """Simulates last night's bug directly: a generator whose output ignores temperature.

    If this returned live=True the whole rail would be theatre, so this is the load-bearing
    test -- it is the guard being shown to catch the real failure, not merely to run.
    """
    def ignores_temperature(_t):
        return ["ACGTACGT", "TTTTGGGG"]

    v = verify_knob_live(ignores_temperature)
    assert v["live"] is False
    assert v["verdict"] == "TEMPERATURE_IS_IGNORED"
    assert v["identical"] == v["n"] == 2


def test_verify_knob_live_PASSES_a_knob_that_responds():
    def responds(t):
        return [f"ACGT{t}", f"TTTT{t}"]

    v = verify_knob_live(responds)
    assert v["live"] is True and v["verdict"] == "TEMPERATURE_IS_LIVE"


def test_zero_sequences_is_INDETERMINATE_not_a_pass():
    """An empty generation must not read as 'the knob works'."""
    v = verify_knob_live(lambda _t: [])
    assert v["live"] is False and v["verdict"] == "INDETERMINATE_NO_SEQUENCES"


def test_the_two_calls_use_DIFFERENT_temperatures():
    """A control that passed the same value twice would always report IGNORED."""
    seen = []

    def record(t):
        seen.append(t)
        return ["A"]

    verify_knob_live(record)
    assert len(seen) == 2 and seen[0] != seen[1]
