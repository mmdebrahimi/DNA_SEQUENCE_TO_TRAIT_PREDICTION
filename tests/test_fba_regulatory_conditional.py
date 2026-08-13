"""Regulatory-constraint experiment -- the CONTROL logic (pure; no cobra, no solver)."""
from __future__ import annotations

import pytest

from dna_decode.fba.conditional_essentiality import (
    CONDITIONS,
    conditionally_essential_genes,
    load_labels,
)
from scripts.fba_regulatory_conditional_test import rate_matched_null


def test_rate_matched_null_is_the_control_that_makes_the_number_mean_something():
    """Forcing a unique route makes far MORE genes look essential, so a per-cell gain could be nothing
    but a better-matched base rate. Measured: calling 96 of 272 cells essential at random scores ~0.517,
    while the pFBA arm scored 0.6157 -- so the gain is not a base-rate artifact."""
    recs = load_labels()
    got = rate_matched_null(recs, 96, n_draws=50)
    assert got["n_draws"] == 50
    assert got["mean"] == pytest.approx(0.517, abs=0.02)
    assert got["max"] < 0.6157          # no random draw reaches the observed value


def test_rate_matched_null_is_deterministic_for_a_given_seed():
    """The control must be reproducible, or 'p < 0.005' is not a number anyone can check."""
    recs = load_labels()
    a = rate_matched_null(recs, 96, n_draws=20, seed0=7)
    b = rate_matched_null(recs, 96, n_draws=20, seed0=7)
    assert a["mean"] == b["mean"] and a["max"] == b["max"]


def test_calling_every_cell_essential_reproduces_the_constant_null():
    """Sanity pin tying the rate-matched control to the constant-predictor null already in the module."""
    recs = load_labels()
    n_cells = len(conditionally_essential_genes(recs)) * len(CONDITIONS)
    got = rate_matched_null(recs, n_cells, n_draws=3)
    assert got["mean"] == pytest.approx(0.4412, abs=1e-3)   # == always_essential


def test_rate_matched_null_refuses_an_impossible_request():
    """More essential calls than there are cells is a caller bug; report it, do not sample-with-error."""
    recs = load_labels()
    assert rate_matched_null(recs, 10**6, n_draws=5)["mean"] is None


# ---- the closed hardcode + the abstention arm ----

def test_rate_matched_null_uses_the_callers_conditions_not_the_four_media():
    """It was hardcoded TWICE: explicitly here, and again by letting switch_accuracy fall through to
    its own 4-media default. Correct for this script, invisible to its tests, and the exact shape of the
    bug that forced the 84.8% retraction."""
    from dna_decode.fba.conditional_essentiality import GeneRecord

    keys = tuple(f"c{i}" for i in range(25))
    recs = [GeneRecord(f"b{i}", f"b{i}", {k: (k == keys[i]) for k in keys}, {}, True)
            for i in range(6)]
    got = rate_matched_null(recs, 10, n_draws=5, conditions=keys)
    assert got["n_draws"] == 5
    # 6 genes x 25 conditions = 150 cells; the 4-media default would have offered only 24 and the
    # request for 10 would still have "worked" -- silently sampling the wrong space.
    assert rate_matched_null(recs, 149, n_draws=1, conditions=keys)["n_draws"] == 1
    assert rate_matched_null(recs, 149, n_draws=1)["mean"] is None      # 4-media default: impossible


def test_abstention_removes_a_cell_from_BOTH_numerator_and_denominator():
    """Deleting a cell from `predicted` is NOT abstention -- switch_accuracy defaults a missing cell to
    False and still counts it, so a dropped cell is silently scored as 'dispensable'. Fail-open."""
    from dna_decode.fba.conditional_essentiality import GeneRecord, switch_accuracy

    keys = ("a", "b")
    r = GeneRecord("b1", "b1", {"a": True, "b": False}, {}, True)
    pred = {"a": {"b1": True}, "b": {"b1": False}}                      # perfect

    full = switch_accuracy([r], pred, conditions=keys)
    assert full["n_cells_scored"] == 2 and full["per_condition_agreement"] == 1.0

    dropped = switch_accuracy([r], {"a": {}, "b": {"b1": False}}, conditions=keys)
    assert dropped["n_cells_scored"] == 2                               # still counted!
    assert dropped["per_condition_agreement"] == 0.5                    # scored as dispensable

    abstained = switch_accuracy([r], pred, conditions=keys, exclude_cells={("b1", "a")})
    assert abstained["n_cells_scored"] == 1 and abstained["n_cells_abstained"] == 1
    assert abstained["per_condition_agreement"] == 1.0


def test_a_gene_with_an_abstained_cell_leaves_the_exact_set_denominator():
    """An incomplete pattern must not be judged on its surviving cells -- that would make exact-set
    EASIER exactly where the solver struggled."""
    from dna_decode.fba.conditional_essentiality import GeneRecord, switch_accuracy

    keys = ("a", "b")
    recs = [GeneRecord("b1", "b1", {"a": True, "b": False}, {}, True),
            GeneRecord("b2", "b2", {"a": False, "b": True}, {}, True)]
    pred = {"a": {"b1": True, "b2": False}, "b": {"b1": False, "b2": True}}

    got = switch_accuracy(recs, pred, conditions=keys, exclude_cells={("b1", "a")})
    assert got["n_scored_exact_set"] == 1          # b1 dropped, b2 kept
    assert got["exact_set_match"] == 1
    assert got["exact_set_match_rate"] == 1.0


def test_the_null_must_be_recomputed_on_the_abstained_denominator():
    """A null over the full cell set compared against a metric over a reduced one is not a control."""
    from dna_decode.fba.conditional_essentiality import constant_baselines

    recs = load_labels()
    subset = conditionally_essential_genes(recs)
    keys = tuple(sorted(CONDITIONS))
    excl = {(subset[i].gene_id, keys[0]) for i in range(10)}

    full = constant_baselines(recs, conditions=keys)
    reduced = constant_baselines(recs, conditions=keys, exclude_cells=excl)
    assert reduced["always_dispensable"]["n_cells_abstained"] == 10
    assert full["always_dispensable"]["n_cells_scored"] != \
        reduced["always_dispensable"]["n_cells_scored"]


def test_exclude_cells_defaults_preserve_every_published_number():
    """Backward-compat pin: no exclude_cells -> identical to the pre-abstention behaviour."""
    from dna_decode.fba.conditional_essentiality import switch_accuracy

    recs = load_labels()
    pred = {c: {r.gene_id: r.paper_fba[c] for r in recs} for c in CONDITIONS}
    got = switch_accuracy(recs, pred)
    assert got["n_cells_abstained"] == 0
    assert got["n_scored_exact_set"] == got["n_conditionally_essential"]
