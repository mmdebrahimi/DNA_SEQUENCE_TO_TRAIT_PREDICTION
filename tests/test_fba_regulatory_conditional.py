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
