"""Every registering route declares a g->p regime, and an undeclared one raises rather than defaults."""
from __future__ import annotations

import pytest

from dna_decode.data.cell_regime import (
    CATALOG_ROUTES, CURATED_CATALOG, LEARNED_ATTEMPT_CLOSED, ROUTE_REGIME, UndeclaredRegime,
    regime_census, regime_for_route, regime_record,
)
from dna_decode.data.cell_registry import cells
from dna_decode.eval.regime import REGIMES

ROUTES = {c.route for c in cells()}


def test_every_registered_route_declares_a_regime():
    """The tripwire. A new decoder route must be classified before it can ship."""
    undeclared = sorted(r for r in ROUTES if r not in ROUTE_REGIME and r not in CATALOG_ROUTES)
    assert not undeclared, (
        f"{len(undeclared)} route(s) register cells with no declared g->p regime: {undeclared}. "
        "Classify each in dna_decode/data/cell_regime.py -- check eval/regime.py first.")


def test_an_unknown_route_raises_rather_than_taking_the_catalog_default():
    """Silently defaulting is the failure this exists to prevent: a learned decoder would inherit
    'curated catalog' and never face the regime screen."""
    with pytest.raises(UndeclaredRegime):
        regime_for_route("dna-somenewthing")


def test_the_regime_column_is_not_constant():
    """If every cell mapped to one regime the mapping would carry no information and this module would
    be ceremony. It was measured before building precisely to check that."""
    seen = {regime_for_route(r) for r in ROUTES}
    assert len(seen) >= 3, seen


def test_the_learned_and_constructed_routes_are_classified_as_such():
    assert regime_for_route("dna-decode-forward") == "constructed_molecular"
    assert regime_for_route("dna-decode-inverse") == "constructed_molecular"
    assert regime_for_route("dna-fba") == "constructed_organism_per_condition"
    assert regime_for_route("dna-amr") == CURATED_CATALOG


def test_every_declared_regime_key_is_a_real_regime():
    keys = {r.key for r in REGIMES}
    for route, key in ROUTE_REGIME.items():
        assert key in keys, f"{route} -> unknown regime {key!r}"
    assert CURATED_CATALOG in keys


def test_a_route_with_a_closed_learned_attempt_says_so():
    """dna-flowering ships as a catalog; the EMBEDDING attempt at that same trait is a closed negative.
    Recording it is what stops the learned version being re-proposed as novel."""
    rec = regime_record("dna-flowering")
    assert "learned_attempt_closed" in rec and "population structure" in rec["learned_attempt_closed"]


def test_a_plain_catalog_route_carries_no_closed_learned_flag():
    """Non-vacuity for the flag above: it must not be present on everything."""
    assert "learned_attempt_closed" not in regime_record("dna-amr")


def test_every_closed_learned_attempt_names_a_route_that_exists():
    for route in LEARNED_ATTEMPT_CLOSED:
        assert route in ROUTES, f"{route} has a closed-learned note but registers no cells"


def test_the_census_raises_on_an_undeclared_route_instead_of_reporting_around_it():
    with pytest.raises(UndeclaredRegime):
        regime_census(ROUTES | {"dna-madeup"})


def test_the_census_never_claims_to_revalidate_a_cell():
    """A WORKS regime is about the regime, not about this cell's accuracy."""
    assert "NOT a re-validation" in regime_census(ROUTES)["note"]
