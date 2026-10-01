"""Gap coverage for `constraints/registry.py` resolution internals (test-epilogue, 2026-10-01).

The plan's own `test_constraints_registry.py` covers the headline contracts: deeper-wins, cross-axis
refusal, unsourced refusal, and `detect_conflicts` on a CLEAN set. Four things it does not reach, and
each is a branch that decides a verdict:

* `_ancestors` is never exercised directly, yet it is what makes the "equal specificity" raise inside
  `constraints_for` unreachable within one axis -- a chain holds exactly ONE bucket per depth.
* precedence is only tested with the universal constraint registered FIRST. `_REGISTRY` iterates in
  insertion order, so the reverse order takes a DIFFERENT branch (`elif axis == UNIVERSAL: continue`).
  An order-dependent precedence rule would be a silent bug.
* `detect_conflicts` is only asserted CLEAN. A conflict detector never seen to fire is the inert-filter
  failure this plan's own `report.py` exists to catch, so its positive branch is pinned here -- together
  with the fact that its taxon sibling-depth loop reports NOTHING by design.
* `register`'s input validation (bad `kind`, malformed scope) and the `Constraint` accessors.
"""

from __future__ import annotations

import pytest

from dna_decode.constraints import registry as R
from dna_decode.constraints.registry import (
    AXIOMATIC,
    INAPPLICABLE,
    REFUTED,
    SATISFIED,
    SOURCED,
    UNIVERSAL,
    Constraint,
    ConstraintConflictError,
    ConstraintVerdict,
)


def _ok(_ctx):
    return ConstraintVerdict(SATISFIED)


@pytest.fixture(autouse=True)
def _isolated_registry():
    """The registry is module-global; every test here owns it outright and leaves it empty."""
    R.clear_registry()
    yield
    R.clear_registry()


def _c(key: str, scope: str, **kw) -> Constraint:
    return Constraint(key, scope, kw.pop("kind", AXIOMATIC), _ok, **kw)


# --------------------------------------------------------------------- _ancestors

def test_ancestors_builds_the_inheritance_chain_and_holds_ONE_bucket_PER_DEPTH():
    """The chain property is load-bearing: it is why two same-depth taxon buckets can never both apply
    to one query, which in turn is why `constraints_for`'s equal-specificity raise cannot fire within
    the taxon axis. Pinned here so a future change to nesting notices it broke that invariant."""
    assert R._ancestors("bacteria.mycoplasma") == ["bacteria", "bacteria.mycoplasma"]
    assert R._ancestors("virus") == ["virus"]
    assert R._ancestors("a.b.c") == ["a", "a.b", "a.b.c"]

    for bucket in R.TAXON_BUCKETS:
        chain = R._ancestors(bucket)
        depths = [b.count(".") + 1 for b in chain]
        assert depths == sorted(set(depths)), bucket
        assert len(depths) == len(set(depths)), f"{bucket} chain has two buckets at one depth"


def test_every_declared_bucket_has_its_ancestors_declared_too():
    """`bacteria.mycoplasma` inherits from `bacteria`, so a nested bucket whose parent is undeclared
    would inherit from a scope no constraint can ever target."""
    for bucket in R.TAXON_BUCKETS:
        for ancestor in R._ancestors(bucket)[:-1]:
            assert ancestor in R.TAXON_BUCKETS, (
                f"{bucket!r} nests under undeclared parent {ancestor!r}")


# --------------------------------------------------------------------- precedence is order-independent

def test_a_specialization_beats_universal_REGARDLESS_of_registration_order():
    """Registration order must not decide a verdict. The two orders take different branches inside
    `constraints_for`, so both are exercised."""
    R.clear_registry()
    R.register(_c("k", UNIVERSAL))
    R.register(_c("k", "taxon:bacteria"))
    assert R.constraints_for(taxon="bacteria")["k"].scope == "taxon:bacteria"

    R.clear_registry()
    R.register(_c("k", "taxon:bacteria"))          # specialization FIRST this time
    R.register(_c("k", UNIVERSAL))
    assert R.constraints_for(taxon="bacteria")["k"].scope == "taxon:bacteria", (
        "precedence is order-dependent: the universal constraint overwrote a specialization")


def test_a_regime_specialization_also_beats_universal_in_either_order():
    regime = R.regime_keys()[0]
    for order in (("u", "r"), ("r", "u")):
        R.clear_registry()
        for which in order:
            R.register(_c("k", UNIVERSAL) if which == "u" else _c("k", f"regime:{regime}"))
        got = R.constraints_for(regime=regime)["k"].scope
        assert got == f"regime:{regime}", f"order {order} resolved to {got}"


def test_a_three_level_chain_resolves_to_the_DEEPEST_bucket():
    """`_ancestors` depth ordering is only meaningful if the deepest of THREE levels wins, not merely
    the deeper of two."""
    R.register(_c("k", UNIVERSAL))
    R.register(_c("k", "taxon:eukaryote"))
    R.register(_c("k", "taxon:eukaryote.ciliate"))
    assert R.constraints_for(taxon="eukaryote.ciliate")["k"].scope == "taxon:eukaryote.ciliate"
    assert R.constraints_for(taxon="eukaryote.nuclear")["k"].scope == "taxon:eukaryote"
    assert R.constraints_for()["k"].scope == UNIVERSAL


def test_a_regime_constraint_is_INVISIBLE_until_a_regime_is_supplied():
    """A query that names no regime must not silently pick up a regime-scoped rule -- otherwise a
    study-design specialization would leak into an unscoped call."""
    regime = R.regime_keys()[0]
    R.register(_c("only_regime", f"regime:{regime}"))
    R.register(_c("always", UNIVERSAL))
    assert sorted(R.constraints_for(taxon="bacteria")) == ["always"]
    assert sorted(R.constraints_for(regime=regime)) == ["always", "only_regime"]


def test_a_taxon_constraint_is_invisible_to_an_unrelated_taxon_and_to_no_taxon():
    R.register(_c("bact_only", "taxon:bacteria"))
    assert R.constraints_for(taxon="fungi") == {}
    assert R.constraints_for() == {}
    assert R.constraints_for(taxon="bacteria.mycoplasma")["bact_only"].scope == "taxon:bacteria"


# --------------------------------------------------------------------- detect_conflicts

def test_detect_conflicts_ACTUALLY_FIRES_on_a_cross_axis_claim():
    """A detector only ever asserted clean is indistinguishable from one that cannot fire. This is the
    positive direction, and it names both offending scopes."""
    regime = R.regime_keys()[0]
    R.register(_c("shared", "taxon:bacteria"))
    R.register(_c("shared", f"regime:{regime}"))
    problems = R.detect_conflicts()
    assert len(problems) == 1, problems
    assert "shared" in problems[0]
    assert "taxon:bacteria" in problems[0] and f"regime:{regime}" in problems[0]
    assert "BOTH axes" in problems[0]


def test_a_universal_constraint_sharing_a_key_with_ONE_axis_is_not_a_conflict():
    """Universal is overridable by either axis, so universal + one specialization is the normal shape
    and must not be reported."""
    R.register(_c("k", UNIVERSAL))
    R.register(_c("k", "taxon:bacteria"))
    assert R.detect_conflicts() == []


def test_detect_conflicts_reports_NOTHING_for_sibling_taxon_branches_BY_DESIGN():
    """`bacteria` and `fungi` may both claim one key: no query resolves both, because a query walks ONE
    ancestor chain. The scan's taxon loop therefore records nothing -- documented here so the empty
    result reads as a deliberate scope limit rather than as a detector that silently does no work."""
    R.register(_c("k", "taxon:bacteria"))
    R.register(_c("k", "taxon:fungi"))
    assert R.detect_conflicts() == []
    assert R.constraints_for(taxon="bacteria")["k"].scope == "taxon:bacteria"
    assert R.constraints_for(taxon="fungi")["k"].scope == "taxon:fungi"


def test_the_cross_axis_conflict_is_refused_at_QUERY_time_too_not_only_scanned():
    """`detect_conflicts` is advisory; `constraints_for` is the enforcement point. Both must agree."""
    regime = R.regime_keys()[0]
    R.register(_c("shared", "taxon:bacteria"))
    R.register(_c("shared", f"regime:{regime}"))
    assert R.detect_conflicts()
    with pytest.raises(ConstraintConflictError, match="orthogonal"):
        R.constraints_for(taxon="bacteria", regime=regime)
    # ...and NOT refused when the query touches only one of the two axes
    assert R.constraints_for(taxon="bacteria")["shared"].scope == "taxon:bacteria"
    assert R.constraints_for(regime=regime)["shared"].scope == f"regime:{regime}"


# --------------------------------------------------------------------- register() validation

def test_register_refuses_an_unrecognised_kind():
    """`kind` decides whether provenance is required, so an unrecognised value must not slip through as
    'not sourced, therefore no citation needed'."""
    with pytest.raises(ValueError, match="kind must be"):
        R.register(Constraint("k", UNIVERSAL, "advisory", _ok))
    assert R.all_constraints() == ()


@pytest.mark.parametrize("scope", ["bogus", "bogus:thing", "TAXON:bacteria", "taxon:"])
def test_register_refuses_a_malformed_scope(scope):
    with pytest.raises(ValueError):
        R.register(Constraint("k", scope, AXIOMATIC, _ok))


@pytest.mark.parametrize("scope", ["taxon", "regime"])
def test_a_COLONLESS_axis_name_raises_IndexError_not_a_clean_ValueError(scope):
    """Rough edge, pinned rather than papered over. `Constraint.bucket` does `split(":", 1)[1]`, so a
    scope that is a bare axis name with no colon raises `IndexError` from the property before
    `register`'s own scope validation can produce its actionable message.

    It is a MESSAGE-QUALITY defect, not a safety one -- nothing is silently accepted, which is what
    matters for a layer whose whole stance is refusal over guessing. Recorded so the behaviour is a
    known limit; if `bucket` is ever hardened, this test should flip to `ValueError`."""
    with pytest.raises(IndexError):
        R.register(Constraint("k", scope, AXIOMATIC, _ok))
    assert R.all_constraints() == (), "a rejected constraint must not land in the registry"


def test_an_unsourced_sourced_constraint_names_WHICH_field_is_missing():
    """The message is the actionable part; a generic refusal makes a caller guess."""
    with pytest.raises(R.UnsourcedConstraintError, match="source_url"):
        R.register(_c("k", UNIVERSAL, kind=SOURCED, verbatim_quote="q"))
    with pytest.raises(R.UnsourcedConstraintError, match="verbatim_quote"):
        R.register(_c("k", UNIVERSAL, kind=SOURCED, source_url="http://x"))


def test_clear_registry_and_all_constraints_round_trip():
    """NOTE: `Constraint` is `frozen=True` but carries a dict `payload`, so instances are UNHASHABLE --
    compare by scope rather than reaching for a set."""
    assert R.all_constraints() == ()
    R.register(_c("a", UNIVERSAL))
    R.register(_c("b", "taxon:bacteria"))
    assert sorted(c.key for c in R.all_constraints()) == ["a", "b"]
    R.clear_registry()
    assert R.all_constraints() == ()


def test_a_frozen_constraint_is_unhashable_so_callers_must_not_set_dedupe_them():
    """Pinned because `frozen=True` ordinarily implies hashable, and a caller that builds a set of
    constraints would crash at runtime rather than at review time."""
    with pytest.raises(TypeError, match="unhashable"):
        {Constraint("k", UNIVERSAL, AXIOMATIC, _ok)}


# --------------------------------------------------------------------- dataclass accessors

@pytest.mark.parametrize("scope,axis,bucket,depth", [
    (UNIVERSAL, UNIVERSAL, "", 0),
    ("taxon:bacteria", "taxon", "bacteria", 1),
    ("taxon:bacteria.mycoplasma", "taxon", "bacteria.mycoplasma", 2),
    ("regime:X", "regime", "X", 1),
])
def test_constraint_axis_bucket_and_depth_accessors(scope, axis, bucket, depth):
    """`depth` is the precedence key, so each shape's value is pinned rather than inferred."""
    c = Constraint("k", scope, AXIOMATIC, _ok)
    assert (c.axis, c.bucket, c.depth) == (axis, bucket, depth)


def test_a_regime_scope_and_a_top_level_taxon_scope_have_EQUAL_depth():
    """Which is exactly why there is no cross-axis ordering: depth cannot break the tie, so the code
    refuses instead of guessing. Pinned so a 'fix' that ranks one axis above the other is visible."""
    taxon = Constraint("k", "taxon:bacteria", AXIOMATIC, _ok)
    regime = Constraint("k", "regime:X", AXIOMATIC, _ok)
    assert taxon.depth == regime.depth == 1


def test_verdict_vocabulary_rejects_anything_outside_the_three_states():
    for bad in ("maybe", "", "SATISFIED", "pass"):
        with pytest.raises(ValueError, match="verdict must be one of"):
            ConstraintVerdict(bad)
    assert ConstraintVerdict(REFUTED).refuted is True
    assert ConstraintVerdict(SATISFIED).refuted is False
    assert ConstraintVerdict(INAPPLICABLE).refuted is False, (
        "inapplicable must never read as refuted -- 'could not check' is not 'found a problem'")


def test_a_constraints_payload_default_is_not_shared_between_instances():
    """A mutable class-level default would make one constraint's payload leak into every other."""
    a, b = Constraint("a", UNIVERSAL, AXIOMATIC, _ok), Constraint("b", UNIVERSAL, AXIOMATIC, _ok)
    a.payload["x"] = 1
    assert b.payload == {}


# --------------------------------------------------------------------- query-side input validation

def test_constraints_for_refuses_an_unknown_taxon_or_regime_rather_than_returning_universal_only():
    """Silently returning the universal set for a typo'd clade would look like a clean, narrow answer."""
    R.register(_c("k", UNIVERSAL))
    with pytest.raises(ValueError, match="unknown taxon bucket"):
        R.constraints_for(taxon="bactria")
    with pytest.raises(ValueError, match="unknown regime key"):
        R.constraints_for(regime="made_up")
    assert sorted(R.constraints_for()) == ["k"]
