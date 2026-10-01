"""Guards for the constraint registry — two axes, precedence, refusals (plan Step 2, 2026-10-01)."""

from __future__ import annotations

import pytest

from dna_decode.constraints import registry as R

OK = lambda ctx: R.ConstraintVerdict(R.SATISFIED)  # noqa: E731


@pytest.fixture(autouse=True)
def _clean_registry():
    R.clear_registry()
    yield
    R.clear_registry()


def _sourced(key, scope, **payload):
    return R.Constraint(key, scope, R.SOURCED, OK, source_url="https://www.ncbi.nlm.nih.gov/x",
                        verbatim_quote="AAs  = FFLL...", payload=payload)


def test_the_regime_axis_is_IMPORTED_from_eval_regime_not_redeclared():
    """A second copy of the regime keys is the drift this repo has hit five times."""
    from dna_decode.eval.regime import REGIMES
    assert R.regime_keys() == tuple(r.key for r in REGIMES)
    assert len(R.regime_keys()) >= 8
    src = (__import__("pathlib").Path(R.__file__)).read_text(encoding="utf-8")
    assert "natural_organism_zeroshot" not in src, "regime keys must not be hardcoded here"


def test_a_deeper_taxon_bucket_overrides_a_shallower_one_and_universal():
    R.register(R.Constraint("genetic_code", R.UNIVERSAL, R.AXIOMATIC, OK, payload={"table": 1}))
    R.register(_sourced("genetic_code", "taxon:bacteria", table=11))
    R.register(_sourced("genetic_code", "taxon:bacteria.mycoplasma", table=4))

    assert R.constraints_for(taxon="bacteria.mycoplasma")["genetic_code"].payload == {"table": 4}
    assert R.constraints_for(taxon="bacteria")["genetic_code"].payload == {"table": 11}
    assert R.constraints_for()["genetic_code"].scope == R.UNIVERSAL


def test_a_sibling_bucket_does_not_inherit_from_another_branch():
    R.register(_sourced("genetic_code", "taxon:bacteria", table=11))
    assert "genetic_code" not in R.constraints_for(taxon="fungi")


def test_an_unsourced_sourced_constraint_is_REFUSED_both_fields_required():
    with pytest.raises(R.UnsourcedConstraintError):
        R.register(R.Constraint("a", "taxon:fungi", R.SOURCED, OK, verbatim_quote="q"))
    with pytest.raises(R.UnsourcedConstraintError):
        R.register(R.Constraint("b", "taxon:fungi", R.SOURCED, OK, source_url="https://x"))
    # non-vacuity: with BOTH present it registers
    assert R.register(_sourced("c", "taxon:fungi")).key == "c"


def test_an_axiomatic_constraint_needs_no_citation():
    """A law is not a citation matter -- mass balance is not someone's finding."""
    assert R.register(R.Constraint("mass_balance", R.UNIVERSAL, R.AXIOMATIC, OK)).kind == R.AXIOMATIC


def test_a_cross_axis_claim_on_one_key_is_REFUSED_not_silently_resolved():
    """The axes are orthogonal and no precedence between them has been validated, so the ambiguity is
    surfaced. Same stance as table_for refusing an unknown clade."""
    rk = R.regime_keys()[0]
    R.register(R.Constraint("collide", "taxon:bacteria", R.AXIOMATIC, OK))
    R.register(R.Constraint("collide", f"regime:{rk}", R.AXIOMATIC, OK))
    with pytest.raises(R.ConstraintConflictError, match="both"):
        R.constraints_for(taxon="bacteria", regime=rk)
    assert R.detect_conflicts(), "detect_conflicts must report the same collision"
    # non-vacuity: querying only one axis is unambiguous and must still work
    assert R.constraints_for(taxon="bacteria")["collide"].axis == "taxon"


def test_a_duplicate_key_scope_registration_is_refused():
    R.register(R.Constraint("dup", "taxon:bacteria", R.AXIOMATIC, OK))
    with pytest.raises(R.ConstraintConflictError):
        R.register(R.Constraint("dup", "taxon:bacteria", R.AXIOMATIC, OK))


def test_unknown_buckets_and_regimes_are_refused_at_registration_and_query():
    with pytest.raises(ValueError, match="taxon bucket"):
        R.register(R.Constraint("x", "taxon:not_a_clade", R.AXIOMATIC, OK))
    with pytest.raises(ValueError, match="regime key"):
        R.register(R.Constraint("y", "regime:not_a_regime", R.AXIOMATIC, OK))
    with pytest.raises(ValueError):
        R.constraints_for(taxon="not_a_clade")


def test_detect_conflicts_is_clean_on_a_well_formed_set():
    R.register(R.Constraint("genetic_code", R.UNIVERSAL, R.AXIOMATIC, OK))
    R.register(_sourced("genetic_code", "taxon:bacteria", table=11))
    R.register(_sourced("genetic_code", "taxon:bacteria.mycoplasma", table=4))
    assert R.detect_conflicts() == []


def test_the_verdict_vocabulary_keeps_inapplicable_distinct_from_satisfied():
    """"could not check" and "checked and clean" are different claims."""
    assert R.INAPPLICABLE != R.SATISFIED
    assert R.ConstraintVerdict(R.REFUTED).refuted is True
    assert R.ConstraintVerdict(R.INAPPLICABLE).refuted is False
    with pytest.raises(ValueError):
        R.ConstraintVerdict("probably_fine")


def test_there_is_no_higher_life_form_axis():
    """Clade nuances SWAP rather than accumulate; depth must not imply "more constrained"."""
    src = (__import__("pathlib").Path(R.__file__)).read_text(encoding="utf-8").lower()
    assert "swap" in src and "higher life form" in src
    assert R.Constraint("k", "taxon:animal.mammal", R.AXIOMATIC, OK).depth == 2
    assert R.Constraint("k", "taxon:bacteria", R.AXIOMATIC, OK).depth == 1
