"""Kill-test: are iML1515's internal reaction bounds BINDING at all?

THE CLAIM UNDER TEST
    The E-Flux null (0 of 1,441 essentiality calls changed) is explained by capacity slack: iML1515's
    internal reaction bounds are so loose that scaling them -- which is exactly what E-Flux does -- is
    absorbed entirely and never touches a binding constraint. Growth is limited by SUBSTRATE UPTAKE
    (the exchange bounds), which E-Flux deliberately leaves alone.

WHY IT MATTERS
    It converts "expression constraints did nothing" from an empirical surprise into a structural fact
    about where the model's binding constraint lives. It also predicts the failure of any future method
    that works by tightening internal bounds.

HOW THIS COULD FAIL
    If internal bounds were near-binding, even a mild (2x) tightening would lower growth. The test
    asserts BOTH directions -- inert when mild, biting when severe -- so it cannot pass vacuously in
    either direction.

DISCOVERED BY A VACUITY GUARD, not by design: an earlier version of the sibling zero-flux test scaled
internal bounds by 0.5 to emulate E-Flux, and its own "did the perturbation actually bite?" assertion
failed, reporting bit-identical growth. That failure is this claim.

Model-only -- no Fitness Browser DB required.
"""
from __future__ import annotations

import pytest

pytest.importorskip("cobra", reason="cobrapy required")

TOL = 1e-6
MILD_SCALES = (0.5, 0.2, 0.1)   # measured INERT on iML1515
SEVERE_SCALE = 0.05             # measured to bite: 0.877 -> 0.715


@pytest.fixture(scope="module")
def model():
    from dna_decode.fba.model import load_model
    return load_model()


def _scaled_growth(model, scale: float) -> float:
    """Growth after scaling every GPR-carrying, non-exchange reaction's bounds by `scale`."""
    from dna_decode.fba.model import wildtype_growth
    exchanges = {r.id for r in model.exchanges}
    with model:
        for rxn in model.reactions:
            if rxn.id in exchanges or not rxn.gene_reaction_rule.strip():
                continue
            if rxn.upper_bound > 0:
                rxn.upper_bound *= scale
            if rxn.lower_bound < 0:
                rxn.lower_bound *= scale
        return float(wildtype_growth(model))


def test_internal_bounds_are_non_binding_under_mild_tightening(model):
    """Tightening every internal bound up to 10x must not move growth at all."""
    from dna_decode.fba.model import wildtype_growth

    base = float(wildtype_growth(model))
    moved = [(s, g) for s in MILD_SCALES
             if abs((g := _scaled_growth(model, s)) - base) > TOL]
    assert not moved, (
        f"an internal-bound tightening DID bind (base {base}): {moved} -- "
        "the capacity-slack explanation for the E-Flux null is false")


def test_the_test_is_not_vacuous_severe_tightening_does_bite(model):
    """The other direction: a severe enough tightening MUST reduce growth.

    Without this, `test_internal_bounds_are_non_binding...` could pass simply because the perturbation
    was never applied correctly.
    """
    from dna_decode.fba.model import wildtype_growth

    base = float(wildtype_growth(model))
    severe = _scaled_growth(model, SEVERE_SCALE)
    assert severe < base - TOL, (
        f"even a {1/SEVERE_SCALE:.0f}x tightening did not reduce growth "
        f"({base} -> {severe}) -- the perturbation is not being applied, so the sibling test proves "
        "nothing")


def test_substrate_uptake_by_contrast_IS_binding(model):
    """The positive control that locates the real constraint: halving carbon uptake must bite.

    This is what makes the finding a statement about WHERE the binding constraint lives, rather than
    a claim that the model is insensitive to everything.
    """
    from dna_decode.fba.model import wildtype_growth

    base = float(wildtype_growth(model))
    with model:
        medium = dict(model.medium)
        assert "EX_glc__D_e" in medium, "glucose uptake not in the default medium"
        medium["EX_glc__D_e"] *= 0.5
        model.medium = medium
        halved = float(wildtype_growth(model))
    assert halved < base - TOL, (
        f"halving the carbon uptake did not reduce growth ({base} -> {halved}) -- "
        "then the model is insensitive to everything and this whole framing is wrong")
