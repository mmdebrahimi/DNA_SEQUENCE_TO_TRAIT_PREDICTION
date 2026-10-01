"""Guards for clade and regime specializations (plan Step 5, 2026-10-01)."""

from __future__ import annotations

import pytest

from dna_decode.constraints import registry as R
from dna_decode.constraints import specializations as S
from dna_decode.constraints import universal as U


@pytest.fixture(autouse=True)
def _registry():
    R.clear_registry()
    U.register_universal()
    S.register_specializations()
    yield
    R.clear_registry()


def test_a_deeper_clade_overrides_a_shallower_one():
    """`bacteria.mycoplasma` must resolve to table 4 while `bacteria` resolves to table 11."""
    myco = R.constraints_for(taxon="bacteria.mycoplasma")[S.GENETIC_CODE_KEY]
    bact = R.constraints_for(taxon="bacteria")[S.GENETIC_CODE_KEY]
    assert myco.payload["table_id"] == 4
    assert bact.payload["table_id"] == 11
    assert myco.scope == "taxon:bacteria.mycoplasma"


def test_each_specialization_reports_the_codons_it_actually_changes():
    """A specialization that changes nothing is decorative; these record exactly what differs."""
    got = {c.payload["clade"]: c.payload["differs_from_standard"]
           for c in S.genetic_code_specializations()}
    assert got["bacteria.mycoplasma"] == ["TGA"]
    assert got["mitochondrion.vertebrate"] == ["AGA", "AGG", "ATA", "TGA"]
    assert got["eukaryote.ciliate"] == ["TAA", "TAG"]
    assert got["bacteria"] == []            # table 11 equals standard -- the control


def test_every_shipped_specialization_carries_both_provenance_fields():
    """The registry refuses a sourced constraint without them; this asserts none slipped through."""
    for c in S.genetic_code_specializations():
        assert c.kind == R.SOURCED, c.scope
        assert c.source_url.startswith("https://www.ncbi.nlm.nih.gov/"), c.scope
        assert "AAs  =" in c.verbatim_quote, c.scope
        assert c.note, c.scope


def test_a_clade_without_a_specialization_still_gets_the_universal_laws():
    cs = R.constraints_for(taxon="fungi")
    assert S.GENETIC_CODE_KEY not in cs
    assert {c.key for c in U.UNIVERSAL_CONSTRAINTS} <= set(cs)


def test_the_shipped_set_has_no_conflicts():
    assert R.detect_conflicts() == []


def test_the_under_sourced_clade_properties_are_DECLARED_but_NOT_registered():
    """THE GATE WORKING ON ITS OWN AUTHOR. Splicing, operons, HGT and ploidy belong on this axis and are
    well established, but shipping them would mean writing citations from recall. They are declared with
    the evidence each needs and deliberately left unregistered."""
    assert set(S.DEFERRED_CLADE_PROPERTIES) == {
        "splicing_applies", "operon_polycistronic", "hgt_acquisition_possible", "ploidy"}
    for key, rec in S.DEFERRED_CLADE_PROPERTIES.items():
        assert rec["evidence_needed"] and rec["why_it_matters"] and rec["buckets"], key
    registered_keys = {c.key for c in R.all_constraints()}
    for key in S.DEFERRED_CLADE_PROPERTIES:
        assert key not in registered_keys, f"{key} must not ship without provenance"


def test_registering_a_deferred_property_without_provenance_is_refused():
    """Non-vacuity for the claim above: attempting it actually fails."""
    with pytest.raises(R.UnsourcedConstraintError):
        R.register(R.Constraint("splicing_applies", "taxon:eukaryote", R.SOURCED,
                                lambda ctx: R.ConstraintVerdict(R.SATISFIED)))


def test_the_regime_axis_is_advisory_metadata_NOT_constraints():
    """A regime describes the study, so it cannot refute a sequence. Registering always-inapplicable
    constraints would manufacture exactly the inert filters report.py exists to catch."""
    assert sum(1 for c in R.all_constraints() if c.axis == "regime") == 0
    table = S.regime_table()
    from dna_decode.eval.regime import REGIMES
    assert set(table) == {r.key for r in REGIMES}
    assert "never refutes a sequence" in table[REGIMES[0].key]["advisory"]


def test_regime_scope_surfaces_the_closed_negative():
    sc = S.regime_scope("natural_organism_zeroshot")
    assert sc["verdict"] == "CLOSED_NEGATIVE" and sc["is_closed"] is True
    assert S.regime_scope("constructed_molecular")["is_closed"] is False
    with pytest.raises(KeyError):
        S.regime_scope("not_a_regime")


def test_no_higher_life_form_axis_is_encoded():
    """Clade nuances swap rather than accumulate; depth must not imply more constraints."""
    src = (__import__("pathlib").Path(S.__file__)).read_text(encoding="utf-8")
    assert "SWAP rather than accumulate" in src
    # a depth-2 bucket does not automatically carry more constraints than a depth-1 one
    shallow = R.constraints_for(taxon="bacteria")
    deep = R.constraints_for(taxon="bacteria.mycoplasma")
    assert set(shallow) == set(deep), "same KEYS; only the resolved payload differs"
