"""Pin the g->p regime boundary, including the SCOPE that keeps getting lost.

The boundary has been mis-stated three times, always by compressing a scoped negative into a general
one. These tests make each of those three compressions fail loudly.

Offline, pure, no fixtures.
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from dna_decode.eval.regime import (CLOSED_NEGATIVE, LOSES_TO_CATALOG, OPEN,  # noqa: E402
                                    REGIMES, REQUIRES_DECONFOUNDING, WORKS, classify_regime,
                                    screen_proposal)


# --- the three compressions that have actually happened ---

def test_organism_level_gp_is_NOT_a_closed_negative_in_general():
    """COMPRESSION 1, made three times. A clean organism-level positive exists (yeast segregant cross,
    12/12 traits at r 0.46-0.80), so a blanket 'organism g->p is closed' is false."""
    assert screen_proposal("constructed", "organism", "supervised").verdict == WORKS
    assert not screen_proposal("constructed", "organism", "supervised").refused


def test_the_negative_is_zero_shot_scoped_not_learning_scoped():
    """COMPRESSION 2. The 0-for-5 is ZERO-SHOT-only; the shipped architecture pairs the deterministic
    catalogue with a SUPERVISED complement. Refusing a supervised proposal would re-commit the error."""
    zero = screen_proposal("natural", "organism", "zero_shot")
    sup = screen_proposal("natural", "organism", "supervised")
    assert zero.verdict == CLOSED_NEGATIVE and zero.refused
    assert sup.verdict == REQUIRES_DECONFOUNDING and not sup.refused, (
        "a supervised natural-population proposal must get CONDITIONS, not a refusal")
    assert any("WITHIN-GROUP" in c for c in sup.conditions)


def test_the_discriminating_variable_is_population_design_not_organism_complexity():
    """COMPRESSION 3. Holding endpoint and method fixed, flipping ONLY the population design flips the
    verdict — which is what 'the discriminator is population design' means operationally."""
    natural = screen_proposal("natural", "organism", "zero_shot").verdict
    constructed = screen_proposal("constructed", "organism", "supervised").verdict
    assert natural == CLOSED_NEGATIVE and constructed == WORKS


# --- the refusal is narrow ---

def test_only_the_measured_dead_regime_is_refused():
    """A boundary that refused broadly would block real work. Exactly one regime refuses."""
    refused = [r for r in REGIMES
               if screen_proposal(r.population, r.endpoint, r.method).refused]
    assert len(refused) == 1, f"expected exactly one refusing regime, got {[r.key for r in refused]}"
    assert refused[0].key == "natural_organism_zeroshot"


def test_an_unscreened_combination_is_open_not_endorsed():
    """Absence of a measurement must not read as promise."""
    res = screen_proposal("constructed", "organism_condition_switch", "zero_shot")
    assert res.verdict == OPEN and not res.refused
    assert "not the same as promising" in res.reason


def test_the_condition_switch_cell_is_open_not_solved():
    """The one genuinely open cell — it must not be reported as working."""
    res = screen_proposal("constructed", "organism_condition_switch", "supervised")
    assert res.verdict == OPEN
    assert "silent" in res.evidence.lower()


# --- the catalog rule, and its inversion ---

def test_a_curated_catalog_beats_a_learned_scorer_wherever_one_exists():
    res = screen_proposal("constructed", "molecular", "supervised", curated_catalog_exists=True)
    assert res.verdict == LOSES_TO_CATALOG
    assert "0.454" in res.evidence, "the BELOW-chance measurement must ship with the verdict"


def test_the_catalog_rule_does_not_fire_on_the_catalog_itself():
    """A deterministic catalog proposal must not be told it loses to a catalog."""
    res = screen_proposal("natural", "molecular", "deterministic_catalog", curated_catalog_exists=True)
    assert res.verdict != LOSES_TO_CATALOG


def test_without_a_catalog_the_working_molecular_regime_still_works():
    assert screen_proposal("constructed", "molecular", "supervised").verdict == WORKS


# --- hygiene ---

def test_scale_is_never_offered_as_the_remedy_for_the_closed_negative():
    """The specific wrong next step. It must be named as wrong in the output, not just omitted."""
    res = screen_proposal("natural", "organism", "zero_shot")
    assert any("do NOT re-run this at larger scale" in c for c in res.conditions)


def test_every_regime_cites_an_artifact_that_exists():
    """A claim whose evidence file is missing is not a claim. Same rule as the doc-citation guard."""
    missing = [r.key for r in REGIMES if not (ROOT / r.artifact).exists()]
    assert not missing, f"regimes cite artifacts that do not exist: {missing}"


def test_bad_axis_values_are_refused_not_guessed():
    for bad in (("eukaryote", "organism", "zero_shot"),
                ("natural", "phenotype", "zero_shot"),
                ("natural", "organism", "finetuned")):
        assert screen_proposal(*bad).verdict == "UNKNOWN"


def test_classify_regime_is_exact_not_fuzzy():
    assert classify_regime("NATURAL", "Organism", "ZERO_SHOT").key == "natural_organism_zeroshot"
    assert classify_regime("natural", "organism", "deterministic_catalog") is None
