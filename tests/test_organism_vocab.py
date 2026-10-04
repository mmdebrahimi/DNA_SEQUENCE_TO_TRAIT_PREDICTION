"""Tests for the canonical organism vocabulary (closed-set identification router, Step 1).

The load-bearing test here is `test_amrfinder_value_and_registry_label_resolve_the_SAME_rule`: the
module's whole simplification -- one routing token instead of two -- rests on `amr_rules`'s
genus-prefix fallback. That fallback was verified empirically before the module was written (the
plan had asserted the opposite). If anyone removes it, this fails loudly rather than letting the
router silently route to a default rule.
"""
from __future__ import annotations

import re
import shutil
import subprocess

import pytest

from dna_decode.data.organism_vocab import (
    BY_CANONICAL,
    UnknownOrganism,
    amrfinder_routable,
    canonical,
    routing_token,
    supported_canonicals,
)

AMRFINDER_IMAGE = "ncbi/amr:4.2.7-2026-03-24.1"


# --------------------------------------------------------------------------- normalization

@pytest.mark.parametrize("spelling", ["Escherichia_coli_Shigella", "Escherichia_coli",
                                      "escherichia_coli", "ESCHERICHIA_COLI", "Escherichia"])
def test_every_observed_ecoli_spelling_collapses_to_one_canonical(spelling):
    """MEASURED in cell_registry: three different E. coli tokens across 10 cells, plus AMRFinder's
    own fourth spelling. All must land on one key."""
    assert canonical(spelling) == "escherichia_coli"


@pytest.mark.parametrize("spelling", ["Klebsiella", "Klebsiella_pneumoniae", "klebsiella"])
def test_both_observed_klebsiella_spellings_collapse_to_one_canonical(spelling):
    """The registry carries `Klebsiella` on 6 AMR cells and `Klebsiella_pneumoniae` on 1 typing cell.
    This inconsistency was found by derivation during execution, not from the plan."""
    assert canonical(spelling) == "klebsiella_pneumoniae"


# ------------------------------------------------- genus-level vs species-level (MEASURED behaviour)

@pytest.mark.parametrize("name,expect", [
    ("Campylobacter_jejuni", "campylobacter"),      # 66 local genomes
    ("Campylobacter_coli", "campylobacter"),        # 34 local genomes
    ("Salmonella_enterica", "salmonella"),          # 60 local genomes
    ("Shigella_flexneri", "escherichia_coli"),      # -O Escherichia covers Shigella
    ("Escherichia coli", "escherichia_coli"),       # space form, as GenBank ORGANISM writes it
])
def test_a_species_in_a_GENUS_level_entry_resolves(name, expect):
    """AMRFinder's `-O Escherichia`/`Campylobacter`/`Salmonella` are GENUS-level values, so any species
    in them routes there. A species-strict reading silently lost 160 real in-set genomes."""
    assert canonical(name) == expect


@pytest.mark.parametrize("other", ["Klebsiella_aerogenes", "Klebsiella_variicola",
                                   "Klebsiella_michiganensis", "Gemmata_obscuriglobus"])
def test_a_different_species_in_a_SPECIES_level_genus_stays_OUT_OF_SET(other):
    """THE DIRECTION THAT MATTERS. AMRFinder lists K. pneumoniae and K. oxytoca separately and has no
    value for aerogenes/variicola/michiganensis, so a Klebsiella genus rule would be WRONG: it would
    answer 'K. pneumoniae' for a different organism and route it to the wrong rule. These 49 local
    genomes are the legitimate out-of-set control for Step 6."""
    with pytest.raises(UnknownOrganism):
        canonical(other)


def test_EXACT_match_beats_the_genus_fallback():
    """Order is load-bearing: K. oxytoca has its own AMRFinder value and must not be captured by any
    broader rule."""
    assert canonical("Klebsiella_oxytoca") == "klebsiella_oxytoca"
    assert routing_token("Klebsiella_oxytoca") == "Klebsiella_oxytoca"
    assert canonical("Klebsiella_pneumoniae") == "klebsiella_pneumoniae"


@pytest.mark.parametrize("bad", ["", "   ", "Homo_sapiens", "not_an_organism",
                                 "Escherchia_coli", "Klebsiella_pneumoniaa"])
def test_an_unsupported_token_RAISES_and_never_defaults(bad):
    """A silent default is how a WRONG organism reaches a rule, and the organism argument selects
    which rule runs.

    Two typo shapes are covered on purpose: `Escherchia_coli` misspells the GENUS (so even the genus
    fallback cannot catch it) and `Klebsiella_pneumoniaa` misspells the species of a SPECIES-level
    entry (whose genus declares no `genus_match`, so there is nothing to fall back to).
    """
    with pytest.raises(UnknownOrganism):
        canonical(bad)


def test_a_typo_in_the_SPECIES_half_of_a_genus_level_entry_DOES_resolve():
    """Recorded as a deliberate consequence, not fuzzy matching.

    This test was originally written the other way round, asserting that `Escherichia_colii` must
    raise. That premise belonged to species-strict matching; once the genus level is honoured --
    because AMRFinder's `-O Escherichia` IS genus-level -- anything in genus Escherichia routes there,
    including a misspelt or novel species. Rewritten rather than deleted so the behaviour is pinned
    and the reversal is visible.

    The cost is bounded and worth naming: an unrecognised Escherichia species resolves to the E. coli
    entry, which is exactly what `-O Escherichia` would do anyway.
    """
    assert canonical("Escherichia_colii") == "escherichia_coli"
    assert canonical("Escherichia_albertii") == "escherichia_coli"


@pytest.mark.parametrize("bad", [None, 42, [], {}])
def test_a_non_string_RAISES_rather_than_crashing_elsewhere(bad):
    with pytest.raises(UnknownOrganism):
        canonical(bad)  # type: ignore[arg-type]


def test_canonical_is_idempotent():
    for c in supported_canonicals():
        assert canonical(c) == c
        assert canonical(canonical(c)) == c


# --------------------------------------------------------------------------- routing token

def test_routing_token_is_the_amrfinder_value():
    assert routing_token("Escherichia_coli_Shigella") == "Escherichia"
    assert routing_token("Klebsiella") == "Klebsiella_pneumoniae"
    assert routing_token("Acinetobacter") == "Acinetobacter_baumannii"


def test_None_and_UNKNOWN_are_DIFFERENT_answers():
    """None = supported but AMRFinder has no `-O` entry (route to that organism's own engine).
    Unknown = not supported at all. Collapsing them would send a fungal genome to AMRFinder."""
    assert routing_token("Candida_auris") is None
    assert routing_token("mycobacterium_tuberculosis") is None
    with pytest.raises(UnknownOrganism):
        routing_token("Saccharomyces_cerevisiae")


def test_the_closed_set_is_non_empty_and_amrfinder_routable_is_a_strict_subset():
    allc, routable = set(supported_canonicals()), set(amrfinder_routable())
    assert routable, "no organism is AMRFinder-routable -- the vocabulary would be useless"
    assert routable < allc, (
        "amrfinder_routable must be a STRICT subset: fungal/TB entries exist precisely to record "
        "organisms AMRFinder cannot route")


# --------------------------------------------------------------------------- non-vacuity

def test_the_vocabulary_is_NOT_a_constant_map():
    """If every organism mapped to one `-O` value the module would be ceremony and should be deleted.
    Mirrors the non-constancy check in tests/test_cell_regime.py."""
    vals = {e.amrfinder_organism for e in BY_CANONICAL.values() if e.amrfinder_organism}
    assert len(vals) >= 2, vals


def test_no_alias_is_claimed_by_two_organisms():
    """Import-time guard; assert it actually held rather than trusting that import succeeded."""
    seen: dict[str, str] = {}
    for c, e in BY_CANONICAL.items():
        for a in {c, *e.registry_labels, *( {e.amrfinder_organism} if e.amrfinder_organism else set())}:
            prev = seen.setdefault(a.lower(), c)
            assert prev == c, f"alias {a!r} shared by {prev!r} and {c!r}"


# --------------------------------------------------------------------------- the load-bearing test

def test_amrfinder_value_and_registry_label_resolve_the_SAME_rule():
    """THE DEPENDENCY THIS MODULE'S DESIGN RESTS ON.

    `amr_rules.calibrated_rule_for` falls back to a GENUS-prefix match, so the AMRFinder `-O` value
    resolves the genus-level registry key. That is why ONE routing token suffices and the plan's
    two-field design was dropped.

    Non-vacuous by construction: it asserts equality for a pair where a rule IS resolved (Klebsiella
    cipro, Acinetobacter meropenem per the live registry keys), so it cannot pass merely because both
    sides returned None.
    """
    from dna_decode.eval.amr_rules import calibrated_rule_for, load_calibrated_registry

    keys = load_calibrated_registry().get("rules", {})
    assert keys, "calibrated registry is empty -- this test would be vacuous"

    resolved_at_least_once = False
    for key in keys:
        reg_org, _, drug = key.partition("|")
        canon = canonical(reg_org)
        amr = BY_CANONICAL[canon].amrfinder_organism
        if amr is None:
            continue
        via_amr = calibrated_rule_for(amr, drug)
        via_reg = calibrated_rule_for(reg_org, drug)
        assert via_amr == via_reg, (
            f"{amr!r} and {reg_org!r} resolve DIFFERENT rules for {drug!r}; the genus-prefix fallback "
            "this module depends on is gone, so one routing token no longer suffices")
        if via_amr is not None:
            resolved_at_least_once = True

    assert resolved_at_least_once, (
        "every comparison returned None -- the test proved nothing. It must exercise at least one "
        "pair where a rule is actually resolved.")


# --------------------------------------------------------------------------- live external surface

def _amrfinder_organism_list() -> set[str]:
    out = subprocess.run(
        ["docker", "run", "--rm", AMRFINDER_IMAGE, "amrfinder", "-l"],
        capture_output=True, text=True, timeout=600)
    m = re.search(r"Available --organism options:\s*(.+)", out.stdout)
    if not m:
        pytest.skip("could not parse `amrfinder -l` output")
    return {t.strip() for t in m.group(1).split(",") if t.strip()}


@pytest.mark.skipif(shutil.which("docker") is None, reason="docker not on PATH")
def test_LIVE_every_amrfinder_value_is_in_the_real_allow_list():
    """The allow-list is PARSED from the pinned image, never hand-copied -- a copied list goes stale
    the moment the AMRFinder DB bumps, and `-O` on an unknown value is a hard error at call time."""
    try:
        allowed = _amrfinder_organism_list()
    except Exception as exc:  # pragma: no cover - environment dependent
        pytest.skip(f"amrfinder image unavailable: {exc}")
    assert len(allowed) > 10, f"suspiciously short allow-list ({len(allowed)}) -- parse likely wrong"
    ours = {e.amrfinder_organism for e in BY_CANONICAL.values() if e.amrfinder_organism}
    missing = sorted(ours - allowed)
    assert not missing, f"not accepted by `amrfinder -O`: {missing}"
