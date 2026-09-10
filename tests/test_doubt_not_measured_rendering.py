"""The doubt layer's four non-clean states must all REACH A HUMAN, and only one silence is honest.

Two defects on 2026-09-10, both on the shipped CLI surface:

  1. NOT-MEASURED rendered as SILENCE. `doubt_one_line` knew `applicable is False` and
     `assessed is False` but not `measured is False`, so every UNMEASURED cell fell through to the
     honest-silence return -- sarscov2-mpro, fungal-fluconazole-erg11 and fungal-voriconazole-erg11
     since 2026-09-02. The completeness signal said "has NOT been measured" in the record while the
     human output said nothing at all, which is the exact failure `doubt_one_line`'s own docstring
     exists to prevent.
  2. An unregistered drug was told "this catalog is position-based" -- a claim about the cell the
     branch cannot verify. For lenacapavir it was FALSE: the CAI catalog is mutant-level, so the CLI
     printed "position-based" two lines above that call's own "MUTANT-LEVEL v0" caveat.

Silence must mean "we checked and found nothing". For these cells nobody checked.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from dna_decode.data.target_site_completeness import UNMEASURED_CELLS  # noqa: E402
from dna_decode.eval.doubt import (  # noqa: E402
    _MUTANT_LEVEL_CELLS, doubt_one_line, target_site_doubt,
)
from dna_decode.eval.position_novelty import KNOWN_CELLS  # noqa: E402

# (drug, gene, catalogued substitution) for each registered mutant-level cell.
REGISTERED = [("efavirenz", "RT", "K103N"), ("nirmatrelvir", "Mpro", "E166V"),
              ("fluconazole", "ERG11", "Y132F"), ("voriconazole", "ERG11", "Y132F"),
              ("lenacapavir", "CA", "M66I")]


def test_a_cell_registered_for_doubt_is_also_registered_for_position_novelty():
    """Registering a cell in ONE map only is a live crash, not a soft miss: routing `hiv-cai` to
    position-novelty before it was known there killed a real lenacapavir call with a bare KeyError."""
    missing = sorted(set(_MUTANT_LEVEL_CELLS.values()) - set(KNOWN_CELLS))
    assert not missing, f"cells in the doubt map but unknown to position_novelty: {missing}"


@pytest.mark.parametrize("drug,gene,sub", REGISTERED)
def test_every_registered_cell_survives_a_real_call(drug, gene, sub):
    """Non-vacuity for the test above: exercising each cell end-to-end is what surfaces a half
    registration, since the parity check alone cannot see a catalog that raises on lookup."""
    block = target_site_doubt(drug, {gene: {sub}}).as_dict()
    assert block["signals"], f"{drug} produced no signals"


@pytest.mark.parametrize("drug,gene,sub", [r for r in REGISTERED if r[0] != "efavirenz"])
def test_an_unmeasured_cell_says_so_out_loud(drug, gene, sub):
    """These cells have no free isolate-level phenotype source (or a TN-starved one). Rendering that
    as silence would tell a reader the call was screened and clean."""
    line = doubt_one_line(target_site_doubt(drug, {gene: {sub}}).as_dict())
    assert line is not None, f"{drug} renders NO line despite being unmeasured"
    assert "NOT MEASURED" in line
    assert "NOT of doubt" in line, "the line must say an absence of measurement is not an absence of doubt"


def test_the_only_honest_silence_is_assessed_and_quiet():
    """efavirenz IS measured, and K103N is catalogued, so 'we checked and found nothing' is exactly
    what silence truthfully means. If this ever returns a line the renderer has become noise."""
    assert doubt_one_line(target_site_doubt("efavirenz", {"RT": {"K103N"}}).as_dict()) is None


def test_a_measured_gap_still_fires_strong():
    """Guards against fixing the silence by muting the signal that matters."""
    line = doubt_one_line(target_site_doubt("efavirenz", {"RT": {"V179F"}}).as_dict())
    assert line and line.startswith("DOUBT [strong]") and "V179F" in line


def test_an_unregistered_drug_is_not_told_its_catalog_is_position_based():
    """The branch cannot see the catalog's shape. Asserting one printed a falsehood for lenacapavir,
    whose catalog is mutant-level; it may only report that the drug is not registered."""
    line = doubt_one_line(target_site_doubt("raltegravir", {"IN": {"Q148H"}}).as_dict())
    assert line is not None, "an unregistered drug must not render as silence either"
    assert "NOT SCREENED" in line
    assert "position-based" not in line, "the renderer asserts a catalog shape it cannot know"


def test_lenacapavir_is_registered_because_its_catalog_is_mutant_level():
    """The concrete instance of defect 2. CAI ships the CAPELLA emergent-substitution set precisely
    BECAUSE capsid polymorphisms defeat a position-based rule, so it is mutant-level by construction."""
    from dna_decode.data.hiv_amr import _HIV_TARGET_CLASSES
    cai = next(c for c in _HIV_TARGET_CLASSES if c.label == "CAI")
    assert cai.major_drms, "CAI is expected to be mutant-level; if it became position-based, unregister it"
    assert _MUTANT_LEVEL_CELLS.get("lenacapavir") == "hiv-cai"
    assert "hiv-cai" in UNMEASURED_CELLS, (
        "hiv-cai has no labelled negative class (CAI validation is 129R/11S), so it must be declared "
        "unmeasured rather than rendering measured-and-empty")


def test_doubt_never_alters_the_call():
    """L2 qualifies, never overrules. Pinned here because this change touched the rendering path that
    sits beside the call in the same output."""
    from dna_decode.data.hiv_amr import call_hiv_observed
    before = call_hiv_observed("lenacapavir", {"CA": {"M66I"}})
    target_site_doubt("lenacapavir", {"CA": {"M66I"}}).as_dict()
    after = call_hiv_observed("lenacapavir", {"CA": {"M66I"}})
    assert before.prediction == after.prediction == "R"
