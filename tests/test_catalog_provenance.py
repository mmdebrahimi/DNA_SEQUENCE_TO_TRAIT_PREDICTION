"""Every shipped catalog names an external authority — the enforceable half of the curation procedure.

WHAT THIS RESOLVES. The curation family carried an open unknown: *whether a curation procedure can be
enforced by test at all, or is irreducibly a review discipline.* Measured 2026-09-01: **partly**.

Per-ENTRY citation is NOT representable. The shipped catalogs are bare collections
(`NNRTI_RT_MAJOR_DRMS: set[str] = {"L100I", ...}`), so there is nowhere to hang a per-entry source
without restructuring the frozen surface -- a large, risky change for no measured benefit.

Per-MODULE authority IS representable, and every catalog already satisfies it. That makes this a
REGRESSION guard, not a migration: it fires when a NEW catalog ships without naming where its facts came
from, which is the fabrication hazard the colour-cell memo names.

Honest limit: naming an authority in a module is much weaker than sourcing each entry to it. This cannot
catch a fabricated entry inside a correctly-cited module. It catches an UNSOURCED CATALOG.
"""
from __future__ import annotations

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]

# Recognised external authorities. Adding one is a deliberate act -- it widens what counts as sourced.
AUTHORITIES = ("Stanford", "HIVDB", "WHO", "CDC", "CLSI", "EUCAST", "NCBI", "AMRFinder", "CARD",
               "CoV-RDB", "Chou", "Rhee", "Margot", "Napier", "OMIA", "CPIC", "PharmVar",
               "VirulenceFinder", "PubMLST", "ResFinder", "EnteroBase", "CAPELLA", "Orth", "BiGG")

# Shipped fact catalogs. A rename fails loudly (test below) rather than silently skipping a module.
CATALOG_MODULES = (
    "dna_decode/data/hiv_amr.py",
    "dna_decode/data/sarscov2_amr.py",
    "dna_decode/data/fungal_amr.py",
    "dna_decode/data/hcmv_amr.py",
    "dna_decode/data/mic_tiers.py",
    "dna_decode/data/tb_who_catalogue.py",
    "dna_decode/eval/amr_rules.py",
    "dna_decode/organism_rules/tb_amr.py",
)

_AUTH = re.compile("|".join(re.escape(a) for a in AUTHORITIES), re.IGNORECASE)


def names_an_authority(text: str) -> bool:
    return bool(_AUTH.search(text))


@pytest.mark.parametrize("rel", CATALOG_MODULES)
def test_every_shipped_catalog_names_an_external_authority(rel):
    """A catalog that cites nothing is indistinguishable from one written from memory."""
    text = (ROOT / rel).read_text(encoding="utf-8")
    assert names_an_authority(text), (
        f"{rel} ships curated biological facts but names no external authority from {AUTHORITIES}. "
        "Cite the source, or do not write the entry.")


@pytest.mark.parametrize("rel", CATALOG_MODULES)
def test_every_listed_catalog_module_exists(rel):
    """The list is hand-maintained; a rename must fail here rather than silently drop coverage."""
    assert (ROOT / rel).is_file(), f"{rel} is listed as a catalog but does not exist"


def test_the_authority_check_is_not_vacuous():
    """A predicate that passes any text proves nothing about the modules it just cleared."""
    assert not names_an_authority(
        "MY_DRMS = {'K103N', 'Y181C'}  # these are the important ones as far as I remember")
    assert names_an_authority("sourced verbatim from the Stanford HIVDB dataset page")


def test_the_procedure_document_exists_and_states_the_measured_limit():
    """The other half of the family's deliverable: what a test CANNOT enforce must be written down,
    or the green suite reads as more assurance than it is."""
    doc = ROOT / "wiki" / "catalog_curation_procedure.md"
    assert doc.is_file(), "wiki/catalog_curation_procedure.md missing"
    text = doc.read_text(encoding="utf-8")
    for required in ("per-entry", "review discipline", "doubt layer", "lock"):
        assert required.lower() in text.lower(), f"procedure does not address {required!r}"
