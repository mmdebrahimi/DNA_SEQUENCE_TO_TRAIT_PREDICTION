"""Guards for the ECOFF catalog and the wild-type classifier.

The catalog's whole job is to REFUSE. An ECOFF is a biological reference value, and a wrong one
produces a fully self-consistent evaluation that is entirely wrong with nothing downstream able to
catch it -- so these tests care much more about the refusal paths than about any value.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from dna_decode.data.ecoff_catalog import (  # noqa: E402
    ECOFFS, SOURCED, EcoffEntry, UnknownDrugError, UnsourcedEcoffError,
    ecoff_for, entry_for, sourced_drugs, unsourced_drugs, validate_catalog,
)
from dna_decode.eval.wildtype_tiering import (  # noqa: E402
    NO_MIC, NWT, WT, classify_wildtype, ecoff_is_resolvable_on_grid, wildtype_counts,
)

SOURCES = ROOT / "wiki" / "ecoff_sources_2026-09-11.json"


# --- the catalog refuses rather than defaulting ---------------------------------------------------

def test_an_unsourced_ecoff_raises_instead_of_returning_a_plausible_default():
    """THE guard. A default here would be indistinguishable from a real value in every downstream
    number, which is exactly how a fabricated cut-off would survive to publication."""
    for drug in unsourced_drugs():
        with pytest.raises(UnsourcedEcoffError):
            ecoff_for(drug)


def test_an_unknown_drug_raises_its_own_error():
    with pytest.raises(UnknownDrugError):
        ecoff_for("not-a-drug")
    with pytest.raises(UnknownDrugError):
        entry_for("not-a-drug")


def test_every_catalog_entry_is_currently_unsourced_and_says_why():
    """Pins the 2026-09-11 state honestly: EUCAST exposes no downloadable table and no stable result
    URL, so nothing could be sourced from this environment. If a value is later added WITH provenance
    this test is expected to be updated -- deliberately, by a human who did the sourcing."""
    assert sourced_drugs() == []
    assert set(unsourced_drugs()) == set(ECOFFS)
    for drug in ECOFFS:
        assert entry_for(drug).note, f"{drug} gives no reason for being unsourced"


def test_a_value_without_provenance_is_a_catalog_violation():
    """NON-VACUITY: this is the shape a remembered number takes when pasted in later. The live catalog
    is clean, so the check is exercised against a deliberately malformed entry."""
    assert validate_catalog() == []
    bad = dict(ECOFFS)
    bad["ciprofloxacin"] = EcoffEntry("ciprofloxacin", "Escherichia coli", value=0.06,
                                      status=SOURCED)          # no url, no quote
    import dna_decode.data.ecoff_catalog as mod
    original = mod.ECOFFS
    try:
        mod.ECOFFS = bad
        problems = mod.validate_catalog()
        assert any("source_url" in p for p in problems)
        assert any("verbatim_quote" in p for p in problems)
    finally:
        mod.ECOFFS = original
    assert validate_catalog() == []


def test_a_value_with_full_provenance_is_usable_and_returns():
    """The positive path, so the refusal is not vacuous -- the catalog must actually work once filled."""
    e = EcoffEntry("gentamicin", "Escherichia coli", value=2.0, status=SOURCED,
                   source_url="https://mic.eucast.org/...", verbatim_quote="ECOFF 2 mg/L",
                   eucast_version="pinned-by-human", retrieved="2026-09-11")
    assert e.is_usable()
    assert not EcoffEntry("gentamicin", "Escherichia coli", value=2.0, status=SOURCED).is_usable()


def test_drug_lookup_normalizes_case_and_whitespace():
    """A caller passing a drug straight out of a cohort TSV must reach the same entry as one passing
    the catalog key. Without normalization a stray space raises `UnknownDrugError` on a drug the
    catalog DOES have -- a refusal that looks identical to the honest one and would be read as
    'unsourced' rather than as a lookup bug."""
    assert entry_for("  Ciprofloxacin ") is ECOFFS["ciprofloxacin"]
    with pytest.raises(UnsourcedEcoffError):
        ecoff_for("GENTAMICIN")          # normalized, found, and THEN refused for lack of provenance
    for empty in ("", "   ", None):
        with pytest.raises(UnknownDrugError):
            entry_for(empty)


def test_validate_catalog_also_catches_a_sourced_hole_and_a_non_positive_value():
    """The other two violation shapes, neither of which the live catalog can exercise. `status=SOURCED`
    with no value is a half-filled entry that `is_usable` would quietly reject rather than report; a
    value of 0 is what a blank or failed parse becomes, and as an ECOFF it makes EVERY isolate
    non-wild-type."""
    bad = dict(ECOFFS)
    bad["gentamicin"] = EcoffEntry("gentamicin", "Escherichia coli", status=SOURCED)
    bad["tetracycline"] = EcoffEntry("tetracycline", "Escherichia coli", value=0.0, status=SOURCED,
                                     source_url="https://example.invalid/", verbatim_quote="blank")
    import dna_decode.data.ecoff_catalog as mod
    original = mod.ECOFFS
    try:
        mod.ECOFFS = bad
        problems = mod.validate_catalog()
        assert any("gentamicin" in p and "no value" in p for p in problems)
        assert any("tetracycline" in p and "non-positive" in p for p in problems)
    finally:
        mod.ECOFFS = original
    assert validate_catalog() == []


def test_a_fully_sourced_entry_is_returned_and_listed_as_sourced():
    """NON-VACUITY for `sourced_drugs() == []`: an implementation that always returned an empty list,
    or an `ecoff_for` that raised unconditionally, would pass every refusal test above. The value here
    is deliberately absurd (1234 mg/L) so it can never be mistaken for a sourced ECOFF."""
    import dna_decode.data.ecoff_catalog as mod
    filled = dict(ECOFFS)
    filled["ciprofloxacin"] = EcoffEntry(
        "ciprofloxacin", "Escherichia coli", value=1234, status=SOURCED,
        source_url="https://example.invalid/not-a-real-ecoff",
        verbatim_quote="synthetic test fixture -- NOT a EUCAST value")
    original = mod.ECOFFS
    try:
        mod.ECOFFS = filled
        got = mod.ecoff_for("ciprofloxacin")
        assert got == 1234.0 and isinstance(got, float)
        assert mod.sourced_drugs() == ["ciprofloxacin"]
        assert "ciprofloxacin" not in mod.unsourced_drugs()
        assert mod.validate_catalog() == []
    finally:
        mod.ECOFFS = original
    assert sourced_drugs() == []


# --- the wild-type classifier ---------------------------------------------------------------------

def test_at_the_ecoff_is_wild_type_and_one_dilution_above_is_not():
    """EUCAST defines the ECOFF as the HIGHEST MIC of the wild-type distribution, so the boundary
    value belongs to the wild type. An off-by-one reclassifies an entire dilution."""
    assert classify_wildtype([2.0], 2.0) == WT
    assert classify_wildtype([4.0], 2.0) == NWT
    assert classify_wildtype([1.0], 2.0) == WT


def test_no_numeric_mic_is_its_own_answer_not_wild_type():
    """Defaulting a missing MIC to WT would silently manufacture susceptible isolates."""
    assert classify_wildtype([], 2.0) == NO_MIC
    assert classify_wildtype([None, float("nan")], 2.0) == NO_MIC


def test_the_median_convention_matches_the_clinical_arm():
    """Both arms must differ ONLY in the anchor, so the aggregation has to be identical."""
    assert classify_wildtype([1.0, 1.0, 32.0], 2.0) == WT      # median 1
    assert classify_wildtype([1.0, 32.0, 32.0], 2.0) == NWT    # median 32


def test_a_nonsense_ecoff_raises():
    for bad in (0, -1, None, float("nan")):
        with pytest.raises(ValueError):
            classify_wildtype([1.0], bad)


def test_counts_roll_up_every_isolate_exactly_once():
    counts = wildtype_counts({"a": [1.0], "b": [32.0], "c": [], "d": [4.0]}, 2.0)
    assert counts == {WT: 1, NWT: 2, NO_MIC: 1}
    assert sum(counts.values()) == 4


def test_an_ecoff_below_the_panel_floor_is_flagged_as_unresolvable():
    """The live risk the plan named in advance: Oxford's ciprofloxacin panel bottoms out at 0.125, so
    a lower ECOFF makes EVERY isolate non-wild-type and the anchor discriminates nothing."""
    oxford_cipro_grid = [0.125, 0.25, 0.5, 8.0]
    assert not ecoff_is_resolvable_on_grid(0.06, oxford_cipro_grid)
    assert not ecoff_is_resolvable_on_grid(0.03, oxford_cipro_grid)
    assert ecoff_is_resolvable_on_grid(0.25, oxford_cipro_grid)
    # gentamicin's grid does leave room
    assert ecoff_is_resolvable_on_grid(2.0, [1.0, 2.0, 4.0, 32.0])


def test_valid_mics_survive_alongside_missing_ones():
    """The all-missing case is pinned above; the MIXED case is the one that can go wrong silently. A
    filter that dropped the whole list on one `None` would return NO_MIC for a measured isolate, and a
    filter that kept the `None`s would shift the median."""
    assert classify_wildtype([None, 4.0, float("nan")], 2.0) == NWT
    assert classify_wildtype([None, 1.0], 2.0) == WT
    assert classify_wildtype([None, 1.0, 1.0, 32.0], 2.0) == WT   # median of the VALID three is 1
    assert classify_wildtype(None, 2.0) == NO_MIC


def test_an_even_length_mic_list_interpolates_to_a_value_off_the_dilution_grid():
    """`statistics.median` averages the middle pair, so two replicates straddling the ECOFF produce a
    midpoint that is not a real doubling-dilution MIC. Pinned rather than assumed, because it decides
    the call for every duplicate-tested isolate."""
    assert classify_wildtype([1.0, 32.0], 2.0) == NWT    # median 16.5
    assert classify_wildtype([1.0, 2.0], 2.0) == WT      # median 1.5


def test_resolvability_fails_closed_on_an_empty_panel_and_ignores_missing_values():
    """No measured MICs means nothing is known about the grid, so the honest answer is 'not
    resolvable' -- returning True there would let a degenerate anchor through on an empty cohort."""
    assert not ecoff_is_resolvable_on_grid(2.0, [])
    assert not ecoff_is_resolvable_on_grid(2.0, [None, float("nan")])
    assert ecoff_is_resolvable_on_grid(2.0, [None, 1.0, float("nan"), 4.0])


def test_an_ecoff_exactly_at_the_panel_floor_still_discriminates():
    """BOUNDARY PIN. At the floor, isolates reported at the floor are WT and everything above is NWT,
    so the anchor still separates; only an ECOFF strictly BELOW the floor makes every isolate NWT.
    That is `ecoff >= min(vals)`, i.e. the degenerate zone is exclusive of the floor -- note the
    function's docstring says 'at or below', which reads as inclusive and does not match."""
    grid = [0.125, 0.25, 8.0]
    assert ecoff_is_resolvable_on_grid(0.125, grid)
    assert not ecoff_is_resolvable_on_grid(0.124, grid)


def test_wildtype_tiering_does_not_import_the_frozen_module():
    """The arm must be evaluable without touching `mic_tiers.py`, which is sha256-pinned in the v2
    prospective lock."""
    src = (ROOT / "dna_decode" / "eval" / "wildtype_tiering.py").read_text(encoding="utf-8")
    assert "mic_tiers" not in src.replace("`mic_tiers.classify_tier`", "").replace(
        "mic_tiers.classify_tier", "").replace("`mic_tiers.py`", "")


# --- the blocked-sourcing record ------------------------------------------------------------------

@pytest.mark.skipif(not SOURCES.exists(), reason="sources artifact not present")
def test_the_sourcing_artifact_records_the_block_and_the_attempts():
    d = json.loads(SOURCES.read_text(encoding="utf-8"))
    assert d["status"].startswith("BLOCKED")
    assert len(d["attempts"]) >= 4, "a one-attempt wall is not a demonstrated wall"
    assert d["frozen_surface_untouched"] is True
    joined = json.dumps(d).lower()
    assert "memory" in joined, "the artifact must state that no value was recalled"
