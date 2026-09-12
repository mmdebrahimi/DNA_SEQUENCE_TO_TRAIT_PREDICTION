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
    number, which is exactly how a fabricated cut-off would survive to publication.

    All four entries are now sourced, so this iterates an empty set on the live catalog -- the
    behaviour is pinned on a synthetic unsourced entry below so the guard cannot go vacuous."""
    for drug in unsourced_drugs():
        with pytest.raises(UnsourcedEcoffError):
            ecoff_for(drug)
    import dna_decode.data.ecoff_catalog as mod
    original = mod.ECOFFS
    try:
        mod.ECOFFS = {**ECOFFS, "ciprofloxacin": EcoffEntry("ciprofloxacin", "Escherichia coli")}
        with pytest.raises(UnsourcedEcoffError):
            mod.ecoff_for("ciprofloxacin")
    finally:
        mod.ECOFFS = original


def test_an_unknown_drug_raises_its_own_error():
    with pytest.raises(UnknownDrugError):
        ecoff_for("not-a-drug")
    with pytest.raises(UnknownDrugError):
        entry_for("not-a-drug")


def test_every_catalog_entry_is_now_sourced_with_real_provenance():
    """SUPERSEDES the earlier 'everything is unsourced' pin, which said in its own docstring that it was
    expected to be updated once values were sourced. They were: the EUCAST search IS addressable by
    query parameter (the earlier failure was a wrong species id, not a JS wall), and each agent's ECOFF
    is rendered INTO a PNG at /search/diagram/<id>, which is why no HTML scrape ever found it.

    Every entry must carry a real EUCAST diagram URL and a verbatim quote naming both the cut-off and
    the observation count -- provenance a remembered number could not supply.
    """
    assert set(sourced_drugs()) == set(ECOFFS)
    assert unsourced_drugs() == []
    for drug in ECOFFS:
        e = entry_for(drug)
        assert e.source_url.startswith("https://mic.eucast.org/search/diagram/"), drug
        assert "Epidemiological cut-off (ECOFF)" in e.verbatim_quote, drug
        assert "observations" in e.verbatim_quote, drug
        assert "2026-09-12" in (e.eucast_version or ""), drug


def test_the_ceftriaxone_anchor_is_flagged_tentative_not_silently_equal_to_the_others():
    """EUCAST prints a TECOFF in PARENTHESES -- set on 3-4 distributions instead of the >=5 an ECOFF
    needs. Ceftriaxone's rests on 908 observations from 4 sources against gentamicin's 78,136 from 82,
    with a CI spanning four doublings. Storing it as an equal peer would let a ceftriaxone result be
    read with the same confidence as a gentamicin one."""
    e = entry_for("ceftriaxone")
    assert "(0.125)" in e.verbatim_quote, "the parenthesised TECOFF marker must survive into the quote"
    assert "TENTATIVE" in e.note and "4 data sources" in (e.note + e.verbatim_quote)
    assert entry_for("gentamicin").note and "TENTATIVE" not in entry_for("gentamicin").note


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
    # Now that every entry is sourced, normalization is shown by the value coming BACK rather than by
    # the refusal that followed it while the catalog was empty.
    assert ecoff_for("GENTAMICIN") == 2.0
    assert ecoff_for("  Tetracycline  ") == 8.0
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


def test_the_sourced_partition_tracks_the_catalog_rather_than_being_hardcoded():
    """NON-VACUITY for the sourced/unsourced split: an implementation that always reported everything
    as sourced would pass the live-state test above, and one that always reported nothing would have
    passed the earlier all-unsourced pin. Both directions are exercised on a patched catalog, with a
    deliberately absurd value (1234 mg/L) that can never be mistaken for a real ECOFF."""
    import dna_decode.data.ecoff_catalog as mod
    original = mod.ECOFFS
    try:
        # one entry knocked back to unsourced -> it must leave `sourced_drugs`
        mod.ECOFFS = {**ECOFFS, "ciprofloxacin": EcoffEntry("ciprofloxacin", "Escherichia coli")}
        assert "ciprofloxacin" in mod.unsourced_drugs()
        assert "ciprofloxacin" not in mod.sourced_drugs()
        assert "gentamicin" in mod.sourced_drugs()
        # and a fully-provenanced synthetic entry must come back
        mod.ECOFFS = {**ECOFFS, "ciprofloxacin": EcoffEntry(
            "ciprofloxacin", "Escherichia coli", value=1234, status=SOURCED,
            source_url="https://example.invalid/not-a-real-ecoff",
            verbatim_quote="synthetic test fixture -- NOT a EUCAST value")}
        got = mod.ecoff_for("ciprofloxacin")
        assert got == 1234.0 and isinstance(got, float)
        assert mod.validate_catalog() == []
    finally:
        mod.ECOFFS = original
    assert ecoff_for("ciprofloxacin") == 0.06, "live catalog must be restored"


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
def test_the_sourcing_artifact_records_the_resolution_and_keeps_the_attempt_log():
    """The artifact began as a BLOCKED record and was RESOLVED the same day. Both halves must survive:
    the attempt log is the evidence the wall was really probed, and the resolution is what corrects it.

    An artifact still claiming BLOCKED after the values were sourced would assert a wall that no longer
    exists -- the stale-claim failure this project keeps paying for.
    """
    d = json.loads(SOURCES.read_text(encoding="utf-8"))
    assert d["status"].startswith("RESOLVED")
    assert len(d["attempts"]) >= 4, "a one-attempt wall is not a demonstrated wall"
    assert d["frozen_surface_untouched"] is True
    joined = json.dumps(d).lower()
    assert "memory" in joined, "the artifact must state that no value was recalled"

    res = d["resolution"]
    assert set(res["values"]) == set(ECOFFS)
    for drug, v in res["values"].items():
        assert v["url"].startswith("https://mic.eucast.org/search/diagram/")
        assert v["ecoff_mg_L"] == ecoff_for(drug), f"{drug}: artifact and catalog disagree"
    assert res["values"]["ceftriaxone"]["tentative"] is True
    # the superseded fields must SAY they are superseded rather than quietly contradicting the outcome
    assert d["consequence"].startswith("[SUPERSEDED")
    assert "SUPERSEDED" in d["attempts"][0]["outcome"]


# --- the feasibility gate (plan Step 2) ----------------------------------------------------------

def test_degeneracy_is_checked_before_power():
    """An anchor that cannot discriminate AT ALL is not merely underpowered. Reporting a degenerate
    drug as UNDERPOWERED would imply a bigger cohort would fix it; nothing about cohort size moves an
    ECOFF that sits below the panel floor."""
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "ecoff_feas", ROOT / "scripts" / "ecoff_tiering_feasibility.py")
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    # below the floor AND a tiny stratum -> degeneracy wins
    assert m.verdict_for(0.06, [0.125, 0.25, 8.0], stratum_n=0) == m.DEGENERATE
    assert m.verdict_for(0.06, [0.125, 0.25, 8.0], stratum_n=5000) == m.DEGENERATE
    # inside the panel, stratum decides
    assert m.verdict_for(2.0, [1.0, 2.0, 4.0, 32.0], stratum_n=25) == m.SCOREABLE
    assert m.verdict_for(2.0, [1.0, 2.0, 4.0, 32.0], stratum_n=19) == m.UNDERPOWERED
    assert m.verdict_for(2.0, [1.0, 2.0, 4.0, 32.0], stratum_n=20) == m.SCOREABLE


def test_the_committed_feasibility_artifact_matches_its_own_verdict_rule():
    """Re-derives every verdict from the artifact's own numbers, so a hand-edited verdict cannot stand."""
    import importlib.util
    path = ROOT / "wiki" / "ecoff_tiering_feasibility_2026-09-11.json"
    if not path.exists():
        pytest.skip("feasibility artifact not generated")
    spec = importlib.util.spec_from_file_location(
        "ecoff_feas", ROOT / "scripts" / "ecoff_tiering_feasibility.py")
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    doc = json.loads(path.read_text(encoding="utf-8"))
    for r in doc["results"]:
        if r["verdict"] in ("NOT_IN_COHORT", m.UNSOURCED):
            continue
        assert r["verdict"] == m.verdict_for(
            r["ecoff_mg_L"], r["measured_grid"], r["discriminating_stratum_n"]), r["drug"]
    assert doc["scoreable_drugs"] == ["gentamicin"]


def test_a_degenerate_drug_really_is_all_non_wildtype():
    """The MEANING of degenerate, pinned: it is not a threshold artefact, it is that every single
    isolate lands above the anchor so the two classes cannot both exist."""
    path = ROOT / "wiki" / "ecoff_tiering_feasibility_2026-09-11.json"
    if not path.exists():
        pytest.skip("feasibility artifact not generated")
    doc = json.loads(path.read_text(encoding="utf-8"))
    degenerate = [r for r in doc["results"] if r["verdict"] == "DEGENERATE_ECOFF_BELOW_PANEL"]
    assert degenerate, "fixture assumption: at least one drug is degenerate on this cohort"
    for r in degenerate:
        assert r["non_wildtype_fraction"] == 1.0, r["drug"]
        assert r["ecoff_mg_L"] < r["panel_floor"], r["drug"]


# --- the evaluation (plan Step 5) ----------------------------------------------------------------

def _eval_mod():
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "ecoff_eval", ROOT / "scripts" / "ecoff_tiering_evaluate.py")
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m


def test_all_four_frozen_verdict_branches():
    m = _eval_mod()
    strong, weak = {"exceeds_perm_max": True}, {"exceeds_perm_max": False}
    assert m.verdict_from_bar(0.30, 0.01, strong, 25) == m.SUPPORTED
    assert m.verdict_from_bar(0.08, 0.005, weak, 25) == m.WEAK
    assert m.verdict_from_bar(0.01, 0.05, strong, 25) == m.FALSIFIED
    assert m.verdict_from_bar(0.30, 0.01, strong, 19) == m.INDETERMINATE


def test_an_equal_rate_is_falsified_not_weak():
    """The predicted ORDERING is strict. Equal carriage means the ECOFF separates nothing, which is a
    failure of the claim rather than a weak version of it."""
    m = _eval_mod()
    assert m.verdict_from_bar(0.05, 0.05, {"exceeds_perm_max": True}, 100) == m.FALSIFIED


def test_a_thin_stratum_is_indeterminate_even_with_a_huge_gap():
    """NON-VACUITY for the floor: an n=3 stratum showing 100%-vs-0% must not read as SUPPORTED."""
    m = _eval_mod()
    assert m.verdict_from_bar(1.0, 0.0, {"exceeds_perm_max": True}, 3) == m.INDETERMINATE


def test_the_permutation_null_holds_stratum_sizes_fixed():
    """The control must shuffle the LABEL, not resample -- otherwise it would not be a null for the
    anchor assignment, which is the only thing under test."""
    import numpy as np
    m = _eval_mod()
    labels = np.concatenate([np.ones(25, int), np.zeros(2652, int)])
    carries = np.concatenate([np.ones(2, int), np.zeros(23, int),
                              np.ones(12, int), np.zeros(2640, int)])
    perm = m.permutation_gap(labels, carries, np.random.default_rng(0), n_perm=200)
    assert perm["perm_max"] > 0 and perm["n_perm"] == 200
    # a label with no association must not clear its own null
    flat = m.permutation_gap(labels, np.zeros(2677, int), np.random.default_rng(0), n_perm=200)
    assert not flat["exceeds_perm_max"]


def test_the_determinant_rule_excludes_streptomycin_determinants():
    """THE definition that decides this result. 3,298 of the cohort's aminoglycoside rows are Subclass
    STREPTOMYCIN; counting them would produce a streptomycin co-carriage finding wearing a gentamicin
    label. Compound subclasses containing GENTAMICIN must still count (substring, as the frozen rule
    does it)."""
    import pandas as pd
    m = _eval_mod()
    amr = pd.DataFrame({
        "Name": ["g1", "g2", "g3", "g4"],
        "Class": ["AMINOGLYCOSIDE"] * 4,
        "Subclass": ["STREPTOMYCIN", "GENTAMICIN", "GENTAMICIN/KANAMYCIN/TOBRAMYCIN", "KANAMYCIN"],
        "Gene symbol": ["aph(6)-Id", "aac(3)-IId", "aac(3)-IIe", "aph(3')-Ia"]})
    got = m.carriers_for_drug(amr, "gentamicin")
    assert got == {"g2", "g3"}, got


@pytest.mark.skipif(not (ROOT / "wiki" / "ecoff_tiering_result_2026-09-12.json").exists(),
                    reason="evaluation artifact not generated")
def test_the_committed_verdict_is_re_derivable_from_its_own_numbers():
    m = _eval_mod()
    d = json.loads((ROOT / "wiki" / "ecoff_tiering_result_2026-09-12.json").read_text(encoding="utf-8"))
    s = d["strata"]
    assert d["verdict"] == m.verdict_from_bar(
        s["clinically_S_and_NWT"]["carriage_rate"], s["clinically_S_and_WT"]["carriage_rate"],
        d["permutation"], s["clinically_S_and_NWT"]["n"])
    assert d["verdict"] == m.WEAK


@pytest.mark.skipif(not (ROOT / "wiki" / "ecoff_tiering_result_2026-09-12.json").exists(),
                    reason="evaluation artifact not generated")
def test_the_deployed_rule_is_sane_on_the_resistant_stratum():
    """A control on the MEASUREMENT itself: if the rule did not recover the clinically-resistant
    isolates, a null result in the susceptible strata would say nothing about the anchor."""
    d = json.loads((ROOT / "wiki" / "ecoff_tiering_result_2026-09-12.json").read_text(encoding="utf-8"))
    assert d["strata"]["clinically_R"]["carriage_rate"] > 0.85
