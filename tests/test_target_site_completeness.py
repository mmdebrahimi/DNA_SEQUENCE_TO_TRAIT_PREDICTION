"""The target-site completeness signal: complementary to position-novelty, and augment-only."""
from __future__ import annotations

import pytest

from dna_decode.data.target_site_completeness import (
    TARGET_SITE_COMPLETENESS, UNMEASURED_CELLS, completeness_units_for, is_measured, matching_units,
)
from dna_decode.eval.doubt import NONE, STRONG, doubt_one_line, target_site_doubt
from dna_decode.eval.position_novelty import flag_for_cell


# --- the two signals cover DIFFERENT blind spots; neither subsumes the other --------------------

def test_position_novelty_is_silent_on_a_non_catalogued_position():
    """The entire reason this signal exists. If position-novelty caught V179F, the index would be
    redundant and should be deleted."""
    assert flag_for_cell(["V179F"], "hiv-nnrti-rt").position_novel is False
    assert flag_for_cell(["K103R"], "hiv-nnrti-rt").position_novel is True


def test_the_completeness_signal_fires_where_position_novelty_cannot():
    sigs = {s.kind: s for s in target_site_doubt("efavirenz", {"RT": {"V179F"}}).signals}
    assert sigs["position_novelty"].tier == NONE
    assert sigs["target_site_completeness"].tier == STRONG


def test_position_novelty_still_fires_on_its_own_shape():
    """Regression: adding a signal must not disturb the incumbent."""
    sigs = {s.kind: s for s in target_site_doubt("efavirenz", {"RT": {"K103R"}}).signals}
    assert sigs["position_novelty"].tier != NONE
    assert sigs["target_site_completeness"].tier == NONE


# --- the human-facing line must name the signal that actually fired ------------------------------

def test_the_printed_line_reports_the_reason_of_the_signal_that_fired():
    """It hardcoded signals[0], so a STRONG tier printed alongside the OTHER signal's 'found nothing'
    prose -- a self-contradicting disclosure, the worst failure mode for a human-facing line."""
    block = target_site_doubt("efavirenz", {"RT": {"V179F"}}).as_dict()
    line = doubt_one_line(block)
    assert block["max_tier"] == STRONG
    assert "V179F" in line and "does not carry" in line
    assert "no uncatalogued substitution" not in line


# --- three states, never collapsed ---------------------------------------------------------------

def test_a_never_measured_cell_says_so_rather_than_reporting_clean():
    sigs = {s.kind: s for s in target_site_doubt("nirmatrelvir", {"Mpro": {"E166V"}}).signals}
    s = sigs["target_site_completeness"]
    assert s.tier == NONE and s.evidence["measured"] is False
    assert "NOT an absence of doubt" in s.reason


def test_measured_and_quiet_is_distinguishable_from_never_measured():
    sigs = {s.kind: s for s in target_site_doubt("efavirenz", {"RT": {"K101P"}}).signals}
    s = sigs["target_site_completeness"]
    assert s.tier == NONE and s.evidence["measured"] is True


def test_an_unsurfaced_path_mentions_both_screens():
    """The not-assessable message named only position-novelty while two screens now exist."""
    sig = target_site_doubt("efavirenz", None).signals[0]
    assert "completeness" in sig.reason and "position-novelty" in sig.reason


# --- the index's entry bar ------------------------------------------------------------------------

def test_every_listed_unit_is_pure_and_survives_familywise_correction():
    for cell, units in TARGET_SITE_COMPLETENESS.items():
        for sub, st in units.items():
            assert st["carriers_labelled_s"] == 0, f"{sub} has a susceptible carrier -- signal ENDS"
            assert st["carriers_labelled_r"] >= 5, f"{sub} underpowered"
            assert st["purity_surprise_p"] <= 0.05 / st["n_units_tested"], f"{sub} fails correction"
            assert st["artifact"].endswith(".json") and st["label"]


def test_a_listed_unit_is_not_already_representable_by_the_catalog():
    """A gap the catalog already covers is not a gap."""
    for cell, units in TARGET_SITE_COMPLETENESS.items():
        for sub in units:
            assert flag_for_cell([sub], cell).position_novel is False


def test_unmeasured_cells_are_declared_and_disjoint_from_measured_ones():
    assert UNMEASURED_CELLS
    assert not (UNMEASURED_CELLS & set(TARGET_SITE_COMPLETENESS))


def test_matching_is_case_insensitive_and_stable():
    assert [s for s, _ in matching_units(["v179f"], "hiv-nnrti-rt")] == ["V179F"]
    assert matching_units([], "hiv-nnrti-rt") == []
    assert is_measured("hiv-nnrti-rt") and not is_measured("nope")
    assert completeness_units_for("nope") == {}


def test_the_block_still_refuses_to_emit_a_call():
    """The layer's load-bearing constraint, re-checked with a second signal present."""
    for obs in ({"RT": {"V179F"}}, {"RT": {"K103R"}}, None):
        d = target_site_doubt("efavirenz", obs).as_dict()   # as_dict() runs assert_no_call
        assert "prediction" not in d and "call" not in d


# ---------------------------------------------------------------------------
# The shipped module vs the artifact it names as its source.
#
# This index is what the CLI PRINTS: `dna-amr --drug efavirenz --observed RT:V179F` renders
# `DOUBT [strong]` from these numbers, not from wiki/. The module's own header says every field
# is "traceable to the probe artifact; nothing is asserted from memory" -- but until now nothing
# enforced it. tests/test_doubt_target_site_probe.py checks the ARTIFACT is self-consistent
# (survives family-wise correction, sits outside the catalogued positions) and the tests above
# check the MODULE behaves; neither compares the two. So a re-run that moved V179F's carrier
# count would update the artifact, keep its own tests green, and leave the CLI printing a number
# no artifact supports.
#
# That is the same drift class as wiki/certification_capstone.json reporting independent_measured=31
# against a registry returning 33 (1ff6020) and the report card guarded in 20fbe1e -- a committed
# number nothing re-derives. Here the pin follows each unit's OWN declared `artifact` path, so a
# newly measured cell is covered automatically rather than needing a hand-maintained mapping
# (cf. the hardcoded-exclusion-list trap this repo has hit five times).
# ---------------------------------------------------------------------------

import json  # noqa: E402
from pathlib import Path  # noqa: E402

_ROOT = Path(__file__).resolve().parents[1]


def _units_by_artifact() -> dict[Path, list[tuple[str, str, dict]]]:
    """Group every measured unit by the artifact IT names, so the pin needs no hand-written map."""
    out: dict[Path, list[tuple[str, str, dict]]] = {}
    for cell, units in TARGET_SITE_COMPLETENESS.items():
        for sub, rec in units.items():
            out.setdefault(_ROOT / rec["artifact"], []).append((cell, sub, rec))
    return out


def test_every_measured_unit_names_an_artifact_that_exists():
    groups = _units_by_artifact()
    assert groups, "no measured units — this pin would be vacuous"
    for art in groups:
        assert art.is_file(), f"module cites a missing artifact: {art}"


def test_module_numbers_match_the_artifact_they_cite():
    """The five run-derived fields must re-derive from the cited probe, not merely look plausible."""
    for art, units in _units_by_artifact().items():
        if not art.is_file():
            continue
        d = json.loads(art.read_text(encoding="utf-8"))
        surviving = {s["substitution"]: s for s in d["surviving"]}
        for cell, sub, rec in units:
            assert sub in surviving, (
                f"{cell}/{sub} is listed as a measured gap but {art.name} does not report it as "
                f"surviving — a gap the artifact does not support"
            )
            a = surviving[sub]
            assert rec["carriers_labelled_r"] == a["carriers"], (
                f"{sub} carriers: module {rec['carriers_labelled_r']} vs artifact {a['carriers']}"
            )
            # The module rounds for readability; pin the rounding so drift cannot hide in it.
            assert rec["purity_surprise_p"] == float(f"{a['p']:.3g}"), (
                f"{sub} p: module {rec['purity_surprise_p']} vs artifact {a['p']:.3g}"
            )
            assert rec["n_units_tested"] == d["n_units_with_min_carriers"]
            assert rec["base_s_rate"] == round(d["base_susceptible_rate"], 4)


def test_no_measured_gap_in_the_artifact_is_dropped_from_the_module():
    """The reverse direction: a surviving unit missing from the index gets NO doubt at call time.

    Scoped PER ARTIFACT, so it also catches a unit attributed to the wrong probe. With one probe
    committed it is subsumed by the on-disk test below (both fire on a renamed unit); it earns its
    keep only once a second probe exists. Kept because the per-artifact attribution is the stricter
    invariant, not because it is independently exercised today.
    """
    for art, units in _units_by_artifact().items():
        if not art.is_file():
            continue
        d = json.loads(art.read_text(encoding="utf-8"))
        listed = {sub for _, sub, _ in units}
        missing = {s["substitution"] for s in d["surviving"]} - listed
        assert not missing, f"{art.name} reports surviving units absent from the index: {sorted(missing)}"


def test_a_purity_sourced_unit_cannot_claim_a_susceptible_carrier():
    """Purity IS the signature that admits a unit; one susceptible carrier ends the signal."""
    for _, units in _units_by_artifact().items():
        for cell, sub, rec in units:
            assert rec["carriers_labelled_s"] == 0, (
                f"{cell}/{sub} claims {rec['carriers_labelled_s']} susceptible carriers, which "
                f"contradicts the purity signature that put it in the index"
            )


def test_each_unit_still_survives_correction_under_the_modules_own_denominator():
    """Checked with the MODULE's numbers — the artifact test does this with the artifact's."""
    for _, units in _units_by_artifact().items():
        for cell, sub, rec in units:
            alpha = 0.05 / max(rec["n_units_tested"], 1)
            assert rec["purity_surprise_p"] <= alpha, (
                f"{cell}/{sub} p={rec['purity_surprise_p']} fails correction over "
                f"{rec['n_units_tested']} units (alpha={alpha:.2e})"
            )


def test_no_probe_artifact_on_disk_has_its_gaps_silently_dropped():
    """Closes a hole in the reverse-direction test above.

    `_units_by_artifact()` is keyed off what the MODULE cites, so if the index dropped every unit
    for an artifact the grouping would be empty and the dropped-gap test would loop over nothing
    and pass. Discovering the probes from disk instead means removing a cell can never hide its
    gaps. Pattern-matched, not hand-listed, so a new probe is covered on arrival.
    """
    probes = sorted((_ROOT / "wiki").glob("doubt_target_site_denominator*probe.json"))
    if not probes:
        return  # gitignored/absent — nothing to compare
    listed = {sub for units in TARGET_SITE_COMPLETENESS.values() for sub in units}
    for art in probes:
        d = json.loads(art.read_text(encoding="utf-8"))
        missing = {s["substitution"] for s in d["surviving"]} - listed
        assert not missing, (
            f"{art.name} reports surviving gaps that appear NOWHERE in the index: {sorted(missing)} "
            f"— a carrier of these gets no doubt at call time"
        )
