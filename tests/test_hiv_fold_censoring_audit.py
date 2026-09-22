"""Guards for the HIV fold-change censoring audit.

The load-bearing properties are (a) the degeneracy bars are IMPORTED from the shipped screen rather
than restated, (b) identifier columns cannot masquerade as measurements, and (c) the audit's SCOPE
rails ship with its numbers -- this audits published numbers, it does not change a decoder.
"""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

ARTIFACT = ROOT / "wiki" / "hiv_fold_censoring_audit_2026-09-22.json"
DATA = ROOT / "data" / "raw" / "hiv"


def _mod():
    spec = importlib.util.spec_from_file_location(
        "hiv_fold_censoring_audit", ROOT / "scripts" / "hiv_fold_censoring_audit.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


# --- the bars are the SHIPPED ones, not a restatement --------------------------------------------

def test_the_degeneracy_bars_are_imported_from_the_shipped_screen():
    """A restated constant drifts from the thing it mirrors, and then the audit silently measures a
    different rule than the one the project enforces elsewhere."""
    m = _mod()
    from scripts.forward_inverse_roundtrip import MAX_MODE_SHARE, MIN_DISTINCT_VALUES
    assert m.MAX_MODE_SHARE is MAX_MODE_SHARE
    assert m.MIN_DISTINCT_VALUES is MIN_DISTINCT_VALUES


# --- identifier columns are not measurements ------------------------------------------------------

def test_an_identifier_column_is_not_mistaken_for_a_fold_column():
    """REGRESSION, found by running this. `SeqID` is all-integer and near-unique with a maximum of
    789,388. Admitted as a fold column it BECOMES the dataset maximum, and the shared-upper-bound
    statistic is then computed against an accession number -- which reported `None` and read as a
    refutation of the ceiling."""
    m = _mod()
    rows = [{"SeqID": str(200000 + i), "EFV": "1.5"} for i in range(60)]
    cols = m.fold_columns(rows)
    assert "SeqID" not in cols


def test_the_identifier_rule_is_not_vacuous():
    """It must reject ids WITHOUT rejecting a real measurement that happens to be integer-valued.
    A fold column repeats values (many isolates share a fold), so it is not near-unique."""
    m = _mod()
    rows = [{"AZT": str(float(1 + i % 5))} for i in range(60)]   # integers, but only 5 distinct
    assert "AZT" in m.fold_columns(rows)


def test_a_flag_or_index_column_is_excluded_by_the_magnitude_floor():
    m = _mod()
    rows = [{"FLAG": str(i % 2)} for i in range(60)]
    assert "FLAG" not in m.fold_columns(rows)


# --- the censoring statistic says what it means ---------------------------------------------------

def test_a_pileup_at_the_maximum_reads_as_right_censored():
    m = _mod()
    vals = [float(i % 40) / 10 + 0.1 for i in range(200)] + [100.0] * 50
    assert m.censoring_profile(vals)["looks_right_censored"] is True


def test_a_well_spread_column_does_not_read_as_censored():
    """NON-VACUITY for the censoring detector: a spread distribution whose max is attained once must
    NOT be flagged, or every column would look censored and the statistic would mean nothing."""
    m = _mod()
    vals = [float(i) / 10 + 0.1 for i in range(300)]
    assert m.censoring_profile(vals)["looks_right_censored"] is False


# --- the committed run ----------------------------------------------------------------------------

@pytest.mark.skipif(not ARTIFACT.exists(), reason="run artifact not generated (HIV data is gitignored)")
def test_one_upper_bound_is_shared_and_nothing_exceeds_it():
    d = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    h = d["headline"]
    assert h["one_upper_bound_across_all_four_datasets"] is True
    assert h["upper_bound"] == 100.0
    assert h["no_value_anywhere_exceeds_it"] is True


@pytest.mark.skipif(not ARTIFACT.exists(), reason="run artifact not generated")
def test_two_columns_legitimately_stay_below_the_bound():
    """Not every column reaches the cap, and that is NOT counter-evidence. Pinned so nobody 'fixes'
    the statistic back to demanding a shared MAXIMUM, which is what reported `None`."""
    d = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    below = {k: v for r in d["results"].values() for k, v in r["columns_below_the_bound"].items()}
    assert set(below) == {"TDF", "BIC"}
    assert all(v < 100.0 for v in below.values())


@pytest.mark.skipif(not ARTIFACT.exists(), reason="run artifact not generated")
def test_the_shipped_bar_is_failed_by_exactly_the_two_heavy_pileups():
    d = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    assert d["headline"]["columns_failing_the_shipped_bar"] == ["NNRTI:NVP", "NRTI:3TC"]


@pytest.mark.skipif(not ARTIFACT.exists(), reason="run artifact not generated")
def test_the_existing_censoring_guard_is_measured_vacuous():
    """The guard handles operator-prefixed folds. If zero exist, it cannot fire -- and the censoring
    that IS present passes through it. A filter that removes nothing is not a control."""
    d = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    g = d["existing_guard_is_vacuous"]
    assert all(n == 0 for n in g["operator_prefixed_values_found"].values())
    assert "VACUOUS" in g["verdict"]


@pytest.mark.skipif(not ARTIFACT.exists(), reason="run artifact not generated")
def test_the_scope_rails_ship_with_the_numbers():
    """Over-claiming here would read as 'the HIV catalogs are wrong'. The artifact must carry, beside
    its numbers, that no v0.1 catalog is deployed and that the curation verdict is untouched."""
    d = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    rails = " ".join(d["what_this_does_not_claim"])
    assert "NOT a decoder change" in rails
    assert "Does NOT reopen the declined NNRTI curation" in rails
    assert "Does NOT show that any published v0.1 gain figure is wrong" in rails
