"""Pins the lru_cache-across-monkeypatch hazard that made a passing test fail in broad selections.

`trust_surface._load` is `@lru_cache(maxsize=None)`. Tests monkeypatch `_card_path` (or `_PKG_CARDS` /
`_WIKI`) so card names resolve to files that do not exist; `_load` then returns None AND MEMOISES IT.
`monkeypatch` restores the attribute at teardown but cannot un-memoise a value computed while the patch was
live, so the poisoned entry outlives the test that made it.

Symptom (root-caused 2026-09-26): test_trust_surface_guard's bacterial-match test passed alone and failed
inside a broad `-k` selection. With `amr_portal_independent_report_card.json` memoised as None,
`lookup_trust` skipped the EBI-portal branch and fell through to the NCBI-PD card, returning
PROVENANCE_DISJOINT instead of INDEPENDENT_MEASURED. The lookup was fine; the card had been cached away.

The FIX is the autouse fixture in tests/conftest.py, which clears the memo around every test. These tests
pin the MECHANISM so the hazard cannot quietly come back -- e.g. if someone drops the cache_clear, or
changes _load's caching in a way that reintroduces cross-test bleed.
"""
from __future__ import annotations

from pathlib import Path

import dna_decode.data.trust_surface as ts

PORTAL = "amr_portal_independent_report_card.json"


def test_load_is_actually_memoised():
    """If this stops being true the hazard is gone -- but so is the reason for the conftest fixture, and
    this file should be revisited rather than silently passing."""
    assert hasattr(ts._load, "cache_clear"), "_load is no longer memoised; revisit the conftest fixture"
    assert hasattr(ts._load, "cache_info")


def test_the_portal_card_resolves_in_a_clean_checkout():
    """Precondition for the rest: the card this hazard hid must genuinely exist, or a passing lookup below
    would prove nothing."""
    assert ts._card_path(PORTAL).exists(), f"{PORTAL} not resolvable in the editable checkout"


def test_a_value_memoised_under_a_patch_outlives_the_patch(monkeypatch, tmp_path):
    """THE MECHANISM, demonstrated in one test rather than inferred across two.

    Patch the path so the card is absent, read it (memoising None), then restore the patch by hand and show
    the None is STILL returned. That is what leaked between tests.
    """
    ts._load.cache_clear()
    real = ts._card_path
    monkeypatch.setattr(ts, "_card_path", lambda name: tmp_path / "absent.json")
    assert ts._load(PORTAL) is None                      # memoised while patched

    monkeypatch.setattr(ts, "_card_path", real)          # patch undone, as teardown would
    assert ts._load(PORTAL) is None, "expected the stale memo to survive -- that IS the hazard"

    ts._load.cache_clear()                               # the only thing that actually fixes it
    assert ts._load(PORTAL) is not None


def test_lookup_recovers_after_a_cache_clear(monkeypatch, tmp_path):
    """End-to-end: the observable symptom was a tier downgrade, so pin the tier, not just the loader."""
    ts._load.cache_clear()
    monkeypatch.setattr(ts, "_card_path", lambda name: tmp_path / "absent.json")
    poisoned = ts.lookup_trust("ciprofloxacin", "Escherichia")
    assert poisoned["tier"] != ts.INDEPENDENT_MEASURED    # no cards at all -> cannot be the portal tier

    monkeypatch.undo()
    ts._load.cache_clear()
    assert ts.lookup_trust("ciprofloxacin", "Escherichia")["tier"] == ts.INDEPENDENT_MEASURED


def test_this_test_module_did_not_leak_a_poisoned_card():
    """Self-check: the tests above deliberately poison the memo. If the conftest fixture were missing or
    not autouse, THIS test would see the damage."""
    assert ts._load(PORTAL) is not None, "a poisoned memo leaked into a later test in this very module"


def test_conftest_registers_the_clearing_fixture_as_autouse():
    """A guard on the guard: the fixture is what contains every current and future poisoner, so its
    autouse-ness is load-bearing and worth pinning in text."""
    src = (Path(__file__).parent / "conftest.py").read_text(encoding="utf-8")
    assert "_clear_trust_surface_card_cache" in src
    assert "autouse=True" in src
    assert "cache_clear" in src
