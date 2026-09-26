"""Shared pytest fixtures."""
from pathlib import Path

import pytest


@pytest.fixture
def fixtures_dir() -> Path:
    """Path to tests/fixtures/. Populated by Step 15 (smoke pipeline + fixtures)."""
    return Path(__file__).parent / "fixtures"


@pytest.fixture
def project_root() -> Path:
    """Path to the project root."""
    return Path(__file__).parent.parent


@pytest.fixture(autouse=True)
def _clear_trust_surface_card_cache():
    """Clear `trust_surface._load`'s memo around EVERY test.

    WHY (a real order-dependent failure, root-caused 2026-09-26)
    -----------------------------------------------------------
    `trust_surface._load` is `@lru_cache(maxsize=None)`. Several tests monkeypatch `_card_path` (or
    `_PKG_CARDS`/`_WIKI`) so that most card names resolve to a file that does not exist, which makes
    `_load` return None -- AND MEMOIZE THAT None. `monkeypatch` faithfully restores the patched attribute
    at teardown, but it cannot un-memoize a value computed while the patch was live. The poisoned entry
    therefore outlives the test that created it.

    The visible symptom was `test_trust_surface_guard.py::test_bacterial_match_unchanged_and_carries_cells`
    passing alone and failing inside a broad `-k` selection: with
    `amr_portal_independent_report_card.json` memoized as None, `lookup_trust` skipped the EBI-portal
    branch and fell through to the NCBI-PD card, returning PROVENANCE_DISJOINT instead of
    INDEPENDENT_MEASURED. Nothing was wrong with the lookup -- the card had been cached away.

    Clearing around every test fixes the CLASS rather than one instance: any current or future test that
    patches a card path is contained, and no test can inherit another's memo. The cost is re-reading a few
    small JSON cards, and correctness beats that. This matters beyond tidiness -- a suite that is red only
    in some selections trains readers to discount red, which is how a genuine new failure gets missed.
    """
    try:
        from dna_decode.data.trust_surface import _load
    except Exception:                      # trust_surface unavailable -> nothing to clear
        yield
        return
    _load.cache_clear()
    yield
    _load.cache_clear()
