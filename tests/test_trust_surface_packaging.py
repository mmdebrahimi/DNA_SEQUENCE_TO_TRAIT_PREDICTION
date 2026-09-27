"""Pins the packaging gate (2026-06-24): the trust-surface report cards ship AS PACKAGE DATA so a built
wheel carries them (without this a wheel install silently degrades every validation badge -- the cards load
from repo-root wiki/ which doesn't exist in site-packages). These tests are FAST (no wheel build): they pin
the pyproject force-include, the drift guard (force-include must cover exactly what _load()s), and the
dual-path resolver. The full artifact-boundary proof (build wheel -> fresh-env install -> badges from the
packaged cards) is scripts/verify_wheel_ships_cards.py (manual; too heavy for CI).
"""
import re
import sys
import tomllib
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from dna_decode.data import trust_surface as ts  # noqa: E402

REPO = Path(__file__).resolve().parent.parent
def _force_included_cards() -> list[str]:
    """What the wheel actually ships -- DERIVED from pyproject, never hand-listed.

    A hand-written list here held 4 entries while pyproject force-included 5, so the HCMV
    card was covered by NO test in this file: it is _load()ed from cell_registry (not
    trust_surface), so the drift guard's regex never saw it either. Deriving removes the
    only way those two can disagree.
    """
    d = tomllib.load(open(REPO / "pyproject.toml", "rb"))
    fi = d["tool"]["hatch"]["build"]["targets"]["wheel"]["force-include"]
    return sorted(
        src.split("/", 1)[1]
        for src, dest in fi.items()
        if src.startswith("wiki/") and dest.startswith("dna_decode/report_cards/")
    )


# Modules that call trust_surface._load(): the cards are consumed from BOTH of these.
_CARD_CONSUMERS = ("dna_decode/data/trust_surface.py", "dna_decode/data/cell_registry.py")

CARDS = _force_included_cards()


def test_pyproject_force_includes_every_card():
    d = tomllib.load(open(REPO / "pyproject.toml", "rb"))
    fi = d["tool"]["hatch"]["build"]["targets"]["wheel"]["force-include"]
    assert CARDS, "derived nothing from pyproject -- the force-include table moved"
    for c in CARDS:
        assert fi.get(f"wiki/{c}") == f"dna_decode/report_cards/{c}", f"{c} not force-included for the wheel"


def test_force_include_covers_exactly_what_is_loaded():
    """Drift guard: if trust_surface starts _load()ing a new card, the force-include list MUST cover it
    (else a wheel install silently degrades that card's badges)."""
    loaded: set[str] = set()
    for mod in _CARD_CONSUMERS:
        src = (REPO / mod).read_text(encoding="utf-8")
        loaded |= set(re.findall(r'_load\("([^"]+)"\)', src))
    assert loaded, f"no _load() call found in any of {_CARD_CONSUMERS} -- the scan broke, not the code"
    # Scans EVERY consumer, not just trust_surface: the HCMV card is loaded from cell_registry,
    # so a trust_surface-only regex reported a clean 4-of-4 while the 5th shipped unguarded.
    assert loaded == set(CARDS), (
        f"cards _load()ed {sorted(loaded)} != cards force-included {CARDS} -- a card loaded but not "
        f"shipped silently degrades its badges in a wheel; one shipped but not loaded is dead weight"
    )


def test_card_path_prefers_packaged_then_falls_back(tmp_path, monkeypatch):
    pkg, wiki = tmp_path / "pkg", tmp_path / "wiki"
    pkg.mkdir(); wiki.mkdir()
    monkeypatch.setattr(ts, "_PKG_CARDS", pkg)
    monkeypatch.setattr(ts, "_WIKI", wiki)
    name = "x.json"
    (wiki / name).write_text("{}", encoding="utf-8")
    assert ts._card_path(name) == wiki / name            # only editable wiki present -> fallback
    (pkg / name).write_text("{}", encoding="utf-8")
    assert ts._card_path(name) == pkg / name             # packaged copy present -> preferred (wheel path)


def test_editable_checkout_resolves_all_cards():
    for c in CARDS:
        assert ts._card_path(c).exists(), f"{c} not resolvable in the editable checkout"
    assert ts.lookup_trust("efavirenz")["tier"] == ts.INDEPENDENT_WETLAB


if __name__ == "__main__":
    import pytest
    raise SystemExit(pytest.main([__file__, "-q"]))


def test_a_source_tree_packaged_copy_must_agree_with_wiki():
    """A STALE packaged copy in a checkout silently SHADOWS the live wiki/ card.

    `_card_path` prefers the packaged copy (correct for a wheel, pinned above). But
    `scripts/build_hcmv_report_card.py` is the one builder of the five that also writes
    dna_decode/report_cards/ in the SOURCE tree, so a checkout can hold both copies --
    and from then on every trust badge reads the packaged one. Planting a divergent
    `tier` into the packaged copy made `lookup_trust` serve it while wiki/ was correct,
    with all four tests above still green: they assert the force-include CONFIG and the
    preference ORDER, never that the two copies AGREE.

    In a real wheel wiki/ does not exist, so the both-exist condition is false and this
    skips -- it constrains only a checkout, which is the only place shadowing can happen.
    """
    for c in CARDS:
        pkg = ts._PKG_CARDS / c
        wiki = ts._WIKI / c
        if not (pkg.exists() and wiki.exists()):
            continue  # wheel (no wiki/) or plain checkout (no packaged copy) — nothing to shadow
        assert pkg.read_bytes() == wiki.read_bytes(), (
            f"{c}: the packaged copy in this source tree DIVERGES from wiki/{c}, and "
            f"_card_path prefers the packaged one — every trust badge is reading the stale "
            f"copy. Rebuild that card, or delete {pkg} to fall back to wiki/."
        )
