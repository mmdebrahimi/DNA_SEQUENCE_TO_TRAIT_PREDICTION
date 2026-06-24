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
CARDS = ["amr_portal_independent_report_card.json", "decoder_validation_report_card.json",
         "hiv_decoder_report_card.json", "tb_report_card.json"]


def test_pyproject_force_includes_the_four_cards():
    d = tomllib.load(open(REPO / "pyproject.toml", "rb"))
    fi = d["tool"]["hatch"]["build"]["targets"]["wheel"]["force-include"]
    for c in CARDS:
        assert fi.get(f"wiki/{c}") == f"dna_decode/report_cards/{c}", f"{c} not force-included for the wheel"


def test_force_include_covers_exactly_what_is_loaded():
    """Drift guard: if trust_surface starts _load()ing a new card, the force-include list MUST cover it
    (else a wheel install silently degrades that card's badges)."""
    src = (REPO / "dna_decode" / "data" / "trust_surface.py").read_text(encoding="utf-8")
    loaded = set(re.findall(r'_load\("([^"]+)"\)', src))
    assert loaded == set(CARDS), f"force-include list (test) drifted from trust_surface._load() calls: {loaded}"


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
