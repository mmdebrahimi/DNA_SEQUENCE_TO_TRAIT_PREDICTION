"""Every disclosure layer on the standing report card must be reachable from a decoder call.

THE GAP THIS CLOSES. Four disclosure layers rendered on `wiki/decoder_validation_report_card.json`
while only two reached a CLI call. `lineage` (10 cells) and `source_concentration` (10 cells) were
card-only, and `prospective` reached the badge ONLY when it CONTRADICTED — so a caller could not tell
"no post-lock data" from "post-lock data that agreed".

It mattered most where it was worst: `escherichia_coli_shigella x gentamicin` reports sens 0.893 from
a cohort that is 95% one BioProject holding no `rmt` carriers, and source-diverse measurements of the
same cell report 0.523. A caller deciding on 0.893 could not see the caveat explaining it.

The layer set is DERIVED from the artifact, never hand-listed — a hand-listed set beside the data that
defines it is the drift bug this repo has now hit four times.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from dna_decode.data import trust_surface as ts  # noqa: E402

CARD = ROOT / "wiki" / "decoder_validation_report_card.json"


def _cells() -> list[dict]:
    if not CARD.exists():
        pytest.skip("report card artifact absent")
    return json.loads(CARD.read_text(encoding="utf-8"))["cells"]


def observed_layers() -> set[str]:
    """Every dict-shaped per-cell disclosure block actually present on the card. DERIVED."""
    return {k for c in _cells() for k, v in c.items() if isinstance(v, dict)}


def test_every_disclosure_layer_on_the_card_reaches_at_least_one_call():
    """The core reachability contract. A layer nobody can see from the tool is not a disclosure."""
    layers = observed_layers()
    assert layers, "no disclosure layers found — this guard would be vacuous"
    unreachable = []
    for layer in sorted(layers):
        if not any(layer in ts.trust_block(c["drug"], c["organism"]) for c in _cells()):
            unreachable.append(layer)
    assert not unreachable, (
        f"disclosure layers render on the report card but reach NO decoder call: {unreachable}. "
        "Attach them in trust_surface.trust_block, augment-only."
    )


def test_the_declared_layer_list_matches_what_the_card_actually_carries():
    """`DISCLOSURE_LAYERS` exists so consumers derive rather than hand-list. It must not drift."""
    missing = observed_layers() - set(ts.DISCLOSURE_LAYERS)
    assert not missing, f"the card carries layers absent from DISCLOSURE_LAYERS: {sorted(missing)}"


# Every function `trust_block` consults to attach a disclosure layer. The "layers off" state is
# simulated by disabling ALL of them: this list was hand-enumerated and went stale the moment a fifth
# layer (`organism_scope`) was added -- the guard then compared a badge WITH that layer against one
# where it had been popped, and failed on a layer that was in fact augment-only. A stale disable-list
# makes the guard cry wolf; a missing entry makes it blind. Keep it in step with DISCLOSURE_LAYERS,
# which `test_every_disclosure_layer_has_a_disable_hook` pins.
_LAYER_SOURCES = {
    "doubt_layer_for": lambda *_a, **_k: None,
    "_cell_layer_for": lambda *_a, **_k: None,
    "overcall_for": lambda *_a, **_k: None,
}


def _disable_all_layer_sources(monkeypatch):
    for name, stub in _LAYER_SOURCES.items():
        monkeypatch.setattr(ts, name, stub)


def test_every_disclosure_layer_has_a_disable_hook():
    """A layer with no entry in _LAYER_SOURCES silently escapes the augment-only guard above."""
    assert len(_LAYER_SOURCES) >= 2
    for name in _LAYER_SOURCES:
        assert hasattr(ts, name), f"{name} no longer exists on trust_surface"
    # organism_scope + doubt_layer are attached by their own named functions; lineage /
    # source_concentration / prospective all route through _cell_layer_for.
    covered = {"doubt_layer", "organism_scope", "lineage", "source_concentration", "prospective"}
    assert set(ts.DISCLOSURE_LAYERS) <= covered, (
        f"a disclosure layer has no disable hook: {set(ts.DISCLOSURE_LAYERS) - covered}")


def test_attaching_the_layers_changes_no_pre_existing_badge_field(monkeypatch):
    """AUGMENT-ONLY, by direct comparison: same badge, layers on vs off."""
    cells = _cells()
    with_layers = {(c["organism"], c["drug"]): ts.trust_block(c["drug"], c["organism"]) for c in cells}
    _disable_all_layer_sources(monkeypatch)
    without = {(c["organism"], c["drug"]): ts.trust_block(c["drug"], c["organism"]) for c in cells}

    for k, a in with_layers.items():
        a, b = dict(a), dict(without[k])
        for layer in ts.DISCLOSURE_LAYERS:
            a.pop(layer, None)
        assert a == b, f"{k}: attaching a disclosure layer altered a pre-existing badge field"


def test_the_augment_only_check_is_not_vacuous():
    """If nothing attaches, the comparison above proves nothing."""
    attached = {layer for c in _cells() for layer in ts.DISCLOSURE_LAYERS
                if layer in ts.trust_block(c["drug"], c["organism"])}
    assert len(attached) >= 3, f"only {sorted(attached)} attach — expected every measured layer"


def test_no_layer_ever_moves_a_cell_tier():
    """A caveat is not a validation tier and must never be read as one."""
    for c in _cells():
        assert ts.trust_block(c["drug"], c["organism"])["tier"] == ts.lookup_trust(
            c["drug"], c["organism"])["tier"], f"{c['organism']}x{c['drug']}: tier moved"


# --- the human-readable half: a JSON-only disclosure is not reachable either ---

def test_single_source_cells_get_a_human_readable_concentration_line():
    """The measured case must actually print, not merely sit in the record."""
    badge = ts.trust_block("gentamicin", "Escherichia_coli_Shigella")
    s = badge.get("source_concentration")
    if not s or not s.get("single_source"):
        pytest.skip("the gentamicin cell is no longer single-source")
    line = ts.concentration_one_line(badge)
    assert line and "SINGLE-SOURCE" in line
    assert "in either direction" in line, "the caveat must not imply the estimate is inflated"


def test_a_source_diverse_cell_prints_no_concentration_line():
    """The line must be decision-relevant, not boilerplate on every call."""
    badge = ts.trust_block("ceftriaxone", "Klebsiella")
    s = badge.get("source_concentration")
    if not s or s.get("single_source"):
        pytest.skip("the klebsiella ceftriaxone cell is not source-diverse on this checkout")
    assert ts.concentration_one_line(badge) is None


def test_concentration_line_is_silent_without_a_measurement():
    """Never-measured must not render as source-diverse."""
    assert ts.concentration_one_line({"tier": "UNKNOWN"}) is None
    assert ts.concentration_one_line({"source_concentration": {"status": "not_measured"}}) is None


# --- lineage renderer: the CI is the result, so it can never be dropped ---

def test_lineage_line_never_prints_a_weighted_point_without_its_interval():
    """Effective lineage N is tiny, so a bare point estimate would mislead. `build_validation_report_
    card._assert_weighted_renderable` refuses to render one at all; a compact renderer that quietly
    dropped the CI would undo that guard."""
    no_ci = {"lineage": {"status": "scored", "raw_N": 60,
                         "effective_lineage_N": {"0.001": {"R": 9, "S": 23}},
                         "cluster_weighted": {"0.001": {"sens": 0.889, "spec": 0.957}},
                         "grade": "clonal"}}
    line = ts.lineage_one_line(no_ci)
    assert line and "0.889" not in line, "a weighted point rendered without its CI"
    assert "effective lineages" in line, "the effective-N disclosure must survive"


def test_lineage_line_uses_the_finest_threshold_available():
    """The coarser rung collapses harder; reporting it alone would overstate the correction."""
    b = {"lineage": {"status": "scored", "raw_N": 60,
                     "effective_lineage_N": {"0.001": {"R": 9, "S": 23}, "0.005": {"R": 2, "S": 18}},
                     "cluster_weighted": {
                         "0.001": {"sens": 0.889, "sens_ci": [0.5, 0.9], "spec": 0.957,
                                   "spec_ci": [0.7, 0.99]},
                         "0.005": {"sens": 0.5, "sens_ci": [0.1, 0.9], "spec": 1.0,
                                   "spec_ci": [0.8, 1.0]}},
                     "grade": "clonal"}}
    line = ts.lineage_one_line(b)
    assert "0.001" in line and "9R/23S" in line


def test_lineage_line_is_silent_when_unscored_or_absent():
    assert ts.lineage_one_line({}) is None
    assert ts.lineage_one_line({"lineage": {"status": "partial"}}) is None
    assert ts.lineage_one_line({"lineage": {"status": "scored"}}) is None      # no weighted block


def test_a_real_clonal_cell_renders_its_correction():
    badge = ts.trust_block("ciprofloxacin", "Klebsiella")
    if not badge.get("lineage"):
        pytest.skip("no lineage block on this checkout")
    line = ts.lineage_one_line(badge)
    assert line and "clonality:" in line and "effective lineages" in line
