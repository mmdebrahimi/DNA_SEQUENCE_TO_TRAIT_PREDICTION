"""Offline tests for the tile-rung selection-value measurement.

The statistics come from `dna_decode/glm/selection.py`, which has its own tests. What is pinned here is
the FROZEN BAR and the two anchors that make the instrument trustworthy — plus the refusal branch,
because a bar that cannot be reached must not emit a headline.
"""
from __future__ import annotations

import gzip
import json
import random
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import glm_tile_selection as sel  # noqa: E402
from glm_tile_predictions import TilePredictionError, table_payload  # noqa: E402
from glm_tile_selection import (  # noqa: E402
    K_VALUES,
    PREREGISTERED,
    cell_for,
    render_md,
    verdict,
)

from dna_decode.glm.selection import (  # noqa: E402
    assert_k1_is_degenerate,
    permutation_null,
    simulate_selection,
)

K = PREREGISTERED["at_k"]
PRIM = PREREGISTERED["primary_scorer"]
BASE = PREREGISTERED["baseline_scorer"]


def _seed(seed: int, *, prim_h: float, base_h: float, sd: float = 0.0) -> dict:
    """One synthetic per-seed record carrying only what the verdict function reads."""
    return {
        "seed": seed,
        "cells": [
            {"scorer": PRIM, "k": K, "headroom_captured": prim_h},
            {"scorer": BASE, "k": K, "headroom_captured": base_h},
        ],
        "permutation_null_by_k": {str(K): {"n_perms": 40, "stdev": sd, "mean": 0.0, "p95": 0.0}},
    }


# --------------------------------------------------------------------------------------------------
# the frozen bar -- all four branches
# --------------------------------------------------------------------------------------------------
def test_verdict_ready_when_both_bars_pass_on_a_majority():
    v, _ = verdict([_seed(i, prim_h=0.80, base_h=0.15) for i in range(10)])
    assert v == "ORACLE_READY_AS_SCORING_HEAD"


def test_verdict_bottleneck_when_bar_a_fails():
    """Bar A failing is the result that JUSTIFIES continued prediction work."""
    v, _ = verdict([_seed(i, prim_h=0.20, base_h=0.02) for i in range(10)])
    assert v == "SELECTION_IS_THE_BOTTLENECK"


def test_verdict_crude_baseline_when_a_passes_but_b_fails():
    """GC doing the work means the biophysical scorer is not what the loop was waiting for."""
    v, _ = verdict([_seed(i, prim_h=0.80, base_h=0.78) for i in range(10)])
    assert v == "CRUDE_BASELINE_SUFFICES"


def test_verdict_REFUSES_when_fewer_seeds_than_the_bar_needs():
    """The bar needs 6 passing seeds, so on 3 seeds it is unreachable BY CONSTRUCTION.

    Non-vacuity matters here: all three supplied seeds pass BOTH bars, so a naive implementation
    would happily print a positive headline from a run whose scope could never satisfy its own bar.
    The sibling designed-grid script earned this branch by printing the opposite headline in the same
    situation.
    """
    v, detail = verdict([_seed(i, prim_h=0.90, base_h=0.05) for i in range(3)])
    assert v == "INDETERMINATE_INSUFFICIENT_SEEDS"
    assert "unreachable by construction" in detail


def test_verdict_refuses_on_no_seeds():
    v, _ = verdict([])
    assert v == "INDETERMINATE_NO_SEEDS_SCORED"


def test_bar_b_margin_is_noise_aware_not_flat():
    """With a large null stdev, a gap that clears a flat 0.10 must NOT clear 2*sd.

    This is the whole reason bar B is not a constant: the measured scorer-realisation spread (0.096)
    is the same order as 0.10, so a flat threshold could be cleared by noise alone.
    """
    wide = [_seed(i, prim_h=0.80, base_h=0.68, sd=0.20) for i in range(10)]   # gap 0.12 < 2*0.20
    v, _ = verdict(wide)
    assert v == "CRUDE_BASELINE_SUFFICES", "a 0.12 gap must not clear a 0.40 noise-aware margin"

    tight = [_seed(i, prim_h=0.80, base_h=0.68, sd=0.0) for i in range(10)]   # gap 0.12 > 0.10
    v2, _ = verdict(tight)
    assert v2 == "ORACLE_READY_AS_SCORING_HEAD"


def test_a_missing_primary_cell_is_skipped_not_counted_as_a_pass():
    seeds = [_seed(i, prim_h=0.90, base_h=0.05) for i in range(6)]
    for s in seeds[:4]:
        s["cells"] = [c for c in s["cells"] if c["scorer"] != PRIM]
    v, _ = verdict(seeds)
    assert v == "INDETERMINATE_INSUFFICIENT_SEEDS", "cells that cannot be scored must not count"


def test_frozen_bar_values_are_what_the_sibling_script_uses():
    """A drifting bar silently redefines the comparison with the designed-grid anchor."""
    assert PREREGISTERED["bar_a_headroom_captured"] == 0.50
    assert PREREGISTERED["bar_b_margin_over_baseline"] == 0.10
    assert PREREGISTERED["at_k"] == 10
    assert PREREGISTERED["primary_scorer"] == "max_tx"
    assert PREREGISTERED["baseline_scorer"] == "gc"


def test_the_comparison_anchor_records_both_pooling_units():
    """0.7937 pooled over AXES; this run pools over SEEDS. Conflating them would be a false compare."""
    a = PREREGISTERED["comparison_anchor"]
    assert a["onehot_mean_headroom_at_k10"] == 0.7937
    assert a["gc_mean_headroom_at_k10"] == 0.1328
    assert "axis" in a["note"].lower() and "seed" in a["note"].lower()


def test_cell_for_finds_and_misses_cleanly():
    s = _seed(0, prim_h=0.5, base_h=0.1)
    assert cell_for(s["cells"], PRIM, K)["headroom_captured"] == 0.5
    assert cell_for(s["cells"], PRIM, 999) is None
    assert cell_for(s["cells"], "nope", K) is None


# --------------------------------------------------------------------------------------------------
# the instrument's two anchors, on synthetic data
# --------------------------------------------------------------------------------------------------
def test_a_perfect_scorer_captures_essentially_all_headroom():
    measured = [float(i) for i in range(200)]
    r = simulate_selection(measured, measured, k=10, scorer="perfect", trials=300, seed=0)
    assert r.headroom_captured is not None and r.headroom_captured > 0.99


def test_a_shuffled_scorer_sits_inside_its_own_null_band():
    """The negative anchor. Without it, a near-zero headroom could not be told from a broken harness."""
    import random
    measured = [float(i) for i in range(200)]
    fake = measured[:]
    random.Random(7).shuffle(fake)
    r = simulate_selection(fake, measured, k=10, scorer="shuffled", trials=400, seed=0)
    nb = permutation_null(measured, k=10, trials=200, n_perms=30, seed=0)
    assert nb["n_perms"] > 0
    assert abs(r.headroom_captured) <= max(0.35, 4.0 * nb["stdev"]), (
        f"shuffled scorer headroom {r.headroom_captured} is outside its own null band "
        f"(stdev {nb['stdev']}) -- that would indicate a harness fault, not a scorer"
    )


def test_k1_headroom_is_none_and_the_self_test_accepts_it():
    measured = [float(i) for i in range(50)]
    r = simulate_selection(measured, measured, k=1, scorer="x", trials=50, seed=0)
    assert r.headroom_captured is None
    assert_k1_is_degenerate(r)          # must not raise


def test_a_fabricated_non_degenerate_k1_cell_raises():
    """k=1 cannot show a lift in a correct implementation, so the self-test localises harness bugs."""
    from dataclasses import replace
    measured = [float(i) for i in range(50)]
    r = simulate_selection(measured, measured, k=1, scorer="x", trials=50, seed=0)
    bogus = replace(r, mean_pct_picked=r.mean_pct_picked + 0.2)
    with pytest.raises(AssertionError, match="k=1 must be degenerate"):
        assert_k1_is_degenerate(bogus)


# --------------------------------------------------------------------------------------------------
# render_md -- the human-facing surface. A number carried only in the JSON is not a disclosure.
# --------------------------------------------------------------------------------------------------
def _cell(scorer: str, k: int, h: float | None) -> dict:
    return {"scorer": scorer, "k": k, "headroom_captured": h}


def _art(per_seed: list[dict] | None = None, **over) -> dict:
    """A minimal artifact carrying exactly what `render_md` reads."""
    if per_seed is None:
        per_seed = []
        for i in range(2):
            per_seed.append({
                "seed": i, "n_test": 300,
                "cells": [_cell(PRIM, 1, None), _cell(PRIM, K, 0.30 + 0.1 * i),
                          _cell(BASE, 1, None), _cell(BASE, K, 0.10)],
                "permutation_null_by_k": {
                    "1": {"n_perms": 0, "mean": 0.0, "stdev": 0.0, "p95": 0.0},
                    str(K): {"n_perms": 40, "mean": 0.0007, "stdev": 0.0961, "p95": 0.1500},
                },
            })
    art = {
        "date": "2026-10-14",
        "verdict": "SELECTION_IS_THE_BOTTLENECK",
        "verdict_detail": "at k=10, 0/2 seeds pass bar A",
        "substrate": "GEO GSE144621 peak tiles, 44106 tiles, leave_peak_out(frac=0.25), ceiling 0.941",
        "k_values": [1, K],
        "scorers": [BASE, PRIM],
        "preregistered": PREREGISTERED,
        "per_seed": per_seed,
        "honest_limits": ["UPPER BOUND on the loop's selection value", "both bars are ASSERTED"],
    }
    art.update(over)
    return art


def test_render_md_carries_the_verdict_its_detail_and_the_substrate():
    md = render_md(_art())
    assert "**VERDICT: `SELECTION_IS_THE_BOTTLENECK`**" in md
    assert "0/2 seeds pass bar A" in md
    assert "GSE144621" in md


def test_render_md_renders_the_k1_column_as_a_dash_not_a_number():
    """`headroom_captured` is None at k=1 BY CONSTRUCTION (the degeneracy the harness self-tests).
    A renderer that coerced that to 0.0000 would publish a fabricated zero for every scorer."""
    md = render_md(_art())
    row = next(ln for ln in md.splitlines() if ln.startswith(f"| `{PRIM}` |"))
    assert row.split("|")[2].strip() == "—", f"k=1 cell should be a dash, got {row!r}"
    assert "0.0000" not in row


def test_render_md_averages_headroom_ACROSS_seeds():
    """0.30 and 0.40 at k=10 must read 0.3500, not either endpoint."""
    assert "0.3500" in render_md(_art())


def test_render_md_counts_a_ZERO_headroom_rather_than_dropping_it_as_falsy():
    """The filter tests `is not None`, not truthiness. A scorer that captured exactly nothing must
    read 0.0000 — dropping it would render an em-dash, i.e. 'never measured', which is a different
    and much kinder claim than 'measured, captured nothing'."""
    seeds = [{
        "seed": 0, "n_test": 300,
        "cells": [_cell(PRIM, K, 0.0), _cell(BASE, K, 0.0)],
        "permutation_null_by_k": {str(K): {"n_perms": 40, "mean": 0.0, "stdev": 0.01, "p95": 0.0}},
    }]
    md = render_md(_art(per_seed=seeds, k_values=[K]))
    assert f"| `{PRIM}` | 0.0000 |" in md


def test_render_md_shows_the_designed_grid_anchor_with_both_numbers():
    """Without the anchor beside it, a tile-rung headroom has nothing to be compared against."""
    md = render_md(_art())
    assert "0.7937" in md and "0.1328" in md
    assert "ORACLE_READY_AS_SCORING_HEAD" in md


def test_render_md_prints_the_null_band_because_a_headroom_without_it_is_uninterpretable():
    md = render_md(_art())
    assert "0.0961" in md, "the per-seed null stdev must reach the rendered page"
    assert "uninterpretable" in md


def test_render_md_omits_a_null_band_that_was_never_computed():
    """`n_perms == 0` means the band is absent, not zero-width. Printing stdev 0.0000 for it would
    claim a precision that was never measured."""
    seeds = [{
        "seed": 0, "n_test": 300,
        "cells": [_cell(PRIM, K, 0.3), _cell(BASE, K, 0.1)],
        "permutation_null_by_k": {str(K): {"n_perms": 0, "mean": 0.0, "stdev": 0.0, "p95": 0.0}},
    }]
    md = render_md(_art(per_seed=seeds, k_values=[K]))
    assert "## Null band" in md
    assert "seed 0 k=10" not in md


def test_render_md_carries_every_honest_limit():
    md = render_md(_art())
    for lim in _art()["honest_limits"]:
        assert lim in md


def test_render_md_renders_a_row_for_every_declared_scorer():
    md = render_md(_art())
    for sc in (PRIM, BASE):
        assert f"| `{sc}` |" in md


# --------------------------------------------------------------------------------------------------
# main() -- the three refusal paths, on synthetic tables. No D:, no model, no network.
# --------------------------------------------------------------------------------------------------
def _write_seed(d: Path, seed: int, *, smoke: bool = False, n: int = 40) -> None:
    d.mkdir(parents=True, exist_ok=True)
    r = random.Random(seed)
    payload = table_payload(
        seed=seed, held=["peakA"], tile_indices=list(range(n)),
        measured=[float(i) for i in range(n)],
        predictions={"gc": [i + r.gauss(0, 8) for i in range(n)],
                     "max_tx": [i + r.gauss(0, 6) for i in range(n)]},
        spearmans={"gc": 0.3000, "max_tx": 0.3968}, smoke=smoke)
    with gzip.open(d / f"seed-{seed:02d}.json.gz", "wt", encoding="utf-8") as fh:
        json.dump(payload, fh)


def test_main_REFUSES_an_empty_input_directory(tmp_path, monkeypatch):
    """No pool is not a 'bottleneck' finding; it must name the producer to run first."""
    monkeypatch.setattr(sys, "argv", ["x", "--in-dir", str(tmp_path)])
    with pytest.raises(SystemExit, match="no seed tables"):
        sel.main()


def test_main_REFUSES_tables_with_no_MANIFEST_beside_them(tmp_path, monkeypatch):
    """Without the manifest there is no record of WHICH table was scored — tile count, split,
    protocol, anchors. The numbers would be perfectly well-formed and unattributable."""
    _write_seed(tmp_path, 0)
    monkeypatch.setattr(sys, "argv", ["x", "--in-dir", str(tmp_path), "--seeds", "1",
                                      "--trials", "20"])
    with pytest.raises(SystemExit, match="manifest.json is absent"):
        sel.main()


def test_main_REFUSES_a_SMOKE_table_at_the_consumer_boundary(tmp_path, monkeypatch):
    """A smoke table skipped the reconcile gate, so the entry point refuses rather than trusting a
    caller to have checked. Ordered BEFORE the manifest check, so this fires even with one present."""
    _write_seed(tmp_path, 0, smoke=True)
    (tmp_path / "manifest.json").write_text(json.dumps({"n_tiles": 40, "split": "x"}),
                                            encoding="utf-8")
    monkeypatch.setattr(sys, "argv", ["x", "--in-dir", str(tmp_path), "--seeds", "1",
                                      "--trials", "20"])
    with pytest.raises(TilePredictionError, match="SMOKE"):
        sel.main()


def test_main_writes_both_the_json_and_the_markdown_sidecar(tmp_path, monkeypatch):
    """End to end on one synthetic seed with a tiny trial count. The MD sidecar is the only surface
    a human reads, so a run that wrote JSON alone would be a silent half-report."""
    tables, out = tmp_path / "tbl", tmp_path / "art.json"
    _write_seed(tables, 0)
    (tables / "manifest.json").write_text(
        json.dumps({"n_tiles": 40, "split": "leave_peak_out(frac=0.25)"}), encoding="utf-8")
    monkeypatch.setattr(sys, "argv", ["x", "--in-dir", str(tables), "--seeds", "1",
                                      "--trials", "60", "--out", str(out)])
    assert sel.main() == 0
    assert out.is_file() and out.with_suffix(".md").is_file()
    art = json.loads(out.read_text(encoding="utf-8"))
    assert art["schema"] == "glm-tile-selection-v1"
    assert art["k_values"] == list(K_VALUES)
    assert art["prediction_table_manifest"]["n_tiles"] == 40
    assert "**VERDICT:" in out.with_suffix(".md").read_text(encoding="utf-8")


def test_main_cannot_emit_a_substantive_verdict_from_an_UNDER_SCOPED_run(tmp_path, monkeypatch):
    """The frozen bar needs 6 passing seeds. A 1-seed run must refuse, not report a headline — the
    sibling designed-grid script once printed the opposite headline in exactly this situation."""
    tables, out = tmp_path / "tbl", tmp_path / "art.json"
    _write_seed(tables, 0)
    (tables / "manifest.json").write_text(json.dumps({"n_tiles": 40, "split": "x"}),
                                          encoding="utf-8")
    monkeypatch.setattr(sys, "argv", ["x", "--in-dir", str(tables), "--seeds", "1",
                                      "--trials", "60", "--out", str(out)])
    sel.main()
    art = json.loads(out.read_text(encoding="utf-8"))
    assert art["verdict"] == "INDETERMINATE_INSUFFICIENT_SEEDS"
    assert art["verdict"] not in PREREGISTERED["verdicts"], (
        "a refusal must not be one of the three substantive verdicts")
