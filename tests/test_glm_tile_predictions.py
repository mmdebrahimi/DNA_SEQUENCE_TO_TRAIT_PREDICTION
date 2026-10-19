"""Offline tests for the shared tile prediction table. No D:, no model, no network.

The script's pure half is separated from its loading half precisely so these can run anywhere. What
they pin is not arithmetic but the REFUSALS -- a table that is silently mis-aligned, partial, or a
smoke run masquerading as real would produce numbers that look entirely plausible downstream, which is
this repo's recorded failure shape.
"""
from __future__ import annotations

import gzip
import json
import sys
from dataclasses import dataclass
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import glm_tile_predictions as pred  # noqa: E402
from glm_tile_predictions import (  # noqa: E402
    RECONCILE,
    SCHEMA,
    SCORERS,
    TilePredictionError,
    load_table,
    manifest_payload,
    max_tx_scalar,
    merge_payload,
    reconcile_problems,
    require_real,
    split_indices,
    table_payload,
    targets,
)


@dataclass
class FakeTile:
    seq: str
    expression: float
    category: str = "tile"


def _tiles(n: int = 6) -> list[FakeTile]:
    return [FakeTile("ACGT" * 10, float(i + 1)) for i in range(n)]


# --------------------------------------------------------------------------------------------------
# split_indices
# --------------------------------------------------------------------------------------------------
def test_split_indices_returns_positions_by_identity():
    tiles = _tiles(6)
    itr, ite = split_indices(tiles, tiles[:4], tiles[4:])
    assert itr == [0, 1, 2, 3]
    assert ite == [4, 5]


def test_split_indices_refuses_an_overlapping_split():
    """A leaking split is the one error that still yields a plausible number."""
    tiles = _tiles(6)
    with pytest.raises(TilePredictionError, match="overlap"):
        split_indices(tiles, tiles[:4], tiles[3:])


def test_split_indices_refuses_a_non_exhaustive_split():
    """A tile in neither half would be silently unscored."""
    tiles = _tiles(6)
    with pytest.raises(TilePredictionError, match="exhaustive"):
        split_indices(tiles, tiles[:2], tiles[4:])


def test_split_indices_refuses_copies_rather_than_mis_joining():
    """If the splitter ever returns copies, identity no longer joins and every row shifts.

    An equal-VALUE copy is the dangerous case: it compares equal but is a different object, so a
    value-based join would succeed at the wrong index.
    """
    tiles = _tiles(4)
    copies = [FakeTile(t.seq, t.expression) for t in tiles[:2]]
    with pytest.raises(TilePredictionError, match="COPIES"):
        split_indices(tiles, copies, tiles[2:])


# --------------------------------------------------------------------------------------------------
# targets -- the asymmetry
# --------------------------------------------------------------------------------------------------
def test_targets_log1p_on_train_and_raw_on_test():
    import math
    tiles = _tiles(4)
    ytr, yte = targets(tiles[:2], tiles[2:])
    assert ytr == [math.log1p(1.0), math.log1p(2.0)]
    assert yte == [3.0, 4.0]


def test_targets_are_not_the_same_transform():
    """Pins the asymmetry itself: if someone 'tidies' this to one transform, the test fails."""
    tiles = _tiles(4)
    ytr, yte = targets(tiles[:2], tiles[:2])
    assert ytr != yte, "train and test targets must use DIFFERENT transforms (log1p vs raw)"


# --------------------------------------------------------------------------------------------------
# max_tx_scalar
# --------------------------------------------------------------------------------------------------
def test_max_tx_scalar_takes_the_strongest_site():
    import math
    v = max_tx_scalar([{"Tx_rate": 2.0}, {"Tx_rate": 9.0}, {"Tx_rate": 1.0}])
    assert v == [math.log1p(9.0)]


def test_max_tx_scalar_on_a_site_free_tile_is_zero_not_an_error():
    """No candidate promoter is a real observation, not missing data."""
    assert max_tx_scalar([]) == [0.0]


def test_max_tx_scalar_clamps_a_negative_rate():
    assert max_tx_scalar([{"Tx_rate": -5.0}]) == [0.0]


# --------------------------------------------------------------------------------------------------
# the reconcile gate -- proven to REFUSE, which is the only thing that makes it evidence
# --------------------------------------------------------------------------------------------------
def test_reconcile_clean_on_the_published_anchors():
    assert reconcile_problems(0, {"gc": 0.3000, "max_tx": 0.3968}) == []


def test_reconcile_tolerates_within_band():
    assert reconcile_problems(0, {"gc": 0.3050, "max_tx": 0.4100}) == []


def test_reconcile_REFUSES_a_drifted_gc():
    probs = reconcile_problems(0, {"gc": 0.2000})
    assert len(probs) == 1 and "gc seed 0" in probs[0]


def test_reconcile_REFUSES_a_drifted_max_tx():
    probs = reconcile_problems(0, {"max_tx": 0.3000})
    assert len(probs) == 1 and "max_tx seed 0" in probs[0]


def test_reconcile_only_gates_scorers_that_actually_ran():
    """`--scorers gc` must not refuse because max_tx was never asked for."""
    assert reconcile_problems(0, {"gc": 0.3000}) == []


def test_reconcile_only_gates_the_anchored_seed():
    """The anchors are seed-0 facts; seed 3's gc legitimately differs."""
    assert reconcile_problems(3, {"gc": 0.9999}) == []


def test_reconcile_anchors_match_the_published_values():
    """A drifting anchor silently redefines what 'faithful' means."""
    assert RECONCILE["gc"]["expected"] == 0.3000
    assert RECONCILE["max_tx"]["expected"] == 0.3968


# --------------------------------------------------------------------------------------------------
# table_payload / load_table round-trip
# --------------------------------------------------------------------------------------------------
def _payload(**over):
    kw = dict(seed=0, held=["peakA"], tile_indices=[4, 5], measured=[3.0, 4.0],
              predictions={"gc": [0.1, 0.2]}, spearmans={"gc": 0.3}, smoke=False)
    kw.update(over)
    return table_payload(**kw)


def test_table_payload_round_trips_through_gzip_json(tmp_path):
    p = tmp_path / "seed-00.json.gz"
    with gzip.open(p, "wt", encoding="utf-8") as fh:
        json.dump(_payload(), fh)
    d = load_table(p)
    assert d["schema"] == SCHEMA
    assert d["tile_indices"] == [4, 5]
    assert d["predictions"]["gc"] == [0.1, 0.2]
    assert d["measured"] == [3.0, 4.0]


def test_table_payload_refuses_a_ragged_prediction_vector():
    with pytest.raises(TilePredictionError, match="predictions for"):
        _payload(predictions={"gc": [0.1]})


def test_table_payload_refuses_ragged_measured():
    with pytest.raises(TilePredictionError, match="measured has"):
        _payload(measured=[3.0])


def test_table_payload_refuses_an_empty_scorer_set():
    """A table with no scorer is not a table; writing one would read as 'nothing to say'."""
    with pytest.raises(TilePredictionError, match="no scorer"):
        _payload(predictions={}, spearmans={})


def test_load_table_refuses_a_foreign_schema(tmp_path):
    p = tmp_path / "x.json.gz"
    bad = _payload()
    bad["schema"] = "something-else-v9"
    with gzip.open(p, "wt", encoding="utf-8") as fh:
        json.dump(bad, fh)
    with pytest.raises(TilePredictionError, match="schema"):
        load_table(p)


def test_load_table_refuses_a_ragged_file_written_by_hand(tmp_path):
    """The write path validates, but a hand-edited or truncated file must be caught on READ too."""
    p = tmp_path / "x.json.gz"
    bad = _payload()
    bad["predictions"]["gc"] = [0.1]
    with gzip.open(p, "wt", encoding="utf-8") as fh:
        json.dump(bad, fh)
    with pytest.raises(TilePredictionError, match="ragged"):
        load_table(p)


# --------------------------------------------------------------------------------------------------
# require_real -- a smoke table must never be scored
# --------------------------------------------------------------------------------------------------
def test_require_real_accepts_real_tables():
    require_real([_payload(), _payload(seed=1)])


def test_require_real_REFUSES_a_smoke_table():
    with pytest.raises(TilePredictionError, match="SMOKE"):
        require_real([_payload(), _payload(seed=1, smoke=True)])


def test_require_real_names_every_smoke_seed():
    with pytest.raises(TilePredictionError) as exc:
        require_real([_payload(seed=2, smoke=True), _payload(seed=5, smoke=True)])
    assert "[2, 5]" in str(exc.value)


# --------------------------------------------------------------------------------------------------
# manifest
# --------------------------------------------------------------------------------------------------
def test_manifest_records_what_identifies_this_table():
    m = manifest_payload(seeds=[0, 1], scorers=["gc", "conv"], n_tiles=44106, coverage=1.0,
                         shards=["decomposition.csv.gz"], conv_params=13217,
                         protocol={"epochs": 1}, smoke=False)
    assert m["n_tiles"] == 44106 and m["expected_tiles"] == 44106
    assert m["conv_parameters"] == 13217
    assert "log1p" in m["targets"] and "raw" in m["targets"]
    assert m["split"].startswith("leave_peak_out")
    assert m["smoke"] is False
    assert m["reconcile_anchors"]["gc"]["expected"] == 0.3000


# --------------------------------------------------------------------------------------------------
# merge_payload -- the staged-run defect found by RUNNING the script
# --------------------------------------------------------------------------------------------------
def test_merge_unions_scorer_columns():
    """The staged workflow: cheap scorers land, then conv is added without destroying them."""
    first = _payload(predictions={"gc": [0.1, 0.2]}, spearmans={"gc": 0.3})
    second = _payload(predictions={"conv": [0.5, 0.6]}, spearmans={"conv": 0.32})
    m = merge_payload(first, second)
    assert sorted(m["predictions"]) == ["conv", "gc"]
    assert m["predictions"]["gc"] == [0.1, 0.2]
    assert m["predictions"]["conv"] == [0.5, 0.6]
    assert sorted(m["spearman_recomputed"]) == ["conv", "gc"]


def test_merge_lets_a_refit_scorer_win():
    """Re-running one scorer should UPDATE it, not keep a stale column beside a fresh one."""
    first = _payload(predictions={"gc": [0.1, 0.2]}, spearmans={"gc": 0.30})
    second = _payload(predictions={"gc": [0.9, 0.8]}, spearmans={"gc": 0.31})
    m = merge_payload(first, second)
    assert m["predictions"]["gc"] == [0.9, 0.8]
    assert m["spearman_recomputed"]["gc"] == 0.31


def test_merge_REFUSES_a_different_split():
    """The dangerous case: uniting columns across splits mis-joins scorer to tile."""
    first = _payload(tile_indices=[4, 5], predictions={"gc": [0.1, 0.2]}, spearmans={"gc": 0.3})
    second = _payload(tile_indices=[7, 9], predictions={"conv": [0.5, 0.6]}, spearmans={"conv": 0.3})
    with pytest.raises(TilePredictionError, match="DIFFERENT splits"):
        merge_payload(first, second)


def test_merge_REFUSES_same_split_but_changed_measured():
    first = _payload(measured=[3.0, 4.0], predictions={"gc": [0.1, 0.2]}, spearmans={"gc": 0.3})
    second = _payload(measured=[3.0, 9.0], predictions={"conv": [0.5, 0.6]}, spearmans={"conv": 0.3})
    with pytest.raises(TilePredictionError, match="different measured"):
        merge_payload(first, second)


def test_merge_REFUSES_mixing_a_smoke_table_with_a_real_one():
    """Merging a smoke column into a real table would launder it past `require_real`."""
    real = _payload(predictions={"gc": [0.1, 0.2]}, spearmans={"gc": 0.3})
    smoky = _payload(smoke=True, predictions={"conv": [0.5, 0.6]}, spearmans={"conv": 0.3})
    with pytest.raises(TilePredictionError, match="smoke differs"):
        merge_payload(real, smoky)


def test_merge_REFUSES_a_different_seed():
    first = _payload(seed=0, predictions={"gc": [0.1, 0.2]}, spearmans={"gc": 0.3})
    second = _payload(seed=1, predictions={"conv": [0.5, 0.6]}, spearmans={"conv": 0.3})
    with pytest.raises(TilePredictionError, match="seed differs"):
        merge_payload(first, second)


def test_merge_leaves_the_scorer_columns_in_SORTED_order():
    """Both dicts are re-sorted after the union, so a late-merged column cannot land out of order.

    Non-vacuous by construction: the merge is fed `z` first and `a` second, so insertion order and
    sorted order disagree.
    """
    first = _payload(predictions={"z": [0.1, 0.2]}, spearmans={"z": 0.3})
    second = _payload(predictions={"a": [0.5, 0.6]}, spearmans={"a": 0.31})
    m = merge_payload(first, second)
    assert list(m["predictions"]) == ["a", "z"]
    assert list(m["spearman_recomputed"]) == ["a", "z"]


# --------------------------------------------------------------------------------------------------
# main() -- the one seam reachable offline: the scorer-name gate fires BEFORE any load.
# --------------------------------------------------------------------------------------------------
def test_main_REFUSES_an_unknown_scorer_before_touching_the_substrate(monkeypatch):
    """A typo in `--scorers` must not be discovered after a multi-hour conv fit.

    This test is also the proof of ORDERING: the gate sits above `sys.path.insert` and the heavy
    imports, so on a machine with no D: and no torch the call still reaches a clean refusal. Move the
    check below the imports and this test stops passing for a different reason.
    """
    monkeypatch.setattr(sys, "argv", ["x", "--scorers", "gc,bogus,conv"])
    with pytest.raises(SystemExit) as exc:
        pred.main()
    assert "bogus" in str(exc.value)
    assert "gc,bogus" not in str(exc.value), "only the unknown names should be reported"


def test_main_names_the_known_scorer_set_in_the_refusal(monkeypatch):
    """The message must be derived from the registry, so adding a scorer cannot leave stale help."""
    monkeypatch.setattr(sys, "argv", ["x", "--scorers", "nope"])
    with pytest.raises(SystemExit) as exc:
        pred.main()
    for name in SCORERS:
        assert name in str(exc.value)


def test_main_strips_whitespace_around_scorer_names(monkeypatch):
    """`--scorers "gc , bogus"` must reject `bogus`, not ` bogus` — otherwise a valid name typed with
    a space would be rejected as unknown."""
    monkeypatch.setattr(sys, "argv", ["x", "--scorers", "gc , bogus "])
    with pytest.raises(SystemExit) as exc:
        pred.main()
    assert "['bogus']" in str(exc.value)


def test_manifest_carries_the_upper_bound_limit():
    """The selection numbers this table feeds are upper bounds; the manifest must say so."""
    m = manifest_payload(seeds=[0], scorers=["gc"], n_tiles=1, coverage=1.0, shards=[],
                         conv_params=None, protocol=None, smoke=False)
    assert any("upper bound" in lim for lim in m["honest_limits"])
