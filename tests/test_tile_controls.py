"""Guards for the non-tile control cohort loader.

The load-bearing test is `test_pos_control_at_150_RAISES_rather_than_returning_empty`: an empty list
and a refusal are the same value to a careless caller, and they mean opposite things. Everything else
pins the parse contract against a synthetic fixture plus the two measured real-data counts.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from dna_decode.glm.tile_controls import (
    EXPECTED_AT_150,
    ControlCohortError,
    control_counts,
    load_control_tiles,
)

HDR = ["variant", "expn_med", "category", "RNA_exp_sum_1", "RNA_exp_sum_2"]


def _write(tmp_path: Path, rows: list[list[str]]) -> Path:
    p = tmp_path / "GSE144621_peak_tile_expression_formatted_std.txt.gz"
    with gzip.open(p, "wt") as fh:
        fh.write("\t".join(HDR) + "\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    return p


def _row(seq: str, cat: str, e: str = "0.5") -> list[str]:
    return [seq, e, cat, "1.0", "2.0"]


def test_filters_on_category_and_length(tmp_path):
    rows = [
        _row("A" * 150, "random"),
        _row("C" * 150, "neg_control"),
        _row("G" * 150, "tile"),          # the thing being controlled FOR
        _row("T" * 100, "random"),        # wrong geometry
    ]
    _write(tmp_path, rows)
    got = load_control_tiles(tmp_path)
    assert control_counts(got) == {"random": 1, "neg_control": 1}


def test_target_column_is_expn_med_and_reps_are_carried(tmp_path):
    _write(tmp_path, [["A" * 150, "7.25", "random", "3.5", "4.5"]])
    # single-category request on purpose: the DEFAULT asks for neg_control too, and this fixture has
    # none, so the default would (correctly) refuse. That refusal is pinned in its own test below.
    (c,) = load_control_tiles(tmp_path, categories=("random",))
    assert c.expression == 7.25 and c.rep1 == 3.5 and c.rep2 == 4.5


def test_non_acgt_is_dropped_like_load_peak_tiles_does(tmp_path):
    _write(tmp_path, [_row("N" + "A" * 149, "random"), _row("A" * 150, "random")])
    assert len(load_control_tiles(tmp_path, categories=("random",))) == 1


def test_pos_control_at_150_RAISES_rather_than_returning_empty(tmp_path):
    """The geometry-incompatibility guard, and the reason this module exists in this shape.

    `pos_control` rows are 48-75 bp in the real substrate. Returning `[]` would read downstream as
    "the encoder could not separate these", a claim about the encoder, when the fact is that the
    category does not exist at this geometry. Padding is refused: it was measured to perturb this
    encoder's pooled embedding by max |delta| 1.232.
    """
    _write(tmp_path, [_row("A" * 60, "pos_control"), _row("A" * 150, "random")])
    with pytest.raises(ControlCohortError) as ei:
        load_control_tiles(tmp_path, categories=("random", "pos_control"))
    assert "pos_control" in str(ei.value)


def test_a_category_absent_entirely_also_raises(tmp_path):
    _write(tmp_path, [_row("A" * 150, "random")])
    with pytest.raises(ControlCohortError):
        load_control_tiles(tmp_path, categories=("random", "neg_control"))


def test_missing_required_column_raises(tmp_path):
    p = tmp_path / "GSE144621_peak_tile_expression_formatted_std.txt.gz"
    with gzip.open(p, "wt") as fh:
        fh.write("variant\tcategory\n" + "A" * 150 + "\trandom\n")
    with pytest.raises(ControlCohortError):
        load_control_tiles(tmp_path)


def test_allow_download_is_refused_not_honoured(tmp_path):
    _write(tmp_path, [_row("A" * 150, "random")])
    with pytest.raises(ControlCohortError):
        load_control_tiles(tmp_path, allow_download=True)


def test_absent_substrate_raises_rather_than_returning_empty(tmp_path):
    with pytest.raises(ControlCohortError):
        load_control_tiles(tmp_path / "nope")


REAL = Path("data/mpra/GSE144621_peak_tile_expression_formatted_std.txt.gz")


@pytest.mark.skipif(not REAL.is_file(), reason="gitignored substrate not cached")
def test_real_substrate_control_counts_are_exactly_the_measured_ones():
    """936 random + 470 neg_control at 150 bp, measured 2026-11-01. These are the 1,406 sequences
    every arm on this rung has discarded via `category == 'tile'`."""
    got = control_counts(load_control_tiles(REAL.parent))
    assert got == EXPECTED_AT_150, got
    assert sum(got.values()) == 1406
