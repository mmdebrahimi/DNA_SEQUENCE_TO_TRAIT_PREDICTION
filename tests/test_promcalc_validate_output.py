"""Guards for the stage-B output validator (`scripts/promcalc_validate_output.py`).

The validator exists for ONE defect: a table of NaN that looks perfect. `_site_field` returns its
default on KeyError / AttributeError / TypeError, so a future field rename yields a ~140 MB table where
every value is `nan` with 0 errors, full coverage and `sites_per_tile_mean` still 73.0 -- every count
reporting success. Only stage C's frozen reconcile gate would catch it ~12 h later, and it would look
indistinguishable from genuine model drift, which is the one diagnosis worth keeping clean.

So the test that matters is NOT "the validator passes a healthy table" -- it is that the validator
FAILS the all-NaN table. A validator only ever observed passing proves nothing, which is the lesson
this session learned the hard way about the smoke gate.

`judge()` is pure, so every threshold is testable from synthetic facts with no gzip on disk.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
VALIDATOR = REPO / "scripts" / "promcalc_validate_output.py"


def _mod():
    spec = importlib.util.spec_from_file_location("_promcalc_validate", VALIDATOR)
    assert spec and spec.loader
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def _facts(*, n_tiles: int, nonfinite_frac: float = 0.0, spt: float = 73.0,
           header: tuple[str, ...] | None = None) -> dict:
    """Synthetic scan facts, shaped exactly as `scan()` returns them."""
    m = _mod()
    n_rows = int(n_tiles * spt)
    per_field = n_rows
    return {
        "shards": ["decomposition.csv.gz"],
        "truncated_shards": [],
        "headers": [header or m.EXPECTED_HEADER],
        "n_rows": n_rows,
        "n_tiles": n_tiles,
        "unparsed_tile_ids": 0,
        "nonfinite": {f: int(per_field * nonfinite_frac) for f in m.FLOAT_FIELDS},
        "present": {f: per_field for f in m.FLOAT_FIELDS},
        "rows_per_tile_hist": [(73, n_tiles)],
    }


# --- THE defect the validator exists for --------------------------------------------------------

def test_an_ALL_NAN_table_is_a_DEFECT_even_though_every_count_looks_perfect():
    """The non-vacuity test. Full coverage, 73.0 sites/tile, correct header, zero errors -- and every
    float is nan. This is precisely what a field-rename version break produces."""
    m = _mod()
    verdict, defects, _ = m.judge(_facts(n_tiles=m.EXPECT_TILES, nonfinite_frac=1.0),
                                  expect_partial=False)
    assert verdict == "DEFECT", "an all-NaN table must NOT be reported scoreable"
    assert any("non-finite" in d for d in defects)
    assert len(defects) >= len(m.FLOAT_FIELDS), (
        "every float field is nan here, so every one of them should be named -- naming only the first "
        "would hide the scale of the break"
    )


def test_a_HEALTHY_full_table_is_SCOREABLE():
    """The other half of non-vacuity: the validator must not reject a good table."""
    m = _mod()
    verdict, defects, _ = m.judge(_facts(n_tiles=m.EXPECT_TILES), expect_partial=False)
    assert verdict == "SCOREABLE", f"healthy table rejected: {defects}"
    assert defects == []


def test_a_FEW_pathological_tiles_are_tolerated_but_a_systemic_rate_is_not():
    """The threshold must separate noise from a break, or it is either useless or a nuisance."""
    m = _mod()
    ok, _, _ = m.judge(_facts(n_tiles=m.EXPECT_TILES, nonfinite_frac=0.001), expect_partial=False)
    bad, _, _ = m.judge(_facts(n_tiles=m.EXPECT_TILES, nonfinite_frac=0.05), expect_partial=False)
    assert ok == "SCOREABLE", "0.1% non-finite is noise, not a break"
    assert bad == "DEFECT", "5% non-finite is systemic and must fail"


# --- coverage, and why --expect-partial exists --------------------------------------------------

def test_a_PARTIAL_table_is_a_defect_by_default():
    """Stage C must never be handed a partial substrate: missing tiles enter as no-promoter
    observations and produce a clean-looking verdict over incomplete data."""
    m = _mod()
    verdict, defects, _ = m.judge(_facts(n_tiles=7000), expect_partial=False)
    assert verdict == "DEFECT"
    assert any("coverage short" in d for d in defects)


def test_expect_partial_downgrades_ONLY_coverage_and_not_the_nan_check():
    """The mid-run flag must not become a blanket pass -- it exists so the NaN and shape checks are
    readable before the run finishes, which is the whole point of running it mid-flight."""
    m = _mod()
    ok, defects, notes = m.judge(_facts(n_tiles=7000), expect_partial=True)
    assert ok == "SCOREABLE" and defects == []
    assert any("expected" in n for n in notes)

    still_bad, bad_defects, _ = m.judge(_facts(n_tiles=7000, nonfinite_frac=1.0), expect_partial=True)
    assert still_bad == "DEFECT", "--expect-partial must NOT suppress the all-NaN defect"
    assert all("coverage" not in d for d in bad_defects)


# --- shape --------------------------------------------------------------------------------------

def test_a_moved_header_is_a_DEFECT():
    m = _mod()
    wrong = tuple(list(m.EXPECTED_HEADER)[:-1] + ["renamed_field"])
    verdict, defects, _ = m.judge(_facts(n_tiles=m.EXPECT_TILES, header=wrong), expect_partial=False)
    assert verdict == "DEFECT"
    assert any("header mismatch" in d for d in defects)


def test_a_moved_sites_per_tile_is_a_DEFECT():
    """~73 forward candidates per 150 bp tile is the published shape; a large move means the model's
    output changed even if nothing is nan."""
    m = _mod()
    for spt in (20.0, 150.0):
        verdict, defects, _ = m.judge(_facts(n_tiles=m.EXPECT_TILES, spt=spt), expect_partial=False)
        assert verdict == "DEFECT", f"sites/tile {spt} should fail"
        assert any("sites/tile" in d for d in defects)


def test_a_missing_field_is_named_rather_than_silently_skipped():
    m = _mod()
    f = _facts(n_tiles=m.EXPECT_TILES)
    f["present"]["dG_total"] = 0
    verdict, defects, _ = m.judge(f, expect_partial=False)
    assert verdict == "DEFECT"
    assert any("dG_total never appeared" in d for d in defects)


def test_no_shards_is_INDETERMINATE_not_scoreable(tmp_path: Path):
    """An absent table is an unknown, never a pass -- the fail direction matters."""
    m = _mod()
    verdict, defects, _ = m.judge(m.scan(tmp_path), expect_partial=False)
    assert verdict == "INDETERMINATE"
    assert any("no decomposition shards" in d for d in defects)


def test_shard_enumeration_never_picks_up_the_SMOKE_output(tmp_path: Path):
    """A bare `decomposition*` glob also matches decomposition_smoke*.csv.gz, and the real run writes
    its smoke output into the SAME directory -- so this separation is live, not hypothetical."""
    m = _mod()
    for name in ("decomposition.csv.gz", "decomposition.part-01.csv.gz",
                 "decomposition_smoke.csv.gz"):
        (tmp_path / name).write_bytes(b"")
    names = [p.name for p in m.shard_paths(tmp_path)]
    assert "decomposition_smoke.csv.gz" not in names, "smoke output must never join a full substrate"
    assert names == ["decomposition.csv.gz", "decomposition.part-01.csv.gz"]
