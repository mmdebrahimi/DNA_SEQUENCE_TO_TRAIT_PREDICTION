"""F3 Step 1: substrate integrity + structure census (offline fixtures only).

`file exists` is not `file intact`. Every `unusable` branch has a fixture, and the truncated
gzip case pins BOTH guards -- the full decompress that catches it and the head-only read that
does NOT -- because the weak guard was written first in this repo and passes exactly the files
that later fail.
"""
from __future__ import annotations

import gzip
import json
import re
import warnings
from pathlib import Path

import numpy as np
import pytest

from scripts.f3_substrate_census import (
    HEAD_ONLY_BYTES,
    PANELS,
    PC1_SELECTION_VERDICT,
    _unpack_loaded,
    build_census,
    head_only_gzip_ok,
    matrix_integrity,
    pc1_excess,
    raw_text_integrity,
    verify_gzip_fully,
)

REPO = Path(__file__).resolve().parent.parent


# --------------------------------------------------------------------------- gzip integrity

def _truncated_gzip(tmp_path: Path) -> Path:
    good = tmp_path / "good.gz"
    with gzip.open(good, "wb") as fh:
        fh.write(b"x" * (1 << 20))
    bad = tmp_path / "trunc.gz"
    bad.write_bytes(good.read_bytes()[: good.stat().st_size // 2])
    return bad


def test_full_decompress_catches_a_truncated_gzip(tmp_path):
    ok, reason, n = verify_gzip_fully(_truncated_gzip(tmp_path))
    assert ok is False
    assert reason


def test_a_head_only_read_PASSES_the_truncated_file_so_it_is_the_wrong_guard(tmp_path):
    """Pinned as a known trap: a truncated gzip decompresses its opening blocks perfectly and
    raises only at the END, so the cheap check passes exactly the files that later fail."""
    bad = _truncated_gzip(tmp_path)
    assert head_only_gzip_ok(bad, HEAD_ONLY_BYTES) is True
    assert verify_gzip_fully(bad)[0] is False


def test_full_decompress_accepts_an_intact_gzip(tmp_path):
    good = tmp_path / "good.gz"
    with gzip.open(good, "wb") as fh:
        fh.write(b"y" * 5000)
    ok, reason, n = verify_gzip_fully(good)
    assert ok is True and reason is None and n == 5000


def test_a_gzipped_text_panel_is_parsed_after_full_decompress(tmp_path):
    p = tmp_path / "g.csv.gz"
    with gzip.open(p, "wt", encoding="utf-8") as fh:
        fh.write("id,a,b\nx,1,2\ny,3,4\n")
    out = raw_text_integrity(p)
    assert out["verdict"] == "intact"
    assert out["gzip_fully_decompressed"] is True


def test_a_truncated_gzip_panel_is_unusable_not_intact(tmp_path):
    out = raw_text_integrity(_truncated_gzip(tmp_path))
    assert out["verdict"].startswith("unusable:truncated_or_corrupt_gzip")


# --------------------------------------------------------------------------- text integrity

def test_absent_path_is_unusable_absent_and_does_not_raise(tmp_path):
    out = raw_text_integrity(tmp_path / "nope.csv")
    assert out["verdict"] == "unusable:absent"
    assert out["exists"] is False


def test_empty_file_is_unusable(tmp_path):
    p = tmp_path / "e.csv"
    p.write_text("")
    assert raw_text_integrity(p)["verdict"] == "unusable:empty_file"


def test_comment_only_file_has_no_data_rows(tmp_path):
    p = tmp_path / "c.csv"
    p.write_text("# just a header comment\n# nrow 0\n")
    out = raw_text_integrity(p)
    assert out["verdict"] == "unusable:no_data_rows"
    assert out["n_comment_lines"] == 2


def test_ragged_rows_are_unusable_and_the_widths_are_named(tmp_path):
    p = tmp_path / "r.csv"
    p.write_text("id,a,b\nx,1,2\ny,1\n")
    out = raw_text_integrity(p)
    assert out["verdict"].startswith("unusable:ragged_rows")
    assert "2" in out["verdict"] and "3" in out["verdict"]
    assert out["field_widths_agree"] is False


def test_the_field_width_histogram_distinguishes_one_stray_row_from_two_shapes(tmp_path):
    p = tmp_path / "r.csv"
    p.write_text("id,a,b\n" + "x,1,2\n" * 50 + "y,1\n")
    out = raw_text_integrity(p)
    assert out["field_width_histogram"] == {"2": 1, "3": 51}


def test_intact_text_panel_records_shape_and_alphabet(tmp_path):
    p = tmp_path / "ok.tsv"
    p.write_text("marker\ts1\ts2\nm1\tB\tR\nm2\tR\tB\n")
    out = raw_text_integrity(p)
    assert out["verdict"] == "intact"
    assert out["delimiter"] == "tab"
    assert out["n_data_lines"] == 3
    assert set(out["value_alphabet_sample"]) == {"B", "R"}


def test_comment_lines_are_excluded_from_the_width_check(tmp_path):
    """R/qtl2 panels here carry '#' preamble lines of a different width than the data."""
    p = tmp_path / "q.csv"
    p.write_text("# Genotypes for something\n# nrow 2\nid,a,b\nx,1,2\ny,3,4\n")
    out = raw_text_integrity(p)
    assert out["verdict"] == "intact"
    assert out["n_comment_lines"] == 2


# --------------------------------------------------------------------------- matrix integrity

def test_matrix_integrity_reports_missingness_and_all_missing_columns():
    X = np.array([[0.0, 1.0, np.nan],
                  [1.0, 0.0, np.nan],
                  [0.0, 0.0, np.nan]])
    out = matrix_integrity(X)
    assert out["verdict"] == "intact"
    assert out["n_individuals"] == 3 and out["n_markers"] == 3
    assert out["n_markers_entirely_missing"] == 1
    assert out["missing_fraction"] == pytest.approx(1 / 3)


def test_matrix_integrity_refuses_a_single_individual():
    assert matrix_integrity(np.zeros((1, 5)))["verdict"] == "unusable:fewer_than_2_individuals"


def test_matrix_integrity_refuses_an_entirely_missing_panel():
    X = np.full((4, 3), np.nan)
    out = matrix_integrity(X)
    assert out["verdict"] == "unusable:every_individual_entirely_missing"
    assert out["n_individuals_entirely_missing"] == 4


def test_matrix_integrity_refuses_a_zero_marker_matrix_by_that_name():
    """The marker check must come BEFORE the all-missing check: with no markers every row is
    vacuously all-missing, so the wrong order would report the wrong cause.

    This is also the reachable source of a non-finite value in the artifact -- `missing_fraction`
    is a mean over nothing -- which is what the strict-JSON guard below exists to catch. The
    warning filter keeps a known empty-slice mean from adding noise to the suite; a noisy suite
    trains readers to discount noise.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        out = matrix_integrity(np.zeros((4, 0)))
    assert out["verdict"] == "unusable:no_markers"
    assert out["n_markers"] == 0
    assert not np.isfinite(out["missing_fraction"])


# --------------------------------------------------------------------------- loader unpacking

# `_unpack_loaded` decides WHICH returned part is the genotype matrix and which are the ids.
# Picking the wrong part, or pairing ids with the wrong axis, misaligns every individual while
# still producing a matrix of a plausible shape.

def test_unpack_takes_the_matrix_and_the_ids_that_match_its_row_count():
    arr = np.arange(6, dtype=float).reshape(3, 2)
    X, ids = _unpack_loaded((arr, ["a", "b", "c"]))
    assert X.shape == (3, 2)
    assert ids == ["a", "b", "c"]


def test_unpack_DISCARDS_ids_whose_count_does_not_match_the_rows():
    """Positional placeholders are the safe answer here: a 2-id list against 3 rows cannot be
    aligned, and using it would label row 2 with row 0's id."""
    arr = np.arange(6, dtype=float).reshape(3, 2)
    X, ids = _unpack_loaded((arr, ["a", "b"]))
    assert ids == ["i0", "i1", "i2"]
    assert len(ids) == X.shape[0]


def test_unpack_only_accepts_ids_that_come_AFTER_the_matrix_in_the_tuple():
    """Pinned as behaviour, not endorsed: the id scan requires the matrix to be found first, so
    an (ids, matrix) ordering silently falls back to placeholders. Harmless today -- the ids are
    positional cluster keys and the placeholders stay row-aligned -- but it is order-dependent,
    so a reorder of some arm's loader return would change this without raising."""
    arr = np.arange(6, dtype=float).reshape(3, 2)
    assert _unpack_loaded((["a", "b", "c"], arr))[1] == ["i0", "i1", "i2"]


def test_unpack_synthesizes_row_aligned_ids_for_a_bare_matrix():
    X, ids = _unpack_loaded(np.zeros((4, 2)))
    assert ids == ["i0", "i1", "i2", "i3"]


def test_unpack_RAISES_when_no_part_is_a_matrix_rather_than_inventing_one():
    with pytest.raises(ValueError, match="no 2-D matrix"):
        _unpack_loaded((["a"], ["b"]))


# --------------------------------------------------------------------------- pc1 against its own null

def _iid_panel(n=60, m=300, seed=0):
    return np.random.default_rng(seed).integers(0, 2, size=(n, m)).astype(float)


def _two_block_panel(n=60, m=300, seed=0):
    rng = np.random.default_rng(seed)
    X = rng.integers(0, 2, size=(n, m)).astype(float)
    a = rng.integers(0, 2, size=m).astype(float)
    b = rng.integers(0, 2, size=m).astype(float)
    for i in range(n // 2):
        X[i] = np.where(rng.random(m) < 0.97, a, X[i])
    for i in range(n // 2, n):
        X[i] = np.where(rng.random(m) < 0.97, b, X[i])
    return X


def test_pc1_excess_is_about_one_on_an_unstructured_panel():
    """The whole reason this wrapper exists: raw pc1 is NOT comparable across panels because its
    no-structure expectation scales as ~1/n, so a small raw value proves nothing on its own.
    Against the panel's OWN within-marker-shuffled null, no structure must read as ~1."""
    out = pc1_excess(_iid_panel(), n_markers=300, n_perm=3)
    assert out["pc1_excess"] == pytest.approx(1.0, abs=0.25)


def test_pc1_excess_is_far_above_one_on_a_two_block_panel():
    """The ORDERING is the claim, asserted against the i.i.d. panel rather than a magic number."""
    iid = pc1_excess(_iid_panel(), n_markers=300, n_perm=3)["pc1_excess"]
    blk = pc1_excess(_two_block_panel(), n_markers=300, n_perm=3)["pc1_excess"]
    assert blk > 1.0
    assert blk > 5 * iid


def test_pc1_excess_subsamples_to_the_requested_marker_count_and_records_it():
    """The marker count must be MATCHED across panels, and the artifact must say which it used:
    a comparison of two panels at different marker counts compares two different quantities."""
    out = pc1_excess(_iid_panel(m=300), n_markers=50, n_perm=2)
    assert out["n_markers_used"] == 50
    assert out["n_perm"] == 2
    # a request above the panel's width is capped, never padded
    assert pc1_excess(_iid_panel(m=40), n_markers=500, n_perm=2)["n_markers_used"] == 40


def test_pc1_excess_is_deterministic_and_does_not_mutate_the_panel():
    """The permutation null shuffles columns, and the census goes on to CLUSTER the same matrix,
    so an in-place shuffle would corrupt the grouping with no error.

    Scope, measured by mutation rather than assumed: removing the explicit `.copy()` alone does
    NOT break this, because the marker subsample is fancy indexing and so already copies. The
    invariant has two independent sources; this pins the caller-visible property, and it does
    fail on the realistic regression (a subsample taken as a VIEW plus an in-place shuffle),
    where determinism breaks first because the second call sees a mutated panel.
    """
    X = _iid_panel(m=120)
    before = X.copy()
    a = pc1_excess(X, n_markers=100, n_perm=2, seed=0)
    b = pc1_excess(X, n_markers=100, n_perm=2, seed=0)
    assert a == b
    assert np.array_equal(X, before)


# --------------------------------------------------------------------------- path provenance

def test_every_arm_loader_panel_takes_its_paths_from_the_arm_module():
    """Not re-hardcoded: drift between the census and the arm would be silent."""
    import importlib

    for spec in PANELS:
        if spec.kind != "arm_loader":
            continue
        mod = importlib.import_module(spec.arm_module)
        assert hasattr(mod, "_GENO"), f"{spec.arm_module} lost _GENO"
        assert hasattr(mod, "_PHENO"), f"{spec.arm_module} lost _PHENO"


def test_the_dgrp_paths_match_its_arm_modules_argparse_defaults():
    """dgrp exposes paths only as argparse defaults, so equality is pinned rather than
    scraped at runtime -- drift is DETECTED instead of silently inherited."""
    spec = next(s for s in PANELS if s.name == "dgrp")
    src = (REPO / "scripts" / "dgrp_learned_decoder.py").read_text(encoding="utf-8")
    cache = re.search(r'"--cache".*?default=Path\("([^"]+)"\)', src, re.S)
    pheno = re.search(r'"--pheno".*?default=Path\("([^"]+)"\)', src, re.S)
    assert cache and pheno, "dgrp arm no longer declares --cache/--pheno defaults"
    assert spec.npz_geno == cache.group(1)
    assert spec.npz_pheno == pheno.group(1)


def test_every_panel_declares_a_design_and_a_basis_for_it():
    """`design` is how the panel was BUILT, not a measurement -- it is what lets the negative
    control claim 'randomized by construction' while the positive control only measures."""
    for spec in PANELS:
        assert spec.design in {"constructed_biparental", "constructed_multiparent",
                               "natural_population"}
        assert spec.design_basis and len(spec.design_basis) > 20


def test_exactly_one_panel_is_a_biparental_cross_and_one_is_a_natural_population():
    """The control pair is selected BY CONSTRUCTION, so each role must be unambiguous."""
    assert sum(1 for s in PANELS if s.design == "constructed_biparental") == 1
    assert sum(1 for s in PANELS if s.design == "natural_population") == 1


# --------------------------------------------------------------------------- the committed artifact

def _committed():
    hits = sorted((REPO / "wiki").glob("f3_substrate_census_*.json"))
    if not hits:
        pytest.skip("census artifact not generated in this checkout")
    return json.loads(hits[-1].read_text(encoding="utf-8"))


def test_the_committed_census_parses_and_carries_its_schema():
    """Hand-edited JSON has shipped invalid in this repo and been cited as evidence."""
    c = _committed()
    assert c["schema"] == "f3-substrate-census-v1"
    assert c["n_panels"] == len(PANELS)


def test_every_panel_in_the_committed_census_carries_an_integrity_verdict():
    for r in _committed()["panels"]:
        v = r.get("integrity")
        assert v == "intact" or v.startswith("unusable:"), f"{r['panel']}: {v!r}"


def test_the_committed_census_names_the_percentile_manufacturing_caveat():
    """A clean cluster shape is not evidence of structure; the artifact must say so."""
    limits = " ".join(_committed()["honest_limits"]).lower()
    assert "manufactures" in limits
    assert "pc1" in limits


def test_the_committed_census_carries_the_pc1_selection_verdict_VERBATIM():
    """The census's load-bearing methodological finding -- that the plan's pre-registered
    'highest pc1' rule ranked the NEGATIVE control first and is therefore unfit for selection.
    Compared against the module constant so artifact/module drift is DETECTED, not inherited:
    an artifact that paraphrased it could read as a weaker claim than the one measured."""
    c = _committed()
    assert c["pc1_selection_verdict"] == PC1_SELECTION_VERDICT
    assert "PC1_UNFIT_FOR_SELECTING_A_STRUCTURED_PANEL" in c["pc1_selection_verdict"]
    assert "BY CONSTRUCTION" in c["pc1_selection_verdict"]


@pytest.mark.parametrize("stem", ["f3_substrate_census", "f3_deconfound_calibration"])
def test_the_committed_f3_artifacts_are_STRICT_json_with_no_nan_or_infinity(stem):
    """`json.dumps` writes a bare `NaN` for a non-finite float, which Python re-reads happily
    and every other JSON parser REJECTS -- so a plain round-trip here would pass while the file
    is unreadable to anything else. Reachable: `matrix_integrity` returns
    `missing_fraction=nan` for a zero-marker matrix.

    This repo has already shipped an invalid JSON artifact, had nothing read it, and cited it
    as evidence."""
    hits = sorted((REPO / "wiki").glob(f"{stem}_*.json"))
    if not hits:
        pytest.skip(f"{stem} artifact not generated in this checkout")
    for p in hits:
        seen: list[str] = []
        # `parse_constant` fires for exactly NaN / Infinity / -Infinity as VALUES. A regex over
        # the raw text was the first attempt and was wrong in the usual direction: it also
        # matched the word inside prose ("NaN-free" in an honest_limits string), so it would
        # have failed on a perfectly valid artifact.
        json.loads(p.read_text(encoding="utf-8"), parse_constant=seen.append)
        assert not seen, f"{p.name} carries non-finite literal(s) {sorted(set(seen))}"
