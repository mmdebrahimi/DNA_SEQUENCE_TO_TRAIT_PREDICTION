"""Guards for `scripts/active_frame_rescore.py`.

The gate that matters is identity-by-reconcile. `data/promoteratlas/` holds two prediction
files with different shapes and one carries a self-declared `rung` label; trusting either the
filename or that label is how a prediction vector gets scored as the wrong arm. So the tests
below pin that a file resolves to the series it actually matches, that a file matching NEITHER
series refuses rather than guesses, and that the tolerance is tight enough to reject a 0.01
perturbation -- a gate loose enough to accept anything is not a gate.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from dna_decode.glm.active_frame import ActiveFrameError
from scripts.active_frame_rescore import (
    RECONCILE_TOLERANCE,
    RescoreError,
    _binned_margin,
    _paired_bootstrap,
    check_split,
    load_committed,
    normalise_prediction_file,
    recompute_pooled,
    resolve_arm_identity,
)

SUBSTRATE = Path("data/mpra/GSE144621_peak_tile_expression_formatted_std.txt.gz")
ASGIVEN = Path("data/promoteratlas/asgiven_arm_pertile_2026-11-02.json")
ATLAS = Path("data/promoteratlas/per_tile_predictions_2026-11-01.json")
COMMITTED = Path("wiki/asgiven_arm_result_2026-11-02.json")

A_NAME = "asgiven-position-resolved-promoteratlas-150bp"
B_NAME = "promoteratlas-base-embeddings-frozen-150bp"


class _Tile:
    def __init__(self, i, expression, active="active"):
        self.seq = f"S{i}"
        self.peak = f"{i}_p_+"
        self.expression = expression
        self.rep1 = expression
        self.rep2 = expression
        self.active = active
        self.category = "tile"


# ----------------------------------------------------------------------------------------
# both committed shapes normalise to one
# ----------------------------------------------------------------------------------------

def test_normalise_accepts_the_per_seed_wrapper_shape():
    obj = {"schema": "x", "rung": "bins10_mean_asgiven",
           "per_seed": {"0": {"tile_indices": [0, 1], "pred": [0.1, 0.2]}}}
    out = normalise_prediction_file(obj)
    assert out == {0: {"tile_indices": [0, 1], "pred": [0.1, 0.2]}}


def test_normalise_accepts_the_flat_top_level_shape():
    obj = {"0": {"tile_indices": [0, 1], "pred": [0.1, 0.2]},
           "1": {"tile_indices": [2], "pred": [0.3]}}
    out = normalise_prediction_file(obj)
    assert sorted(out) == [0, 1]


def test_normalise_discards_the_self_declared_rung_label():
    """Identity comes from the reconcile, so the label must not survive into the result."""
    obj = {"rung": "bins10_mean_asgiven",
           "per_seed": {"0": {"tile_indices": [0], "pred": [0.1]}}}
    out = normalise_prediction_file(obj)
    assert "rung" not in json.dumps(out)


def test_normalise_refuses_a_file_with_no_usable_blocks():
    with pytest.raises(RescoreError, match="no per-seed"):
        normalise_prediction_file({"schema": "x", "rung": "y"})


# ----------------------------------------------------------------------------------------
# identity by reconcile, and the tolerance is not vacuous
# ----------------------------------------------------------------------------------------

def _committed(series_a, series_b, n_test=None):
    return {
        "series": {A_NAME: series_a, B_NAME: series_b},
        "medians": {A_NAME: 0.0, B_NAME: 0.0},
        "split": {s: {"n_test": (n_test or {}).get(s, 10), "n_train": 1, "held_out_peaks": 1}
                  for s in series_a},
        "comparator_provenance": {},
    }


def test_a_file_resolves_to_the_series_it_matches():
    a = {0: 0.5000, 1: 0.5100}
    b = {0: 0.4000, 1: 0.4100}
    out = resolve_arm_identity({0: 0.5000, 1: 0.5100}, _committed(a, b))
    assert out["resolved_to"] == A_NAME
    assert out["n_matching_series"] == 1
    assert out["candidates"][B_NAME]["matches"] is False


def test_a_file_matching_neither_series_refuses_rather_than_guessing():
    a = {0: 0.5000, 1: 0.5100}
    b = {0: 0.4000, 1: 0.4100}
    out = resolve_arm_identity({0: 0.1234, 1: 0.9876}, _committed(a, b))
    assert out["resolved_to"] is None
    assert out["n_matching_series"] == 0


def test_the_tolerance_rejects_a_0_01_perturbation():
    """If a 0.01 shift still matched, the gate would accept a different arm's predictions."""
    a = {0: 0.5000, 1: 0.5100}
    b = {0: 0.4000, 1: 0.4100}
    out = resolve_arm_identity({0: 0.5100, 1: 0.5200}, _committed(a, b))
    assert out["resolved_to"] is None, (
        f"a uniform +0.01 shift resolved anyway at tolerance {RECONCILE_TOLERANCE}"
    )


def test_one_seed_off_by_more_than_tolerance_breaks_the_match():
    """All seeds must match: a single outlier is the signature of a mis-aligned vector."""
    a = {0: 0.5000, 1: 0.5100, 2: 0.5200}
    b = {0: 0.1, 1: 0.2, 2: 0.3}
    out = resolve_arm_identity({0: 0.5000, 1: 0.5100, 2: 0.5900}, _committed(a, b))
    assert out["resolved_to"] is None


# ----------------------------------------------------------------------------------------
# the split gate catches what Spearman agreement cannot
# ----------------------------------------------------------------------------------------

def test_n_test_disagreement_fails_even_when_spearman_agrees():
    """Two different splits can agree on a correlation; they cannot agree on row counts."""
    tiles = [_Tile(i, float(i)) for i in range(40)]
    norm = {0: {"tile_indices": list(range(10)), "pred": [float(i) for i in range(10)]}}
    committed = _committed({0: 1.0}, {0: 0.0}, n_test={0: 999})
    out = check_split(norm, tiles, committed)
    assert out["passed"] is False
    assert out["per_seed"][0]["size_matches_committed"] is False


def test_the_split_gate_also_checks_the_exact_index_set():
    """Equal sizes do not imply the same rows, which is the stronger failure to catch."""
    tiles = [_Tile(i, float(i)) for i in range(40)]
    norm = {0: {"tile_indices": list(range(10)), "pred": [float(i) for i in range(10)]}}
    committed = _committed({0: 1.0}, {0: 0.0}, n_test={0: 10})
    out = check_split(norm, tiles, committed)
    # size matches by construction; the index set almost certainly does not
    assert out["per_seed"][0]["size_matches_committed"] is True
    assert out["per_seed"][0]["exact_index_set_matches_splitter"] is False
    assert out["passed"] is False


# ----------------------------------------------------------------------------------------
# the binned margin table cannot silently drop rows
# ----------------------------------------------------------------------------------------

def test_binned_margin_has_equal_bins_and_conserves_rows():
    n, n_bins = 97, 8
    v = list(range(n))
    out = _binned_margin(v, v, [x * 0.5 for x in v], v, n_bins)
    assert out["n_bins"] == n_bins
    assert out["n_total"] == n, "row counts must sum back to the input"
    sizes = [r["n"] for r in out["bins"]]
    assert max(sizes) - min(sizes) <= 1, f"bins are not equal-sized: {sizes}"


def test_paired_bootstrap_reports_a_ci_that_excludes_zero_on_a_clear_margin():
    out = _paired_bootstrap([0.30, 0.31, 0.32, 0.29, 0.33])
    assert out["excludes_zero"] is True
    assert out["ci_lo"] > 0


def test_paired_bootstrap_does_not_claim_significance_on_a_straddling_margin():
    out = _paired_bootstrap([-0.2, 0.25, -0.15, 0.3, -0.1, 0.2])
    assert out["excludes_zero"] is False


# ----------------------------------------------------------------------------------------
# real data, gated
# ----------------------------------------------------------------------------------------

_REAL = SUBSTRATE.exists() and ASGIVEN.exists() and ATLAS.exists() and COMMITTED.exists()


@pytest.mark.skipif(not _REAL, reason="absent substrate or committed prediction files")
def test_the_committed_medians_are_what_the_plan_records():
    c = load_committed(COMMITTED)
    assert c["medians"][A_NAME] == pytest.approx(0.5086, abs=1e-4)
    assert c["medians"][B_NAME] == pytest.approx(0.5076, abs=1e-4)


@pytest.mark.skipif(not _REAL, reason="absent substrate or committed prediction files")
def test_each_real_file_resolves_to_its_own_series_and_discriminates_against_the_other():
    from dna_decode.glm.genomewide import load_peak_tiles

    loaded, _st = load_peak_tiles("data/mpra", allow_download=False)
    tiles = [t for t in loaded if t.category == "tile" and len(t.seq) == 150]
    committed = load_committed(COMMITTED)

    for path, expect in ((ASGIVEN, A_NAME), (ATLAS, B_NAME)):
        norm = normalise_prediction_file(json.loads(path.read_text(encoding="utf-8")))
        pooled = recompute_pooled(norm, tiles)
        ident = resolve_arm_identity(pooled, committed)
        assert ident["resolved_to"] == expect, f"{path.name} resolved to {ident['resolved_to']}"
        assert ident["candidates"][expect]["max_abs_delta"] == pytest.approx(0.0, abs=1e-4)
        other = B_NAME if expect == A_NAME else A_NAME
        assert ident["candidates"][other]["matches"] is False
        # the gate must DISCRIMINATE, not merely accept: the wrong series must be clearly off
        assert ident["candidates"][other]["max_abs_delta"] > 10 * RECONCILE_TOLERANCE
        assert check_split(norm, tiles, committed)["passed"] is True


@pytest.mark.skipif(not _REAL, reason="absent substrate or committed prediction files")
def test_perturbing_one_stored_prediction_vector_makes_the_gate_fire():
    """Non-vacuity on real data: a tampered vector must not resolve."""
    from dna_decode.glm.genomewide import load_peak_tiles

    loaded, _st = load_peak_tiles("data/mpra", allow_download=False)
    tiles = [t for t in loaded if t.category == "tile" and len(t.seq) == 150]
    committed = load_committed(COMMITTED)

    norm = normalise_prediction_file(json.loads(ASGIVEN.read_text(encoding="utf-8")))
    assert resolve_arm_identity(recompute_pooled(norm, tiles), committed)["resolved_to"] == A_NAME
    # reverse one seed's predictions: same values, different alignment
    norm[0]["pred"] = list(reversed(norm[0]["pred"]))
    assert resolve_arm_identity(recompute_pooled(norm, tiles), committed)["resolved_to"] is None


@pytest.mark.skipif(not _REAL, reason="absent substrate or committed prediction files")
def test_gc_seed_zero_reproduces_the_committed_0_3000_anchor():
    """A free validation: the published anchor falls out of this script's own GC fit.

    `scripts/glm_encoder_gate.py::PREREGISTERED` pins gc seed 0 at 0.3000 +/- 0.01, sourced to
    `wiki/glm_genomewide_oracle_2026-10-07.json`. Reproducing it here proves this script reads
    the same data under the same split before any frame is changed.
    """
    from dna_decode.glm.genomewide import leave_peak_out, load_peak_tiles
    from scripts.glm_encoder_gate import ridge_gc

    loaded, _st = load_peak_tiles("data/mpra", allow_download=False)
    tiles = [t for t in loaded if t.category == "tile" and len(t.seq) == 150]
    train, test, _held = leave_peak_out(tiles, frac=0.25, seed=0)
    r, _pred = ridge_gc(train, test)
    assert r == pytest.approx(0.3000, abs=0.01)


@pytest.mark.skipif(not _REAL, reason="absent substrate or committed prediction files")
def test_restrict_predictions_refuses_an_index_from_a_different_tile_load():
    """The committed vectors index into a 44,106-tile load; a shorter load must not score."""
    from dna_decode.glm.genomewide import load_peak_tiles

    loaded, _st = load_peak_tiles("data/mpra", allow_download=False)
    tiles = [t for t in loaded if t.category == "tile" and len(t.seq) == 150]
    norm = normalise_prediction_file(json.loads(ASGIVEN.read_text(encoding="utf-8")))
    blk = norm[0]
    with pytest.raises(ActiveFrameError, match="outside the"):
        recompute_pooled({0: blk}, tiles[:100])
