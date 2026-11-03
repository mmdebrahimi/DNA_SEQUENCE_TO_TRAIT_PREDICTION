"""Tests for the expression-heads scorer.

SCOPE, STATED SO IT IS NOT MISREAD AS COVERAGE IT IS NOT. The head arithmetic -- the windowing,
the asymmetric `padding='same'` split, the edge-contamination closed form -- is covered by
`tests/test_atlas_heads.py` (50 tests) and is deliberately NOT re-tested here. What this file
covers is everything BETWEEN that arithmetic and the verdict: the streaming block's resume and
PARTIAL stamping, the gate wiring, and the re-derivation of the verdict from the artifact. Those
are the parts that can silently turn 4.6 hours of correct compute into an unreadable or
over-claimed result.

THE REFUSALS ARE THE POINT. A capped `--limit` run stamps the block PARTIAL, and the scorer must
REFUSE to read it as the arm's result; a block whose shape disagrees with its own sidecar must
refuse too. Both are cheap to get wrong in a way no number would reveal -- a PARTIAL block scores
perfectly well, it just answers a different question than the one the arm asked.
"""
from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

import scripts.expression_heads_score as S
from dna_decode.glm import atlas_heads as AH
from scripts.freeze_expression_heads_bar import (
    CHECK_ORDER,
    NULL_SIGNAL_RATIO,
    PRIMARY_REDUCTION,
    SECONDARY_REDUCTION,
)
from scripts.positional_arm_score import ScoreError

HEADS = ["trspred-hossain2020", "trspred-lafleur2022", "trspred-urtecho2018",
         "trspred-yu2021", "tslpred-kosuri2013"]


# --------------------------------------------------------------------- column layout


def test_the_layout_covers_every_head_crossed_with_every_reduction():
    lay = S.column_layout(HEADS)
    assert lay["_total"] == len(HEADS) * len(S.REDUCTIONS) == 10
    for h in HEADS:
        for r in S.REDUCTIONS:
            assert "%s::%s" % (h, r) in lay


def test_the_layout_columns_are_unique_and_contiguous():
    """A duplicated column index would have two heads writing over each other, and the only
    symptom would be two heads reporting the same Spearman."""
    lay = S.column_layout(HEADS)
    idx = sorted(v for k, v in lay.items() if k != "_total")
    assert idx == list(range(lay["_total"]))


def test_the_layout_is_order_stable_under_input_shuffling():
    """It is keyed by SORTED name, so the block written by one invocation is readable by the next
    regardless of the order the heads were loaded in. The resume path compares layouts for
    equality, so an unstable layout would silently discard a part-finished block."""
    a = S.column_layout(HEADS)
    b = S.column_layout(list(reversed(HEADS)))
    assert a == b


def test_the_layout_sorts_trspred_before_tslpred():
    """Alphabetical, which puts the four transcription heads before the translation one. Pinned
    because the natural reading ('transcription first, then translation') happens to agree with
    sort order here and would stop agreeing if a head were renamed."""
    lay = S.column_layout(HEADS)
    assert lay["trspred-hossain2020::%s" % PRIMARY_REDUCTION] == 0
    assert lay["tslpred-kosuri2013::%s" % PRIMARY_REDUCTION] > 0


# --------------------------------------------------------------------- sidecar + block loading


def _write_block(d: Path, *, status="COMPLETE", n_rows=8, heads=None, n_done=None,
                 shape=None) -> dict:
    heads = heads or HEADS
    lay = S.column_layout(heads)
    d.mkdir(parents=True, exist_ok=True)
    arr = np.lib.format.open_memmap(d / S.BLOCK_NAME, mode="w+", dtype="float64",
                                    shape=shape or (n_rows, lay["_total"]))
    arr[:] = np.random.default_rng(0).normal(size=arr.shape)
    arr.flush()
    side = {"schema": "expression-heads-scalars-v1", "status": status, "n_rows": n_rows,
            "n_done": n_rows if n_done is None else n_done, "layout": lay,
            "heads": sorted(heads), "reductions": list(S.REDUCTIONS), "route": S.ROUTE_TRUE}
    (d / S.SIDECAR_NAME).write_text(json.dumps(side), encoding="utf-8")
    return side


def test_a_missing_sidecar_reads_as_absent(tmp_path):
    assert S.read_sidecar(tmp_path) is None


def test_a_malformed_sidecar_reads_as_absent(tmp_path):
    """Absent rather than raising, so a half-written sidecar restarts the extraction instead of
    wedging it. The extraction is resumable; a crash during the sidecar write costs one batch."""
    (tmp_path / S.SIDECAR_NAME).write_text('{"status": "COMP', encoding="utf-8")
    assert S.read_sidecar(tmp_path) is None


def test_load_scalars_refuses_when_there_is_no_sidecar(tmp_path):
    with pytest.raises(ScoreError, match="no sidecar"):
        S.load_scalars(tmp_path)


def test_load_scalars_REFUSES_a_partial_block(tmp_path):
    """The decisive refusal. A `--limit` run scores perfectly well and answers a DIFFERENT
    question than the arm asked, so reading it as the arm's result is the failure mode with no
    symptom."""
    _write_block(tmp_path, status="PARTIAL", n_rows=8, n_done=3)
    with pytest.raises(ScoreError, match="not COMPLETE"):
        S.load_scalars(tmp_path)


def test_load_scalars_allows_a_partial_block_only_when_asked(tmp_path):
    _write_block(tmp_path, status="PARTIAL", n_rows=8, n_done=3)
    arr, side = S.load_scalars(tmp_path, allow_partial=True)
    assert side["status"] == "PARTIAL"
    assert arr.shape[0] == 8


def test_load_scalars_refuses_a_block_whose_shape_contradicts_its_sidecar(tmp_path):
    """A sidecar is a claim about a file. If the two disagree, every column index is suspect."""
    _write_block(tmp_path, n_rows=8, shape=(8, 4))
    with pytest.raises(ScoreError, match="sidecar says"):
        S.load_scalars(tmp_path)


def test_load_scalars_accepts_a_complete_block(tmp_path):
    side = _write_block(tmp_path, n_rows=8)
    arr, got = S.load_scalars(tmp_path)
    assert arr.shape == (8, side["layout"]["_total"])
    assert got["status"] == "COMPLETE"


# --------------------------------------------------------------------- streaming extraction


class _FakeHead:
    """Only `.name` and `.state_dict` are read by `extract`."""

    def __init__(self, name):
        self.name = name
        self.state_dict = {"fc1.weight": np.zeros(AH.FC1_SHAPE)}


@pytest.fixture
def stub_head_math(monkeypatch):
    """Replace the head arithmetic with a deterministic stand-in.

    The arithmetic is covered by test_atlas_heads.py. Stubbing it here is what makes the
    streaming, resume and stamping logic testable at all -- the real path needs the PromoterAtlas
    weights, torch, and minutes of compute per batch.
    """
    import torch

    def per_window(z, _sd):
        # one row per sequence, N_WINDOWS columns, values keyed to the batch content so a
        # misaligned write is visible rather than averaged away
        n = z.shape[0]
        base = z.reshape(n, -1)[:, :1]
        return base.expand(n, AH.n_windows(150)).clone()

    def reduce_windows(pw, red):
        return pw.mean(dim=1) if red == PRIMARY_REDUCTION else pw.max(dim=1).values

    monkeypatch.setattr(AH, "head_forward_per_window", per_window)
    monkeypatch.setattr(AH, "reduce_windows", reduce_windows)
    return torch


def _tiles(n, length=150):
    rng = np.random.default_rng(7)
    return [SimpleNamespace(seq="".join(rng.choice(list("ACGT"), size=length)),
                            expression=float(i), peak="p%d" % (i // 3)) for i in range(n)]


def _encode(seqs):
    """A latent whose first element encodes the tile's own sequence, so misalignment shows."""
    out = np.zeros((len(seqs), AH.LATENT_CHANNELS, 150), dtype=np.float32)
    for i, s in enumerate(seqs):
        out[i, 0, 0] = float(sum(ord(c) for c in s) % 10_000)
    return out


def test_extract_writes_a_complete_block_and_stamps_it(tmp_path, stub_head_math):
    tiles = _tiles(10)
    heads = [_FakeHead(h) for h in HEADS[:2]]
    side = S.extract(tiles=tiles, heads=heads, block_dir=tmp_path, batch=4,
                     encode=_encode, verbose=False)
    assert side["status"] == "COMPLETE"
    assert side["n_done"] == side["n_rows"] == 10
    arr, _ = S.load_scalars(tmp_path)
    assert np.isfinite(arr).all()


def test_extract_stamps_PARTIAL_under_a_limit(tmp_path, stub_head_math):
    """`--limit` is the capped diagnostic path, and the stamp is what stops it being read as the
    arm's result."""
    side = S.extract(tiles=_tiles(10), heads=[_FakeHead(HEADS[0])], block_dir=tmp_path,
                     limit=4, batch=2, encode=_encode, verbose=False)
    assert side["status"] == "PARTIAL"
    assert side["n_done"] == 4 and side["n_rows"] == 10
    with pytest.raises(ScoreError, match="not COMPLETE"):
        S.load_scalars(tmp_path)


def test_extract_RESUMES_from_a_partial_block_without_recomputing(tmp_path, stub_head_math):
    """The property a 4.6-hour run depends on: a death mid-flight must cost one batch, not the
    whole run. Verified by counting encoder calls on the second pass.
    """
    tiles, heads = _tiles(10), [_FakeHead(HEADS[0])]
    S.extract(tiles=tiles, heads=heads, block_dir=tmp_path, limit=6, batch=2,
              encode=_encode, verbose=False)

    seen = []

    def counting_encode(seqs):
        seen.extend(seqs)
        return _encode(seqs)

    side = S.extract(tiles=tiles, heads=heads, block_dir=tmp_path, batch=2,
                     encode=counting_encode, verbose=False)
    assert side["status"] == "COMPLETE"
    assert len(seen) == 4, "resume recomputed %d sequences; only the last 4 were outstanding" \
                           % len(seen)


def test_extract_resume_preserves_the_rows_already_written(tmp_path, stub_head_math):
    """Resuming must not disturb finished rows. Checked by value, not by the sidecar's count."""
    tiles, heads = _tiles(10), [_FakeHead(HEADS[0])]
    S.extract(tiles=tiles, heads=heads, block_dir=tmp_path, limit=6, batch=2,
              encode=_encode, verbose=False)
    first, _ = S.load_scalars(tmp_path, allow_partial=True)
    before = np.array(first[:6, :])
    S.extract(tiles=tiles, heads=heads, block_dir=tmp_path, batch=2, encode=_encode,
              verbose=False)
    after, _ = S.load_scalars(tmp_path)
    assert np.allclose(before, after[:6, :])


def test_extract_RESTARTS_when_the_head_set_changed(tmp_path, stub_head_math):
    """A block written for one head set has a different layout, so its columns mean different
    things. Resuming across that change would read one head's scalars as another's."""
    tiles = _tiles(6)
    S.extract(tiles=tiles, heads=[_FakeHead(HEADS[0])], block_dir=tmp_path, limit=4, batch=2,
              encode=_encode, verbose=False)
    side = S.extract(tiles=tiles, heads=[_FakeHead(h) for h in HEADS[:3]], block_dir=tmp_path,
                     batch=2, encode=_encode, verbose=False)
    assert side["layout"]["_total"] == 3 * len(S.REDUCTIONS)
    assert side["n_done"] == 6


def test_extract_RESTARTS_when_the_tile_count_changed(tmp_path, stub_head_math):
    """Row i in a 10-row block is not row i in an 11-row cohort."""
    S.extract(tiles=_tiles(10), heads=[_FakeHead(HEADS[0])], block_dir=tmp_path, limit=6,
              batch=2, encode=_encode, verbose=False)
    side = S.extract(tiles=_tiles(12), heads=[_FakeHead(HEADS[0])], block_dir=tmp_path,
                     batch=4, encode=_encode, verbose=False)
    assert side["n_rows"] == 12 and side["n_done"] == 12


def test_extract_refuses_a_non_finite_latent(tmp_path, stub_head_math):
    """A NaN latent would propagate into a scalar that still has a Spearman, so it must stop the
    run rather than be reduced away."""
    def bad_encode(seqs):
        z = _encode(seqs)
        z[0, 0, 0] = np.nan
        return z

    with pytest.raises(ScoreError, match="non-finite"):
        S.extract(tiles=_tiles(4), heads=[_FakeHead(HEADS[0])], block_dir=tmp_path, batch=2,
                  encode=bad_encode, verbose=False)


def test_extract_drives_a_heartbeat_when_one_is_supplied(tmp_path, stub_head_math):
    """Without this the status surface reports NOT_STARTED for a job that has been running for
    hours -- the exact conflation the launcher's state set exists to prevent, and the surface an
    absent reader would be checking."""
    from scripts.detached_launch import Heartbeat, classify, STATE_RUNNING

    hb = Heartbeat("hb-probe", units_total=10, state_dir=tmp_path / "state",
                   cadence_seconds=3600)
    S.extract(tiles=_tiles(10), heads=[_FakeHead(HEADS[0])], block_dir=tmp_path / "blk",
              batch=5, encode=_encode, verbose=False, heartbeat=hb)
    st = classify("hb-probe", state_dir=tmp_path / "state")
    assert st.state == STATE_RUNNING, st.reason
    assert st.units_done == 10 and st.units_total == 10


def test_the_final_heartbeat_is_FORCED_past_the_cadence(tmp_path, stub_head_math):
    """A long cadence must not swallow the last beat: the final unit is the one a returning
    reader most wants confirmed, and a stale count beside a finished job reads as a stall."""
    from scripts.detached_launch import Heartbeat, read_heartbeat

    hb = Heartbeat("hb-final", units_total=10, state_dir=tmp_path / "state",
                   cadence_seconds=99_999)
    S.extract(tiles=_tiles(10), heads=[_FakeHead(HEADS[0])], block_dir=tmp_path / "blk",
              batch=5, encode=_encode, verbose=False, heartbeat=hb)
    assert read_heartbeat(hb.path)["units_done"] == 10


def test_extract_works_with_NO_heartbeat(tmp_path, stub_head_math):
    """The heartbeat is optional, so the library path and the tests above stay unaffected."""
    side = S.extract(tiles=_tiles(6), heads=[_FakeHead(HEADS[0])], block_dir=tmp_path,
                     batch=3, encode=_encode, verbose=False, heartbeat=None)
    assert side["status"] == "COMPLETE"


def test_extract_records_the_faithful_route_in_the_sidecar(tmp_path, stub_head_math):
    """The bar rejected the cheap route; the block must record which route produced it, or a
    later reader cannot tell whether the frozen choice was honoured."""
    side = S.extract(tiles=_tiles(4), heads=[_FakeHead(HEADS[0])], block_dir=tmp_path, batch=2,
                     encode=_encode, verbose=False)
    assert side["route"] == S.ROUTE_TRUE


# --------------------------------------------------------------------- gate wiring


def test_the_gate_sequence_is_derived_from_the_frozen_bar():
    """Derived, not restated: the scorer must run exactly the gates the bar froze, in that order.
    A restated copy drifts from the bar it claims to implement."""
    expected = tuple((g, o) for g, o in CHECK_ORDER if isinstance(o, str))
    assert S.GATE_SEQUENCE == expected
    assert [g for g, _ in S.GATE_SEQUENCE] == \
        ["split_equality", "geometry", "head_validity", "permutation_null"]


def test_every_gate_in_the_sequence_has_an_implementation():
    """A gate named in the order with no function raises KeyError inside `run` -- AFTER the hours
    of extraction. This is the check that makes that a test failure instead of a wasted run."""
    for gate, _ in S.GATE_SEQUENCE:
        assert gate in S.GATE_FUNCS, "gate %r is in the order with no implementation" % gate


def test_no_implemented_gate_is_missing_from_the_sequence():
    """The other direction: a gate with a function but no place in the order never runs, and its
    silence would read as a pass."""
    assert set(S.GATE_FUNCS) == {g for g, _ in S.GATE_SEQUENCE}


# --------------------------------------------------------------------- individual gates


def _bar_stub(margin=0.02, n_windows=65):
    return {"inherited_not_chosen": {"margin": {"value": margin}},
            "windowing": {"n_windows": n_windows}}


def test_split_equality_passes_when_the_sizes_reproduce():
    ctx = {"committed_n_test": {"0": 2, "1": 3}, "seeds": [0, 1],
           "splits": {0: ([], [1, 2], None), 1: ([], [1, 2, 3], None)}}
    assert S.gate_split_equality(ctx)["passed"]


def test_split_equality_fails_on_a_single_differing_seed():
    """One differing seed is enough: the comparison is against committed per-seed numbers, so a
    single mismatch means the gain for that seed is computed against another split's atlas."""
    ctx = {"committed_n_test": {"0": 2, "1": 3}, "seeds": [0, 1],
           "splits": {0: ([], [1, 2], None), 1: ([], [1, 2], None)}}
    g = S.gate_split_equality(ctx)
    assert not g["passed"] and "differ" in g["detail"]


def test_geometry_passes_on_the_pinned_fc1_and_window_count():
    ctx = {"heads": [_FakeHead(HEADS[0])], "length": 150, "bar": _bar_stub()}
    g = S.gate_geometry(ctx)
    assert g["passed"] and g["n_windows"] == 65


def test_geometry_fails_on_a_drifted_fc1():
    h = _FakeHead(HEADS[0])
    h.state_dict["fc1.weight"] = np.zeros((64, 9999))
    g = S.gate_geometry({"heads": [h], "length": 150, "bar": _bar_stub()})
    assert not g["passed"] and HEADS[0] in g["heads_with_wrong_fc1"]


def test_geometry_fails_when_the_window_count_contradicts_the_bar():
    """150 - 86 + 1 = 65 is frozen in the bar. A substrate of another length is a different arm."""
    ctx = {"heads": [_FakeHead(HEADS[0])], "length": 200, "bar": _bar_stub()}
    assert not S.gate_geometry(ctx)["passed"]


def _validity_ctx(values_by_head):
    heads = [_FakeHead(h) for h in values_by_head]
    lay = S.column_layout(list(values_by_head))
    n = len(next(iter(values_by_head.values())))
    arr = np.zeros((n, lay["_total"]))
    for h, vals in values_by_head.items():
        arr[:, lay["%s::%s" % (h, PRIMARY_REDUCTION)]] = vals
    return {"heads": heads, "layout": lay, "arr": arr}


def test_head_validity_passes_on_a_varying_finite_head():
    ctx = _validity_ctx({HEADS[0]: [1.0, 2.0, 3.0, 4.0]})
    assert S.gate_head_validity(ctx)["passed"]


def test_head_validity_FAILS_on_a_collapsed_head():
    """A head emitting one value has Spearman 0 against anything, which reads as a clean null
    rather than as a broken head."""
    ctx = _validity_ctx({HEADS[0]: [2.5, 2.5, 2.5, 2.5]})
    g = S.gate_head_validity(ctx)
    assert not g["passed"]
    assert g["per_head"][HEADS[0]]["std"] == 0.0


def test_head_validity_fails_on_a_non_finite_head():
    ctx = _validity_ctx({HEADS[0]: [1.0, 2.0, np.nan, 4.0]})
    assert not S.gate_head_validity(ctx)["passed"]


def test_head_validity_names_every_bad_head_not_just_the_first():
    ctx = _validity_ctx({HEADS[0]: [1.0, 1.0, 1.0, 1.0], HEADS[1]: [1.0, 2.0, 3.0, 4.0],
                         HEADS[2]: [0.0, 0.0, 0.0, 0.0]})
    g = S.gate_head_validity(ctx)
    assert not g["passed"]
    assert g["per_head"][HEADS[1]]["std"] > 0


# --------------------------------------------------------------------- the permutation null


def _null_ctx(n=400, *, signal=True, margin=0.02):
    rng = np.random.default_rng(3)
    y = rng.normal(size=n)
    pred = y + rng.normal(scale=0.2, size=n) if signal else rng.normal(size=n)
    lay = S.column_layout([HEADS[0]])
    arr = np.zeros((n, lay["_total"]))
    arr[:, lay["%s::%s" % (HEADS[0], PRIMARY_REDUCTION)]] = pred
    idx = list(range(n))
    from dna_decode.glm.genomewide import spearman
    real = float(spearman(list(pred), list(y)))
    return {
        "arr": arr, "layout": lay, "seeds": [0], "splits": {0: ([], idx, None)},
        "bar": _bar_stub(margin=margin), "best_head": HEADS[0], "best_median": real,
        "ys": lambda i: [float(y[k]) for k in i],
    }, real


def test_the_null_is_the_PREDICTOR_SIDE_form_the_bar_amended_to():
    """Fixed predictor, permuted TEST labels, no refit. The inherited shuffled-TRAINING-label form
    was measured vacuous on a one-feature predictor (it returned +-the entire signal)."""
    ctx, _ = _null_ctx()
    g = S.gate_permutation_null(ctx)
    assert g["form"] == "fixed predictor, permuted TEST labels, NO refit"


def test_the_null_passes_against_a_genuinely_predictive_head():
    ctx, real = _null_ctx(signal=True)
    g = S.gate_permutation_null(ctx)
    assert g["passed"], g["detail"]
    assert abs(g["median"]) < 0.02
    assert g["max_abs_over_real"] < NULL_SIGNAL_RATIO
    assert g["real_median_spearman"] == round(real, 4)


def test_the_null_reports_BOTH_conditions_separately():
    """Median near zero and magnitude small relative to the arm's own signal are different
    claims; a single boolean would hide which one failed."""
    ctx, _ = _null_ctx()
    g = S.gate_permutation_null(ctx)
    assert set(["condition_median_passed", "condition_magnitude_passed"]) <= set(g)
    assert g["passed"] == (g["condition_median_passed"] and g["condition_magnitude_passed"])


def test_the_null_FAILS_the_magnitude_condition_on_a_signal_free_head():
    """The rail that catches the vacuity. With no real signal the ratio max|null|/real blows up,
    so the gate fires even though each individual null is small."""
    ctx, real = _null_ctx(signal=False)
    g = S.gate_permutation_null(ctx)
    assert abs(real) < 0.2, "fixture precondition: this head must carry no real signal"
    assert not g["condition_magnitude_passed"]
    assert not g["passed"]


def test_the_null_reports_the_sampling_scale_it_is_judged_against():
    """1/sqrt(n_test) is the spread a permutation Spearman has by construction. Reporting it is
    what stops a null of ordinary size being read as a failure -- the mistake that sank the
    first version of this rail."""
    ctx, _ = _null_ctx(n=400)
    g = S.gate_permutation_null(ctx)
    assert g["sampling_scale_1_over_sqrt_n_test"] == round(1 / np.sqrt(400), 4)


def test_the_null_is_deterministic():
    """Fixed seeds, so a rerun of the artifact reproduces the same nulls."""
    ctx, _ = _null_ctx()
    assert S.gate_permutation_null(ctx)["per_seed"] == S.gate_permutation_null(ctx)["per_seed"]


# --------------------------------------------------------------------- score_head


def _score_ctx(n=200, *, seeds=(0, 1), atlas=0.30):
    rng = np.random.default_rng(11)
    y = rng.normal(size=n)
    pred = y + rng.normal(scale=0.5, size=n)
    lay = S.column_layout([HEADS[0]])
    arr = np.zeros((n, lay["_total"]))
    for red in S.REDUCTIONS:
        arr[:, lay["%s::%s" % (HEADS[0], red)]] = pred
    idx = list(range(n))
    return {
        "arr": arr, "layout": lay, "seeds": list(seeds),
        "splits": {s: ([], idx, None) for s in seeds},
        "committed_atlas": {str(s): atlas for s in seeds},
        "pct": [2.5, 97.5], "n_boot": 200,
        "ys": lambda i: [float(y[k]) for k in i],
    }


def test_score_head_reports_the_primary_as_eligible():
    row = S.score_head(_score_ctx(), HEADS[0], PRIMARY_REDUCTION)
    assert row["eligible_as_the_arms_result"] is True


def test_score_head_reports_the_SECONDARY_as_INELIGIBLE():
    """Pre-declared in the bar. On the positional arm the highest-scoring rung was an ineligible
    secondary, and that declaration is the only reason it could not be re-read as a win."""
    row = S.score_head(_score_ctx(), HEADS[0], SECONDARY_REDUCTION)
    assert row["eligible_as_the_arms_result"] is False


def test_score_head_computes_the_gain_against_the_COMMITTED_atlas_per_seed():
    ctx = _score_ctx(atlas=0.30)
    row = S.score_head(ctx, HEADS[0], PRIMARY_REDUCTION)
    for s, v in row["per_seed"].items():
        assert v["atlas_committed"] == 0.30
        assert v["gain"] == round(v["spearman"] - 0.30, 4)


def test_score_head_counts_seed_wins():
    ctx = _score_ctx(atlas=-1.0)     # every seed must beat an impossible-to-lose comparator
    row = S.score_head(ctx, HEADS[0], PRIMARY_REDUCTION)
    assert row["seed_wins_over_atlas"] == row["n_seeds"]


def test_score_head_does_no_fitting():
    """DIRECT Spearman, sign as given. Scoring through a one-feature ridge would borrow the sign
    from the training labels and the arm would stop being zero-shot, which is the claim its own
    name makes. Checked by anti-correlating the predictor: a fitted scorer would sign-flip and
    report a POSITIVE rho."""
    ctx = _score_ctx()
    col = ctx["layout"]["%s::%s" % (HEADS[0], PRIMARY_REDUCTION)]
    ctx["arr"][:, col] *= -1.0
    row = S.score_head(ctx, HEADS[0], PRIMARY_REDUCTION)
    assert row["median_spearman"] < 0, "the sign was not preserved as given"


# --------------------------------------------------------------------- verdict re-derivation


def _artifact(gates_pass=True, *, ci=(0.05, 0.09), missing=None, failing=None):
    gates = {}
    for gate, _ in S.GATE_SEQUENCE:
        if gate == missing:
            continue
        ok = gates_pass and gate != failing
        gates[gate] = {"passed": ok, "detail": "%s %s" % (gate, "ok" if ok else "bad")}
    rec = {"gates": gates,
           "best_primary": {"gain_bootstrap": {"ci_lo": ci[0], "ci_hi": ci[1]}}}
    return rec


BAR = {"inherited_not_chosen": {"margin": {"value": 0.02}}}


def test_the_verdict_is_rederivable_from_the_artifact_alone():
    v, _ = S.verdict_from_artifact(_artifact(ci=(0.05, 0.09)), BAR)
    assert v == "HEADS_ADDS"


def test_a_gain_below_the_margin_rederives_as_TIES():
    v, _ = S.verdict_from_artifact(_artifact(ci=(0.001, 0.019)), BAR)
    assert v == "HEADS_TIES"


def test_an_interval_below_zero_rederives_as_NO_GAIN():
    v, _ = S.verdict_from_artifact(_artifact(ci=(-0.2, -0.05)), BAR)
    assert v == "HEADS_NO_GAIN"


@pytest.mark.parametrize("gate,outcome", list(S.GATE_SEQUENCE))
def test_each_failing_gate_rederives_its_own_outcome(gate, outcome):
    v, _ = S.verdict_from_artifact(_artifact(failing=gate), BAR)
    assert v == outcome


@pytest.mark.parametrize("gate,outcome", list(S.GATE_SEQUENCE))
def test_an_UNEVALUATED_gate_rederives_as_indeterminate_not_as_a_pass(gate, outcome):
    """A gate absent from the artifact must not read as a pass. Silence is not evidence."""
    v, why = S.verdict_from_artifact(_artifact(missing=gate), BAR)
    assert v == outcome
    assert "never evaluated" in why


def test_the_gates_are_rederived_in_the_frozen_ORDER():
    """With two gates failing, the EARLIER one must be reported -- validity before decision, and
    split equality before everything, because every later number is meaningless if the test sets
    differ."""
    rec = _artifact()
    rec["gates"]["split_equality"]["passed"] = False
    rec["gates"]["permutation_null"]["passed"] = False
    v, _ = S.verdict_from_artifact(rec, BAR)
    assert v == "INDETERMINATE_SPLIT_MISMATCH"


def test_a_missing_primary_rederives_as_TIES_not_as_a_win():
    rec = _artifact()
    rec.pop("best_primary")
    v, why = S.verdict_from_artifact(rec, BAR)
    assert v == "HEADS_TIES"
    assert "no primary" in why


def test_every_rederivable_verdict_is_a_declared_outcome():
    for ci in [(-0.5, -0.1), (-0.01, 0.01), (0.0, 0.019), (0.03, 0.2)]:
        v, _ = S.verdict_from_artifact(_artifact(ci=ci), BAR)
        assert v in S.OUTCOMES
    for gate, _ in S.GATE_SEQUENCE:
        v, _ = S.verdict_from_artifact(_artifact(failing=gate), BAR)
        assert v in S.OUTCOMES


# --------------------------------------------------------------------- the frozen bar is present


def test_the_scorer_points_at_the_frozen_bar_that_exists():
    assert S.BAR_PATH.is_file(), "the scorer refuses without this; it must be frozen first"
    bar = json.loads(S.BAR_PATH.read_text(encoding="utf-8"))
    assert bar["arm"] == S.ARM


def test_the_frozen_bar_supplies_every_field_the_run_reads():
    """Read the bar through the exact key paths `run` uses, so a bar that is missing one fails
    HERE rather than after the extraction."""
    bar = json.loads(S.BAR_PATH.read_text(encoding="utf-8"))
    inh = bar["inherited_not_chosen"]
    assert inh["split"]["n_seeds"] > 0
    assert 0 < inh["split"]["frac"] < 1
    assert inh["split"]["n_test_per_seed"]
    assert inh["comparator"]["atlas_per_seed"]
    assert inh["comparator"]["atlas_median"]
    assert inh["margin"]["value"] == 0.02
    assert bar["multiplicity"]["ci_percentiles"]
    assert bar["windowing"]["n_windows"] == 65
    for key in ("scoring", "metric_rule", "honest_limits_fixed_in_advance", "forecast"):
        assert key in bar, "run() reads bar[%r]" % key


def test_the_committed_split_sizes_cover_every_seed_the_run_will_use():
    """An absent seed would raise a KeyError inside the split gate after the extraction."""
    bar = json.loads(S.BAR_PATH.read_text(encoding="utf-8"))
    inh = bar["inherited_not_chosen"]
    for s in range(inh["split"]["n_seeds"]):
        assert str(s) in inh["split"]["n_test_per_seed"]
        assert str(s) in inh["comparator"]["atlas_per_seed"]
