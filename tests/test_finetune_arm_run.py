"""Tests for the long-horizon fine-tune runner.

THE CENTRAL PROPERTY IS TESTED, NOT ASSUMED. The whole design rests on one claim: the frozen
budget's verdict read off a long run's PREFIX equals the verdict of a run that stopped at that
budget. If that is false, the longer horizons are not free diagnostics -- they are a different
experiment, and the arm would need three separate runs at 7x the cost. So it is checked by
training both ways on a tiny synthetic cohort and comparing.

A REAL ENCODER IS NOT USED. `build` is injectable, and these tests pass a tiny conv+linear pair.
The PromoterAtlas weights are 1.8M parameters and a full-cohort epoch is 20 minutes; nothing here
is about the model's capacity, it is about the run's bookkeeping -- resume, budget prefix,
wall-clock bound, and the refusal of a wrong target transform.
"""
from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

import scripts.finetune_arm_run as R

torch = pytest.importorskip("torch")


# --------------------------------------------------------------------- fixtures


def _tiles(n_peaks=12, per_peak=5, length=24, seed=0):
    """A tiny cohort whose expression is a learnable function of the sequence.

    Learnable on purpose: a curve that cannot move would make every reproduction test pass
    trivially, and the resume check would prove nothing.
    """
    rng = np.random.default_rng(seed)
    tiles = []
    for p in range(n_peaks):
        for _ in range(per_peak):
            s = "".join(rng.choice(list("ACGT"), size=length))
            expr = float(s.count("G") + s.count("C"))
            tiles.append(SimpleNamespace(seq=s, expression=expr, peak="peak%d" % p,
                                         category="tile"))
    return tiles


class _TinyEnc(torch.nn.Module):
    def __init__(self):
        super().__init__()
        self.c = torch.nn.Conv1d(4, 8, 3, padding=1)

    def forward(self, x, return_latent=False):
        z = torch.relu(self.c(x))
        return (None, None, z) if return_latent else z


def _build():
    return _TinyEnc(), torch.nn.Linear(8, 1)


BAR = {
    "frozen_hyperparameters": {"n_epochs_budget": 4, "patience": 2, "learning_rate": 0.01,
                               "batch_size": 8},
    "inherited_not_chosen": {"comparator": {"n_seeds": 2, "split_frac": 0.25}},
    "ceiling": {"ceiling_hours": 1.0},
    "horizon": {"minutes_per_epoch": 0.001, "diagnostic_horizon_epochs": 10},
    "route": {"chosen": "local_cpu"},
}


def _run(ckpt_dir, **kw):
    kw.setdefault("n_seeds", 2)
    kw.setdefault("horizon_epochs", 8)
    kw.setdefault("n_null_seeds", 1)
    return R.run(bar=BAR, tiles=_tiles(), ckpt_dir=Path(ckpt_dir), build=_build,
                 verbose=False, **kw)


# --------------------------------------------------------------------- the target transform


def test_the_log1p_transform_is_what_the_runner_uses():
    got = R._targets([0.0, 1.0, 3.0])
    assert torch.allclose(got, torch.log1p(torch.tensor([0.0, 1.0, 3.0])))


@pytest.mark.parametrize("bad", ["raw", "none", "log", "identity", "log2"])
def test_training_on_anything_but_log1p_is_REFUSED(bad):
    """The C3 correction, pinned. Training on RAW expression would confound the objective
    transform with the representation and add a third reading on top of the
    budget-versus-representation pair the bar names."""
    with pytest.raises(R.RunRefused, match="refused"):
        R._targets([1.0, 2.0], bad)


def test_the_refusal_explains_WHY_not_just_THAT():
    with pytest.raises(R.RunRefused) as e:
        R._targets([1.0], "raw")
    msg = str(e.value)
    assert "ridge_spearman" in msg
    assert "OBJECTIVE TRANSFORM" in msg and "REPRESENTATION" in msg


def test_a_seed_trained_with_a_bad_transform_cannot_be_started(tmp_path):
    with pytest.raises(R.RunRefused):
        R.train_one_seed(_tiles(), seed=0, horizon_epochs=2, lr=0.01, batch_size=8,
                         split_frac=0.25, ckpt_dir=tmp_path, build=_build,
                         transform="raw", verbose=False)


# --------------------------------------------------------------------- the frozen prefix rule


def test_the_frozen_budget_reads_only_the_prefix():
    curve = [0.1, 0.2, 0.3, 0.25, 0.9, 0.95]
    got = R.frozen_budget_value(curve, frozen_epochs=4, patience=2)
    assert got["value"] == 0.3
    assert got["best_epoch"] == 2
    assert got["n_epochs_read"] == 4


def test_a_SPECTACULAR_LATE_EPOCH_CANNOT_REACH_THE_VERDICT():
    """The property that makes the longer horizons diagnostics rather than a bigger arm. Epoch 9
    here would win outright and must be invisible to the frozen rule."""
    curve = [0.01] * 4 + [0.99] * 6
    got = R.frozen_budget_value(curve, frozen_epochs=4, patience=2)
    assert got["value"] == 0.01
    assert got["value"] != max(curve)


def test_an_empty_prefix_yields_no_value_rather_than_zero():
    """None, not 0.0 -- a missing measurement must not read as a measured absence of signal,
    which is the same trap as the scorer's NaN-to-0.0 mapping."""
    got = R.frozen_budget_value([], frozen_epochs=4, patience=2)
    assert got["value"] is None


def test_the_frozen_rule_detects_early_stopping_inside_the_budget():
    curve = [0.1, 0.5, 0.4, 0.3]
    got = R.frozen_budget_value(curve, frozen_epochs=4, patience=2)
    assert got["best_epoch"] == 1 and got["stopped_early"] is True


def test_the_frozen_rule_reports_no_early_stop_when_still_improving():
    got = R.frozen_budget_value([0.1, 0.2, 0.3, 0.4], frozen_epochs=4, patience=2)
    assert got["stopped_early"] is False


def test_a_nonpositive_budget_is_refused():
    with pytest.raises(R.RunRefused):
        R.frozen_budget_value([0.1], frozen_epochs=0, patience=2)


def test_the_diagnostic_tail_is_marked_INELIGIBLE():
    d = R.diagnostic_tail([0.1] * 4 + [0.9] * 6, frozen_epochs=4)
    assert d["eligible_as_the_arms_result"] is False
    assert "INELIGIBLE" in d["status"]
    assert d["n_epochs"] == 6
    assert d["best_beyond_budget"] == 0.9
    assert d["best_inside_budget"] == 0.1
    assert d["improvement_beyond_budget"] == 0.8


def test_the_diagnostic_tail_is_empty_when_the_run_stopped_at_the_budget():
    d = R.diagnostic_tail([0.1, 0.2, 0.3, 0.4], frozen_epochs=4)
    assert d["n_epochs"] == 0 and d["best_beyond_budget"] is None
    assert d["improvement_beyond_budget"] is None


# --------------------------------------------------------------------- planning


def test_plan_run_reserves_the_nulls_BEFORE_sizing_the_horizon():
    """Sizing against the whole ceiling and then adding the null would overrun the window or drop
    the null, and the bar requires the null."""
    with_null = R.plan_run(ceiling_hours=100.0, seconds_per_epoch=3600.0, n_seeds=10,
                           frozen_epochs=7, n_null_seeds=3)
    without = R.plan_run(ceiling_hours=100.0, seconds_per_epoch=3600.0, n_seeds=10,
                         frozen_epochs=7, n_null_seeds=0)
    assert with_null["horizon_epochs"] < without["horizon_epochs"]
    assert with_null["null_epochs_reserved"] == 21


def test_the_plan_fits_inside_its_ceiling():
    p = R.plan_run(ceiling_hours=129.6, seconds_per_epoch=1225.2, n_seeds=10,
                   frozen_epochs=7, n_null_seeds=3)
    assert p["fits"] is True
    assert p["projected_hours"] <= p["ceiling_hours"]
    assert p["horizon_epochs"] > p["frozen_budget_epochs"]


def test_a_ceiling_too_small_for_the_frozen_budget_REFUSES_to_start():
    """Better to refuse than to start a multi-day run that cannot reach its own verdict."""
    with pytest.raises(R.RunRefused, match="must not start"):
        R.plan_run(ceiling_hours=2.0, seconds_per_epoch=1225.2, n_seeds=10,
                   frozen_epochs=7, n_null_seeds=3)


def test_the_horizon_shrinks_with_the_ceiling():
    big = R.plan_run(ceiling_hours=129.6, seconds_per_epoch=1225.2, n_seeds=10, frozen_epochs=7)
    small = R.plan_run(ceiling_hours=60.0, seconds_per_epoch=1225.2, n_seeds=10, frozen_epochs=7)
    assert small["horizon_epochs"] < big["horizon_epochs"]


@pytest.mark.parametrize("kw", [
    {"ceiling_hours": 0}, {"seconds_per_epoch": 0}, {"n_seeds": 0}, {"frozen_epochs": 0}])
def test_plan_run_refuses_nonpositive_inputs(kw):
    base = dict(ceiling_hours=100.0, seconds_per_epoch=100.0, n_seeds=10, frozen_epochs=7)
    base.update(kw)
    with pytest.raises(R.RunRefused):
        R.plan_run(**base)


# --------------------------------------------------------------------- the resume contract


def test_the_epoch_shuffle_seed_depends_only_on_seed_and_epoch():
    """The whole resume contract. Epoch 20's data order must not depend on whether epochs 0-19
    ran in this process or in one that died days ago."""
    assert R.epoch_generator_seed(3, 20) == R.epoch_generator_seed(3, 20)
    assert R.epoch_generator_seed(3, 20) != R.epoch_generator_seed(3, 21)
    assert R.epoch_generator_seed(3, 20) != R.epoch_generator_seed(4, 20)


def test_the_epoch_shuffle_seeds_are_wide_apart_across_seeds_and_epochs():
    """A collision would silently give two different (seed, epoch) pairs the same data order."""
    vals = {R.epoch_generator_seed(s, e) for s in range(12) for e in range(60)}
    assert len(vals) == 12 * 60


def test_a_RESUMED_seed_reproduces_an_uninterrupted_one_BIT_FOR_BIT(tmp_path):
    """The property a ten-day run depends on, checked by value rather than by argument."""
    tiles = _tiles()
    kw = dict(lr=0.01, batch_size=8, split_frac=0.25, build=_build, verbose=False)

    whole = R.train_one_seed(tiles, seed=0, horizon_epochs=6, ckpt_dir=tmp_path / "a", **kw)

    part = R.train_one_seed(tiles, seed=0, horizon_epochs=3, ckpt_dir=tmp_path / "b", **kw)
    assert len(part["val_spearman_by_epoch"]) == 3
    resumed = R.train_one_seed(tiles, seed=0, horizon_epochs=6, ckpt_dir=tmp_path / "b", **kw)

    assert resumed["val_spearman_by_epoch"] == whole["val_spearman_by_epoch"], (
        "resume diverged:\n  whole   %r\n  resumed %r"
        % (whole["val_spearman_by_epoch"], resumed["val_spearman_by_epoch"]))


def test_a_resumed_seed_does_not_recompute_the_epochs_it_already_has(tmp_path):
    tiles = _tiles()
    kw = dict(lr=0.01, batch_size=8, split_frac=0.25, build=_build, verbose=False)
    R.train_one_seed(tiles, seed=0, horizon_epochs=3, ckpt_dir=tmp_path, **kw)
    before = R.read_curve(tmp_path, 0)["val_spearman_by_epoch"]
    resumed = R.train_one_seed(tiles, seed=0, horizon_epochs=5, ckpt_dir=tmp_path, **kw)
    assert resumed["val_spearman_by_epoch"][:3] == before


def test_a_FINISHED_seed_is_read_not_retrained(tmp_path):
    """Measured by making retraining IMPOSSIBLE: if the run touched the finished seed it would
    raise, so a clean return proves it was read from the committed curve."""
    _run(tmp_path, horizon_epochs=3, n_seeds=1, n_null_seeds=0)
    first = R.read_curve(tmp_path, 0)["val_spearman_by_epoch"]

    def exploding_build():
        raise AssertionError("a finished seed was retrained")

    again = R.run(bar=BAR, tiles=_tiles(), ckpt_dir=Path(tmp_path), build=exploding_build,
                  n_seeds=1, horizon_epochs=3, n_null_seeds=0, verbose=False)
    assert again["per_seed"]["0"]["val_spearman_by_epoch"] == first


def test_only_the_UNFINISHED_seeds_are_retrained(tmp_path):
    """Seed 0 finished, seed 1 absent: the resumed run must train exactly one of them."""
    R.train_one_seed(_tiles(), seed=0, horizon_epochs=3, lr=0.01, batch_size=8,
                     split_frac=0.25, ckpt_dir=tmp_path, build=_build, verbose=False)
    trained = []

    def counting_build():
        trained.append(1)
        return _build()

    R.run(bar=BAR, tiles=_tiles(), ckpt_dir=Path(tmp_path), build=counting_build,
          n_seeds=2, horizon_epochs=3, n_null_seeds=0, verbose=False)
    assert sum(trained) == 1, "trained %d seeds; only seed 1 was outstanding" % sum(trained)


def test_a_malformed_curve_file_reads_as_absent(tmp_path):
    tmp_path.mkdir(parents=True, exist_ok=True)
    (tmp_path / "seed_0_curve.json").write_text('{"seed": 0, "val', encoding="utf-8")
    assert R.read_curve(tmp_path, 0) is None


def test_the_curve_json_is_readable_without_loading_the_checkpoint(tmp_path):
    """A small sidecar beside the heavy tensor file, so progress can be inspected -- and the
    artifact assembled -- without torch or a 20 MB load."""
    R.train_one_seed(_tiles(), seed=0, horizon_epochs=2, lr=0.01, batch_size=8,
                     split_frac=0.25, ckpt_dir=tmp_path, build=_build, verbose=False)
    raw = json.loads((tmp_path / "seed_0_curve.json").read_text(encoding="utf-8"))
    assert len(raw["val_spearman_by_epoch"]) == 2
    assert raw["target_transform"] == "log1p"


def test_the_checkpoint_leaves_no_tmp_file_behind(tmp_path):
    R.train_one_seed(_tiles(), seed=0, horizon_epochs=2, lr=0.01, batch_size=8,
                     split_frac=0.25, ckpt_dir=tmp_path, build=_build, verbose=False)
    assert list(tmp_path.glob("*.tmp")) == []


# --------------------------------------------------------------------- prefix == stopped run


def test_THE_FROZEN_VERDICT_OFF_A_LONG_PREFIX_EQUALS_A_RUN_STOPPED_AT_THE_BUDGET(tmp_path):
    """The claim the entire design rests on.

    If this is false the longer horizons are not free diagnostics but a different experiment,
    and the arm would need three runs at sevenfold the cost. Trained both ways from the same
    seed, so only the horizon differs.
    """
    tiles = _tiles()
    kw = dict(lr=0.01, batch_size=8, split_frac=0.25, build=_build, verbose=False)
    budget = 4

    stopped = R.train_one_seed(tiles, seed=1, horizon_epochs=budget,
                               ckpt_dir=tmp_path / "stopped", **kw)
    long_run = R.train_one_seed(tiles, seed=1, horizon_epochs=12,
                                ckpt_dir=tmp_path / "long", **kw)

    a = R.frozen_budget_value(stopped["val_spearman_by_epoch"], frozen_epochs=budget, patience=2)
    b = R.frozen_budget_value(long_run["val_spearman_by_epoch"], frozen_epochs=budget, patience=2)

    assert long_run["val_spearman_by_epoch"][:budget] == stopped["val_spearman_by_epoch"], (
        "the long run's prefix diverged from the stopped run; the horizons are not a shared "
        "trajectory")

    # The DECISION-BEARING fields only. `n_epochs_available` and `detail` legitimately differ --
    # the long run really does have 12 epochs and 8 beyond the budget, and reporting that is
    # correct. Comparing the whole dict conflated a descriptive difference with a verdict
    # difference, which is what the first version of this test did.
    decision = ("value", "best_epoch", "n_epochs_read", "stopped_early")
    assert {k: a[k] for k in decision} == {k: b[k] for k in decision}, (
        "the frozen verdict differs between the prefix and the stopped run: %r vs %r"
        % ({k: a[k] for k in decision}, {k: b[k] for k in decision}))
    assert a["n_epochs_available"] != b["n_epochs_available"], (
        "non-vacuity: the two runs must differ in length, or this test compares a run with "
        "itself")


def test_the_long_run_records_more_epochs_than_the_budget(tmp_path):
    """Non-vacuity for the test above: if the long run had only `budget` epochs the equality
    would hold trivially and prove nothing."""
    long_run = R.train_one_seed(_tiles(), seed=1, horizon_epochs=12, lr=0.01, batch_size=8,
                                split_frac=0.25, ckpt_dir=tmp_path, build=_build, verbose=False)
    assert len(long_run["val_spearman_by_epoch"]) == 12 > 4


def test_the_curves_actually_MOVE(tmp_path):
    """Non-vacuity for every reproduction test: a flat curve would make them all pass."""
    rec = R.train_one_seed(_tiles(), seed=0, horizon_epochs=6, lr=0.01, batch_size=8,
                           split_frac=0.25, ckpt_dir=tmp_path, build=_build, verbose=False)
    curve = rec["val_spearman_by_epoch"]
    assert len(set(curve)) > 1, "the synthetic curve is constant: %r" % curve


# --------------------------------------------------------------------- the wall clock


def test_the_wall_clock_stop_FIRES_and_stamps_the_run_BOUNDED(tmp_path):
    """A run that stopped early must never be readable as a complete one: overrunning the window
    leaves every seed half-done, which is worse than fewer seeds finished cleanly."""
    rec = R.run(bar=BAR, tiles=_tiles(), ckpt_dir=Path(tmp_path), build=_build,
                n_seeds=2, horizon_epochs=8, n_null_seeds=1, verbose=False,
                wall_ceiling_hours=1e-9, now=0.0)
    assert rec["status"] == R.STATUS_BOUNDED
    assert rec["status"] != R.STATUS_COMPLETE
    assert "NOT readable as a complete one" in rec["why_status"]


def test_an_unbounded_run_stamps_COMPLETE(tmp_path):
    rec = _run(tmp_path, horizon_epochs=3, n_seeds=2, n_null_seeds=1)
    assert rec["status"] == R.STATUS_COMPLETE
    assert rec["n_seeds_done"] == 2


def test_a_bounded_run_still_reports_the_seeds_it_DID_finish(tmp_path):
    """Partial progress is kept: the checkpoints are the point of a resumable design."""
    R.train_one_seed(_tiles(), seed=0, horizon_epochs=3, lr=0.01, batch_size=8,
                     split_frac=0.25, ckpt_dir=tmp_path, build=_build, verbose=False)
    rec = R.run(bar=BAR, tiles=_tiles(), ckpt_dir=Path(tmp_path), build=_build,
                n_seeds=2, horizon_epochs=3, n_null_seeds=0, verbose=False,
                wall_ceiling_hours=1e-9, now=0.0)
    assert rec["status"] == R.STATUS_BOUNDED
    assert "0" in rec["per_seed"]


def test_the_wall_clock_stops_at_an_EPOCH_boundary_not_mid_epoch(tmp_path):
    """Mid-epoch truncation would leave the checkpoint and the curve disagreeing about how much
    training the weights have had."""
    rec = R.train_one_seed(_tiles(), seed=0, horizon_epochs=8, lr=0.01, batch_size=8,
                           split_frac=0.25, ckpt_dir=tmp_path, build=_build, verbose=False,
                           deadline=0.0)
    assert rec["val_spearman_by_epoch"] == []
    assert rec["bounded_by_wall_clock"] is True


# --------------------------------------------------------------------- the null


def test_the_null_shuffles_only_the_TRAINING_labels(tmp_path):
    """Everything else -- split, data order, hyperparameters -- is untouched, so the null differs
    from the real run in exactly one respect."""
    tiles = _tiles()
    kw = dict(seed=0, horizon_epochs=4, lr=0.01, batch_size=8, split_frac=0.25,
              build=_build, verbose=False)
    real = R.train_one_seed(tiles, ckpt_dir=tmp_path / "r", **kw)
    null = R.train_one_seed(tiles, ckpt_dir=tmp_path / "n", shuffle_labels=True, kind="null",
                            **kw)
    assert real["n_train"] == null["n_train"] and real["n_val"] == null["n_val"]
    assert null["shuffled_labels"] is True and real["shuffled_labels"] is False


def test_the_shuffle_genuinely_permutes_the_training_labels(tmp_path):
    """The mechanism, which IS checkable here: the permutation must actually reorder. Measured
    rather than assumed, because a null whose shuffle silently no-ops is the worst case -- it
    reports a high value and looks like a finding."""
    tiles = _tiles()
    tr, _va, _rest = R.leave_peak_out_indices(tiles, frac=0.25, seed=0)
    raw = [float(tiles[i].expression) for i in tr]
    rng = np.random.default_rng(770_000 + 0)
    shuffled = list(np.array(raw)[rng.permutation(len(raw))])
    assert raw != shuffled
    assert sorted(raw) == sorted(shuffled), "a permutation must preserve the multiset"


def test_THIS_FIXTURE_CANNOT_SETTLE_WHETHER_THE_NULL_IS_VACUOUS(tmp_path):
    """A deliberate negative result about the TEST, recorded rather than worked around.

    The bar states the shuffled-label null is non-vacuous here "because the fine-tune genuinely
    fits". This fixture cannot check that claim, and the reason is worth pinning: its UNTRAINED
    model already scores ~0.76, because mean-pooling makes the representation essentially base
    composition and the synthetic target IS base composition. A null that fails to collapse
    against that baseline says nothing about the fine-tune at all.

    So the assertion is the precondition, not the conclusion. Whether the REAL arm's null
    collapses is an empirical question the real run answers -- which is why the run now records
    an untrained baseline beside every null, so a high null is diagnosable instead of ambiguous.
    """
    tiles = _tiles()
    rec = R.train_one_seed(tiles, seed=0, horizon_epochs=2, lr=0.01, batch_size=8,
                           split_frac=0.25, ckpt_dir=tmp_path, build=_build, verbose=False)
    assert rec["untrained_val_spearman"] is not None
    assert rec["untrained_val_spearman"] > 0.5, (
        "this fixture's untrained baseline is %r; if it ever drops near zero the fixture COULD "
        "test vacuity and this test should be replaced by one that does"
        % rec["untrained_val_spearman"])


def test_the_untrained_baseline_is_recorded_for_every_seed(tmp_path):
    """It is the field that distinguishes a fine-tune that learned unjustified structure from a
    random readout of an already-informative representation. Both produce a high null."""
    rec = _run(tmp_path, horizon_epochs=3, n_seeds=2, n_null_seeds=1)
    for s in rec["per_seed"]:
        assert rec["per_seed"][s]["untrained_val_spearman"] is not None
    assert rec["null"]["untrained_baseline_per_replicate"]
    assert "TWO causes" in rec["null"]["how_to_read_a_HIGH_null"]


def test_the_untrained_baseline_survives_a_resume(tmp_path):
    """It is measured once, before the first gradient step, so a resumed seed must carry the
    ORIGINAL value rather than re-measuring it against partly-trained weights."""
    kw = dict(seed=0, lr=0.01, batch_size=8, split_frac=0.25, build=_build, verbose=False,
              ckpt_dir=tmp_path)
    first = R.train_one_seed(_tiles(), horizon_epochs=2, **kw)
    resumed = R.train_one_seed(_tiles(), horizon_epochs=4, **kw)
    assert resumed["untrained_val_spearman"] == first["untrained_val_spearman"]


def test_the_null_is_stored_under_its_own_kind_so_it_cannot_overwrite_a_real_seed(tmp_path):
    """Shared-key collision: a null written over seed 0's checkpoint would silently replace the
    arm's own result with its control."""
    kw = dict(seed=0, horizon_epochs=2, lr=0.01, batch_size=8, split_frac=0.25,
              build=_build, verbose=False, ckpt_dir=tmp_path)
    R.train_one_seed(_tiles(), **kw)
    R.train_one_seed(_tiles(), shuffle_labels=True, kind="null", **kw)
    assert (tmp_path / "seed_0_curve.json").is_file()
    assert (tmp_path / "null_0_curve.json").is_file()
    assert R.read_curve(tmp_path, 0)["shuffled_labels"] is False
    assert R.read_curve(tmp_path, 0, kind="null")["shuffled_labels"] is True


def test_the_run_reports_the_null_scope_limit(tmp_path):
    rec = _run(tmp_path, horizon_epochs=3, n_null_seeds=1)
    assert "NOT enough to estimate a null distribution" in rec["null"]["scope_limit"]
    assert rec["null"]["kind"] == "shuffled_training_labels"


# --------------------------------------------------------------------- the assembled record


def test_the_record_separates_the_frozen_budget_from_the_diagnostic_tail(tmp_path):
    rec = _run(tmp_path, horizon_epochs=8, n_seeds=1, n_null_seeds=1)
    row = rec["per_seed"]["0"]
    assert row["frozen_budget"]["n_epochs_read"] == 4
    assert row["diagnostic_tail"]["eligible_as_the_arms_result"] is False
    assert row["diagnostic_tail"]["n_epochs"] == 4


def test_the_record_states_that_late_epochs_cannot_reach_the_verdict(tmp_path):
    rec = _run(tmp_path, horizon_epochs=6, n_seeds=1, n_null_seeds=1)
    assert "cannot reach" in rec["ineligible_diagnostics_note"]
    assert "INELIGIBLE" in rec["ineligible_diagnostics_note"]


def test_the_record_names_the_target_transform(tmp_path):
    rec = _run(tmp_path, horizon_epochs=3, n_seeds=1, n_null_seeds=0)
    assert rec["target_transform"] == "log1p"


def test_the_record_reports_peak_blocked_splits_with_no_shared_peaks(tmp_path):
    """Peak-blocked, or a tile from the same peak in both halves leaks."""
    rec = _run(tmp_path, horizon_epochs=2, n_seeds=1, n_null_seeds=0)
    assert rec["per_seed"]["0"]["peaks_shared_between_splits"] == 0


def test_a_non_finite_validation_prediction_is_REFUSED(tmp_path):
    """The scorer maps NaN to 0.0, so divergence would otherwise be recorded as 'no signal'."""
    class _NanEnc(torch.nn.Module):
        def forward(self, x, return_latent=False):
            z = torch.full((x.shape[0], 8, x.shape[2]), float("nan"))
            return (None, None, z) if return_latent else z

    with pytest.raises(R.RunRefused, match="non-finite"):
        R.train_one_seed(_tiles(), seed=0, horizon_epochs=1, lr=0.01, batch_size=8,
                         split_frac=0.25, ckpt_dir=tmp_path,
                         build=lambda: (_NanEnc(), torch.nn.Linear(8, 1)), verbose=False)


def test_the_heartbeat_unit_count_is_monotone_across_seeds():
    """So the progress fraction rises rather than resetting at each seed boundary."""
    assert R.heartbeat_units(0, 5, 10) < R.heartbeat_units(1, 0, 10)
    assert R.heartbeat_units(1, 3, 10) < R.heartbeat_units(1, 4, 10)
