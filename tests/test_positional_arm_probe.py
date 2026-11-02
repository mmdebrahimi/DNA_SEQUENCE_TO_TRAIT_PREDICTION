"""Tests for the positional arm's feasibility probe.

Everything load-bearing here is a PURE helper, so the file is green on a host with no checkpoint,
no D: and no substrate. The model-touching path is exercised only by the real run and is skipped
cleanly when the clone does not resolve -- see the gated tests at the bottom.

The two properties worth the most are `fits_in_ram` returning None for UNKNOWN (never True, because
a probe that could not measure RAM has not shown a rung is affordable) and `choose_lag_set`
REFUSING when nothing fits (rather than quietly scoring a lag set nobody chose).
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pytest

from scripts.positional_arm_probe import (
    EXPECTED_CHANNELS,
    EXPECTED_TILES,
    LAG_CANDIDATES,
    PARITY_ATOL,
    RAM_HEADROOM,
    ProbeFailure,
    choose_lag_set,
    design_matrix_bytes,
    extrapolate_wall_seconds,
    fits_in_ram,
    free_ram_bytes,
    main,
    parity_report,
    rung_feature_counts,
    stored_block_bytes,
    stride_sample,
)

REPO = Path(__file__).resolve().parents[1]
ARTIFACT = REPO / "wiki/positional_arm_feasibility_2026-11-01.json"
PA_REPO = Path(os.environ.get("PA_REPO", r"C:/Users/LAPTOP~1/AppData/Local/Temp/pa_probe/PromoterAtlas"))
CKPT = PA_REPO / "trained_weights/base_model/promoteratlas-base.pt"


# --------------------------------------------------------------------------------------------------
# byte arithmetic
# --------------------------------------------------------------------------------------------------
def test_design_matrix_bytes_is_rows_times_features_times_itemsize_times_copies():
    assert design_matrix_bytes(100, 10, itemsize=8, copies=2) == 100 * 10 * 8 * 2
    assert design_matrix_bytes(1, 1, itemsize=1, copies=1) == 1


def test_design_matrix_bytes_reproduces_the_measured_largest_rung():
    """44,106 rows x 2,560 autocorr features, float64, two copies -> the 1.807 GB the probe printed."""
    got = design_matrix_bytes(EXPECTED_TILES, EXPECTED_CHANNELS * 20)
    assert round(got / 1e9, 3) == 1.807


def test_stored_block_bytes_reproduces_the_measured_block():
    """5,120 float32 per tile over 44,106 tiles -> the ~0.90 GB the plan predicted."""
    assert round(stored_block_bytes(EXPECTED_TILES, 5120) / 1e9, 2) == 0.90


@pytest.mark.parametrize("kwargs", [
    {"n_rows": 0, "n_features": 10}, {"n_rows": 10, "n_features": 0},
    {"n_rows": -1, "n_features": 10},
])
def test_design_matrix_bytes_refuses_degenerate_shapes(kwargs):
    with pytest.raises(ValueError):
        design_matrix_bytes(**kwargs)


def test_stored_block_bytes_refuses_degenerate_shapes():
    with pytest.raises(ValueError):
        stored_block_bytes(0, 10)
    with pytest.raises(ValueError):
        stored_block_bytes(10, 0)


# --------------------------------------------------------------------------------------------------
# fits_in_ram -- UNKNOWN must never read as affordable
# --------------------------------------------------------------------------------------------------
def test_fits_in_ram_applies_the_headroom_factor():
    # need 100 bytes, headroom 1.5 -> requires 150 free
    assert fits_in_ram(100, 150, headroom=1.5) is True
    assert fits_in_ram(100, 149, headroom=1.5) is False
    assert fits_in_ram(100, 100, headroom=1.0) is True


def test_fits_in_ram_returns_None_when_free_ram_is_unknown_and_never_True():
    """The whole point: an unmeasured host must not be treated as a host with room."""
    for need in (1, 10**6, 10**12):
        got = fits_in_ram(need, None)
        assert got is None
        assert got is not True


def test_fits_in_ram_refuses_a_headroom_below_one():
    with pytest.raises(ValueError, match="headroom"):
        fits_in_ram(100, 1000, headroom=0.9)


def test_free_ram_bytes_is_either_a_positive_int_or_None():
    got = free_ram_bytes()
    assert got is None or (isinstance(got, int) and got > 0)


# --------------------------------------------------------------------------------------------------
# extrapolation
# --------------------------------------------------------------------------------------------------
def test_extrapolate_wall_seconds_is_linear_in_the_cohort_size():
    assert extrapolate_wall_seconds(1000, 100, 10.0) == pytest.approx(100.0)
    assert extrapolate_wall_seconds(100, 100, 7.5) == pytest.approx(7.5)


def test_extrapolate_wall_seconds_refuses_a_zero_sample_or_zero_elapsed():
    with pytest.raises(ValueError, match="n_sample"):
        extrapolate_wall_seconds(100, 0, 5.0)
    with pytest.raises(ValueError, match="elapsed_seconds"):
        extrapolate_wall_seconds(100, 10, 0.0)


# --------------------------------------------------------------------------------------------------
# rung feature counts -- the figures the bar inherits
# --------------------------------------------------------------------------------------------------
def test_rung_feature_counts_match_the_geometry_the_plan_names():
    c = rung_feature_counts(20)
    assert c["bins10_mean_oriented"] == 1280
    assert c["bins10_mean_asgiven"] == 1280
    assert c["bins10_max_oriented"] == 1280
    assert c["bins5_mean_oriented"] == 640
    assert c["bins2_mean_oriented"] == 256
    assert c["autocorr"] == 2560


def test_rung_feature_counts_scale_with_the_lag_set_only_on_the_autocorr_rung():
    a, b = rung_feature_counts(20), rung_feature_counts(9)
    assert a["autocorr"] == 2560 and b["autocorr"] == 1152
    for k in a:
        if k != "autocorr":
            assert a[k] == b[k], "only the autocorr rung may depend on the lag set"


def test_rung_feature_counts_refuses_degenerate_inputs():
    with pytest.raises(ValueError):
        rung_feature_counts(0)
    with pytest.raises(ValueError):
        rung_feature_counts(5, channels=0)


# --------------------------------------------------------------------------------------------------
# choose_lag_set -- the pre-declared rule, and its refusal
# --------------------------------------------------------------------------------------------------
def test_choose_lag_set_takes_the_largest_candidate_when_ram_is_plentiful():
    name, lags, dec = choose_lag_set(
        LAG_CANDIDATES, n_rows=EXPECTED_TILES, channels=EXPECTED_CHANNELS, free_bytes=10**12,
    )
    assert name == "dense_1_20"
    assert lags == list(range(1, 21))
    assert len(dec) == len(LAG_CANDIDATES), "every candidate must be recorded, not just the winner"
    assert [d["n_lags"] for d in dec] == sorted([d["n_lags"] for d in dec], reverse=True)


def test_choose_lag_set_steps_down_when_the_largest_does_not_fit():
    """Only ever gives up resolution to a resource limit -- never for a nicer number."""
    need_log9 = design_matrix_bytes(EXPECTED_TILES, EXPECTED_CHANNELS * 9)
    free = int(need_log9 * RAM_HEADROOM)  # exactly enough for log_spaced_9, not for dense
    name, lags, dec = choose_lag_set(
        LAG_CANDIDATES, n_rows=EXPECTED_TILES, channels=EXPECTED_CHANNELS, free_bytes=free,
    )
    assert name == "log_spaced_9"
    assert len(lags) == 9
    assert [d for d in dec if d["candidate"] == "dense_1_20"][0]["fits_in_ram"] is False


def test_choose_lag_set_refuses_rather_than_picking_a_set_nobody_chose():
    with pytest.raises(ProbeFailure, match="no candidate"):
        choose_lag_set(LAG_CANDIDATES, n_rows=EXPECTED_TILES, channels=EXPECTED_CHANNELS,
                       free_bytes=1)


def test_choose_lag_set_refuses_when_ram_is_unmeasured():
    """`None` is UNKNOWN; it must not be promoted to a fit."""
    with pytest.raises(ProbeFailure, match="no candidate"):
        choose_lag_set(LAG_CANDIDATES, n_rows=EXPECTED_TILES, channels=EXPECTED_CHANNELS,
                       free_bytes=None)


def test_every_declared_lag_candidate_is_legal_for_the_real_geometry():
    """A candidate that the reduction itself would refuse must never be selectable."""
    from dna_decode.glm.positional_readout import MIN_OVERLAP, channel_autocorrelation

    z = np.zeros((1, 2, 150))
    z[0, 0, :] = np.arange(150)
    for name, lags in LAG_CANDIDATES.items():
        assert len(set(lags)) == len(lags), "%s has duplicate lags" % name
        assert min(lags) >= 1, "%s includes lag 0 or negative" % name
        assert 150 - max(lags) >= MIN_OVERLAP, "%s violates MIN_OVERLAP at L=150" % name
        assert channel_autocorrelation(z, lags).shape == (1, 2 * len(lags))


# --------------------------------------------------------------------------------------------------
# stride_sample
# --------------------------------------------------------------------------------------------------
def test_stride_sample_is_sorted_unique_in_range_and_spans_the_cohort():
    idx = stride_sample(44106, 512)
    assert idx == sorted(set(idx))
    assert len(idx) == 512
    assert 0 <= min(idx) and max(idx) <= 44105
    assert min(idx) < 100, "a stride must start near the beginning"
    assert max(idx) > 44106 - 200, "a stride must reach the end -- a head slice would not"


def test_stride_sample_is_deterministic():
    assert stride_sample(1000, 37) == stride_sample(1000, 37)


def test_stride_sample_saturates_at_the_cohort_size():
    assert stride_sample(10, 50) == list(range(10))


def test_stride_sample_refuses_degenerate_inputs():
    with pytest.raises(ValueError):
        stride_sample(0, 10)
    with pytest.raises(ValueError):
        stride_sample(10, 0)


# --------------------------------------------------------------------------------------------------
# parity_report -- the gate that makes every rung comparable to 0.5076
# --------------------------------------------------------------------------------------------------
def test_parity_report_passes_on_identical_blocks():
    rng = np.random.default_rng(1)
    a = rng.normal(size=(50, 128))
    rep = parity_report(a, a.copy())
    assert rep["passed"] is True
    assert rep["max_abs_deviation"] == 0.0
    assert rep["n_rows"] == 50


def test_parity_report_fails_on_a_perturbed_block_so_the_gate_is_not_vacuous():
    rng = np.random.default_rng(2)
    a = rng.normal(size=(50, 128))
    b = a.copy()
    b[7, 11] += 1e-2
    rep = parity_report(a, b, atol=PARITY_ATOL)
    assert rep["passed"] is False
    assert rep["max_abs_deviation"] == pytest.approx(1e-2)


def test_parity_report_tolerance_is_a_boundary_not_a_suggestion():
    a = np.zeros((4, 4))
    b = np.zeros((4, 4))
    b[0, 0] = 1e-4
    assert parity_report(a, b, atol=1e-4)["passed"] is True
    b[0, 0] = 1.01e-4
    assert parity_report(a, b, atol=1e-4)["passed"] is False


def test_parity_report_refuses_mismatched_shapes_rather_than_broadcasting():
    with pytest.raises(ProbeFailure, match="shapes differ"):
        parity_report(np.zeros((10, 128)), np.zeros((9, 128)))
    with pytest.raises(ProbeFailure, match="shapes differ"):
        parity_report(np.zeros((10, 128)), np.zeros((10, 1)))


@pytest.mark.parametrize("bad", [np.nan, np.inf])
def test_parity_report_refuses_non_finite_inputs(bad):
    a = np.zeros((5, 5))
    b = np.zeros((5, 5))
    b[2, 2] = bad
    with pytest.raises(ProbeFailure, match="non-finite"):
        parity_report(a, b)


# --------------------------------------------------------------------------------------------------
# exit codes -- a broken probe must NOT read as a clean feasibility
# --------------------------------------------------------------------------------------------------
def test_main_exits_non_zero_when_a_parity_deviation_exceeds_tolerance(monkeypatch, capsys):
    """The composition that gives the end-to-end property: `run` raises on a failed parity report
    (asserted above), and `main` maps any ProbeFailure to a non-zero exit."""
    import scripts.positional_arm_probe as probe

    def boom(**_):
        rep = parity_report(np.zeros((3, 3)), np.full((3, 3), 1e-2))
        assert rep["passed"] is False
        raise ProbeFailure("parity deviation %.3e exceeds atol" % rep["max_abs_deviation"])

    monkeypatch.setattr(probe, "run", boom)
    assert probe.main([]) == 1
    assert "PROBE FAILED" in capsys.readouterr().out


def test_main_exits_non_zero_on_any_failure_not_just_probe_failure(monkeypatch):
    import scripts.positional_arm_probe as probe

    monkeypatch.setattr(probe, "run", lambda **_: (_ for _ in ()).throw(RuntimeError("boom")))
    assert probe.main([]) == 1


def test_main_does_not_write_the_artifact_on_failure(monkeypatch, tmp_path):
    """A feasibility artifact left behind by a failed probe would later read as a clean run."""
    import scripts.positional_arm_probe as probe

    target = tmp_path / "should_not_exist.json"
    monkeypatch.setattr(probe, "OUT", target)
    monkeypatch.setattr(probe, "run", lambda **_: (_ for _ in ()).throw(ProbeFailure("no")))
    assert probe.main([]) == 1
    assert not target.exists()


# --------------------------------------------------------------------------------------------------
# the committed artifact -- skipped cleanly where the probe has not run
# --------------------------------------------------------------------------------------------------
@pytest.mark.skipif(not ARTIFACT.is_file(), reason="feasibility artifact absent; run the probe")
def test_the_committed_feasibility_artifact_is_internally_consistent():
    rec = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    assert rec["schema"] == "positional-arm-feasibility-v1"
    assert rec["verdict"] == "FEASIBLE"
    assert rec["parity"]["passed"] is True
    assert rec["model"]["n_params"] == rec["model"]["expected_params"] == 1_847_428
    assert rec["substrate"]["n_cohort"] == EXPECTED_TILES
    assert rec["forward_contract"]["latent_shape"][1:] == [EXPECTED_CHANNELS, 150]

    # the recorded rung bytes must be re-derivable from the recorded feature counts
    for rr in rec["resources"]["rungs"]:
        assert rr["peak_design_bytes"] == design_matrix_bytes(
            rec["substrate"]["n_cohort"], rr["n_features"])

    # the chosen lag set must be one of the declared candidates, and the largest that fit
    chosen = rec["lag_set"]["chosen"]
    assert chosen in LAG_CANDIDATES
    assert rec["lag_set"]["lags"] == LAG_CANDIDATES[chosen]
    fitting = [d for d in rec["lag_set"]["candidates_considered"] if d["fits_in_ram"] is True]
    assert rec["lag_set"]["n_lags"] == max(d["n_lags"] for d in fitting)

    # the stored block must be the sum of the three stored feature blocks
    counts = rung_feature_counts(rec["lag_set"]["n_lags"])
    assert rec["stored_block"]["n_features_total"] == (
        counts["bins10_mean_oriented"] + counts["bins10_max_oriented"] + counts["autocorr"])
    assert rec["stored_block"]["bytes"] == stored_block_bytes(
        rec["substrate"]["n_cohort"], rec["stored_block"]["n_features_total"])


@pytest.mark.skipif(not ARTIFACT.is_file(), reason="feasibility artifact absent; run the probe")
def test_the_artifact_records_the_substrate_facts_the_plan_asserted():
    rec = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    assert rec["substrate"]["strand_counts"] == {"+": 21821, "-": 22285}
    assert rec["substrate"]["n_peaks"] == 2433
    assert rec["substrate"]["n_duplicate_sequences"] == 689, (
        "the plan names 689 known duplicate sequences; a different count means a different cohort")


@pytest.mark.skipif(not ARTIFACT.is_file(), reason="feasibility artifact absent; run the probe")
def test_the_artifact_does_not_conflate_the_checkpoint_file_with_the_clone():
    """A 686 MB 'checkpoint' is the clone. Pinned because the conflation was live in the notes."""
    ck = json.loads(ARTIFACT.read_text(encoding="utf-8"))["checkpoint"]
    assert ck["bytes"] < ck["clone_bytes"]
    assert 15 < ck["mb"] < 40, "the base checkpoint file is ~21 MB"
    assert ck["clone_mb"] > 500, "the clone is ~686 MB"
    assert len(ck["sha256"]) == 64


@pytest.mark.skipif(not CKPT.is_file(), reason="PromoterAtlas clone absent (volatile temp dir)")
@pytest.mark.skipif(not ARTIFACT.is_file(), reason="feasibility artifact absent; run the probe")
def test_the_recorded_checkpoint_hash_still_matches_the_file_on_disk():
    """The clone lives in a volatile temp directory; this is what makes a swap detectable."""
    from scripts.positional_arm_probe import sha256_of

    rec = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    assert sha256_of(CKPT) == rec["checkpoint"]["sha256"]
    assert CKPT.stat().st_size == rec["checkpoint"]["bytes"]
