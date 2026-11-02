"""Tests for the expression-heads feasibility probe.

The probe's job is to make four things measurable before the arm runs: the geometry lock, the base
byte-identity, the cost of each windowing route, and the noise floor that separates a real edge
effect from conv-blocking arithmetic. These tests pin the refusals and the artifact's completeness;
the head arithmetic itself is pinned in `test_atlas_heads.py` and is not re-tested here.

The one claim worth stating plainly: **this probe must not choose the windowing route.** That is an
acceptance-bar decision in Step 4, and a probe that quietly picked one would be freezing a bar
before the bar exists.
"""

from __future__ import annotations

import json
import pathlib

import numpy as np
import pytest

torch = pytest.importorskip("torch", reason="the probe's route measurement needs torch")

from dna_decode.glm import atlas_heads as AH
from dna_decode.glm.atlas_embed import DEFAULT_PA_REPO
from scripts import expression_heads_probe as P

ARTIFACT = pathlib.Path(__file__).resolve().parents[1] / "wiki/expression_heads_feasibility_2026-11-02.json"
BASE_CK = pathlib.Path(DEFAULT_PA_REPO) / "trained_weights/base_model/promoteratlas-base.pt"
EXPR_DIR = pathlib.Path(DEFAULT_PA_REPO) / "trained_weights/expression_prediction"
_CLONE = BASE_CK.is_file() and EXPR_DIR.is_dir()
needs_clone = pytest.mark.skipif(not _CLONE, reason="PromoterAtlas clone absent (volatile temp dir)")


# --------------------------------------------------------------------------------------------------
# pure arithmetic -- known answers, no clone, no model
# --------------------------------------------------------------------------------------------------
def test_route_cost_projection_scales_linearly_and_then_by_the_head_count():
    """10 of 100 tiles in 2 s -> 20 s per head -> 100 s over five heads."""
    p = P.route_cost_projection(100, 10, 2.0, 5)
    assert p["seconds_per_head"] == pytest.approx(20.0)
    assert p["n_heads"] == 5
    assert p["seconds_all_heads"] == pytest.approx(100.0)
    # hours is a 2-decimal DISPLAY field; seconds is the figure to reason from
    assert p["hours_all_heads"] == round(100.0 / 3600.0, 2)


def test_route_cost_projection_refuses_a_nonsense_head_count():
    with pytest.raises(ValueError, match="n_heads must be"):
        P.route_cost_projection(100, 10, 2.0, 0)


def test_the_head_count_multiplies_the_cost_because_each_head_runs_its_own_trunk():
    one = P.route_cost_projection(1000, 10, 1.0, 1)["seconds_all_heads"]
    five = P.route_cost_projection(1000, 10, 1.0, 5)["seconds_all_heads"]
    assert five == pytest.approx(5 * one)


def test_uncontaminated_deviation_reads_ONLY_the_columns_the_closed_form_calls_identical():
    """Planted: a big difference inside a contaminated column, a small one outside it."""
    w, c = AH.HEAD_SEQ_LENGTH, 3
    cheap = np.zeros((c, w))
    true_ = np.zeros((c, w))
    bad = sorted(AH.contaminated_columns(30, 150, w))
    good = [i for i in range(w) if i not in set(bad)]

    cheap[:, bad[0]] = 99.0          # inside the edge effect -- must be IGNORED
    cheap[:, good[0]] = 1e-7         # outside it -- this is the floor being measured

    out = P.uncontaminated_deviation(cheap, true_, 30, 150, w)
    assert out["n_columns"] == len(good)
    assert out["max_abs_dev"] == pytest.approx(1e-7)
    assert out["max_abs_dev"] < 1.0, "a contaminated column leaked into the noise floor"


def test_uncontaminated_deviation_declines_to_report_a_floor_when_every_column_is_contaminated():
    """At a width small enough that the two edges cover everything there is no floor to read."""
    w = 4
    out = P.uncontaminated_deviation(np.zeros((2, w)), np.zeros((2, w)), 1, 10, w)
    assert out["n_columns"] == 0
    assert out["max_abs_dev"] is None
    assert "no floor" in out["why"]


def test_spearman_between_is_rank_based_and_so_is_blind_to_a_monotone_rescale():
    """The property that makes the heads' denormalisation irrelevant to a Spearman verdict."""
    a = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    assert P.spearman_between(a, 3.7 * a + 11.0) == pytest.approx(1.0)
    assert P.spearman_between(a, np.exp(a)) == pytest.approx(1.0)
    assert P.spearman_between(a, -a) == pytest.approx(-1.0)


# --------------------------------------------------------------------------------------------------
# the refusals
# --------------------------------------------------------------------------------------------------
@needs_clone
def test_verify_checkpoint_REFUSES_a_moved_hash(monkeypatch):
    """The clone lives in a volatile temp dir, so a moved hash voids every inherited figure."""
    monkeypatch.setattr(P, "EXPECTED_SHA256", "0" * 64)
    with pytest.raises(P.ProbeRefused, match="checkpoint MOVED"):
        P.verify_checkpoint()


@needs_clone
def test_verify_checkpoint_REFUSES_a_moved_byte_count(monkeypatch):
    monkeypatch.setattr(P, "EXPECTED_BYTES", 1)
    with pytest.raises(P.ProbeRefused, match="checkpoint MOVED"):
        P.verify_checkpoint()


def test_verify_checkpoint_refuses_an_absent_clone(tmp_path):
    with pytest.raises(P.ProbeRefused, match="absent"):
        P.verify_checkpoint(tmp_path)


def test_main_returns_nonzero_when_the_clone_is_absent(tmp_path):
    rc = P.main(["--pa-repo", str(tmp_path), "--out", str(tmp_path / "x.json")])
    assert rc == 2, "a refusal must exit non-zero, never 0 with a thin artifact"
    assert not (tmp_path / "x.json").exists(), "a refused probe must not write an artifact"


@needs_clone
def test_main_returns_nonzero_on_a_moved_hash_and_writes_nothing(monkeypatch, tmp_path):
    monkeypatch.setattr(P, "EXPECTED_SHA256", "f" * 64)
    out = tmp_path / "x.json"
    assert P.main(["--out", str(out)]) == 2
    assert not out.exists()


def test_verify_heads_propagates_a_base_drift_refusal(tmp_path):
    """A head whose base weights differ makes this a FINE-TUNED-encoder arm -- a different arm."""
    base = {"a.weight": torch.zeros(2, 2)}
    sd = {
        "conv1.weight": torch.zeros(AH.LATENT_CHANNELS, AH.LATENT_CHANNELS, AH.CONV_KERNEL),
        "conv1.bias": torch.zeros(AH.LATENT_CHANNELS),
        "conv2.weight": torch.zeros(AH.LATENT_CHANNELS, AH.LATENT_CHANNELS, AH.CONV_KERNEL),
        "conv2.bias": torch.zeros(AH.LATENT_CHANNELS),
        "fc1.weight": torch.zeros(*AH.FC1_SHAPE), "fc1.bias": torch.zeros(AH.FC1_SHAPE[0]),
        "fc2.weight": torch.zeros(1, AH.FC1_SHAPE[0]), "fc2.bias": torch.zeros(1),
        "base_model.a.weight": torch.ones(2, 2),            # DRIFTED
    }
    d = tmp_path / "trained_weights" / "expression_prediction"
    d.mkdir(parents=True)
    torch.save({"model_state_dict": sd}, d / "promoteratlas-trspred-x2020.pt")
    with pytest.raises(AH.AtlasHeadsError, match="NOT byte-identical"):
        AH.load_heads(d, base_state_dict=base)


# --------------------------------------------------------------------------------------------------
# the verdict function -- closed set, and it must NOT pick a route
# --------------------------------------------------------------------------------------------------
def _route_row(*, const=False):
    return {"per_reduction": {"mean": {"cheap_is_constant": const, "true_is_constant": False}}}


def _head_row(*, identical=124, checked=124):
    return {"base_tensors_identical": identical, "base_tensors_checked": checked}


def _proj(cheap=0.24, true_=4.64):
    return {"cheap_route": {"hours_all_heads": cheap}, "true_route": {"hours_all_heads": true_}}


def test_a_collapsed_head_is_INDETERMINATE_not_a_clean_null():
    v, why = P.classify_feasibility({"h": _route_row(const=True)}, _proj(), {"h": _head_row()})
    assert v == "INDETERMINATE_COLLAPSED_HEAD"
    assert "clean null" in why


def test_no_measured_head_is_INDETERMINATE():
    v, _ = P.classify_feasibility({}, _proj(), {})
    assert v == "INDETERMINATE_NO_HEAD_MEASURED"


def test_base_drift_is_INDETERMINATE():
    v, _ = P.classify_feasibility({"h": _route_row()}, _proj(), {"h": _head_row(identical=123)})
    assert v == "INDETERMINATE_BASE_DRIFT"


def test_the_feasible_verdict_REFUSES_to_choose_a_route_and_says_so():
    v, why = P.classify_feasibility({"h": _route_row()}, _proj(), {"h": _head_row()})
    assert v == "FEASIBLE_ROUTE_CHOICE_IS_THE_BARS"
    assert "deliberately NOT made here" in why
    assert "4.64" in why and "0.24" in why, "both route costs must be in the verdict text"


def test_every_verdict_the_classifier_can_return_names_itself_as_either_feasible_or_indeterminate():
    seen = {
        P.classify_feasibility({}, _proj(), {})[0],
        P.classify_feasibility({"h": _route_row(const=True)}, _proj(), {"h": _head_row()})[0],
        P.classify_feasibility({"h": _route_row()}, _proj(), {"h": _head_row(identical=1)})[0],
        P.classify_feasibility({"h": _route_row()}, _proj(), {"h": _head_row()})[0],
    }
    assert len(seen) == 4, "the four branches must be distinguishable"
    assert all(v.startswith("INDETERMINATE_") or v.startswith("FEASIBLE_") for v in seen)


# --------------------------------------------------------------------------------------------------
# the committed artifact -- a head must not be able to go missing quietly
# --------------------------------------------------------------------------------------------------
@pytest.mark.skipif(not ARTIFACT.is_file(), reason="probe has not been run")
def test_the_artifact_names_all_FIVE_heads_each_with_its_own_deviation_figure():
    rec = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    assert rec["n_heads"] == 5
    assert set(rec["heads"]) == set(rec["routes"]), "every loaded head must be measured"
    assert len(rec["routes"]) == 5
    for name, row in rec["routes"].items():
        dev = row["per_reduction"]["mean"]["deviation"]
        assert dev["max_abs_dev"] is not None, "%s has no deviation figure" % name
        assert row["seconds_cheap_route"] > 0 and row["seconds_true_route"] > 0
        assert row["n_windows"] == 65


@pytest.mark.skipif(not ARTIFACT.is_file(), reason="probe has not been run")
def test_the_artifact_records_the_geometry_lock_and_its_arithmetic():
    g = json.loads(ARTIFACT.read_text(encoding="utf-8"))["geometry"]
    assert g["fc1_in_features"] == 11008
    assert g["head_seq_length"] == 86 and g["substrate_length"] == 150
    assert g["substrate_would_flatten_to"] == 19200
    assert g["direct_feed_is_impossible"] is True
    assert g["n_windows"] == 65
    assert g["edge_columns_left"] == 9 and g["edge_columns_right"] == 11


@pytest.mark.skipif(not ARTIFACT.is_file(), reason="probe has not been run")
def test_the_artifact_records_every_heads_base_identity_as_124_of_124_at_zero_deviation():
    heads = json.loads(ARTIFACT.read_text(encoding="utf-8"))["heads"]
    for name, row in heads.items():
        assert row["base_tensors_checked"] == row["base_tensors_identical"] == 124, name
        assert row["base_max_abs_dev"] == 0.0, name
        assert row["fc1_shape"] == [64, 11008] and row["implied_seq_length"] == 86


@pytest.mark.skipif(not ARTIFACT.is_file(), reason="probe has not been run")
def test_the_artifact_states_the_multiplicity_consequence_of_FIVE_heads():
    rec = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    m = rec["multiplicity_consequence"]
    assert "alpha/5" in m and "divisor" in m
    assert "HARDER" in m, "the honest direction of an inherited miscount must be stated"


@pytest.mark.skipif(not ARTIFACT.is_file(), reason="probe has not been run")
def test_the_artifact_scopes_the_circularity_to_the_rung_where_it_actually_applies():
    c = json.loads(ARTIFACT.read_text(encoding="utf-8"))["circularity_scope"]
    assert c["circular_head"] == "promoteratlas-trspred-urtecho2018"
    assert c["circular_on_rung"] == "GSE108535"
    assert c["this_substrate"] == "GSE144621"
    assert c["is_circular_here"] is False
    assert "the rung where it IS circular" in c["recorded_because"]


@pytest.mark.skipif(not ARTIFACT.is_file(), reason="probe has not been run")
def test_the_artifact_records_the_tslpred_endpoint_mismatch_as_a_LIMIT_not_an_exclusion():
    e = json.loads(ARTIFACT.read_text(encoding="utf-8"))["endpoint_scope"]
    assert e["translation_heads"] == ["promoteratlas-tslpred-kosuri2013"]
    assert len(e["transcription_heads"]) == 4
    assert "TRANSLATION" in e["limit"]
    assert "rather than resolved by exclusion" in e["limit"]


@pytest.mark.skipif(not ARTIFACT.is_file(), reason="probe has not been run")
def test_the_artifact_records_the_count_trap_beside_both_counts():
    ck = json.loads(ARTIFACT.read_text(encoding="utf-8"))["checkpoint"]
    assert ck["n_parameters_via_model_parameters"] == 1_847_428
    assert ck["n_elements_via_state_dict"] == 1_847_492
    assert ck["buffer_elements_excluded_by_parameters"] == 64
    assert "BUFFERS" in ck["count_trap"]
    assert ck["sha256"] == P.EXPECTED_SHA256 and ck["bytes"] == P.EXPECTED_BYTES


@pytest.mark.skipif(not ARTIFACT.is_file(), reason="probe has not been run")
def test_the_artifact_names_wall_clock_and_NOT_memory_as_the_binding_constraint():
    r = json.loads(ARTIFACT.read_text(encoding="utf-8"))["resources"]
    assert r["n_tiles"] == 44106 and r["n_windows"] == 65 and r["n_heads"] == 5
    assert r["ram_is_not_the_binding_constraint"] is True
    assert "WALL CLOCK" in r["binding_constraint"]
    assert "UNMEASURED" in r["fits_in_ram_semantics"]


@pytest.mark.skipif(not ARTIFACT.is_file(), reason="probe has not been run")
def test_the_artifact_frames_a_null_as_bounding_TRANSFER_rather_than_the_heads():
    z = json.loads(ARTIFACT.read_text(encoding="utf-8"))["zero_shot_framing"]
    assert "TRANSFER" in z
    assert "does NOT bound the heads" in z


@pytest.mark.skipif(not ARTIFACT.is_file(), reason="probe has not been run")
def test_the_artifact_records_both_route_costs_so_the_bar_can_choose_from_a_number():
    p = json.loads(ARTIFACT.read_text(encoding="utf-8"))["cost_projection"]
    cheap = p["cheap_route"]["hours_all_heads"]
    true_ = p["true_route"]["hours_all_heads"]
    assert cheap > 0 and true_ > cheap, "the true route must be measured as the costlier one"
    assert p["base_forward_pass"]["shared_across_heads"] is True
    assert "byte-identical" in p["base_forward_pass"]["why_shared"]


@pytest.mark.skipif(not ARTIFACT.is_file(), reason="probe has not been run")
def test_the_artifact_shows_the_cheap_routes_error_is_FAR_above_the_arithmetic_noise_floor():
    """The measurement that makes the route choice a real decision rather than a formality.

    If the edge effect sat at the float32 floor there would be nothing to trade off. It does not:
    the reduced-value deviation is orders of magnitude larger than the uncontaminated-column floor.
    """
    rec = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    for name, row in rec["routes"].items():
        floor_key = [k for k in row if k.startswith("noise_floor_at_offset_")][0]
        floor = row[floor_key]["max_abs_dev"]
        edge = row["per_reduction"]["mean"]["deviation"]["max_abs_dev"]
        assert floor is not None, name
        assert edge > 100 * floor, (
            "%s: edge effect %.3e is not clearly above the noise floor %.3e, so the route choice "
            "would be measuring arithmetic" % (name, edge, floor)
        )


@pytest.mark.skipif(not ARTIFACT.is_file(), reason="probe has not been run")
def test_the_artifact_shows_the_two_routes_agree_on_RANK_which_is_what_a_spearman_arm_scores():
    """Reported, not buried: the cheap route moves VALUES a lot and RANKS very little."""
    rec = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    rhos = [row["per_reduction"]["mean"]["spearman_cheap_vs_true"]
            for row in rec["routes"].values()]
    assert len(rhos) == 5
    assert min(rhos) > 0.95, "rank agreement %r is lower than the bar's premise assumes" % rhos


@pytest.mark.skipif(not ARTIFACT.is_file(), reason="probe has not been run")
def test_the_committed_artifact_is_feasible_and_did_not_choose_a_route():
    rec = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    assert rec["verdict"] == "FEASIBLE_ROUTE_CHOICE_IS_THE_BARS"
    assert rec["schema"] == "expression-heads-feasibility-v1"
    assert "chosen_route" not in rec and "route" not in rec, (
        "the probe must not record a route CHOICE; that is the bar's job in Step 4"
    )
