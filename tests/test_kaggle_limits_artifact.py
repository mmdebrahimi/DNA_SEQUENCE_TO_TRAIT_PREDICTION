"""Pins the Kaggle accelerator-limits artifact that Step 4 of the bounded-chain plan produced.

WHY A TEST FOR A FACT-RECORDING FILE. The artifact exists because the plan's highest-severity risk
was a scheduling question nobody had answered: whether the gradeable prefix fits inside one week's
GPU allowance. Its figures then decide how the chain is run. A recorded figure that nothing checks
drifts from its own arithmetic, and this file caught exactly that: an earlier draft paired the
401-epoch full-run count with the 380-epoch cost, because the smoke's `projected_hours_at_horizon`
is training-only and excludes the required nulls. The plan body makes the same pairing.

So the reconciliation below is the point of this file. It does NOT re-measure anything and it cannot
tell you whether the weekly quota is still 30 h -- re-run `kaggle quota` for that. What it can do is
refuse an artifact whose own numbers disagree with each other or with the committed bar.

A DELIBERATE DIVERGENCE FROM THE PLAN'S TEST STRATEGY, recorded rather than quietly dropped. Step 4
said "every number carries a source URL and a read date". The weekly quota has no URL: it came from
`kaggle quota`, the authenticated first-party account endpoint, after six separate web routes to
kaggle.com returned title-only bodies. Requiring a URL would have forced either a fabricated locator
or a worse source, so the assertion here is `source` + `read_date` + a verbatim quote of the command
output. That is stronger provenance than a doc URL, not weaker.

Pure structural + arithmetic checks against committed files: no network, no torch, no Kaggle.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "wiki" / "kaggle_session_limits_2026-11-15.json"
BAR = ROOT / "wiki" / "finetune_arm_acceptance_bar.json"

#: The runner's own null-replicate count. Read from source rather than restated, because the whole
#: null-phase cost derives from it and a drift here would silently move every hour figure below.
RUNNER = ROOT / "scripts" / "finetune_arm_run.py"


@pytest.fixture(scope="module")
def art() -> dict:
    return json.loads(ARTIFACT.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def bar() -> dict:
    return json.loads(BAR.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def n_null_seeds() -> int:
    for line in RUNNER.read_text(encoding="utf-8").splitlines():
        if line.startswith("N_NULL_SEEDS"):
            return int(line.split("=")[1].split("#")[0].strip())
    raise AssertionError("N_NULL_SEEDS not found in the runner")


def test_the_artifact_exists_and_declares_its_schema(art):
    assert art["schema"] == "kaggle-accelerator-limits-v1"
    assert art["measured_at"] == "2026-11-15"


# ---------------------------------------------------------------------------------------------------
# provenance
# ---------------------------------------------------------------------------------------------------
def test_the_resolved_quota_carries_its_source_and_read_date(art):
    q = art["weekly_gpu_quota"]
    assert q["status"] == "RESOLVED"
    assert q["source"] and q["read_date"], q
    assert q["verbatim_output"].strip(), "a resolved figure must ship the output it was read from"
    # The quoted output has to actually contain the number the artifact claims.
    assert "%.2f" % q["total_hours"] in q["verbatim_output"]
    assert "%.2f" % q["remaining_hours_at_read"] in q["verbatim_output"]


def test_the_quota_numbers_are_internally_consistent(art):
    q = art["weekly_gpu_quota"]
    assert abs(q["used_hours_at_read"] + q["remaining_hours_at_read"] - q["total_hours"]) < 0.011, q


def test_an_unresolved_status_carries_what_was_tried(art):
    cap = art["gpu_session_cap"]
    assert cap["status"] == "UNRESOLVED"
    tried = cap["what_was_tried"]
    assert isinstance(tried, list) and len(tried) >= 3, tried
    assert all(isinstance(t, str) and len(t) > 40 for t in tried), (
        "each logged route must say what was tried AND how it failed, not just name a URL")
    # An unresolved figure must say why it does not block, or it reads as an open blocker.
    assert cap["why_it_does_not_block"]
    assert cap["how_it_resolves_empirically"]


def test_no_figure_is_asserted_for_the_unresolved_cap(art):
    """The failure mode this guards is adopting one of the disagreeing community numbers.

    9 h, 12 h and 24 h all circulate, and the repo's own pusher comment asserts 12 h. An UNRESOLVED
    status that nonetheless carried a number would be a choice dressed as a measurement.
    """
    cap = art["gpu_session_cap"]
    numeric = {k: v for k, v in cap.items() if isinstance(v, (int, float))}
    assert not numeric, "the unresolved cap carries numeric fields: %r" % numeric


# ---------------------------------------------------------------------------------------------------
# the arithmetic -- the check that caught a real error
# ---------------------------------------------------------------------------------------------------
def test_every_hour_figure_reconciles_against_the_measured_rate(art, bar, n_null_seeds):
    """Derive each cost from spe and the COMMITTED bar; the artifact may only agree.

    This is what found the defect: `projected_hours_at_horizon` 134.69 h in the smoke artifact is
    spe * n_seeds * diagnostic_horizon_epochs / 3600 -- 380 epochs, TRAINING ONLY, with the 21
    required null epochs excluded. An earlier draft paired that cost with the 401-epoch count.
    """
    c = art["cost_arithmetic_at_the_measured_rate"]
    spe = c["seconds_per_epoch_measured"]
    h = bar["horizon"]
    n_seeds = int(h["n_seeds"])
    frozen = int(bar["frozen_hyperparameters"]["n_epochs_budget"])
    horizon = int(h["diagnostic_horizon_epochs"])

    def hours(epochs: int) -> float:
        return epochs * spe / 3600.0

    budget_ep = n_seeds * frozen
    null_ep = n_null_seeds * frozen
    prefix_ep = budget_ep + null_ep
    train_ep = n_seeds * horizon
    full_ep = train_ep + null_ep

    assert c["budget_phase_epochs"] == budget_ep
    assert c["null_phase_epochs"] == null_ep
    assert c["gradeable_prefix_epochs"] == prefix_ep
    assert c["training_only_epochs"] == train_ep
    assert c["full_run_including_nulls_epochs"] == full_ep

    for key, ep in (("budget_phase_hours", budget_ep),
                    ("null_phase_hours", null_ep),
                    ("gradeable_prefix_hours", prefix_ep),
                    ("training_only_hours", train_ep),
                    ("full_run_including_nulls_hours", full_ep)):
        assert abs(c[key] - hours(ep)) < 0.02, (
            "%s is %r but %d epochs at %.4f s cost %.4f h" % (key, c[key], ep, spe, hours(ep)))

    # And the two totals must NOT be equal -- that equality was the bug.
    assert c["training_only_hours"] != c["full_run_including_nulls_hours"], (
        "the training-only and full-run costs are the same number again; the nulls are missing from "
        "one of them")


def test_the_reconciliation_is_NON_VACUOUS(art, bar, n_null_seeds):
    """Reintroduce the original defect in a throwaway copy; the check must reject it."""
    import copy
    broken = copy.deepcopy(art)
    c = broken["cost_arithmetic_at_the_measured_rate"]
    c["full_run_including_nulls_hours"] = c["training_only_hours"]  # the exact original error

    spe = c["seconds_per_epoch_measured"]
    full_ep = c["full_run_including_nulls_epochs"]
    assert abs(c["full_run_including_nulls_hours"] - full_ep * spe / 3600.0) >= 0.02, (
        "pairing the full-run epoch count with the training-only cost no longer registers as an "
        "error, so the reconciliation above cannot be catching anything")


def test_the_epoch_counts_come_from_the_committed_bar_not_from_prose(bar, n_null_seeds):
    """Anchors the arithmetic's inputs, so a bar edit cannot leave the artifact quietly stale."""
    assert int(bar["horizon"]["n_seeds"]) == 10
    assert int(bar["frozen_hyperparameters"]["n_epochs_budget"]) == 7
    assert int(bar["horizon"]["diagnostic_horizon_epochs"]) == 38
    assert n_null_seeds == 3
    # `run()` reads the frozen budget from frozen_hyperparameters, NOT from horizon. The bar happens
    # to agree today; if it ever stops agreeing, the arithmetic must follow the field run() reads.
    assert int(bar["horizon"]["frozen_budget_epochs"]) == int(
        bar["frozen_hyperparameters"]["n_epochs_budget"]), (
        "the bar's two statements of the frozen budget have diverged; run() reads "
        "frozen_hyperparameters.n_epochs_budget and the cost arithmetic must use that one")


# ---------------------------------------------------------------------------------------------------
# the verdict must follow from the artifact's own figures
# ---------------------------------------------------------------------------------------------------
def test_the_binding_constraint_verdict_follows_from_the_numbers(art):
    """A prose verdict that does not follow from its own artifact is the thing to catch."""
    c = art["cost_arithmetic_at_the_measured_rate"]
    q = art["weekly_gpu_quota"]
    f = art["the_outcome_determinative_finding"]
    # The claim is that the prefix does not fit in the remaining quota.
    assert c["gradeable_prefix_hours"] > q["remaining_hours_at_read"], (
        "the artifact says the prefix does not fit the quota, but its own figures say it does")
    assert "WEEKLY QUOTA IS THE BINDING CONSTRAINT" in f["verdict"]
    assert "DOES NOT FIT IN ONE WEEK" in f["verdict"]
    # And it must record that this costs no frozen pre-registration, which is why it needs no
    # authority call -- the distinction the plan's Open Question 1 collapsed.
    assert f["crucially_it_changes_no_frozen_pre_registration"]


def test_the_reachable_epochs_this_week_follow_from_the_remaining_quota(art):
    c = art["cost_arithmetic_at_the_measured_rate"]
    q = art["weekly_gpu_quota"]
    spe = c["seconds_per_epoch_measured"]
    derived = int(q["remaining_hours_at_read"] * 3600.0 // spe)
    assert abs(c["epochs_reachable_with_remaining_quota_this_week"] - derived) <= 1, (
        "claimed %r reachable epochs, but %.2f h remaining at %.2f s/epoch gives %d"
        % (c["epochs_reachable_with_remaining_quota_this_week"],
           q["remaining_hours_at_read"], spe, derived))


def test_the_artifact_states_its_own_limits(art):
    lim = art["honest_limits"]
    assert isinstance(lim, list) and len(lim) >= 3
    joined = " ".join(lim).lower()
    assert "floating" in joined, "the quota is a point reading of a figure Kaggle floats; say so"
    assert "six routes" in joined or "routes" in joined


def test_the_per_session_overhead_finding_is_recorded_with_its_decision(art):
    """The measurement epoch costs quota on every chained session and was left in place.

    Recorded because it is a real cost against the binding constraint, and because leaving it was a
    judgement -- not an oversight -- that a later reader is entitled to re-examine.
    """
    o = art["per_session_overhead_finding"]
    assert o["decision"].startswith("LEFT IN PLACE")
    assert o["why_not_fixed_now"]
    assert o["what_a_fix_would_look_like_if_it_is_ever_worth_it"]
    assert o["what_was_considered_and_rejected"], (
        "the rejected option -- promoting the smoke epoch to seed 0's epoch 1 -- is the tempting "
        "one, so it has to be on the record as rejected rather than unconsidered")


# ---------------------------------------------------------------------------------------------------
# the supersession, added after session 1 measured the real steady-state rate
# ---------------------------------------------------------------------------------------------------
def test_the_stale_verdict_cannot_be_quoted_WITHOUT_its_supersession(art):
    """The failure this guards: a reader lifts the headline and misses that it was overturned.

    The original verdict is preserved verbatim on purpose -- it was the correct reading of the only
    rate measured at the time -- so the protection has to be that it carries a pointer, not that it
    was edited away.
    """
    f = art["the_outcome_determinative_finding"]
    assert "SUPERSEDED_BY" in f, (
        "the superseded verdict has lost its pointer; anyone quoting it would state that the prefix "
        "does not fit, which measurement has since inverted")
    assert "SUPERSEDED_LATER_THE_SAME_DAY" in f["SUPERSEDED_BY"]
    assert "never alone" in f["SUPERSEDED_BY"]
    # And the original wording must still be there, unedited.
    assert "DOES NOT FIT IN ONE WEEK" in f["verdict"]


def test_the_supersession_states_the_measurement_and_the_inversion(art):
    s = art["SUPERSEDED_LATER_THE_SAME_DAY"]
    # The measured numbers, not a summary of them.
    for n in ("1276.04", "1137.0", "938.7", "886.7", "945.5", "923.6"):
        assert n in s["what_changed"], "the supersession omits the measured value %s" % n
    assert "INVERT" in s["the_conclusion_INVERTS"].upper() or "FITS" in s["the_conclusion_INVERTS"]
    assert "24.18" in s["the_conclusion_INVERTS"] and "28.14" in s["the_conclusion_INVERTS"]
    # It must say what did NOT change, or the reader over-reads the inversion.
    assert "still the binding constraint" in s["what_does_NOT_change"]
    assert "UNRESOLVED" in s["what_does_NOT_change"]


def test_the_supersession_does_not_claim_the_original_was_WRONG(art):
    """A correction that mislabels its predecessor teaches the wrong lesson.

    The original arithmetic was right at the rate then known. What changed is the rate, which makes
    the multi-week schedule unnecessary rather than mistaken -- and that distinction is what keeps
    the multi-week option available if the margin disappears.
    """
    s = art["SUPERSEDED_LATER_THE_SAME_DAY"]
    txt = s["what_that_does_to_the_finding_below"]
    assert "does not make it wrong" in txt
    assert "UNNECESSARY rather than mistaken" in txt
    assert "remains available" in txt


def test_the_supersession_bounds_its_own_evidence(art):
    """Three epochs is thin, and the inversion is load-bearing on them."""
    s = art["SUPERSEDED_LATER_THE_SAME_DAY"]
    trust = s["how_much_to_trust_it"]
    assert "Three steady-state epochs" in trust
    assert "would NOT fit" in trust, (
        "the supersession must state what would happen if the small sample is unrepresentative")
    assert "self-correcting" in trust



# --------------------------------------------------------------------------------------------------
# Session two re-measured the same rate and the first supersession was itself ~7% low. The protection
# is the same shape as the one above: the stale layer keeps its wording and gains a pointer, and the
# live figure lives in its own layer. These tests exist because the FIRST layer's own text told a
# reader to re-read the logs rather than quote 923.6 s, and nothing enforced that.
# --------------------------------------------------------------------------------------------------

LAYER_2 = "SUPERSEDED_AGAIN_BY_SESSION_TWO_2026_11_16"


def test_the_first_supersession_now_carries_a_pointer_to_the_second(art):
    """Without this, a reader quotes 923.6 s and a 4 h margin as the live schedule.

    The first layer is still correct AT the rate it measured, which is why it is preserved unedited
    rather than rewritten -- so the only available protection is that it points forward.
    """
    s = art["SUPERSEDED_LATER_THE_SAME_DAY"]
    assert "SUPERSEDED_BY" in s, (
        "the first supersession has no forward pointer; anyone quoting its 923.6 s rate or its "
        "'roughly 4 h of margin' would be stating a superseded schedule as live")
    ptr = s["SUPERSEDED_BY"]
    assert LAYER_2 in ptr
    assert "never quote" in ptr
    assert "923.6" in ptr and "990.0" in ptr, (
        "the pointer must name BOTH rates, or a reader cannot tell which figure it is retiring")
    # The superseded layer's own wording must survive the pointer being added.
    assert "923.6" in s["what_changed"]
    assert "24.18" in s["the_conclusion_INVERTS"] and "28.14" in s["the_conclusion_INVERTS"]


def test_the_second_supersession_states_the_measurement_and_its_basis(art):
    s = art[LAYER_2]
    for n in ("923.6", "990.0", "955.0", "1029.4"):
        assert n in s["what_changed"], "the second supersession omits the measured value %s" % n
    # The BASIS is what makes 990 preferable to 923.6; a bare pair of numbers does not say which wins.
    assert "25 steady-state epochs" in s["what_changed"]
    assert "THREE epochs of ONE seed" in s["what_changed"]


def test_the_second_supersession_narrows_rather_than_claiming_another_inversion(art):
    """An inversion and a narrowing teach different things, and only one of them happened.

    Session one INVERTED the verdict (does not fit -> fits). Session two did not invert anything; it
    moved the margin. Reporting it as a second inversion would imply the prefix stopped fitting.
    """
    s = art[LAYER_2]
    txt = s["the_conclusion_NARROWS_rather_than_inverting_again"]
    assert "STILL FITS" in txt
    assert "2.4 h" in txt and "4 h" in txt, "the narrowing must name both margins"
    # Two independent derivations, because one of them assumes an epoch model and the other does not.
    assert "0.2955" in txt and "17.14" in txt and "19.57" in txt


def test_the_second_supersession_does_not_relabel_the_first_as_a_shape_error(art):
    """The single-epoch smoke was a measurement-SHAPE error. This one is ordinary sampling.

    Collapsing the two would discard the first diagnosis, which is the reusable one: the smoke's
    figure was wrong because of how it was taken, not because of how many samples it had.
    """
    s = art[LAYER_2]
    txt = s["why_the_session_one_figure_was_low"]
    assert "ordinary sampling" in txt
    assert "NOT a new measurement-shape error" in txt
    assert "that diagnosis stands" in txt


def test_the_second_supersession_names_the_first_cut_if_the_margin_binds(art):
    """A thinning margin is only actionable if the artifact says what to spend it on."""
    s = art[LAYER_2]
    txt = s["the_first_thing_to_cut_if_the_margin_becomes_binding"]
    assert "rate-measurement epoch" in txt
    assert "changes no behaviour" in txt
    assert "that reason has expired" in txt, (
        "the cut was deferred for a stated reason; the artifact must record that the reason lapsed")


def test_the_session_cap_floor_moved_but_is_still_unresolved(art):
    """Session two ran 8.57 h un-killed, which is evidence about the cap without resolving it."""
    s = art[LAYER_2]
    txt = s["what_does_NOT_change"]
    assert "STILL UNRESOLVED" in txt
    assert "8.57" in txt and "1.45" in txt, (
        "the empirical floor on the session cap moved from session one's 1.45 h to 8.57 h; both "
        "numbers belong here or the reader cannot see that the floor rose")


def test_the_second_supersession_bounds_its_own_evidence(art):
    s = art[LAYER_2]
    trust = s["how_much_to_trust_990"]
    assert "ONE session on ONE T4" in trust
    assert "Session three re-measures" in trust
    assert "live risk" in trust


def test_three_disagreeing_measurements_of_one_quantity_are_all_recorded(art):
    """The rate has been measured three times and no two agree; a reader needs all three.

    Quoting any one of them as 'the' rate is the error this artifact keeps having to correct.
    """
    txt = art[LAYER_2]["why_the_rate_keeps_needing_re_measurement"]
    for n in ("1276.7", "923.6", "990.0", "1323.3"):
        assert n in txt, "the three-measurement record omits %s" % n
    assert "no two agree" in txt
    assert "recurs every link" in txt


# ---------------------------------------------------------------------------
# Layer THREE -- session three, 2026-11-16.
#
# The artifact's house pattern is to supersede in LAYERS with forward pointers and
# never by editing, because each layer's arithmetic is correct at the rate IT
# measured. These pins exist so a layer cannot go stale SILENTLY: the hazard is not
# a wrong number, it is a reader quoting a dead margin as live.
# ---------------------------------------------------------------------------

LAYER_3 = "SUPERSEDED_AGAIN_BY_SESSION_THREE_2026_11_16"


def test_layer_two_carries_a_forward_pointer_to_layer_three(art):
    """A superseded layer that does not say so is the whole failure mode this guards."""
    assert LAYER_3 in art, "the third supersession layer is missing"
    ptr = art[LAYER_2].get("SUPERSEDED_BY")
    assert ptr, "layer two states a ~2.4 h margin and must point at what superseded it"
    assert LAYER_3 in ptr, "the pointer must NAME the superseding layer, not merely exist"
    assert "990.0" in ptr and "960.0" in ptr, (
        "the pointer must name BOTH rates so a reader can see which direction it moved")
    assert "preserved unedited" in ptr
    assert "NOT the live margin" in ptr, (
        "the stale margin is the specific thing a reader would otherwise quote")


def test_layer_three_records_the_third_measurement_and_the_pooled_estimate(art):
    s = art[LAYER_3]
    txt = s["what_changed"]
    for n in ("960.0", "990.0", "923.6", "974.7"):
        assert n in txt, "layer three omits %s" % n
    assert "26" in txt and "51" in txt, (
        "the epoch counts are what make 960.0 weightier than 923.6; both belong here")


def test_layer_three_refuses_the_marching_corrections_reading(art):
    """Three measurements that do not agree are a BAND, not a trend.

    Reading 923.6 -> 990.0 -> 960.0 as successive corrections would invite a fourth
    session to 'continue' it. The artifact must say the quantity has a band.
    """
    s = art[LAYER_3]
    txt = s["the_honest_reading_of_three_measurements"]
    assert "NOT a sequence of corrections" in txt
    assert "BAND" in txt
    assert "974.7" in txt
    assert "925.2-988.9" in txt, (
        "the within-session spread is the evidence that this is variance, not a trend")
    assert "not to continue a trend" in txt


def test_layer_three_preserves_the_measurement_shape_diagnosis(art):
    """Layer one's diagnosis is a DIFFERENT error class and must not be relabelled.

    Collapsing it into 'small sample' would retire a diagnosis that still applies --
    the in-session rate epoch over-reads its own steady state every single link.
    """
    s = art[LAYER_3]
    kept = " ".join(s["what_does_NOT_change"])
    assert "measurement-SHAPE diagnosis from layer one STANDS" in kept
    assert "different kind of error" in kept
    assert "1299.1" in kept, "the recurrence needs this session's own measurement"
    assert "30.00" in kept, "the weekly quota is a layer-one fact that still holds"


def test_layer_three_loosens_rather_than_narrows_the_margin(art):
    s = art[LAYER_3]
    txt = s["the_conclusion_does_not_invert_or_narrow_it_LOOSENS_SLIGHTLY"]
    assert "OVER-shot" in txt, (
        "layer two's narrowing was slightly too aggressive and the artifact must own that")
    assert "2.87" in txt and "11.00" in txt and "28" in txt
    assert "ONE more session" in txt


def test_layer_three_records_a_second_identical_session_cap_observation(art):
    """8.57 h twice is evidence about the cap; it is still not a resolution."""
    s = art[LAYER_3]
    txt = s["the_session_cap_gains_a_SECOND_identical_observation"]
    assert "8.57" in txt
    assert "UNRESOLVED" in txt, "two observations of a floor do not resolve a cap"
    assert "TWICE" in txt
    assert "9.5" in txt, "session four's ceiling is the first probe above the floor"


def test_layer_three_confirms_the_overshoot_bound_independently(art):
    s = art[LAYER_3]
    txt = s["the_overshoot_bound_is_CONFIRMED_on_an_independent_session"]
    assert "2028" in txt and "711" in txt, (
        "both sessions' overshoots belong here: one crossed a phase boundary and one did not")
    assert "phase boundary" in txt.lower()
    assert "LESS than one steady epoch" in txt


def test_layer_three_keeps_the_rate_epoch_and_gives_the_new_reason(art):
    """The deferral reason lapsed in layer two; retaining it now needs a FRESH reason.

    Carrying a lapsed justification forward unchanged is how a stale deferral survives.
    """
    s = art[LAYER_3]
    txt = s["the_first_thing_to_cut_was_NOT_cut_and_the_reason_is_now_different"]
    assert "1299.1" in txt
    assert "NEW reason, not the lapsed one" in txt
    assert "caught layer two" in txt, (
        "the epoch that looks redundant is what caught the over-correction; that is the reason")


def test_layer_three_bounds_its_own_evidence_including_the_unmeasured_null(art):
    """The pooled rate is the best estimate yet AND governs only a quarter of what is left."""
    s = art[LAYER_3]
    trust = s["how_much_to_trust_974_7"]
    assert "51" in trust and "9 independent seeds" in trust
    assert "NEVER been measured" in trust, (
        "a null epoch's cost governs 75% of the remaining work and has never been timed; an "
        "artifact that bounds only the measured rate would hide the dominant uncertainty")
    assert "75%" in trust


def test_layer_three_points_at_its_own_detail_artifact(art):
    s = art[LAYER_3]
    detail = ROOT / s["detail"]
    assert detail.exists(), "layer three cites %s, which does not exist" % s["detail"]
    body = json.loads(detail.read_text(encoding="utf-8"))
    assert body["session"] == 3
    assert body["reconcile_verdict"] == "LINK_OK"
