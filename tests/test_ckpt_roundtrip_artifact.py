"""Pins the checkpoint-dataset round-trip artifact, and specifically the distinction the plan named.

THE PLAN'S OWN TEST STRATEGY asked for one thing here: that the artifact distinguish
`mount_preserved_directory_name` from `guard_fired`, "so a later reader cannot collapse them". They
call for opposite responses -- a fired guard means the protection worked, a failed mount assumption
means code must change -- and in the event BOTH happened at once, which is exactly the reading a
collapsed record would lose.

The artifact then recorded a THIRD finding the plan did not anticipate (the mount path is not
`/kaggle/input/<slug>` at all), so these tests also check that the two falsifications stay separate
from each other rather than merging into a single "the mount did not work".

Pure structural checks against a committed JSON file: no network, no Kaggle, no torch.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "wiki" / "ckpt_dataset_roundtrip_2026-11-15.json"
RUNNER = ROOT / "scripts" / "finetune_arm_run.py"


@pytest.fixture(scope="module")
def art() -> dict:
    return json.loads(ARTIFACT.read_text(encoding="utf-8"))


def test_the_artifact_declares_its_schema_and_date(art):
    assert art["schema"] == "finetune-arm-ckpt-roundtrip-v1"
    assert art["measured_at"] == "2026-11-15"


def test_it_distinguishes_a_fired_guard_from_a_failed_mount_assumption(art):
    """The distinction the plan asked for, as two SEPARATE booleans rather than one verdict."""
    pair = art["the_plan_s_two_named_outcomes"]
    assert set(pair) >= {"guard_fired", "mount_preserved_directory_name"}
    assert isinstance(pair["guard_fired"], bool)
    assert isinstance(pair["mount_preserved_directory_name"], bool)
    # Both happened, which is the designed-for pair and the one a collapsed record would lose.
    assert pair["guard_fired"] is True
    assert pair["mount_preserved_directory_name"] is False
    assert pair["how_to_read_that_pair"], (
        "the pair needs a reading, or a later reader treats a fired guard as a failure")


def test_the_two_falsifications_stay_separate(art):
    """One assumption failed in two independent ways and they need different fixes.

    Collapsing them into "the mount did not work" would hide that the FILENAME rule and the SEARCH
    ROOT are different changes, either of which alone leaves the round trip broken.
    """
    f1 = art["finding_1_the_directory_level_does_not_survive_the_upload"]
    f2 = art["finding_2_the_mount_path_is_not_predictable"]
    assert f1["status"].startswith("FALSIFIED")
    assert f2["status"].startswith("FALSIFIED")
    assert f1["fix"] != f2["fix"], "the two findings were given the same fix; they are not the same"
    assert "FILENAME" in f1["fix"]
    assert "HINT" in f2["fix"] or "fall back" in f2["fix"]
    # Each must carry what was MEASURED, not just what was concluded.
    assert isinstance(f1["what_was_measured"], list) and len(f1["what_was_measured"]) >= 3
    assert f2["what_was_measured"]


def test_the_smoke_exclusion_is_recorded_as_load_bearing(art):
    """The new rule is stricter than the old one, and the reason has to be on the record.

    The old parent-name rule excluded smoke output by ACCIDENT (different directory). A filename rule
    must exclude it deliberately, because `advance()` reads any curve that reaches the requested
    depth and a smoke curve is shaped exactly like a real one.
    """
    f1 = art["finding_1_the_directory_level_does_not_survive_the_upload"]
    key = "the_fix_is_STRICTER_not_looser"
    assert key in f1, "the artifact must record that the replacement rule is stricter, not looser"
    txt = f1[key]
    assert "smoke" in txt and "ACCIDENT" in txt.upper()
    # And the pattern quoted in the artifact must be the one actually shipped.
    shipped = None
    for line in RUNNER.read_text(encoding="utf-8").splitlines():
        if line.startswith("CKPT_NAME_RE = re.compile("):
            shipped = line.split('r"', 1)[1].rsplit('"', 1)[0]
    assert shipped, "could not read the shipped pattern from the runner"
    assert shipped in txt, (
        "the artifact quotes a pattern that is not the shipped one: shipped %r" % shipped)


def test_the_unexcluded_hypothesis_is_marked_as_unexcluded(art):
    """A readiness race was neither demonstrated nor excluded, and must not read as measured."""
    f3 = art["finding_3_a_dataset_must_be_attachable_before_the_push"]
    assert f3["status"] == "OPEN_BUT_BOUNDED"
    assert "do_not_over_read" in f3
    assert "neither demonstrated nor excluded" in f3["observation"]


def test_the_end_to_end_confirmation_is_not_predicted(art):
    """The two findings were each the opposite of what was expected, so the third is left open."""
    c = art["end_to_end_confirmation"]
    assert c["status"] in {"IN_FLIGHT", "CONFIRMED", "FAILED"}
    assert c["probe"], "the confirmation must name the probe it depends on"
    assert c["how_to_read_the_outcome"], (
        "both outcomes must be given a reading BEFORE the result lands, or the reading is chosen "
        "after seeing the number")


def test_the_gpu_cost_of_this_step_is_recorded_as_zero_with_its_evidence(art):
    """GPU quota is the binding constraint, so a step that spent none should say how it managed it."""
    c = art["cost"]
    assert c["gpu_quota_consumed_by_this_step"] == "0.00 h"
    assert "CPU session" in c["how"]
    assert "0.41" in c["how"] and "0.88" in c["how"], (
        "the zero-cost claim needs the quota readings that support it, not just an assertion")


def test_it_states_what_has_NOT_been_round_tripped(art):
    """Every file in the probe was synthetic; claiming the round trip is proven would overstate it."""
    lim = art["honest_limits"]
    assert isinstance(lim, list) and len(lim) >= 3
    joined = " ".join(lim)
    assert "synthetic" in joined
    assert "No real checkpoint has made the round trip yet" in joined


def test_the_confirmation_records_the_real_mount_path_AND_why_it_is_not_baked(art):
    """The finding a later session is most likely to undo, so it is pinned with its reasoning.

    Now that `/kaggle/input/datasets/<owner>/<slug>` is known, baking it looks like an obvious
    simplification. It is the opposite: it re-introduces a prediction about Kaggle's layout, which is
    the thing that failed twice in one afternoon. The recursive fallback is layout-INDEPENDENT and
    survives a further change; a baked literal would not.
    """
    c = art["end_to_end_confirmation"]
    assert c["status"] == "CONFIRMED"
    assert c["seeded_count"] >= 1
    assert c["THE_REAL_MOUNT_PATH"] == "/kaggle/input/datasets/<owner>/<slug>"
    why = c["why_the_baked_path_is_STILL_not_changed_to_it"]
    assert "prediction" in why and "layout-independent" in why.lower()
    # The verbatim log is the evidence; a confirmation without it is an assertion.
    assert any("found in: /kaggle/input/datasets/" in line for line in c["verbatim_log"])


def test_the_confirmation_scopes_what_it_does_NOT_show(art):
    """The probe is a --smoke build and never calls run(), so it cannot show resume SKIPS epochs."""
    c = art["end_to_end_confirmation"]
    scope = c["what_this_does_and_does_not_establish"]
    assert "does NOT establish" in scope
    assert "smoke" in scope and "run()" in scope
    assert "epochs_added_by_this_link" in scope, (
        "the artifact must name what WILL test the unestablished half, or the gap reads as closed")


def test_the_incidental_cpu_rate_is_recorded_without_being_acted_on(art):
    """A route argument found while testing something else. Recorded, deliberately not acted on."""
    inc = art["incidental_measurement_worth_keeping"]
    assert "1758.7" in inc["what"] and "1276.04" in inc["what"]
    assert "authority call" in inc["do_not_act_on_this_here"]
    assert inc["caveat"], "a single measurement on each of two devices needs its caveat"
    assert "stable ratio" in inc["caveat"]

