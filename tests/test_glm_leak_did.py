"""Tests for the leak DiD: the geometry gate, the delta orientation, the bar, and the role split.

Four things here can fail silently and each would produce a confident wrong number, so each is
pinned:

1. **The geometry gate.** Window selection is `random.Random(seed).sample(pool, n)` and `sample`
   has no prefix property, so artifacts scored at different `n` describe DIFFERENT windows. Equal
   length does not detect it -- two different 120-window draws are both length 120. Pairing across
   draws must REFUSE, not compute.
2. **The delta orientation.** `committed - assessed` with positive meaning the assessed model is
   better. A flipped sign would invert every verdict while every interval still looked healthy.
3. **The bar.** All four outcomes must be reachable, and an unresolved interval must never count as
   a win however large its point estimate.
4. **The role split.** `run2_reseed_seed1` is the CONTROL, and nothing in that label says so --
   only its `ckpt_dir` does. The first version of the synthesis read it as a second treatment and
   published a sign disagreement that does not exist, so the derivation is pinned here.
"""
from __future__ import annotations

import json

import pytest

from scripts.glm_leak_did import (
    PREDICTED_DID_NATS,
    DidError,
    check_geometry,
    deltas,
    did_with_ci,
    verdict_for,
)
from scripts.glm_leak_did_synthesis import CONTROL, TREATMENT, SynthError, apply_gate, role_for


def _art(eval_set: str, *, n=120, positions=30, seed=0, selection_seed=0, vals=None):
    rows = [{"left_only": v, "bidirectional": v} for v in (vals or [0.0] * n)]
    return {"run_args": {"n_windows": n, "positions": positions, "seed": seed,
                         "selection_seed": selection_seed, "eval_set": eval_set},
            "per_window": rows}


def _arts(**over):
    base = {
        "committed_bearing": (None, _art("ortholog_bearing")),
        "committed_clean": (None, _art("withheld_orthology_free")),
        "assessed_bearing": (None, _art("ortholog_bearing")),
        "assessed_clean": (None, _art("withheld_orthology_free")),
    }
    base.update({k: (None, v) for k, v in over.items()})
    return base


# ------------------------------------------------------------------------------- geometry gate


def test_geometry_gate_accepts_a_consistent_set():
    geoms = check_geometry(_arts())
    assert set(geoms) == {"committed_bearing", "committed_clean",
                          "assessed_bearing", "assessed_clean"}


def test_geometry_gate_REFUSES_a_120_against_a_500_draw():
    """The advisor's case and the real hazard: a different n is a different window set."""
    with pytest.raises(DidError, match="draw geometry"):
        check_geometry(_arts(assessed_bearing=_art("ortholog_bearing", n=500),
                             assessed_clean=_art("withheld_orthology_free", n=500)))


@pytest.mark.parametrize("field,bad", [("positions", 10), ("seed", 7), ("selection_seed", 3)])
def test_geometry_gate_refuses_any_differing_draw_parameter(field, bad):
    """Not just n: positions, the scoring seed and the selection seed each change the draw."""
    with pytest.raises(DidError, match="draw geometry"):
        check_geometry(_arts(assessed_clean=_art("withheld_orthology_free", **{field: bad})))


def test_geometry_gate_refuses_a_set_on_the_wrong_eval_set():
    """A clean-labelled artifact actually scored on the bearing half must not pass."""
    with pytest.raises(DidError, match="expected"):
        check_geometry(_arts(assessed_clean=_art("ortholog_bearing")))


def test_geometry_gate_refuses_when_rows_disagree_with_declared_n():
    art = _art("ortholog_bearing", n=120)
    art["per_window"] = art["per_window"][:100]
    with pytest.raises(DidError, match="rows vs n_windows"):
        check_geometry(_arts(assessed_bearing=art))


# ------------------------------------------------------------------------------- delta orientation


def test_delta_is_committed_minus_assessed_so_positive_means_assessed_better():
    committed = _art("ortholog_bearing", n=2, vals=[1.5, 1.5])
    assessed = _art("ortholog_bearing", n=2, vals=[1.4, 1.6])
    d = deltas(committed, assessed, arm="left_only")
    assert d == pytest.approx([0.1, -0.1])
    # a LOWER assessed cross-entropy (1.4 < 1.5) is the assessed model being better -> positive
    assert d[0] > 0


def test_deltas_refuse_mismatched_lengths():
    with pytest.raises(DidError, match="window counts differ"):
        deltas(_art("ortholog_bearing", n=3), _art("ortholog_bearing", n=2), arm="left_only")


# ------------------------------------------------------------------------------- the bar


def _stat(did, lo, hi):
    return {"did_nats": did, "ci95": [lo, hi],
            "excludes_zero": bool(lo > 0 or hi < 0),
            "excludes_predicted": bool(lo > PREDICTED_DID_NATS or hi < PREDICTED_DID_NATS)}


def test_every_named_outcome_is_reachable():
    from scripts.glm_leak_did import PREREGISTERED

    reached = {
        verdict_for(_stat(-0.0054, -0.0070, -0.0040))["VERDICT"],   # excl zero, contains predicted
        verdict_for(_stat(-0.0020, -0.0035, -0.0008))["VERDICT"],   # excl zero AND excl predicted
        verdict_for(_stat(-0.0008, -0.0031, 0.0013))["VERDICT"],    # contains zero, excl predicted
        verdict_for(_stat(-0.0030, -0.0090, 0.0030))["VERDICT"],    # contains both
    }
    assert reached == set(PREREGISTERED["outcomes"])


def test_the_500_window_left_only_cell_reads_leak_excluded():
    """The real measured cell, so the bar is pinned against a number that actually occurred."""
    v = verdict_for(_stat(-0.000844, -0.003106, 0.001353))
    assert v["VERDICT"] == "LEAK_EXCLUDED_AT_PREDICTED_SIZE"


def test_an_interval_spanning_zero_is_never_a_resolved_leak():
    v = verdict_for(_stat(-0.0050, -0.0120, 0.0020))
    assert v["VERDICT"] == "STILL_UNRESOLVED"


def test_did_ci_is_deterministic_and_orientation_preserving():
    a = did_with_ci([0.002] * 40, [0.0] * 40, seed=1, boot=200)
    b = did_with_ci([0.002] * 40, [0.0] * 40, seed=1, boot=200)
    assert a == b
    assert a["did_nats"] == pytest.approx(0.002)


def test_did_refuses_an_empty_half():
    with pytest.raises(DidError):
        did_with_ci([], [0.1], seed=1, boot=10)


# ------------------------------------------------------------------------------- the role split


_MANIFEST = {"arms": {
    "run1_representative": {"ckpt_dir": "data/glm_ckpt_representative"},
    "run2_reseed_seed1": {"ckpt_dir": "data/glm_ckpt_reseed_control"},
}}


def test_roles_are_derived_from_ckpt_dir():
    assert role_for("run1_representative", _MANIFEST) == TREATMENT
    assert role_for("run2_reseed_seed1", _MANIFEST) == CONTROL


def test_the_role_comes_from_the_ckpt_dir_and_NOT_from_the_label():
    """THE defect that shipped: nothing in `run2_reseed_seed1` says `control`.

    A label engineered to look like a treatment while pointing at the control directory must still
    resolve to CONTROL, and vice versa -- otherwise the derivation is really trusting the name.
    """
    m = {"arms": {"looks_like_a_treatment": {"ckpt_dir": "data/glm_ckpt_reseed_control"},
                  "reseed_control_sounding_name": {"ckpt_dir": "data/glm_ckpt_representative"}}}
    assert role_for("looks_like_a_treatment", m) == CONTROL
    assert role_for("reseed_control_sounding_name", m) == TREATMENT


def test_role_refuses_an_unknown_checkpoint_directory():
    m = {"arms": {"mystery": {"ckpt_dir": "data/glm_ckpt_something_else"}}}
    with pytest.raises(SynthError, match="neither a treatment nor a control"):
        role_for("mystery", m)


def test_role_refuses_a_label_absent_from_the_manifest():
    with pytest.raises(SynthError, match="not in the kernel manifest"):
        role_for("never_ran", _MANIFEST)


def test_the_live_manifest_still_assigns_both_roles():
    """Pinned against the real file: if a future kernel renames a directory, this fires."""
    from scripts.glm_leak_did_synthesis import MANIFEST

    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    roles = {lab: role_for(lab, manifest) for lab in
             ("run1_representative", "run2_reseed_seed1")}
    assert roles == {"run1_representative": TREATMENT, "run2_reseed_seed1": CONTROL}


# ------------------------------------------------------------------------------- the gate


def _cell(did, lo, hi, *, arm="left_only", label="x"):
    return {"arm": arm, "assessed_label": label, "result": _stat(did, lo, hi)}


def test_gate_requires_BOTH_halves_exceeding_the_floor_and_excluding_zero():
    treat = _cell(-0.0200, -0.0300, -0.0100)   # big and excludes zero
    ctrl = _cell(0.0008, -0.0016, 0.0033)
    assert apply_gate(treat, ctrl)["VERDICT"] == "LEAK_IDENTIFIED"


def test_gate_fails_when_the_interval_includes_zero_even_though_it_out_sizes_the_floor():
    """The real 500-window shape: exceeds the floor, CI still spans zero -> not identifiable."""
    g = apply_gate(_cell(-0.001970, -0.004237, 0.000411), _cell(0.000922, -0.001581, 0.003514))
    assert g["exceeds_floor"] is True
    assert g["treatment_ci_excludes_zero"] is False
    assert g["VERDICT"] == "LEAK_NOT_IDENTIFIABLE"


def test_gate_fails_when_the_effect_does_not_out_size_its_floor():
    g = apply_gate(_cell(-0.0005, -0.0009, -0.0001), _cell(0.0028, -0.0018, 0.0078))
    assert g["exceeds_floor"] is False
    assert g["VERDICT"] == "LEAK_NOT_IDENTIFIABLE"
