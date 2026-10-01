"""Guards for the one-hot unseen-position headroom probe (2026-10-01).

This probe's job was to KILL a build before it was paid for: the orthogonal-modality channel (T2) needed
one-hot to have a fillable gap, and it measured that the gap is absorbed. These tests pin the parts that
make that negative trustworthy.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
_spec = importlib.util.spec_from_file_location(
    "hiv_onehot_unseen_position_headroom", REPO / "scripts" / "hiv_onehot_unseen_position_headroom.py")
P = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(P)

ART = REPO / "wiki" / "hiv_onehot_unseen_headroom_2026-10-01.json"


def _art():
    if not ART.exists():
        pytest.skip("artifact absent (gitignored Stanford data missing at build time)")
    return json.loads(ART.read_text(encoding="utf-8"))


def test_the_rule_was_frozen_before_any_number_and_reaches_all_three_verdicts():
    """A rule that can only return one verdict is not a test. self_check drives all three."""
    P.self_check()
    assert P.FROZEN_RULE["frozen_before_any_number_was_seen"] is True
    assert P.FROZEN_RULE["frozen"] == "2026-10-01"


def test_the_floor_anchor_reproduces_on_every_probed_gene():
    """Strata measured against a floor that does not reproduce are not interpretable."""
    d = _art()
    assert d["anchor"], "no anchor recorded"
    for cls, an in d["anchor"].items():
        assert an["reproduces"] is True, f"{cls}: {an}"
        assert abs(an["here"] - an["committed"]) <= P.ANCHOR_TOL


def test_the_blindness_is_MEASURED_not_assumed():
    """THE LOAD-BEARING CHECK. The whole hypothesis rests on one-hot giving zero weight to a column with
    no training signal. That is measured directly: max |coef| over those columns must be ~0."""
    d = _art()
    for cls, c in d["cells"].items():
        z = c["zero_training_weight_absmax"]
        assert z is not None, f"{cls}: the blindness was never measured"
        assert z == 0.0 or z < 1e-9, f"{cls}: zero-signal columns carry weight {z} -- hypothesis wrong"


def test_the_verdict_is_the_frozen_rule_recomputed_not_a_stored_string():
    d = _art()
    for cls, c in d["cells"].items():
        assert P.verdict(c)["verdict"] == c["verdict_block"]["verdict"], cls


def test_the_measured_outcome_is_NO_HEADROOM_on_the_powered_gene():
    """RT is the powered cell; it says a position-level prior has nothing to add."""
    d = _art()
    rt = d["cells"]["NNRTI"]
    assert rt["verdict_block"]["verdict"] == "NO_HEADROOM_ONEHOT_ALREADY_COVERS_IT"
    # exposure below the bar AND the gap in the wrong direction -- either alone would suffice
    assert rt["exposure"]["frac_blindspot_with_unseen"] < P.MIN_EXPOSED_FRAC
    assert rt["auroc_gap_unexposed_minus_exposed"] < P.AUROC_GAP


def test_the_exposed_stratum_does_not_score_worse_which_is_the_actual_finding():
    """Not 'we found no signal' -- the model scores the same or BETTER where it is structurally blind."""
    d = _art()
    rt = d["cells"]["NNRTI"]
    exposed = rt["strata"]["exposed"]["auroc"]
    unexposed = rt["strata"]["unexposed"]["auroc"]
    assert exposed is not None and unexposed is not None
    assert exposed >= unexposed - 0.02, (
        f"exposed {exposed} vs unexposed {unexposed} -- if the exposed stratum IS materially worse the "
        "verdict must be revisited")


def test_underpowered_strata_are_WITHHELD_not_reported():
    """INSTI's exposed stratum holds 8 R; an AUROC there would be noise dressed as a number."""
    d = _art()
    ins = d["cells"]["INSTI"]
    st = ins["strata"]["exposed"]
    if st["n_R"] < P.MIN_STRATUM or st["n_S"] < P.MIN_STRATUM:
        assert st["auroc"] is None, "an underpowered stratum must withhold its AUROC"
        assert "refused" in st
        assert ins["verdict_block"]["verdict"] == "UNDERPOWERED"


def test_the_redundancy_MECHANISM_is_measured_and_its_exception_is_recorded():
    """The cause must not be asserted -- this repo had three asserted causes refuted by measurement.
    Redundancy predicts a SMALL typical blind fraction; the RT maximum of 1.0 (one entirely-blind isolate)
    is a real exception that must stay visible rather than be smoothed away."""
    d = _art()
    for cls, c in d["cells"].items():
        r = c["redundancy"]
        assert r["mean_active_columns_per_isolate"] > 1, f"{cls}: redundancy needs several active columns"
        med = r["among_exposed_median_blind_fraction"]
        assert med is not None and med < 0.25, (
            f"{cls}: median blind fraction {med} is too high for redundancy to be the explanation")
        assert "among_exposed_MAX_blind_fraction" in r, "the exception must be reported, not hidden"
        assert "n_isolates_majority_blind" in r


def test_the_probe_bounds_ROOM_and_does_not_claim_a_prior_would_capture_it():
    d = _art()
    scope = " ".join(d["honest_scope"]).lower()
    assert "does not say such a prior would" in scope or "whether any room exists" in scope
    assert "per fold" in scope, "exposure must be stated as per-fold, the only honest denominator"


def test_the_feature_set_subtlety_is_disclosed():
    """A zero-training-signal column still EXISTS (featurization is global); what is missing is signal.
    Omitting that would make the measurement look like something it is not."""
    d = _art()
    scope = " ".join(d["honest_scope"]).lower()
    assert "globally" in scope and "still exists" in scope
