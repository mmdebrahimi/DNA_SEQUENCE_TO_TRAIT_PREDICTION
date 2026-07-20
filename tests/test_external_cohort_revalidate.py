"""Unit tests for the external-cohort scorer (no network/Docker — fake predict)."""
from __future__ import annotations

import scripts.external_cohort_revalidate as ecr


# --------------------------------------------------------------------------- #
# Pinned organism triple (must match the frozen E. coli provdisjoint cells)
# --------------------------------------------------------------------------- #
def test_pinned_organism_triple():
    assert ecr.AMRFINDER_ORGANISM == "Escherichia"
    assert ecr.REGISTRY_ORGANISM == "Escherichia_coli_Shigella"


# --------------------------------------------------------------------------- #
# score_label_set
# --------------------------------------------------------------------------- #
def test_score_perfect():
    free = {"SAMN_r": "GCA_r.1", "SAMN_s": "GCA_s.1"}
    labels = {"SAMN_r": "R", "SAMN_s": "S"}
    preds = {"GCA_r.1": "R", "GCA_s.1": "S"}
    conf = ecr.score_label_set(free, labels, lambda g: preds[g])
    assert conf["n_scored"] == 2
    assert conf["sens"] == 1.0 and conf["spec"] == 1.0
    assert conf["n_excluded_no_assembly"] == 0


def test_score_excludes_assembly_required():
    free = {"SAMN_r": "GCA_r.1"}                       # SAMN_s has no assembly
    labels = {"SAMN_r": "R", "SAMN_s": "S"}
    conf = ecr.score_label_set(free, labels, lambda g: "R")
    assert conf["n_scored"] == 1               # only the FREE one scored
    assert conf["n_excluded_no_assembly"] == 1


def test_score_abstain_excluded_from_n():
    free = {"SAMN_r": "GCA_r.1", "SAMN_s": "GCA_s.1"}
    labels = {"SAMN_r": "R", "SAMN_s": "S"}
    preds = {"GCA_r.1": "R", "GCA_s.1": "ABSTAIN"}
    conf = ecr.score_label_set(free, labels, lambda g: preds[g])
    assert conf["abstain"] == 1
    assert conf["n_scored"] == 1               # abstain not in n_scored


def test_score_indeterminate_uncounted():
    # INDETERMINATE (absent main.tsv) is neither scored nor an abstain (mirrors provdisjoint).
    free = {"SAMN_r": "GCA_r.1"}
    labels = {"SAMN_r": "R"}
    conf = ecr.score_label_set(free, labels, lambda g: "INDETERMINATE")
    assert conf["n_scored"] == 0
    assert conf["abstain"] == 0


# --------------------------------------------------------------------------- #
# predict_records / conf_from_records (Step-5 substrate)
# --------------------------------------------------------------------------- #
def test_predict_records_one_call_per_free_strain():
    free = {"SAMN_r": "GCA_r.1", "SAMN_s": "GCA_s.1"}
    labels = {"SAMN_r": "R", "SAMN_s": "S", "SAMN_noasm": "R"}
    calls = []

    def predict(g):
        calls.append(g)
        return "R" if g == "GCA_r.1" else "S"

    records, n_excl = ecr.predict_records(free, labels, predict)
    assert len(records) == 2 and n_excl == 1          # SAMN_noasm excluded
    assert len(calls) == 2                            # one predict per FREE strain
    rec_r = next(r for r in records if r["biosample"] == "SAMN_r")
    assert rec_r["gca"] == "GCA_r.1" and rec_r["prediction"] == "R" and rec_r["y"] == 1


def test_conf_from_records_matches_score_label_set():
    free = {"SAMN_r": "GCA_r.1", "SAMN_s": "GCA_s.1"}
    labels = {"SAMN_r": "R", "SAMN_s": "S"}
    preds = {"GCA_r.1": "R", "GCA_s.1": "S"}
    records, n_excl = ecr.predict_records(free, labels, lambda g: preds[g])
    conf_direct = ecr.conf_from_records(records, n_excl)
    conf_wrapper = ecr.score_label_set(free, labels, lambda g: preds[g])
    assert conf_direct == conf_wrapper


# --------------------------------------------------------------------------- #
# build_artifact
# --------------------------------------------------------------------------- #
def test_assert_manifest_alignment_ok():
    ecr.assert_manifest_alignment({"SAMEA1", "SAMEA2"}, {"SAMEA1", "SAMEA2", "SAMEA3"})  # subset -> no raise


def test_assert_manifest_alignment_drift_raises():
    import pytest
    with pytest.raises(ecr.ManifestDriftError):
        ecr.assert_manifest_alignment({"SAMEA1", "SAMEA_EXTRA"}, {"SAMEA1"})


def test_build_artifact_carries_run_id():
    art = ecr.build_artifact("oxford", "ciprofloxacin", strict={}, relaxed={}, buckets={},
                             leakage_control="x", run_id="run_42")
    assert art["run_id"] == "run_42"


def test_build_artifact_schema():
    art = ecr.build_artifact("spain_probac", "ciprofloxacin",
                             strict={"n_scored": 10, "acc": 0.9, "sens": 0.9, "spec": 0.9},
                             relaxed={"n_scored": 14, "acc": 0.85, "sens": 0.86, "spec": 0.84},
                             buckets={"buckets": {"HIGH_R": 5}}, leakage_control="x")
    assert art["_schema"] == "external-validation-v1"
    assert art["evidence_tier"] == "external_clinical"
    assert art["amrfinder_organism"] == "Escherichia"
    assert art["registry_organism"] == "Escherichia_coli_Shigella"
    assert art["primary_metric"] == "strict"
    assert art["strict"]["n_scored"] == 10 and art["relaxed"]["n_scored"] == 14


def test_build_artifact_organism_override():
    """The organism triple is parameterizable (Klebsiella etc.) but E. coli stays the default."""
    art = ecr.build_artifact("ar_bank_kleb", "ceftriaxone", strict={}, relaxed={}, buckets={},
                             leakage_control="x", amrfinder_organism="Klebsiella_pneumoniae",
                             registry_organism="Klebsiella")
    assert art["amrfinder_organism"] == "Klebsiella_pneumoniae"
    assert art["registry_organism"] == "Klebsiella"
    assert art["organism"] == "Klebsiella"
    # default is still E. coli (backward-compatible)
    d = ecr.build_artifact("oxford", "ciprofloxacin", strict={}, relaxed={}, buckets={},
                           leakage_control="x")
    assert d["amrfinder_organism"] == "Escherichia" and d["registry_organism"] == "Escherichia_coli_Shigella"


# --------------------------------------------------------------------------- #
# smoke_predict (fail-fast before the full loop)
# --------------------------------------------------------------------------- #
def test_smoke_ok():
    sm = ecr.smoke_predict({"SAMN_a": "GCA_a.1"}, lambda g: "R")
    assert sm["ok"] is True and sm["gca"] == "GCA_a.1"


def test_smoke_fail_on_indeterminate():
    sm = ecr.smoke_predict({"SAMN_a": "GCA_a.1"}, lambda g: "INDETERMINATE")
    assert sm["ok"] is False and "broken" in sm["reason"]


def test_smoke_empty_genomes():
    sm = ecr.smoke_predict({}, lambda g: "R")
    assert sm["ok"] is False


def test_smoke_picks_deterministic_first():
    seen = []
    ecr.smoke_predict({"SAMN_b": "GCA_z.1", "SAMN_a": "GCA_a.1"}, lambda g: seen.append(g) or "S")
    assert seen == ["GCA_a.1"]                # sorted-first, deterministic


# --------------------------------------------------------------------------- #
# powering_gate (fail-open guard)
# --------------------------------------------------------------------------- #
def _powered_conf(nR=12, nS=12):
    # all correct: tp=nR, fn=0, tn=nS, fp=0
    return {"n_scored": nR + nS, "tp": nR, "fn": 0, "tn": nS, "fp": 0}


def test_powering_pass():
    pg = ecr.powering_gate(_powered_conf(12, 12), n_attempted_free=24, n_indeterminate=0)
    assert pg["hard_fail"] is False and pg["degraded"] is False
    assert pg["scored_R"] == 12 and pg["scored_S"] == 12


def test_powering_hard_fail_zero_scored():
    pg = ecr.powering_gate({"n_scored": 0, "tp": 0, "fn": 0, "tn": 0, "fp": 0},
                           n_attempted_free=20, n_indeterminate=20)
    assert pg["hard_fail"] is True
    assert any("n_scored == 0" in r for r in pg["reasons"])


def test_powering_hard_fail_below_class_floor():
    # 4R/4S scored, no indeterminates -> non-empty but UNDERPOWERED -> hard fail
    pg = ecr.powering_gate(_powered_conf(4, 4), n_attempted_free=8, n_indeterminate=0)
    assert pg["hard_fail"] is True
    assert any("underpowered" in r for r in pg["reasons"])


def test_powering_pilot_threshold_override():
    # lowering min_per_class to a documented pilot threshold clears the floor
    pg = ecr.powering_gate(_powered_conf(4, 4), n_attempted_free=8, n_indeterminate=0, min_per_class=4)
    assert pg["hard_fail"] is False


def test_powering_degraded_high_indeterminate():
    # 12R/12S scored (powered) but 10 of 34 attempted-FREE were indeterminate -> ~0.29 > 0.20
    pg = ecr.powering_gate(_powered_conf(12, 12), n_attempted_free=34, n_indeterminate=10)
    assert pg["hard_fail"] is False
    assert pg["degraded"] is True
    assert pg["indeterminate_fraction"] > 0.20


def test_powering_reads_only_exclusions_dont_count():
    # attempted-FREE denominator excludes reads-only; 0 indeterminate of 24 attempted -> not degraded
    pg = ecr.powering_gate(_powered_conf(12, 12), n_attempted_free=24, n_indeterminate=0)
    assert pg["degraded"] is False and pg["indeterminate_fraction"] == 0.0


def test_build_artifact_carries_powering():
    pg = ecr.powering_gate(_powered_conf(12, 12), n_attempted_free=24, n_indeterminate=0)
    art = ecr.build_artifact("oxford", "ciprofloxacin", strict=_powered_conf(12, 12),
                             relaxed={}, buckets={}, leakage_control="x",
                             powering=pg, run_degraded=pg["degraded"])
    assert art["powering"]["scored_R"] == 12
    assert art["run_degraded"] is False


# --------------------------------------------------------------------------- #
# gate_ok (fail-closed)
# --------------------------------------------------------------------------- #
def test_gate_pass():
    ok, _ = ecr.gate_ok({"verdict": "PASS"}, allow_degraded=False)
    assert ok is True


def test_gate_fail_closed_on_fail():
    ok, reason = ecr.gate_ok({"verdict": "FAIL", "reasons": ["leak"]}, allow_degraded=False)
    assert ok is False and "fail-closed" in reason


def test_gate_fail_closed_on_missing():
    ok, _ = ecr.gate_ok(None, allow_degraded=False)
    assert ok is False


def test_gate_degraded_override():
    ok, reason = ecr.gate_ok({"verdict": "FAIL", "reasons": []}, allow_degraded=True)
    assert ok is True and "DEGRADED" in reason


def test_gate_missing_preflight_degraded_override():
    # No preflight artifact AT ALL + --allow-degraded -> proceed (distinct branch from FAIL+degraded).
    ok, reason = ecr.gate_ok(None, allow_degraded=True)
    assert ok is True and "DEGRADED" in reason


# --------------------------------------------------------------------------- #
# build_artifact degraded flag
# --------------------------------------------------------------------------- #
def test_build_artifact_degraded_flag():
    art = ecr.build_artifact("oxford", "ciprofloxacin", strict={}, relaxed={},
                             buckets={}, leakage_control="x", degraded=True)
    assert art["independence_degraded"] is True
    # default is non-degraded
    art2 = ecr.build_artifact("oxford", "ciprofloxacin", strict={}, relaxed={},
                              buckets={}, leakage_control="x")
    assert art2["independence_degraded"] is False


def test_predict_records_label_case_insensitive():
    # Lowercase "r"/"s" labels normalize via .upper() -> correct y + uppercased label.
    free = {"SAMN_r": "GCA_r.1", "SAMN_s": "GCA_s.1"}
    labels = {"SAMN_r": "r", "SAMN_s": " s "}
    records, _ = ecr.predict_records(free, labels, lambda g: "R")
    by_bs = {r["biosample"]: r for r in records}
    assert by_bs["SAMN_r"]["label"] == "R" and by_bs["SAMN_r"]["y"] == 1
    assert by_bs["SAMN_s"]["label"] == "S" and by_bs["SAMN_s"]["y"] == 0


def test_assert_manifest_alignment_empty_selected_ok():
    # An empty scored set is trivially a subset -> no drift error.
    ecr.assert_manifest_alignment(set(), {"SAMEA1", "SAMEA2"})


def test_powering_degraded_zero_attempted_no_zero_division():
    # n_attempted_free == 0 -> indeterminate_fraction is 0.0 (guarded), not a crash.
    pg = ecr.powering_gate({"n_scored": 0, "tp": 0, "fn": 0, "tn": 0, "fp": 0},
                           n_attempted_free=0, n_indeterminate=0)
    assert pg["indeterminate_fraction"] == 0.0
    assert pg["degraded"] is False


def test_predict_records_empty_labels():
    records, n_excl = ecr.predict_records({"SAMN_a": "GCA_a.1"}, {}, lambda g: "R")
    assert records == [] and n_excl == 0
    conf = ecr.conf_from_records(records, n_excl)
    assert conf["n_scored"] == 0 and conf["acc"] is None
