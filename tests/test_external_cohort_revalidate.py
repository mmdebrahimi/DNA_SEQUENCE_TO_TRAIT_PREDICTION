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
