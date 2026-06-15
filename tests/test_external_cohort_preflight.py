"""Offline unit tests for the Gate-0 preflight verdict logic.

No network: pure verdict functions are tested directly; `preflight()` is tested
with an injected resolver backed by a fake fetch.
"""
from __future__ import annotations

import json

import scripts.external_cohort_preflight as pf
from dna_decode.eval.biosample_resolver import BioSampleResolver


# --------------------------------------------------------------------------- #
# classify_availability
# --------------------------------------------------------------------------- #
def test_classify_availability():
    out = pf.classify_availability({"SAMN1": ["GCA_1.1"], "SAMN2": [], "SAMN3": ["GCA_3.1"]})
    assert out["n_free"] == 2
    assert out["n_assembly_required"] == 1
    assert out["assembly_required"] == ["SAMN2"]


# --------------------------------------------------------------------------- #
# leakage_verdict
# --------------------------------------------------------------------------- #
def test_leakage_pass_disjoint_resolved():
    v = pf.leakage_verdict({"GCA_a": "SAMN_a", "GCA_b": "SAMN_b"}, cohort_biosamples={"SAMN_x"})
    assert v["passed"] is True
    assert v["overlap_biosamples"] == []
    assert v["n_unresolved"] == 0


def test_leakage_fail_on_overlap():
    v = pf.leakage_verdict({"GCA_a": "SAMN_shared"}, cohort_biosamples={"SAMN_shared", "SAMN_y"})
    assert v["passed"] is False
    assert v["fail_overlap"] is True
    assert v["overlap_biosamples"] == ["SAMN_shared"]


def test_leakage_fail_closed_on_unresolved_over_threshold():
    # 2 of 10 unresolved = 20% > 5% -> fail-closed
    m = {f"GCA_{i}": (None if i < 2 else f"SAMN_{i}") for i in range(10)}
    v = pf.leakage_verdict(m, cohort_biosamples={"SAMN_X"})
    assert v["fail_unresolved"] is True
    assert v["passed"] is False
    assert v["unresolved_fraction"] == 0.2


def test_leakage_within_threshold_passes():
    # 0 unresolved, disjoint
    m = {f"GCA_{i}": f"SAMN_{i}" for i in range(30)}
    v = pf.leakage_verdict(m, cohort_biosamples={"SAMN_OUT"})
    assert v["passed"] is True


def test_leakage_exactly_at_threshold_passes():
    # 1 of 20 unresolved = exactly 5.0% — the gate is strictly `>` so this PASSES.
    m = {f"GCA_{i}": (None if i == 0 else f"SAMN_{i}") for i in range(20)}
    v = pf.leakage_verdict(m, cohort_biosamples={"SAMN_OUT"})
    assert v["unresolved_fraction"] == 0.05
    assert v["fail_unresolved"] is False
    assert v["passed"] is True


def test_leakage_just_over_threshold_fails():
    # 2 of 20 unresolved = 10% > 5% -> fail-closed (the adjacent boundary case).
    m = {f"GCA_{i}": (None if i < 2 else f"SAMN_{i}") for i in range(20)}
    v = pf.leakage_verdict(m, cohort_biosamples={"SAMN_OUT"})
    assert v["fail_unresolved"] is True
    assert v["passed"] is False


def test_leakage_empty_tuning_set():
    # No tuning accessions -> 0/0 fraction is 0.0, not a ZeroDivisionError; passes.
    v = pf.leakage_verdict({}, cohort_biosamples={"SAMN_OUT"})
    assert v["unresolved_fraction"] == 0.0
    assert v["passed"] is True


# --------------------------------------------------------------------------- #
# overall_verdict
# --------------------------------------------------------------------------- #
def _good_leakage():
    return pf.leakage_verdict({"GCA_a": "SAMN_a"}, cohort_biosamples={"SAMN_x"})


def _avail(n_free):
    return pf.classify_availability({f"SAMN{i}": ["GCA"] for i in range(n_free)} or {"SAMN0": []})


def test_overall_pass():
    out = pf.overall_verdict(_good_leakage(), _avail(3), mic_open=True)
    assert out["verdict"] == "PASS"
    assert out["reasons"] == []


def test_overall_fail_mic_gated():
    out = pf.overall_verdict(_good_leakage(), _avail(3), mic_open=False)
    assert out["verdict"] == "FAIL"
    assert any("MTA-gated" in r or "not openly" in r for r in out["reasons"])


def test_overall_fail_mic_unconfirmed():
    out = pf.overall_verdict(_good_leakage(), _avail(3), mic_open=None)
    assert out["verdict"] == "FAIL"


def test_overall_fail_zero_free():
    out = pf.overall_verdict(_good_leakage(), _avail(0), mic_open=True)
    assert out["verdict"] == "FAIL"
    assert any("free pilot N is zero" in r for r in out["reasons"])


def test_overall_fail_closed_on_incomplete_manifest():
    out = pf.overall_verdict(_good_leakage(), _avail(3), mic_open=True,
                             manifest_incomplete=True, allow_degraded=False)
    assert out["verdict"] == "FAIL"
    assert any("INCOMPLETE leakage manifest" in r for r in out["reasons"])


def test_overall_incomplete_manifest_degraded_override():
    out = pf.overall_verdict(_good_leakage(), _avail(3), mic_open=True,
                             manifest_incomplete=True, allow_degraded=True)
    assert out["verdict"] == "PASS"                      # override clears the block...
    assert out["manifest_degraded_override"] is True     # ...but is stamped degraded


def test_overall_accumulates_multiple_reasons():
    # mic-gated AND zero-free AND an overlap leak all surface (no short-circuit).
    leak = pf.leakage_verdict({"GCA_a": "SAMN_shared"}, cohort_biosamples={"SAMN_shared"})
    out = pf.overall_verdict(leak, _avail(0), mic_open=False)
    assert out["verdict"] == "FAIL"
    assert any("overlap" in r for r in out["reasons"])
    assert any("MTA-gated" in r for r in out["reasons"])
    assert any("free pilot N is zero" in r for r in out["reasons"])
    assert len(out["reasons"]) >= 3


# --------------------------------------------------------------------------- #
# preflight() orchestration with injected resolver (no network, no manifest deps)
# --------------------------------------------------------------------------- #
def _uid_aware_fetch(read_run_tsv: str, esearch_term_to_uid: dict[str, str], uid_to_record: dict[str, dict]):
    """Fake fetch that routes esearch by term substring and esummary by uid, so the
    cohort-biosample lookups and the tuning-accession lookup get DISTINCT records."""
    def fetch(url: str) -> str:
        if "result=read_run" in url:
            return read_run_tsv
        if "esearch.fcgi" in url:
            for term, uid in esearch_term_to_uid.items():
                if term in url:
                    return json.dumps({"esearchresult": {"idlist": [uid]}})
            return json.dumps({"esearchresult": {"idlist": []}})
        if "esummary.fcgi" in url:
            for uid, rec in uid_to_record.items():
                if f"id={uid}" in url:
                    return json.dumps({"result": {"uids": [uid], uid: rec}})
            return json.dumps({"result": {"uids": []}})
        if "result=assembly" in url:  # ENA fallback unused in these fixtures
            return "assembly_accession\tsample_accession\n"
        raise AssertionError(url)
    return fetch


def test_preflight_pass(tmp_path, monkeypatch):
    fetch = _uid_aware_fetch(
        read_run_tsv="run_accession\tsample_accession\nERR1\tSAMEA_C1\nERR2\tSAMEA_C2\n",
        esearch_term_to_uid={"SAMEA_C1": "1", "SAMEA_C2": "2", "GCA_TUNING.1": "999"},
        uid_to_record={
            "1": {"assemblyaccession": "GCA_C1.1", "biosampleaccn": "SAMEA_C1"},
            "2": {"assemblyaccession": "GCA_C2.1", "biosampleaccn": "SAMEA_C2"},
            "999": {"assemblyaccession": "GCA_TUNING.1", "biosampleaccn": "SAMN_TUNING_DISJOINT"},
        },
    )
    resolver = BioSampleResolver(cache_path=tmp_path / "cache.json", fetch=fetch)

    # tuning manifest -> a single tuning accession that resolves to a DISJOINT biosample
    class _M:
        incomplete = False
        cohorts = [1, 2, 3]
    monkeypatch.setattr(pf, "build_manifest", lambda: _M())
    monkeypatch.setattr(pf, "prior_accessions", lambda m, exclude_cohort: {"GCA_TUNING.1"})

    art = pf.preflight("PRJEB_TEST", "spain_probac", mic_open=True,
                       resolver=resolver, wiki_dir=tmp_path, write=True)
    assert art["verdict"] == "PASS"
    assert art["assembly_availability"]["n_free"] >= 1
    assert art["leakage"]["overlap_biosamples"] == []
    assert (tmp_path / f"external_preflight_spain_probac_{art['date']}.json").exists()


def test_preflight_fail_on_leak(tmp_path, monkeypatch):
    # The tuning accession resolves to a biosample that IS in the cohort -> overlap -> FAIL.
    fetch = _uid_aware_fetch(
        read_run_tsv="run_accession\tsample_accession\nERR1\tSAMEA_SHARED\n",
        esearch_term_to_uid={"SAMEA_SHARED": "1", "GCA_TUNING.1": "1"},  # tuning maps to same uid -> SAMEA_SHARED
        uid_to_record={"1": {"assemblyaccession": "GCA_x.1", "biosampleaccn": "SAMEA_SHARED"}},
    )
    resolver = BioSampleResolver(cache_path=tmp_path / "cache.json", fetch=fetch)

    class _M:
        incomplete = False
        cohorts = [1]
    monkeypatch.setattr(pf, "build_manifest", lambda: _M())
    monkeypatch.setattr(pf, "prior_accessions", lambda m, exclude_cohort: {"GCA_TUNING.1"})

    art = pf.preflight("PRJEB_TEST", "oxford", mic_open=True,
                       resolver=resolver, wiki_dir=tmp_path, write=False)
    assert art["verdict"] == "FAIL"
    assert "SAMEA_SHARED" in art["leakage"]["overlap_biosamples"]


def test_preflight_fail_closed_on_incomplete_manifest(tmp_path, monkeypatch):
    fetch = _uid_aware_fetch(
        read_run_tsv="run_accession\tsample_accession\nERR1\tSAMEA_C1\n",
        esearch_term_to_uid={"SAMEA_C1": "1", "GCA_TUNING.1": "999"},
        uid_to_record={
            "1": {"assemblyaccession": "GCA_C1.1", "biosampleaccn": "SAMEA_C1"},
            "999": {"assemblyaccession": "GCA_TUNING.1", "biosampleaccn": "SAMN_DISJOINT"},
        },
    )
    resolver = BioSampleResolver(cache_path=tmp_path / "cache.json", fetch=fetch)

    class _M:                       # manifest that failed to load a source
        incomplete = True
        cohorts = [1]
        warnings = ["parquet load failed: x"]
    monkeypatch.setattr(pf, "build_manifest", lambda: _M())
    monkeypatch.setattr(pf, "prior_accessions", lambda m, exclude_cohort: {"GCA_TUNING.1"})

    art = pf.preflight("PRJEB_TEST", "oxford", mic_open=True, resolver=resolver,
                       wiki_dir=tmp_path, write=False, allow_degraded=False)
    assert art["verdict"] == "FAIL"
    assert art["manifest_complete"] is False
    # with override -> proceeds, stamped degraded
    art2 = pf.preflight("PRJEB_TEST", "oxford", mic_open=True, resolver=resolver,
                        wiki_dir=tmp_path, write=False, allow_degraded=True)
    assert art2["verdict"] == "PASS"
    assert art2["manifest_degraded"] is True
