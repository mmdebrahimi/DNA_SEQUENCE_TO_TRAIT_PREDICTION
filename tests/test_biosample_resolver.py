"""Offline unit tests for the bidirectional BioSample resolver.

No network: every test injects a fake `fetch(url) -> str` returning fixture text,
or exercises the pure parsers directly.
"""
from __future__ import annotations

import json

from dna_decode.eval.biosample_resolver import (
    BioSampleResolver,
    _reconcile_biosample,
    parse_ena_assembly,
    parse_ena_read_run,
    parse_ena_read_run_records,
    parse_entrez_assembly_summary,
    parse_entrez_esearch,
)


# --------------------------------------------------------------------------- #
# parse_ena_read_run_records (multi-field; additive, Step 1)
# --------------------------------------------------------------------------- #
def test_parse_read_run_records_all_fields():
    tsv = ("run_accession\tsample_accession\tsample_alias\tsecondary_sample_accession\n"
           "ERR1\tSAMEA1\taliasA\tERS1\n"
           "ERR2\tSAMEA2\taliasB\tERS2\n")
    fields = ("run_accession", "sample_accession", "sample_alias", "secondary_sample_accession")
    recs = parse_ena_read_run_records(tsv, fields)
    assert recs[0] == {"run_accession": "ERR1", "sample_accession": "SAMEA1",
                       "sample_alias": "aliasA", "secondary_sample_accession": "ERS1"}
    assert len(recs) == 2


def test_parse_read_run_records_missing_optional_column():
    # only run+sample present; requested alias/secondary absent -> omitted, not crash
    tsv = "run_accession\tsample_accession\nERR1\tSAMEA1\n"
    recs = parse_ena_read_run_records(tsv, ("run_accession", "sample_accession",
                                           "sample_alias", "secondary_sample_accession"))
    assert recs == [{"run_accession": "ERR1", "sample_accession": "SAMEA1"}]


def test_parse_read_run_records_empty():
    assert parse_ena_read_run_records("", ("run_accession",)) == []
    assert parse_ena_read_run_records("nomatch\tcols\nx\ty\n", ("run_accession",)) == []


def test_read_run_records_for_project(tmp_path):
    tsv = ("run_accession\tsample_accession\tsample_alias\tsecondary_sample_accession\n"
           "ERR1\tSAMEA1\taliasA\tERS1\n")

    def fetch(url):
        assert "result=read_run" in url and "sample_alias" in url
        return tsv

    r = BioSampleResolver(cache_path=tmp_path / "c.json", fetch=fetch)
    recs = r.read_run_records_for_project("PRJEB1")
    assert recs[0]["sample_alias"] == "aliasA"


# --------------------------------------------------------------------------- #
# Pure parsers
# --------------------------------------------------------------------------- #
def test_parse_ena_read_run():
    tsv = ("run_accession\tsample_accession\n"
           "ERR111\tSAMEA1\n"
           "ERR222\tSAMEA2\n"
           "\t\n"          # blank-ish row -> skipped
           "ERR333\t\n")   # missing biosample -> skipped
    assert parse_ena_read_run(tsv) == [("ERR111", "SAMEA1"), ("ERR222", "SAMEA2")]


def test_parse_ena_read_run_no_header():
    assert parse_ena_read_run("") == []
    assert parse_ena_read_run("foo\tbar\nx\ty\n") == []  # no required columns


def test_parse_ena_assembly():
    tsv = "assembly_accession\tsample_accession\nGCA_1.1\tSAMEA1\nGCA_2.1\tSAMEA2\n"
    assert parse_ena_assembly(tsv) == ["GCA_1.1", "GCA_2.1"]


def test_parse_entrez_esearch():
    js = json.dumps({"esearchresult": {"idlist": ["100", "200"]}})
    assert parse_entrez_esearch(js) == ["100", "200"]
    assert parse_entrez_esearch("not json") == []


def test_parse_entrez_assembly_summary():
    js = json.dumps({"result": {
        "uids": ["100"],
        "100": {"assemblyaccession": "GCA_9.1", "biosampleaccn": "SAMN9"},
    }})
    assert parse_entrez_assembly_summary(js) == [{"assembly": "GCA_9.1", "biosample": "SAMN9"}]


def test_parse_ena_assembly_missing_column_is_empty():
    # No assembly_accession column -> [] (mirrors read_run's no-required-column case).
    assert parse_ena_assembly("foo\tbar\nx\ty\n") == []
    assert parse_ena_assembly("") == []


def test_parse_entrez_assembly_summary_pascalcase_keys():
    # NCBI sometimes returns PascalCase; the `or` fallback branch must pick them up.
    js = json.dumps({"result": {
        "uids": ["7"],
        "7": {"AssemblyAccession": "GCA_7.2", "BioSampleAccn": "SAMN7"},
    }})
    assert parse_entrez_assembly_summary(js) == [{"assembly": "GCA_7.2", "biosample": "SAMN7"}]


def test_parse_entrez_assembly_summary_bad_json():
    assert parse_entrez_assembly_summary("not json") == []


# --------------------------------------------------------------------------- #
# _reconcile_biosample (C4 disagreement -> None)
# --------------------------------------------------------------------------- #
def test_reconcile_agree():
    assert _reconcile_biosample("SAMN1", "SAMN1") == ("SAMN1", "entrez+ena")


def test_reconcile_disagree_is_none():
    assert _reconcile_biosample("SAMN1", "SAMN2") == (None, "disagreement")


def test_reconcile_single_source():
    assert _reconcile_biosample("SAMN1", None) == ("SAMN1", "entrez")
    assert _reconcile_biosample(None, "SAMEA9") == ("SAMEA9", "ena")


def test_reconcile_neither():
    assert _reconcile_biosample(None, None) == (None, "unresolved")


# --------------------------------------------------------------------------- #
# Resolver composition with a fake fetch + cache
# --------------------------------------------------------------------------- #
def _fake_fetch_factory(routes: dict[str, str]):
    """Return a fetch(url) that matches the first substring key present in the URL."""
    def fetch(url: str) -> str:
        for key, body in routes.items():
            if key in url:
                return body
        raise AssertionError(f"unexpected url: {url}")
    return fetch


def test_runs_for_project(tmp_path):
    routes = {"result=read_run": "run_accession\tsample_accession\nERR1\tSAMEA1\nERR2\tSAMEA2\n"}
    r = BioSampleResolver(cache_path=tmp_path / "cache.json", fetch=_fake_fetch_factory(routes))
    assert r.runs_for_project("PRJEB1") == [("ERR1", "SAMEA1"), ("ERR2", "SAMEA2")]


def test_biosample_to_assemblies_entrez_primary(tmp_path):
    routes = {
        "esearch.fcgi": json.dumps({"esearchresult": {"idlist": ["55"]}}),
        "esummary.fcgi": json.dumps({"result": {"uids": ["55"],
                                                 "55": {"assemblyaccession": "GCA_5.1", "biosampleaccn": "SAMN5"}}}),
    }
    r = BioSampleResolver(cache_path=tmp_path / "cache.json", fetch=_fake_fetch_factory(routes))
    assert r.biosample_to_assemblies("SAMN5") == ["GCA_5.1"]


def test_biosample_to_assemblies_ena_fallback(tmp_path):
    # Entrez returns no UIDs -> ENA fallback supplies the assembly.
    routes = {
        "esearch.fcgi": json.dumps({"esearchresult": {"idlist": []}}),
        "result=assembly": "assembly_accession\tsample_accession\nGCA_7.1\tSAMEA7\n",
    }
    r = BioSampleResolver(cache_path=tmp_path / "cache.json", fetch=_fake_fetch_factory(routes))
    assert r.biosample_to_assemblies("SAMEA7") == ["GCA_7.1"]


def test_biosample_to_assemblies_empty_is_valid(tmp_path):
    routes = {
        "esearch.fcgi": json.dumps({"esearchresult": {"idlist": []}}),
        "result=assembly": "assembly_accession\tsample_accession\n",  # header only
    }
    r = BioSampleResolver(cache_path=tmp_path / "cache.json", fetch=_fake_fetch_factory(routes))
    assert r.biosample_to_assemblies("SAMEA_noasm") == []  # reads-only, not an error


def test_cache_hit_avoids_refetch(tmp_path):
    calls = {"n": 0}

    def fetch(url):
        calls["n"] += 1
        if "esearch" in url:
            return json.dumps({"esearchresult": {"idlist": ["1"]}})
        return json.dumps({"result": {"uids": ["1"], "1": {"assemblyaccession": "GCA_1.1", "biosampleaccn": "SAMN1"}}})

    r = BioSampleResolver(cache_path=tmp_path / "cache.json", fetch=fetch)
    r.biosample_to_assemblies("SAMN1")
    n_after_first = calls["n"]
    r.biosample_to_assemblies("SAMN1")  # cache hit
    assert calls["n"] == n_after_first


def test_cache_persists_and_reloads(tmp_path):
    routes = {
        "esearch.fcgi": json.dumps({"esearchresult": {"idlist": ["1"]}}),
        "esummary.fcgi": json.dumps({"result": {"uids": ["1"],
                                                 "1": {"assemblyaccession": "GCA_1.1", "biosampleaccn": "SAMN1"}}}),
    }
    path = tmp_path / "cache.json"
    r1 = BioSampleResolver(cache_path=path, fetch=_fake_fetch_factory(routes))
    r1.biosample_to_assemblies("SAMN1")
    r1.save_cache()

    def boom(url):
        raise AssertionError("should not fetch — cache should serve this")

    r2 = BioSampleResolver(cache_path=path, fetch=boom)
    assert r2.biosample_to_assemblies("SAMN1") == ["GCA_1.1"]


def test_assembly_to_biosample_disagreement_none(tmp_path):
    routes = {
        "esearch.fcgi": json.dumps({"esearchresult": {"idlist": ["1"]}}),
        "esummary.fcgi": json.dumps({"result": {"uids": ["1"],
                                                 "1": {"assemblyaccession": "GCA_1.1", "biosampleaccn": "SAMN_E"}}}),
        "result=assembly": "assembly_accession\tsample_accession\nGCA_1.1\tSAMN_DIFFERENT\n",
    }
    r = BioSampleResolver(cache_path=tmp_path / "cache.json", fetch=_fake_fetch_factory(routes))
    assert r.assembly_to_biosample("GCA_1.1") is None  # disagreement -> unresolved


def test_assembly_to_biosample_agreement_returns_value(tmp_path):
    # Both Entrez + ENA agree -> the value resolves (the leakage-positive direction).
    routes = {
        "esearch.fcgi": json.dumps({"esearchresult": {"idlist": ["1"]}}),
        "esummary.fcgi": json.dumps({"result": {"uids": ["1"],
                                                 "1": {"assemblyaccession": "GCA_1.1", "biosampleaccn": "SAMN_SAME"}}}),
        "result=assembly": "assembly_accession\tsample_accession\nGCA_1.1\tSAMN_SAME\n",
    }
    r = BioSampleResolver(cache_path=tmp_path / "cache.json", fetch=_fake_fetch_factory(routes))
    assert r.assembly_to_biosample("GCA_1.1") == "SAMN_SAME"


def test_biosample_to_assemblies_entrez_raises_falls_back_to_ena(tmp_path):
    # Entrez fetch RAISES (network error) -> the except branch falls through to ENA.
    def fetch(url):
        if "eutils.ncbi" in url:
            raise RuntimeError("entrez down")
        if "result=assembly" in url:
            return "assembly_accession\tsample_accession\nGCA_8.1\tSAMEA8\n"
        raise AssertionError(url)

    r = BioSampleResolver(cache_path=tmp_path / "cache.json", fetch=fetch)
    assert r.biosample_to_assemblies("SAMEA8") == ["GCA_8.1"]


def test_assembly_to_biosample_both_sources_raise_is_none(tmp_path):
    # Both Entrez and ENA raise -> reconcile(None, None) -> None (unresolved), no crash.
    def boom(url):
        raise RuntimeError("network down")

    r = BioSampleResolver(cache_path=tmp_path / "cache.json", fetch=boom)
    assert r.assembly_to_biosample("GCA_x.1") is None


def test_assembly_to_biosample_caches_none(tmp_path):
    # A resolved-None result is cached so a second call does not re-fetch.
    calls = {"n": 0}

    def fetch(url):
        calls["n"] += 1
        if "esearch" in url:
            return json.dumps({"esearchresult": {"idlist": []}})
        return "assembly_accession\tsample_accession\n"   # ENA: header only -> None

    r = BioSampleResolver(cache_path=tmp_path / "cache.json", fetch=fetch)
    assert r.assembly_to_biosample("GCA_none.1") is None
    n_after_first = calls["n"]
    assert r.assembly_to_biosample("GCA_none.1") is None   # cache hit
    assert calls["n"] == n_after_first


def test_biosample_to_assemblies_multiple_deduped_sorted(tmp_path):
    # Two UIDs -> two assemblies (one duplicated) -> sorted unique list.
    routes = {
        "esearch.fcgi": json.dumps({"esearchresult": {"idlist": ["1", "2"]}}),
        "esummary.fcgi": json.dumps({"result": {"uids": ["1", "2"],
                                                 "1": {"assemblyaccession": "GCF_9.1", "biosampleaccn": "SAMN5"},
                                                 "2": {"assemblyaccession": "GCF_2.1", "biosampleaccn": "SAMN5"}}}),
    }
    r = BioSampleResolver(cache_path=tmp_path / "cache.json", fetch=_fake_fetch_factory(routes))
    assert r.biosample_to_assemblies("SAMN5") == ["GCF_2.1", "GCF_9.1"]
