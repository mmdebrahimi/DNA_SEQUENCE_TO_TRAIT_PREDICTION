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
    parse_entrez_assembly_summary,
    parse_entrez_esearch,
)


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
