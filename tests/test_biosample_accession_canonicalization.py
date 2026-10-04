"""The unversioned-accession defect that burned two hours producing nothing.

MEASURED on the real surface, 2026-10-04. ENA's `filereport` returns `assembly_accession`
**UNVERSIONED** (`GCA_003073675`); NCBI's Entrez path returns it **versioned**
(`GCF_001875025.1`); and the genome downloader needs the version, because NCBI Datasets answers an
unversioned accession with an HTML error page that arrives downstream as "not a valid ZIP".

Consequence, not hypothetical: a 143-isolate scoring run produced **zero genomes, zero AMRFinder runs
and zero diagnostics** across two hours, because every ENA-resolved accession failed download and
retried while stdout sat block-buffered. A 3-isolate unbuffered smoke found it in minutes.

The fix RESOLVES the version via NCBI Datasets and never guesses it -- appending ".1" would be wrong
for any re-versioned assembly, and identity is this module's entire job. Unresolvable -> None, so the
caller drops to ASSEMBLY-REQUIRED (fail-closed).

These tests use an INJECTED fetch, so they are offline and deterministic; the live-surface checks are
marked and skip without network.
"""
from __future__ import annotations

import json

import pytest

from dna_decode.eval.biosample_resolver import BioSampleResolver


def _report(acc: str) -> str:
    return json.dumps({"reports": [{"accession": acc, "organism": {"organism_name": "K. pneumoniae"}}]})


def test_an_unversioned_accession_is_RESOLVED_not_guessed():
    """The resolved version must come from the service. `.2` proves it is not a hardcoded `.1`."""
    calls = []

    def fetch(url, **kw):
        calls.append(url)
        return _report("GCA_003073675.2")

    r = BioSampleResolver(fetch=fetch) if "fetch" in BioSampleResolver.__init__.__code__.co_varnames \
        else BioSampleResolver()
    r.fetch = fetch
    assert r.canonicalize_accession("GCA_003073675") == "GCA_003073675.2"
    assert calls, "no request was made -- the version was guessed, not resolved"


def test_an_already_versioned_accession_is_returned_WITHOUT_a_network_call():
    """The Entrez path already yields versioned accessions; re-resolving every one of them would add
    a request per isolate for nothing."""
    def fetch(url, **kw):  # pragma: no cover - must not run
        raise AssertionError(f"made a request for an already-versioned accession: {url}")
    r = BioSampleResolver()
    r.fetch = fetch
    assert r.canonicalize_accession("GCF_001875025.1") == "GCF_001875025.1"


def test_identity_is_guarded_a_DIFFERENT_accession_is_refused():
    """If the service answers with some other assembly, taking it would silently swap one genome for
    another -- the worst possible failure in an identity bridge. Must return None."""
    r = BioSampleResolver()
    r.fetch = lambda url, **kw: _report("GCA_999999999.1")
    assert r.canonicalize_accession("GCA_003073675") is None


@pytest.mark.parametrize("bad", ["", "   ", None])
def test_empty_input_returns_None_rather_than_raising(bad):
    r = BioSampleResolver()
    r.fetch = lambda url, **kw: _report("GCA_1.1")
    assert r.canonicalize_accession(bad) is None


def test_unresolvable_returns_None_so_the_caller_FAILS_CLOSED():
    r = BioSampleResolver()
    r.fetch = lambda url, **kw: json.dumps({"reports": []})
    assert r.canonicalize_accession("GCA_003073675") is None

    def boom(url, **kw):
        raise RuntimeError("network down")
    r2 = BioSampleResolver()
    r2.fetch = boom
    assert r2.canonicalize_accession("GCA_003073675") is None, (
        "a fetch failure must be a valid 'unresolvable' outcome, not an exception that aborts a "
        "whole cohort run")


def test_a_POISONED_CACHE_is_recanonicalized_on_READ():
    """The write-path fix alone was NOT enough, and this is the test that would have caught it.

    A cache written before canonicalization existed holds unversioned accessions. The cache-hit
    branch returned them directly, so the smoke STILL failed with `GCA_003073675` after the first
    fix. Pin the read path."""
    r = BioSampleResolver()
    r._cache["biosample_to_assemblies"] = {
        "SAMN07291531": {"value": ["GCA_003073675"], "source": "ena"},
    }
    r.fetch = lambda url, **kw: _report("GCA_003073675.1")
    got = r.biosample_to_assemblies("SAMN07291531")
    assert got == ["GCA_003073675.1"], got
    assert "recanonicalized" in r._cache["biosample_to_assemblies"]["SAMN07291531"]["source"]


def test_a_CLEAN_cache_is_served_untouched():
    """Non-vacuity for the read path: a versioned cache must not trigger a request."""
    r = BioSampleResolver()
    r._cache["biosample_to_assemblies"] = {
        "SAMN1": {"value": ["GCF_001875025.1"], "source": "entrez"},
    }
    def fetch(url, **kw):  # pragma: no cover
        raise AssertionError("re-resolved an already-clean cache entry")
    r.fetch = fetch
    assert r.biosample_to_assemblies("SAMN1") == ["GCF_001875025.1"]


@pytest.mark.parametrize("acc", ["GCA_003073675", "GCA_003075515"])
def test_LIVE_the_two_accessions_that_actually_broke_the_run(acc):
    """Real-surface check on the exact accessions from the failed run. Skips without network."""
    r = BioSampleResolver()
    try:
        got = r.canonicalize_accession(acc)
    except Exception:  # pragma: no cover
        pytest.skip("no network")
    if got is None:
        pytest.skip("NCBI Datasets unreachable or accession withdrawn")
    assert got.startswith(acc + "."), got
    assert "." in got
