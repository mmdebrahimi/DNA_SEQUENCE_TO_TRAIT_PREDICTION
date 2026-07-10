"""Offline tests for the prospective-cohort fetch script (pure parsers + funnel; no network)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.fetch_prospective_cohort import (  # noqa: E402
    build_cohort_rows, cells_by_taxon, parse_bvbrc_amr, parse_bvbrc_genomes,
    parse_datasets_report, phenotype_to_label, write_cohort_tsv,
)


def test_phenotype_to_label():
    assert phenotype_to_label("Resistant") == "R"
    assert phenotype_to_label("Susceptible") == "S"
    assert phenotype_to_label("Intermediate") is None
    assert phenotype_to_label("Non-susceptible") is None
    assert phenotype_to_label("") is None and phenotype_to_label(None) is None


def test_cells_by_taxon_covers_scored_grid():
    c = cells_by_taxon()
    assert c[562] >= {"ciprofloxacin", "ceftriaxone", "gentamicin", "tetracycline"}   # E. coli
    assert c[573] >= {"ciprofloxacin", "meropenem", "tetracycline"}                   # Klebsiella
    assert "ciprofloxacin" in c[194]                                                  # Campylobacter


def test_parse_bvbrc_genomes_keeps_only_public_with_assembly():
    rows = [
        {"genome_id": "562.1", "assembly_accession": "GCA_001.1", "public": True},
        {"genome_id": "562.2", "assembly_accession": "GCF_002.1", "public": True},
        {"genome_id": "562.3", "assembly_accession": "", "public": True},          # no assembly -> drop
        {"genome_id": "562.4", "assembly_accession": "GCA_004.1", "public": False}, # not public -> drop
    ]
    assert parse_bvbrc_genomes(rows) == {"562.1": "GCA_001.1", "562.2": "GCF_002.1"}


def test_parse_bvbrc_amr_excludes_computational_and_offtarget_and_ambiguous():
    rows = [
        {"genome_id": "562.1", "antibiotic": "ciprofloxacin", "resistant_phenotype": "Resistant",
         "laboratory_typing_method": "Disk diffusion"},                              # keep -> R
        {"genome_id": "562.2", "antibiotic": "ciprofloxacin", "resistant_phenotype": "Susceptible",
         "laboratory_typing_method": "Computational Prediction"},                    # circular -> drop
        {"genome_id": "562.3", "antibiotic": "colistin", "resistant_phenotype": "Resistant",
         "laboratory_typing_method": "MIC"},                                         # off-target drug -> drop
        {"genome_id": "562.4", "antibiotic": "ciprofloxacin", "resistant_phenotype": "Intermediate",
         "laboratory_typing_method": "MIC"},                                         # ambiguous -> drop
    ]
    got = parse_bvbrc_amr(rows, {"ciprofloxacin"})
    assert got == [{"genome_id": "562.1", "drug": "ciprofloxacin", "label": "R", "method": "Disk diffusion"}]


def test_parse_datasets_report():
    rep = {"reports": [{"assembly_info": {"release_date": "2026-07-01T00:00:00Z",
                                          "biosample": {"accession": "SAMN123"}, "assembly_status": "current"}}]}
    out = parse_datasets_report(rep)
    assert out == {"release_date": "2026-07-01", "biosample": "SAMN123", "status": "current"}
    assert parse_datasets_report({}) == {"release_date": "", "biosample": "", "status": ""}


def test_build_cohort_rows_funnel():
    amr = [
        {"genome_id": "g1", "drug": "ciprofloxacin", "label": "R"},   # post-lock -> eligible
        {"genome_id": "g2", "drug": "ceftriaxone", "label": "S"},     # pre-lock -> excluded
        {"genome_id": "g3", "drug": "tetracycline", "label": "R"},    # undatable -> excluded (fail-closed)
        {"genome_id": "g4", "drug": "gentamicin", "label": "R"},      # no GCA -> dropped before resolve
    ]
    gid_to_gca = {"g1": "GCA_1.1", "g2": "GCA_2.1", "g3": "GCA_3.1"}
    gca_release = {
        "GCA_1.1": {"release_date": "2026-07-01", "biosample": "SAMN1"},
        "GCA_2.1": {"release_date": "2026-05-01", "biosample": "SAMN2"},
        "GCA_3.1": {"release_date": "", "biosample": ""},
    }
    rows, stats = build_cohort_rows(amr, gid_to_gca, gca_release, "2026-06-13")
    assert rows == [{"biosample": "SAMN1", "first_public_date": "2026-07-01", "gca": "GCA_1.1",
                     "drug": "ciprofloxacin", "label": "R"}]
    assert stats["amr_records"] == 4 and stats["with_gca"] == 3 and stats["eligible"] == 1
    assert stats["excluded_pre_or_undatable"] == 1   # g2 (pre-lock); g3 fails at resolve (no date), not here


def test_write_cohort_tsv_format(tmp_path):
    rows = [{"biosample": "SAMN1", "first_public_date": "2026-07-01", "gca": "GCA_1.1",
             "drug": "ciprofloxacin", "label": "R"}]
    p = tmp_path / "prospective_cohort.tsv"
    write_cohort_tsv(rows, p)
    lines = p.read_text(encoding="utf-8").splitlines()
    assert lines[0] == "biosample\tfirst_public_date\tgca\tdrug\tlabel"
    assert lines[1] == "SAMN1\t2026-07-01\tGCA_1.1\tciprofloxacin\tR"


def test_write_cohort_tsv_header_only_when_empty(tmp_path):
    p = tmp_path / "empty.tsv"
    write_cohort_tsv([], p)
    assert p.read_text(encoding="utf-8").splitlines() == ["biosample\tfirst_public_date\tgca\tdrug\tlabel"]


if __name__ == "__main__":
    import pytest
    sys.exit(pytest.main([__file__, "-q"]))


# ============================================================================================
# Outage-vs-genuine-zero invariant (added 2026-07-10, after BV-BRC returned HTTP 200 + a 503
# error envelope and the sweep reported "0 recent genomes").
#
# THE LOAD-BEARING RULE: an accrual pipeline's honest signal is "0 eligible isolates have accrued
# yet". A DEAD SOURCE must never be able to emit that same signal, or an outage silently reads as
# a real zero and the prospective number is quietly never earned.
# ============================================================================================

import json  # noqa: E402
import pytest  # noqa: E402

from scripts import fetch_prospective_cohort as F  # noqa: E402


class _Resp:
    def __init__(self, body: bytes):
        self._body = body

    def read(self):
        return self._body

    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False


def test_get_json_raises_on_bvbrc_error_envelope(monkeypatch):
    """BV-BRC answers an outage with HTTP 200 + {"status":500,...}. Must raise, not reach the parsers."""
    envelope = json.dumps({"status": 500, "message": "... 503 Service Unavailable ..."}).encode()
    monkeypatch.setattr(F.urllib.request, "urlopen", lambda *a, **k: _Resp(envelope))
    monkeypatch.setattr(F.time, "sleep", lambda *_: None)
    with pytest.raises(F.SourceUnavailable):
        F._get_json("https://example.invalid", retries=2, backoff=0)


def test_get_json_retries_then_raises_on_timeout(monkeypatch):
    calls = []

    def _boom(*a, **k):
        calls.append(1)
        raise TimeoutError("read timed out")

    monkeypatch.setattr(F.urllib.request, "urlopen", _boom)
    monkeypatch.setattr(F.time, "sleep", lambda *_: None)
    with pytest.raises(F.SourceUnavailable):
        F._get_json("https://example.invalid", retries=3, backoff=0)
    assert len(calls) == 3  # retried, not one-shot


def test_get_json_returns_payload_on_success(monkeypatch):
    monkeypatch.setattr(F.urllib.request, "urlopen", lambda *a, **k: _Resp(b'[{"genome_id": "1.1"}]'))
    assert F._get_json("https://example.invalid") == [{"genome_id": "1.1"}]


def test_overall_status_ok_only_when_every_taxon_succeeded():
    assert F.overall_status({562: "ok", 573: "ok"}) == "OK"


def test_overall_status_partial_when_some_taxon_failed():
    assert F.overall_status({562: "ok", 573: "source_unavailable: 503"}) == "PARTIAL_SOURCE_UNAVAILABLE"


def test_overall_status_unavailable_when_all_failed():
    assert F.overall_status({562: "source_unavailable: 503"}) == "SOURCE_UNAVAILABLE"


def test_overall_status_unavailable_on_empty_status_map():
    """No taxon was even queried -> cannot claim a zero."""
    assert F.overall_status({}) == "SOURCE_UNAVAILABLE"


def test_fetch_live_records_source_failure_and_continues(monkeypatch):
    def _recent(taxon, lock_date, limit):
        if taxon == 562:
            raise F.SourceUnavailable("503 Service Unavailable")
        return [], "taxon_lineage_ids"

    monkeypatch.setattr(F, "cells_by_taxon", lambda: {562: {"ciprofloxacin"}, 194: {"ciprofloxacin"}})
    monkeypatch.setattr(F, "_bvbrc_recent_genomes", _recent)

    rows, _agg, taxon_status, _filt = F.fetch_live("2026-06-13", 10, 10, 0.0)
    assert rows == []
    assert taxon_status[194] == "ok"                            # a GENUINE zero
    assert taxon_status[562].startswith("source_unavailable")   # an OUTAGE — distinguishable
    assert F.overall_status(taxon_status) == "PARTIAL_SOURCE_UNAVAILABLE"


def test_fetch_live_genuine_zero_is_ok(monkeypatch):
    monkeypatch.setattr(F, "cells_by_taxon", lambda: {194: {"ciprofloxacin"}})
    monkeypatch.setattr(F, "_bvbrc_recent_genomes", lambda *a: ([], "taxon_lineage_ids"))
    rows, _agg, taxon_status, _f = F.fetch_live("2026-06-13", 10, 10, 0.0)
    assert rows == [] and F.overall_status(taxon_status) == "OK"


def _patch_fetch(monkeypatch, rows, taxon_status, filt=None):
    """Patch BOTH source entry points.

    `main`'s default source is `ncbi_pd`; patching only `fetch_live` (the bvbrc path) let these tests fall
    through to a REAL network stream of the PD metadata and hang.
    """
    def _stub(*a, **k):
        return (rows, {"eligible": len(rows)}, taxon_status, filt or {})

    monkeypatch.setattr(F, "fetch_live", _stub)
    monkeypatch.setattr(F, "fetch_live_ncbi_pd", _stub)


def test_main_refuses_cohort_tsv_when_source_unavailable(monkeypatch, tmp_path):
    _patch_fetch(monkeypatch, [], {562: "source_unavailable: 503"})
    rc = F.main(["--out-dir", str(tmp_path)])
    assert rc == 2
    assert not (tmp_path / "prospective_cohort.tsv").exists()   # <-- the invariant
    st = json.loads((tmp_path / "prospective_cohort_status.json").read_text())
    assert st["overall_status"] == "SOURCE_UNAVAILABLE"
    assert st["cohort_tsv_written"] is False


def test_main_partial_outage_exits_1_and_writes_no_tsv(monkeypatch, tmp_path):
    _patch_fetch(monkeypatch, [], {562: "ok", 573: "source_unavailable: timeout"})
    assert F.main(["--out-dir", str(tmp_path)]) == 1
    assert not (tmp_path / "prospective_cohort.tsv").exists()


def test_main_genuine_zero_writes_header_only_tsv_and_exits_0(monkeypatch, tmp_path):
    _patch_fetch(monkeypatch, [], {562: "ok", 573: "ok"})
    assert F.main(["--out-dir", str(tmp_path)]) == 0
    assert (tmp_path / "prospective_cohort.tsv").exists()       # ACCRUING is a real, writable state
    st = json.loads((tmp_path / "prospective_cohort_status.json").read_text())
    assert st["overall_status"] == "OK" and st["n_eligible_rows"] == 0


def test_status_artifact_records_taxon_filter_coverage_caveat(monkeypatch, tmp_path):
    _patch_fetch(monkeypatch, [], {562: "ok"}, filt={562: "taxon_id_only"})
    F.main(["--out-dir", str(tmp_path)])
    st = json.loads((tmp_path / "prospective_cohort_status.json").read_text())
    assert st["taxon_filter_used"]["562"] == "taxon_id_only"
    assert "UNDER-COUNT" in st["coverage_caveat"]


def test_recent_genomes_falls_back_to_taxon_id_when_lineage_field_rejected(monkeypatch):
    seen = []

    def _get(url, **kw):
        seen.append(url)
        if "taxon_lineage_ids" in url:
            raise F.SourceUnavailable("400 unsupported field")
        return [{"genome_id": "562.1", "assembly_accession": "GCA_1", "public": True}]

    monkeypatch.setattr(F, "_get_json", _get)
    rows, filt = F._bvbrc_recent_genomes(562, "2026-06-13", 10)
    assert filt == "taxon_id_only" and len(rows) == 1
    assert any("taxon_lineage_ids" in u for u in seen) and any("eq(taxon_id,562)" in u for u in seen)


def test_recent_genomes_prefers_lineage_filter_when_supported(monkeypatch):
    monkeypatch.setattr(F, "_get_json", lambda url, **kw: [{"genome_id": "562.1", "assembly_accession": "GCA_1"}])
    _rows, filt = F._bvbrc_recent_genomes(562, "2026-06-13", 10)
    assert filt == "taxon_lineage_ids"


def test_recent_genomes_raises_when_narrow_query_returns_non_list(monkeypatch):
    def _get(url, **kw):
        if "taxon_lineage_ids" in url:
            raise F.SourceUnavailable("outage")
        return {"status": 500}

    monkeypatch.setattr(F, "_get_json", _get)
    with pytest.raises(F.SourceUnavailable):
        F._bvbrc_recent_genomes(562, "2026-06-13", 10)


# ---------------------------------------------------------------- NCBI-PD source (added 2026-07-10)

def test_parse_ast_phenotypes_keeps_only_wanted_RS():
    got = F.parse_ast_phenotypes(
        "ciprofloxacin=R,gentamicin=S,tetracycline=I,amoxicillin=R", {"ciprofloxacin", "gentamicin", "tetracycline"})
    assert got == {"ciprofloxacin": "R", "gentamicin": "S"}   # I dropped, off-panel drug dropped


def test_parse_ast_phenotypes_real_pd_format_is_comma_separated_and_quoted():
    """The live PD field is comma-separated AND wrapped in literal double quotes."""
    field = '"ampicillin=ND,cefazolin=ND,ceftriaxone=R,ciprofloxacin=S,gentamicin=ND"'
    got = F.parse_ast_phenotypes(field, {"ceftriaxone", "ciprofloxacin", "gentamicin"})
    assert got == {"ceftriaxone": "R", "ciprofloxacin": "S"}   # ND dropped


def test_parse_ast_phenotypes_matches_drug_at_first_and_last_position():
    """Regression: the surrounding quotes ride on the FIRST and LAST tokens.

    A naive `f"{drug}=R" in ast.split(",")` never matches an end-position drug, silently under-counting.
    """
    assert F.parse_ast_phenotypes('"ciprofloxacin=R"', {"ciprofloxacin"}) == {"ciprofloxacin": "R"}
    assert F.parse_ast_phenotypes('"ciprofloxacin=R,gentamicin=S"',
                                  {"ciprofloxacin", "gentamicin"}) == {"ciprofloxacin": "R", "gentamicin": "S"}


def test_parse_ast_phenotypes_drops_nd_and_null():
    assert F.parse_ast_phenotypes('"ciprofloxacin=ND"', {"ciprofloxacin"}) == {}
    assert F.parse_ast_phenotypes("NULL", {"ciprofloxacin"}) == {}


def test_parse_ast_phenotypes_empty_and_malformed():
    assert F.parse_ast_phenotypes(None, {"ciprofloxacin"}) == {}
    assert F.parse_ast_phenotypes("garbage,=,cipro", {"ciprofloxacin"}) == {}


def test_earliest_public_date_takes_the_earliest_and_is_fail_closed():
    assert F.earliest_public_date("2026-07-01", "2026-06-20") == "2026-06-20"
    assert F.earliest_public_date(None, "2026-06-20") == "2026-06-20"
    assert F.earliest_public_date(None, None) is None      # -> INELIGIBLE downstream
    assert F.earliest_public_date("", "2026") is None      # year-only is not a usable date


def test_earliest_public_date_feeds_fail_closed_eligibility():
    from dna_decode.eval.prospective_lock import is_prospective_eligible
    assert not is_prospective_eligible(F.earliest_public_date(None, None), "2026-06-13").eligible


def test_cells_by_group_covers_the_pd_mappable_scored_grid():
    groups = F.cells_by_group()
    assert set(groups) == set(F.PD_ORG_GROUPS)
    assert "meropenem" in groups["Klebsiella"]
    assert "ciprofloxacin" in groups["Escherichia_coli_Shigella"]


def test_fetch_live_ncbi_pd_records_source_failure(monkeypatch):
    monkeypatch.setattr(F, "cells_by_group", lambda: {"Klebsiella": {"ciprofloxacin"}})

    def _boom(group, timeout=600):
        raise F.SourceUnavailable("ftp 500")
        yield  # pragma: no cover

    monkeypatch.setattr(F, "_pd_metadata_lines", _boom)
    rows, _agg, status, _f = F.fetch_live_ncbi_pd("2026-06-13", 0.0)
    assert rows == []
    assert status["Klebsiella"].startswith("source_unavailable")
    assert F.overall_status(status) == "SOURCE_UNAVAILABLE"


def test_fetch_live_ncbi_pd_end_to_end_prefilter(monkeypatch):
    """Pre-lock SRA release is dropped by the cheap filter; post-lock survives to the datasets check."""
    header = "asm_acc\tAST_phenotypes\tsra_release_date\tbiosample_acc"
    rows_tsv = [
        header,
        "GCA_OLD\tciprofloxacin=R\t2026-01-01\tSAMN1",   # pre-lock -> dropped
        "GCA_NEW\tciprofloxacin=R\t2026-07-01\tSAMN2",   # post-lock -> kept
        "\tciprofloxacin=R\t2026-07-01\tSAMN3",          # no assembly -> dropped
        "GCA_NOAST\t\t2026-07-01\tSAMN4",                # no AST -> dropped
    ]
    monkeypatch.setattr(F, "cells_by_group", lambda: {"Klebsiella": {"ciprofloxacin"}})
    monkeypatch.setattr(F, "_pd_metadata_lines", lambda g, timeout=600: iter(rows_tsv))
    monkeypatch.setattr(F, "_ncbi_release", lambda acc: {"release_date": "2026-07-02", "biosample": "SAMN2"})
    monkeypatch.setattr(F.time, "sleep", lambda *_: None)

    rows, agg, status, filt = F.fetch_live_ncbi_pd("2026-06-13", 0.0)
    assert status["Klebsiella"] == "ok" and filt["Klebsiella"] == "ncbi_pd_metadata"
    assert agg["with_ast_and_asm"] == 2 and agg["prefilter_post_lock"] == 1
    assert len(rows) == 1
    assert rows[0]["gca"] == "GCA_NEW" and rows[0]["label"] == "R"
    assert rows[0]["first_public_date"] == "2026-07-01"   # earliest of (sra 07-01, asm 07-02)


def test_fetch_live_ncbi_pd_excludes_when_assembly_predates_lock(monkeypatch):
    """SRA says post-lock but the assembly was public BEFORE the lock -> earliest date wins -> excluded."""
    rows_tsv = ["asm_acc\tAST_phenotypes\tsra_release_date\tbiosample_acc",
                "GCA_X\tciprofloxacin=R\t2026-07-01\tSAMN9"]
    monkeypatch.setattr(F, "cells_by_group", lambda: {"Klebsiella": {"ciprofloxacin"}})
    monkeypatch.setattr(F, "_pd_metadata_lines", lambda g, timeout=600: iter(rows_tsv))
    monkeypatch.setattr(F, "_ncbi_release", lambda acc: {"release_date": "2026-05-01", "biosample": "SAMN9"})
    monkeypatch.setattr(F.time, "sleep", lambda *_: None)

    rows, agg, _s, _f = F.fetch_live_ncbi_pd("2026-06-13", 0.0)
    assert rows == [] and agg["excluded_pre_or_undatable"] == 1


def test_main_row_cap_cannot_claim_accruing(monkeypatch, tmp_path):
    monkeypatch.setattr(F, "fetch_live_ncbi_pd", lambda *a, **k: ([], {"eligible": 0}, {"Klebsiella": "ok"}, {}))
    rc = F.main(["--out-dir", str(tmp_path), "--source", "ncbi_pd", "--row-cap", "10"])
    assert rc == 1
    st = json.loads((tmp_path / "prospective_cohort_status.json").read_text())
    assert st["overall_status"] == "TRUNCATED_SMOKE_RUN"
    assert not (tmp_path / "prospective_cohort.tsv").exists()
