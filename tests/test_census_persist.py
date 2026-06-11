"""Pin census self-persistence (scripts/ncbi_pd_provenance_census.py M1 amendment).

The normalizer maps group->organism and refuses degraded (error/capped) rows; the upsert is idempotent and
never overwrites a prior good powering verdict with degraded data.
"""
from __future__ import annotations

import json

from scripts.ncbi_pd_provenance_census import census_result_to_sidecar_row, upsert_census_result


def test_normalizer_maps_group_to_organism():
    row = census_result_to_sidecar_row(
        {"group": "Klebsiella", "drug": "ceftriaxone", "other_R": 505, "other_S": 410, "powered": True},
        "2026-06-10", 20)
    assert row["organism"] == "Klebsiella" and row["drug"] == "ceftriaxone"
    assert row["other_R"] == 505 and row["powered"] is True and "group" not in row


def test_normalizer_refuses_error_and_capped():
    assert census_result_to_sidecar_row({"group": "X", "drug": "y", "error": "boom"}, "d", 20) is None
    assert census_result_to_sidecar_row(
        {"group": "X", "drug": "y", "other_R": 5, "other_S": 5, "capped": True}, "d", 20) is None
    assert census_result_to_sidecar_row({"group": "X", "drug": "y"}, "d", 20) is None  # missing counts


def test_upsert_idempotent_one_row_per_pair(tmp_path):
    p = tmp_path / "census.json"
    r1 = {"group": "Klebsiella", "drug": "gentamicin", "other_R": 100, "other_S": 100, "powered": True}
    assert upsert_census_result(r1, "2026-06-10", 20, path=p) is True
    # re-run same pair with new counts -> replace, not duplicate
    r2 = {"group": "Klebsiella", "drug": "gentamicin", "other_R": 317, "other_S": 339, "powered": True}
    assert upsert_census_result(r2, "2026-06-11", 20, path=p) is True
    doc = json.loads(p.read_text())
    rows = [r for r in doc["results"] if (r["organism"], r["drug"]) == ("Klebsiella", "gentamicin")]
    assert len(rows) == 1 and rows[0]["other_R"] == 317 and rows[0]["date"] == "2026-06-11"


def test_upsert_skips_degraded_keeps_prior(tmp_path):
    p = tmp_path / "census.json"
    upsert_census_result({"group": "Kleb", "drug": "ceft", "other_R": 9, "other_S": 9, "powered": False},
                         "2026-06-10", 20, path=p)
    before = p.read_text()
    # a capped run for the SAME pair must NOT overwrite the prior row
    assert upsert_census_result({"group": "Kleb", "drug": "ceft", "other_R": 1, "other_S": 1, "capped": True},
                                "2026-06-11", 20, path=p) is False
    assert p.read_text() == before  # untouched


def test_upsert_preserves_header(tmp_path):
    p = tmp_path / "census.json"
    upsert_census_result({"group": "A", "drug": "b", "other_R": 50, "other_S": 50, "powered": True},
                         "2026-06-10", 20, path=p)
    doc = json.loads(p.read_text())
    assert doc["_schema"] == "provdisjoint-census-results-v1" and doc["min_per_class"] == 20
    assert "ecosystem_excluded" in doc
