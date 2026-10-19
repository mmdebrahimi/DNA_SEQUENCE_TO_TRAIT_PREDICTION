"""Offline regression suite for the CABBAGE G8 lineage sidecar (plan Steps 1 + 4).

No network, no Docker, no D: dependency. Real-data assertions skip when the
gitignored artifacts are absent, following the established pattern.

The collision guard is the load-bearing test here and it is written to be
NON-VACUOUS: the two cells carry deliberately DIFFERENT values, so a merge bug
cannot pass by coincidence, and one test simulates the bug to prove the guard
would actually catch it.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from dna_decode.eval import cabbage_lineage as cl
from dna_decode.eval.rejection_gates import MIN_EFFECTIVE_LINEAGES, g8_dedup_collapses_balance


# ---------------------------------------------------------------------------
# Namespace separation. Two of the four CABBAGE cells share a key with the
# provdisjoint lineage sidecar (Campylobacter x ciprofloxacin,
# Escherichia_coli_Shigella x gentamicin). That collision is measured, not
# hypothetical, which is why the sidecar is a separate file.
# ---------------------------------------------------------------------------

def test_sidecar_is_not_the_provdisjoint_file(tmp_path: Path) -> None:
    p = cl.sidecar_path(tmp_path)
    assert p.name == cl.SIDECAR_NAME
    assert p.name != cl.FORBIDDEN_SIDECAR_NAME


def test_the_forbidden_path_is_refused_by_name() -> None:
    with pytest.raises(cl.SidecarPathRefused):
        cl.assert_not_forbidden("wiki/" + cl.FORBIDDEN_SIDECAR_NAME)


def test_the_loader_constant_still_names_the_file_we_guard_against() -> None:
    """If the report card renames its sidecar, this guard must be re-pointed.

    The provdisjoint sidecar is read by an EXACT filename constant, not a glob,
    so path inequality is the real protection and this pins the other half.
    """
    src = Path("scripts/build_validation_report_card.py")
    if not src.exists():
        pytest.skip("report-card builder not present")
    needle = 'LINEAGE_SIDECAR = "' + cl.FORBIDDEN_SIDECAR_NAME + '"'
    assert needle in src.read_text(encoding="utf-8")


def test_collision_guard_is_non_vacuous(tmp_path: Path) -> None:
    """Same cell key, DIFFERENT values, two files, neither clobbers the other.

    The values are deliberately distinct so a merge bug cannot pass by coincidence.
    """
    key = cl.cabbage_cell_key("Campylobacter jejuni", "ciprofloxacin")
    cl.write_cells([{"key": key, "effective_lineage_N_R": 11}], root=tmp_path)

    neighbour = tmp_path / "wiki" / cl.FORBIDDEN_SIDECAR_NAME
    neighbour.write_text(
        json.dumps({"cells": [{"key": key, "effective_lineage_N_R": 99}]}), encoding="utf-8"
    )

    ours = cl.load_cells(tmp_path)
    theirs = json.loads(neighbour.read_text(encoding="utf-8"))["cells"]
    assert [c["effective_lineage_N_R"] for c in ours] == [11]
    assert [c["effective_lineage_N_R"] for c in theirs] == [99]


def test_the_collision_guard_would_actually_catch_the_bug() -> None:
    """Simulate the merge bug: one shared file means one number vanishes."""
    key = cl.cabbage_cell_key("Campylobacter jejuni", "ciprofloxacin")
    merged = cl.merge_cells(
        [{"key": key, "effective_lineage_N_R": 99}],
        [{"key": key, "effective_lineage_N_R": 11}],
    )
    assert len(merged) == 1
    assert merged[0]["effective_lineage_N_R"] == 11, (
        "a same-key merge silently overwrites, which is why the files are separate"
    )


# ---------------------------------------------------------------------------
# Merge-or-refuse. This sidecar accrues one cell per run, so a partial write is
# indistinguishable from a complete one unless dropping a cell raises.
# ---------------------------------------------------------------------------

def test_merge_is_additive_and_idempotent(tmp_path: Path) -> None:
    a = {"key": "Escherichia coli x gentamicin", "effective_lineage_N_R": 3}
    b = {"key": "Salmonella enterica x ceftriaxone", "effective_lineage_N_R": 7}
    cl.write_cells([a], root=tmp_path)
    cl.write_cells([b], root=tmp_path)
    cl.write_cells([b], root=tmp_path)
    assert {c["key"] for c in cl.load_cells(tmp_path)} == {a["key"], b["key"]}


def test_a_one_cell_run_cannot_replace_a_four_cell_artifact(tmp_path: Path) -> None:
    prior = [{"key": "cell-" + str(i), "n": i} for i in range(4)]
    cl.write_cells(prior, root=tmp_path)
    with pytest.raises(cl.CellDropRefused) as excinfo:
        cl.assert_no_cell_dropped(cl.load_cells(tmp_path), [{"key": "cell-0", "n": 0}])
    assert "cell-1" in str(excinfo.value)
    assert len(cl.load_cells(tmp_path)) == 4, "a refusal must not mutate disk"


# ---------------------------------------------------------------------------
# The bar is imported, never restated, and the three G8 states stay distinct.
# ---------------------------------------------------------------------------

def test_the_bar_is_imported_from_the_gate(tmp_path: Path) -> None:
    cl.write_cells([{"key": "k"}], root=tmp_path)
    doc = json.loads(cl.sidecar_path(tmp_path).read_text(encoding="utf-8"))
    assert doc["min_effective_lineages_bar"] == MIN_EFFECTIVE_LINEAGES
    assert doc["bar_source"].endswith("MIN_EFFECTIVE_LINEAGES")
    assert doc["schema"] == cl.SCHEMA


def test_partial_emits_no_lineage_number_so_g8_stays_insufficient() -> None:
    """A partial cohort leaves G8 INSUFFICIENT_DATA, never a number from a subset."""
    assert g8_dedup_collapses_balance({}).verdict == "insufficient_data"
    empty = cl.cluster_structure({}, {})
    assert empty["effective_lineage_N_R"] == 0
    assert empty["largest_cluster_fraction"] is None


def test_both_trip_and_pass_render() -> None:
    clusters = {"r1": 0, "r2": 1, "s1": 2, "s2": 3, "s3": 4}
    labels = {"r1": "R", "r2": "R", "s1": "S", "s2": "S", "s3": "S"}
    n, which = cl.min_effective_lineages(clusters, labels)
    assert (n, which) == (2, "R")
    assert g8_dedup_collapses_balance({"min_effective_lineages": n}).verdict == "trip"

    n2, _ = cl.min_effective_lineages(dict(clusters, r3=5), dict(labels, r3="R"))
    assert n2 == MIN_EFFECTIVE_LINEAGES
    assert g8_dedup_collapses_balance({"min_effective_lineages": n2}).verdict == "pass"


# ---------------------------------------------------------------------------
# Structure, not just counts. The TB arm reported 43 discordant clusters while
# those clusters hid 97% of the isolates, and the count alone looked healthy.
# ---------------------------------------------------------------------------

def test_a_dominant_clone_is_surfaced_behind_a_healthy_cluster_count() -> None:
    clusters = {"big" + str(i): 0 for i in range(77)}
    labels = {k: "R" for k in clusters}
    for i in range(23):
        clusters["s" + str(i)] = i + 1
        labels["s" + str(i)] = "S"
    st = cl.cluster_structure(clusters, labels)
    assert st["n_clusters"] == 24, "the cluster COUNT looks healthy"
    assert st["largest_cluster_fraction"] == pytest.approx(0.77)
    assert st["effective_lineage_N_R"] == 1, "77 isolates collapse to ONE lineage"


def test_discordant_clusters_are_excluded_from_both_classes_and_reported() -> None:
    clusters = {"a": 0, "b": 0, "c": 1, "d": 2}
    labels = {"a": "R", "b": "S", "c": "R", "d": "S"}
    st = cl.cluster_structure(clusters, labels)
    assert st["n_discordant_clusters"] == 1
    assert st["discordant_isolate_fraction"] == pytest.approx(0.5)
    assert st["effective_lineage_N_R"] == 1
    assert st["effective_lineage_N_S"] == 1


def test_a_weighted_point_without_its_interval_is_refused() -> None:
    with pytest.raises(cl.WeightedPointNotRenderable):
        cl.assert_weighted_renderable(cl.WeightedPoint(label="sens", value=0.5, ci=None))
    good = cl.weighted_point("sens", 3, 4)
    assert good.ci is not None
    assert good.value == pytest.approx(0.75)
    cl.assert_weighted_renderable(cl.weighted_point("empty", 0, 0))


# ---------------------------------------------------------------------------
# Real artifacts, skipped when absent.
# ---------------------------------------------------------------------------

def test_the_four_scored_cells_still_carry_honest_limits() -> None:
    cells = sorted(Path("wiki").glob("cabbage_scored_cell_*.json"))
    if not cells:
        pytest.skip("CABBAGE scored-cell artifacts not present")
    for f in cells:
        doc = json.loads(f.read_text(encoding="utf-8"))
        assert "honest_limits" in doc
