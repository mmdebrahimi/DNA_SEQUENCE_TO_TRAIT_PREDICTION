"""Offline tests for the finer Mash lineage collapse of the independent TB cohort.

Pure helpers only — no Docker, no D: cache, no network. The Mash invocation itself is the one
non-pure seam and is exercised on the real surface by the script's own smoke path.
"""
from __future__ import annotations

import json

import numpy as np
import pytest

from scripts.tb_independent_mash_collapse import (
    BARCODE_REFERENCE,
    SKETCH_SIZE,
    THRESHOLDS,
    cluster_structure,
    collapse_at,
    discordant_isolate_fraction,
    granularity_matched_rung,
    load_results,
    parse_phylip_lower_triangle,
    rung_is_degenerate,
)

# `mash triangle` relaxed-Phylip lower triangle, exactly as the container emits it
# (tab-separated, container-absolute paths as names).
TRIANGLE = "\n".join([
    "4",
    "/asm/A.fna",
    "/asm/B.fna\t0.001",
    "/asm/C.fna\t0.010\t0.011",
    "/asm/D.fna\t0.0005\t0.0012\t0.012",
]) + "\n"


def test_parse_triangle_names_are_stripped_basenames():
    names, _ = parse_phylip_lower_triangle(TRIANGLE)
    assert names == ["A", "B", "C", "D"]


def test_parse_triangle_matrix_is_symmetric_with_zero_diagonal():
    names, m = parse_phylip_lower_triangle(TRIANGLE)
    assert m.shape == (4, 4)
    assert (m.diagonal() == 0).all()
    assert np.allclose(m, m.T)
    i = names.index
    assert m[i("A"), i("B")] == pytest.approx(0.001)
    assert m[i("C"), i("B")] == pytest.approx(0.011)
    assert m[i("D"), i("A")] == pytest.approx(0.0005)


def test_parse_triangle_rejects_empty():
    with pytest.raises(ValueError):
        parse_phylip_lower_triangle("")


def test_parse_triangle_rejects_truncated_body():
    """A declared count larger than the rows present must raise, never silently under-cluster."""
    with pytest.raises(ValueError):
        parse_phylip_lower_triangle("3\n/asm/A.fna\n/asm/B.fna\t0.01\n")


def test_parse_triangle_tolerates_crlf():
    """Regression: the list-file CRLF bug. Parsing must not depend on line-ending flavor."""
    names, m = parse_phylip_lower_triangle(TRIANGLE.replace("\n", "\r\n"))
    assert names == ["A", "B", "C", "D"]
    assert np.allclose(m, m.T)


def test_sketch_size_is_fine_enough_for_monomorphic_tb():
    """TB pairwise Mash distances live near 1e-5..1e-3; the default s=1000 quantizes them to 0."""
    assert SKETCH_SIZE >= 10000


def test_sweep_spans_raw_limit_to_barcode_granularity():
    """The sweep must bracket BOTH degenerate ends, else the curve hides one of them.

    Grounded in the real cohort: at 1e-5 Mash yields ~2501 clusters over 2845 isolates (~raw), at 1e-3 it
    yields 43 (coarser than the barcode's ~110). A sweep that failed to bracket those would report a
    threshold-dependent number as if it were threshold-independent.
    """
    assert min(THRESHOLDS) <= 1e-5
    assert max(THRESHOLDS) >= 1e-3


def test_granularity_matched_rung_picks_closest_cluster_count_to_barcode():
    target = BARCODE_REFERENCE["rifampicin"]["n_clusters_total"]  # 110

    def rung(thr, n_r, n_s, disc):
        return {"threshold": thr, "n_clusters_total": n_r + n_s + disc,
                "drugs": {"rifampicin": {"n_clusters_R": n_r, "n_clusters_S": n_s, "n_discordant": disc}}}

    sweep = [rung(1e-5, 900, 900, 5), rung(3e-4, 40, 60, 8), rung(1e-3, 10, 20, 3)]
    picked = granularity_matched_rung(sweep, "rifampicin")
    assert picked["threshold"] == 3e-4  # 108 is nearest to 110
    assert target == 110


def _results(**per_strain):
    return {
        sid: {"strain_id": sid, "rif_label": lab, "rif_pred": pred,
              "inh_label": lab, "inh_pred": pred}
        for sid, (lab, pred) in per_strain.items()
    }


def test_collapse_at_collapses_a_clone_to_one_vote():
    """3 identical R isolates + 1 S: raw would count 3 R votes; collapsed counts 1."""
    names = ["A", "B", "C", "S1"]
    m = np.array([
        [0.0, 0.00001, 0.00001, 0.02],
        [0.00001, 0.0, 0.00001, 0.02],
        [0.00001, 0.00001, 0.0, 0.02],
        [0.02, 0.02, 0.02, 0.0],
    ])
    res = _results(A=("R", "R"), B=("R", "R"), C=("R", "R"), S1=("S", "S"))
    out = collapse_at(m, names, res, threshold=0.001)

    rif = out["drugs"]["rifampicin"]
    assert out["n_clusters_total"] == 2
    assert rif["tp"] == 1 and rif["tn"] == 1
    assert rif["n_clusters_R"] == 1 and rif["n_clusters_S"] == 1
    assert rif["sens"] == 1.0 and rif["spec"] == 1.0
    assert rif["effective_lineage_n_R"] == 1


def test_collapse_at_excludes_mixed_label_lineage_as_discordant():
    """A clone carrying both R and S labels is DISCORDANT — never majority-voted into sens/spec."""
    names = ["A", "B"]
    m = np.array([[0.0, 0.00001], [0.00001, 0.0]])
    res = _results(A=("R", "R"), B=("S", "R"))
    out = collapse_at(m, names, res, threshold=0.001)

    rif = out["drugs"]["rifampicin"]
    assert out["n_clusters_total"] == 1
    assert rif["n_discordant"] == 1
    assert rif["n_scored"] == 0
    assert rif["sens"] is None and rif["spec"] is None


def test_collapse_at_finer_threshold_splits_the_clone_back_apart():
    """Threshold -> 0 returns the clonality-INFLATED raw view (each isolate its own lineage)."""
    names = ["A", "B", "S1"]
    m = np.array([
        [0.0, 0.00005, 0.02],
        [0.00005, 0.0, 0.02],
        [0.02, 0.02, 0.0],
    ])
    res = _results(A=("R", "R"), B=("R", "R"), S1=("S", "S"))

    coarse = collapse_at(m, names, res, threshold=0.001)
    fine = collapse_at(m, names, res, threshold=0.00001)

    assert coarse["n_clusters_total"] == 2
    assert coarse["drugs"]["rifampicin"]["n_clusters_R"] == 1
    # finer: the two R isolates separate -> 2 R lineages == the raw isolate count
    assert fine["n_clusters_total"] == 3
    assert fine["drugs"]["rifampicin"]["n_clusters_R"] == 2


def test_collapse_at_ignores_isolates_without_an_R_S_label():
    names = ["A", "B"]
    m = np.array([[0.0, 0.02], [0.02, 0.0]])
    res = _results(A=("R", "R"), B=("U", "S"))  # B has an unusable label
    out = collapse_at(m, names, res, threshold=0.001)
    assert out["drugs"]["rifampicin"]["n_isolates"] == 1


def test_collapse_at_reports_wilson_ci_alongside_every_point():
    names = ["A", "S1"]
    m = np.array([[0.0, 0.02], [0.02, 0.0]])
    res = _results(A=("R", "R"), S1=("S", "S"))
    rif = collapse_at(m, names, res, threshold=0.001)["drugs"]["rifampicin"]
    lo, hi = rif["sens_ci95"]
    assert 0.0 <= lo <= 1.0 and 0.0 <= hi <= 1.0
    assert lo < hi  # n=1 -> a wide, non-degenerate interval


def test_cluster_structure_detects_a_blob():
    """One representative absorbing the cohort must be visible as a structure metric."""
    clusters = {f"s{i}": 0 for i in range(90)} | {f"t{i}": i + 1 for i in range(10)}
    st = cluster_structure(clusters, n_isolates=100)
    assert st["n_clusters"] == 11
    assert st["largest_cluster_fraction"] == pytest.approx(0.90)
    assert st["n_singletons"] == 10


def test_discordant_isolate_fraction_counts_isolates_not_clusters():
    """43 discordant CLUSTERS can hide 97% of ISOLATES — the isolate share is the honest metric."""
    clusters = {"a": 0, "b": 0, "c": 0, "d": 1}
    labels = {"a": "R", "b": "S", "c": "R", "d": "S"}  # cluster 0 is mixed -> 3 of 4 isolates excluded
    assert discordant_isolate_fraction(clusters, labels) == pytest.approx(0.75)


def test_rung_is_degenerate_on_blob():
    rung = {
        "n_clusters_total": 11,
        "structure": {"largest_cluster_fraction": 0.77},
        "drugs": {"rifampicin": {"discordant_isolate_fraction": 0.02}},
    }
    assert rung_is_degenerate(rung, n_isolates=100)


def test_rung_is_degenerate_on_mass_discordance():
    rung = {
        "n_clusters_total": 11,
        "structure": {"largest_cluster_fraction": 0.05},
        "drugs": {"rifampicin": {"discordant_isolate_fraction": 0.97}},
    }
    assert rung_is_degenerate(rung, n_isolates=100)


def test_rung_is_degenerate_when_barely_collapsed():
    """~1 cluster per isolate IS the clonality-inflated raw view, not a lineage collapse."""
    rung = {
        "n_clusters_total": 88,
        "structure": {"largest_cluster_fraction": 0.02},
        "drugs": {"rifampicin": {"discordant_isolate_fraction": 0.01}},
    }
    assert rung_is_degenerate(rung, n_isolates=100)


def test_rung_not_degenerate_when_balanced_and_collapsed():
    rung = {
        "n_clusters_total": 20,
        "structure": {"largest_cluster_fraction": 0.10},
        "drugs": {"rifampicin": {"discordant_isolate_fraction": 0.05}},
    }
    assert not rung_is_degenerate(rung, n_isolates=100)


def test_load_results_parses_jsonl(tmp_path):
    p = tmp_path / "results.jsonl"
    p.write_text(
        json.dumps({"strain_id": "X", "rif_label": "R", "inh_label": "S",
                    "rif_pred": "R", "inh_pred": "S"}) + "\n\n", encoding="utf-8")
    out = load_results(tmp_path)
    assert out["X"]["rif_label"] == "R"
