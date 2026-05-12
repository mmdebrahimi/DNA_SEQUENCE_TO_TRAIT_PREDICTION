"""Tests for Step 10 — Mash/ANI phylogeny clustering."""
from __future__ import annotations

import subprocess
from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pytest

from dna_decode.eval.phylogeny import (
    DEFAULT_ANI_THRESHOLD,
    DistanceMatrix,
    MashError,
    _parse_mash_dist_line,
    cluster_by_ani,
    compute_mash_distances,
)


# ---- DistanceMatrix dataclass ----


def test_distance_matrix_shape_validation():
    with pytest.raises(ValueError, match="shape"):
        DistanceMatrix(strain_ids=["a", "b", "c"], matrix=np.zeros((2, 2)))


def test_distance_matrix_empty_ok():
    dm = DistanceMatrix(strain_ids=[], matrix=np.zeros((0, 0)))
    assert dm.matrix.shape == (0, 0)


# ---- _parse_mash_dist_line ----


def test_parse_mash_dist_line_basic():
    line = "ref.fna\tquery.fna\t0.05\t0.001\t100/1000\n"
    assert _parse_mash_dist_line(line) == 0.05


def test_parse_mash_dist_line_empty_raises():
    with pytest.raises(MashError, match="empty"):
        _parse_mash_dist_line("")


def test_parse_mash_dist_line_few_fields_raises():
    with pytest.raises(MashError, match="unexpected"):
        _parse_mash_dist_line("only\ttwo\n")


def test_parse_mash_dist_line_non_numeric_raises():
    with pytest.raises(MashError, match="non-numeric"):
        _parse_mash_dist_line("ref\tquery\tNOT_A_FLOAT\t0.5\t10/100\n")


# ---- compute_mash_distances (mocked subprocess) ----


def _mock_mash_result(distance: float) -> MagicMock:
    """Build a CompletedProcess-like with one tab-formatted dist line."""
    r = MagicMock(spec=subprocess.CompletedProcess)
    r.returncode = 0
    r.stdout = f"a\tb\t{distance}\t0.001\t100/1000\n"
    r.stderr = ""
    return r


def test_compute_mash_distances_two_strains(tmp_path: Path):
    """Mocked mash dist returns 0.01 for the only pair."""
    genomes = {"a": tmp_path / "a.fna", "b": tmp_path / "b.fna"}
    for p in genomes.values():
        p.write_text(">x\nAGCT\n")

    with patch("subprocess.run", return_value=_mock_mash_result(0.01)):
        dm = compute_mash_distances(genomes)

    assert dm.strain_ids == ["a", "b"]
    assert dm.matrix.shape == (2, 2)
    assert dm.matrix[0, 1] == pytest.approx(0.01)
    assert dm.matrix[1, 0] == pytest.approx(0.01)  # symmetric
    assert dm.matrix[0, 0] == 0  # diagonal
    assert dm.matrix[1, 1] == 0


def test_compute_mash_distances_empty(tmp_path: Path):
    dm = compute_mash_distances({})
    assert dm.strain_ids == []
    assert dm.matrix.shape == (0, 0)


def test_compute_mash_distances_missing_binary_raises():
    with patch("subprocess.run", side_effect=FileNotFoundError("mash not found")):
        with pytest.raises(MashError, match="not found on PATH"):
            compute_mash_distances({"a": Path("/tmp/a.fna"), "b": Path("/tmp/b.fna")})


def test_compute_mash_distances_nonzero_exit_raises():
    bad = MagicMock(spec=subprocess.CompletedProcess)
    bad.returncode = 1
    bad.stdout = ""
    bad.stderr = "mash error: bad file"
    with patch("subprocess.run", return_value=bad):
        with pytest.raises(MashError, match="exit 1"):
            compute_mash_distances({"a": Path("/tmp/a.fna"), "b": Path("/tmp/b.fna")})


# ---- cluster_by_ani ----


def test_cluster_by_ani_two_close_pairs_into_one_cluster():
    """Distance below threshold → strains land in the same cluster."""
    dm = DistanceMatrix(
        strain_ids=["a", "b", "c"],
        matrix=np.array(
            [[0.0, 0.005, 0.5], [0.005, 0.0, 0.5], [0.5, 0.5, 0.0]], dtype=np.float32
        ),
    )
    clusters = cluster_by_ani(dm, threshold=0.02)
    # a + b cluster together (distance 0.005 < 0.02); c alone
    assert clusters["a"] == clusters["b"]
    assert clusters["c"] != clusters["a"]


def test_cluster_by_ani_all_far_apart_separate_clusters():
    dm = DistanceMatrix(
        strain_ids=["a", "b", "c"],
        matrix=np.array(
            [[0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]], dtype=np.float32
        ),
    )
    clusters = cluster_by_ani(dm, threshold=0.02)
    assert len(set(clusters.values())) == 3


def test_cluster_by_ani_all_close_one_cluster():
    dm = DistanceMatrix(
        strain_ids=["a", "b", "c"],
        matrix=np.array(
            [[0.0, 0.005, 0.005], [0.005, 0.0, 0.005], [0.005, 0.005, 0.0]],
            dtype=np.float32,
        ),
    )
    clusters = cluster_by_ani(dm, threshold=0.02)
    assert len(set(clusters.values())) == 1


def test_cluster_by_ani_empty():
    dm = DistanceMatrix(strain_ids=[], matrix=np.zeros((0, 0)))
    assert cluster_by_ani(dm) == {}


def test_cluster_compact_ids_zero_indexed():
    dm = DistanceMatrix(
        strain_ids=["a", "b"], matrix=np.array([[0.0, 0.5], [0.5, 0.0]], dtype=np.float32)
    )
    clusters = cluster_by_ani(dm, threshold=0.02)
    assert set(clusters.values()) == {0, 1}
