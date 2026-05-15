"""Tests for Step 10 — Mash/ANI phylogeny clustering.

Batched-call refactor (2026-05-15): `compute_mash_distances` now issues
ONE `mash sketch` + ONE `mash dist sketch.msh sketch.msh` instead of
N*(N-1)/2 pairwise `mash dist` calls. Tests mock the all-pairs output
shape, not the per-pair shape.
"""
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


# ---- _parse_mash_dist_line (single-line legacy helper) ----


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


# ---- compute_mash_distances (batched, mocked subprocess) ----


def _all_pairs_stdout(strain_distances: dict[tuple[str, str], float]) -> str:
    """Build a mash-dist all-pairs stdout string.

    Each entry maps (ref_name, query_name) -> distance. Names should match
    the staged "<strain_id>.fna" filename (with absolute path prefix is
    fine — the parser uses Path(...).name).
    """
    lines = []
    for (ref, query), d in strain_distances.items():
        lines.append(f"{ref}\t{query}\t{d}\t0.001\t100/1000")
    return "\n".join(lines) + "\n"


def _mock_subprocess_factory(dist_stdout: str):
    """Return a side-effect callable yielding empty stdout for `mash sketch`
    then `dist_stdout` for `mash dist`. Matches the batched 2-call pattern."""
    sketch_proc = MagicMock(spec=subprocess.CompletedProcess)
    sketch_proc.returncode = 0
    sketch_proc.stdout = ""
    sketch_proc.stderr = ""

    dist_proc = MagicMock(spec=subprocess.CompletedProcess)
    dist_proc.returncode = 0
    dist_proc.stdout = dist_stdout
    dist_proc.stderr = ""

    return [sketch_proc, dist_proc]


def test_compute_mash_distances_two_strains(tmp_path: Path):
    """Batched 2-call shape: 1 sketch + 1 dist all-pairs."""
    genomes = {"a": tmp_path / "a.fna", "b": tmp_path / "b.fna"}
    for p in genomes.values():
        p.write_text(">x\nAGCT\n")

    pairs_stdout = _all_pairs_stdout({
        ("a.fna", "a.fna"): 0.0,
        ("a.fna", "b.fna"): 0.01,
        ("b.fna", "a.fna"): 0.01,
        ("b.fna", "b.fna"): 0.0,
    })

    with patch(
        "dna_decode.eval.phylogeny.subprocess.run",
        side_effect=_mock_subprocess_factory(pairs_stdout),
    ) as m:
        dm = compute_mash_distances(genomes)

    assert m.call_count == 2  # batched: sketch + dist
    assert dm.strain_ids == ["a", "b"]
    assert dm.matrix.shape == (2, 2)
    assert dm.matrix[0, 1] == pytest.approx(0.01)
    assert dm.matrix[1, 0] == pytest.approx(0.01)
    assert dm.matrix[0, 0] == 0
    assert dm.matrix[1, 1] == 0


def test_compute_mash_distances_three_strains_symmetric_matrix(tmp_path: Path):
    genomes = {sid: tmp_path / f"{sid}.fna" for sid in ("a", "b", "c")}
    for p in genomes.values():
        p.write_text(">x\nAGCT\n")

    pairs_stdout = _all_pairs_stdout({
        ("a.fna", "a.fna"): 0.0,
        ("a.fna", "b.fna"): 0.02,
        ("a.fna", "c.fna"): 0.05,
        ("b.fna", "a.fna"): 0.02,
        ("b.fna", "b.fna"): 0.0,
        ("b.fna", "c.fna"): 0.04,
        ("c.fna", "a.fna"): 0.05,
        ("c.fna", "b.fna"): 0.04,
        ("c.fna", "c.fna"): 0.0,
    })

    with patch(
        "dna_decode.eval.phylogeny.subprocess.run",
        side_effect=_mock_subprocess_factory(pairs_stdout),
    ):
        dm = compute_mash_distances(genomes)

    assert dm.matrix[0, 1] == pytest.approx(0.02)
    assert dm.matrix[0, 2] == pytest.approx(0.05)
    assert dm.matrix[1, 2] == pytest.approx(0.04)
    # symmetric (ref-major output covers both halves)
    np.testing.assert_array_equal(dm.matrix, dm.matrix.T)
    # diagonal zero
    for i in range(3):
        assert dm.matrix[i, i] == 0


def test_compute_mash_distances_batched_call_count_independent_of_n(tmp_path: Path):
    """Regression pin: N strains -> 2 subprocess calls (NOT N*(N-1)/2)."""
    genomes = {sid: tmp_path / f"{sid}.fna" for sid in ("a", "b", "c", "d", "e")}
    for p in genomes.values():
        p.write_text(">x\nAGCT\n")

    pairs = {}
    for i in "abcde":
        for j in "abcde":
            pairs[(f"{i}.fna", f"{j}.fna")] = 0.0 if i == j else 0.03

    with patch(
        "dna_decode.eval.phylogeny.subprocess.run",
        side_effect=_mock_subprocess_factory(_all_pairs_stdout(pairs)),
    ) as m:
        compute_mash_distances(genomes)

    assert m.call_count == 2, "batched mash should be 1 sketch + 1 dist"


def test_compute_mash_distances_empty(tmp_path: Path):
    dm = compute_mash_distances({})
    assert dm.strain_ids == []
    assert dm.matrix.shape == (0, 0)


def test_compute_mash_distances_missing_binary_raises(tmp_path: Path):
    genomes = {"a": tmp_path / "a.fna", "b": tmp_path / "b.fna"}
    for p in genomes.values():
        p.write_text(">x\nAGCT\n")
    with patch(
        "dna_decode.eval.phylogeny.subprocess.run",
        side_effect=FileNotFoundError("mash not found"),
    ):
        with pytest.raises(MashError, match="not found on PATH"):
            compute_mash_distances(genomes)


def test_compute_mash_distances_nonzero_exit_raises(tmp_path: Path):
    genomes = {"a": tmp_path / "a.fna", "b": tmp_path / "b.fna"}
    for p in genomes.values():
        p.write_text(">x\nAGCT\n")
    bad = MagicMock(spec=subprocess.CompletedProcess)
    bad.returncode = 1
    bad.stdout = ""
    bad.stderr = "mash error: bad file"
    with patch("dna_decode.eval.phylogeny.subprocess.run", return_value=bad):
        with pytest.raises(MashError, match="exit 1"):
            compute_mash_distances(genomes)


def test_compute_mash_distances_unknown_strain_in_dist_output_raises(tmp_path: Path):
    """Defensive: if mash dist references a file not in the cohort, raise."""
    genomes = {"a": tmp_path / "a.fna", "b": tmp_path / "b.fna"}
    for p in genomes.values():
        p.write_text(">x\nAGCT\n")
    bogus_stdout = _all_pairs_stdout({
        ("a.fna", "b.fna"): 0.01,
        ("ghost.fna", "a.fna"): 0.5,  # not in cohort
    })
    with patch(
        "dna_decode.eval.phylogeny.subprocess.run",
        side_effect=_mock_subprocess_factory(bogus_stdout),
    ):
        with pytest.raises(MashError, match="unknown staged file"):
            compute_mash_distances(genomes)


def test_compute_mash_distances_docker_path_calls_docker_runner(tmp_path: Path):
    """use_docker=True routes through tools.docker_runner.run, not subprocess."""
    genomes = {"a": tmp_path / "a.fna", "b": tmp_path / "b.fna"}
    for p in genomes.values():
        p.write_text(">x\nAGCT\n")

    sketch_proc = MagicMock(spec=subprocess.CompletedProcess)
    sketch_proc.returncode = 0
    sketch_proc.stdout = ""
    sketch_proc.stderr = ""

    dist_proc = MagicMock(spec=subprocess.CompletedProcess)
    dist_proc.returncode = 0
    dist_proc.stdout = _all_pairs_stdout({
        ("a.fna", "a.fna"): 0.0,
        ("a.fna", "b.fna"): 0.05,
        ("b.fna", "a.fna"): 0.05,
        ("b.fna", "b.fna"): 0.0,
    })
    dist_proc.stderr = ""

    with patch(
        "tools.docker_runner.run",
        side_effect=[sketch_proc, dist_proc],
    ) as m:
        dm = compute_mash_distances(genomes, use_docker=True)

    assert m.call_count == 2
    # First positional arg of docker_runner.run is the image
    assert "mash" in m.call_args_list[0][0][0]
    assert dm.matrix[0, 1] == pytest.approx(0.05)


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
