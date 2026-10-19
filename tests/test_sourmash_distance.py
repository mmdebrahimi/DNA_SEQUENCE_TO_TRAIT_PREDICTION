"""Tests for the sourmash FracMinHash distance producer.

Collectable WITHOUT sourmash installed -- the module must not import it at module
level. The real-sourmash cases skip via `importlib.util.find_spec`.

Synthetic fixtures use ~20 kb sequences at `scaled=100`. That is deliberate: at the
default `scaled=1000` a 20 kb sequence sketches to ~20 hashes, sourmash sets
`size_is_inaccurate=True`, and `.ani` comes back `None` -- so the fixture would trip
`SketchTooSmall` rather than exercise the distance path. Measured boundary: ~207
hashes is already enough for a float ANI.
"""
from __future__ import annotations

import importlib.util
import random
from pathlib import Path

import numpy as np
import pytest

from dna_decode.eval.clonality import greedy_representative_clusters_from_matrix
from dna_decode.eval.phylogeny import DistanceMatrix
from dna_decode.eval.sourmash_distance import (
    EmptySketch,
    SketchTooSmall,
    SourmashUnavailable,
    ani_to_distance,
    compute_sourmash_distances,
)

REPO_ROOT = Path(__file__).resolve().parent.parent

HAVE_SOURMASH = importlib.util.find_spec("sourmash") is not None
needs_sourmash = pytest.mark.skipif(
    not HAVE_SOURMASH, reason="sourmash not installed (optional `lineage` extra)"
)


# --------------------------------------------------------------------------- #
# ani_to_distance -- pure, always runnable
# --------------------------------------------------------------------------- #

def test_ani_to_distance_endpoints():
    assert ani_to_distance(1.0) == 0.0
    assert ani_to_distance(0.0) == 1.0
    assert ani_to_distance(0.999) == pytest.approx(0.001)
    assert ani_to_distance(0.995) == pytest.approx(0.005)


def test_ani_to_distance_is_monotone_decreasing():
    anis = [0.0, 0.25, 0.5, 0.9, 0.99, 0.999, 1.0]
    dists = [ani_to_distance(a) for a in anis]
    assert dists == sorted(dists, reverse=True)
    assert len(set(dists)) == len(dists)


@pytest.mark.parametrize("bad", [-0.001, 1.001, 2.0, -1.0])
def test_ani_to_distance_rejects_out_of_range(bad):
    with pytest.raises(ValueError):
        ani_to_distance(bad)


@pytest.mark.parametrize("bad", [None, "0.99", object()])
def test_ani_to_distance_rejects_non_numbers(bad):
    """`None` is the real hazard: it is what sourmash returns on a too-small sketch."""
    with pytest.raises(ValueError):
        ani_to_distance(bad)


# --------------------------------------------------------------------------- #
# availability + argument validation
# --------------------------------------------------------------------------- #

def test_module_does_not_import_sourmash_at_module_level():
    """The test file above imported the module; that must not have needed sourmash."""
    src = Path(
        "dna_decode/eval/sourmash_distance.py"
    ).read_text(encoding="utf-8")
    body = "\n".join(
        line for line in src.splitlines()
        if not line.startswith((" ", "\t", "#"))
    )
    assert "import sourmash" not in body
    assert "from sourmash" not in body


def test_sourmash_unavailable_is_named_not_bare(monkeypatch, tmp_path):
    """A missing optional extra must raise SourmashUnavailable, not ImportError."""
    import dna_decode.eval.sourmash_distance as mod

    def boom():
        raise SourmashUnavailable("simulated absence")

    monkeypatch.setattr(mod, "_minhash_cls", boom)
    fasta = tmp_path / "a.fna"
    fasta.write_text(">a\nACGTACGTAC\n", encoding="utf-8")
    with pytest.raises(SourmashUnavailable):
        mod.compute_sourmash_distances({"a": fasta})


def test_sourmash_unavailable_subclasses_runtime_error_not_import_error():
    assert issubclass(SourmashUnavailable, RuntimeError)
    assert not issubclass(SourmashUnavailable, ImportError)


def test_empty_input_returns_empty_matrix():
    dm = compute_sourmash_distances({})
    assert isinstance(dm, DistanceMatrix)
    assert dm.strain_ids == []
    assert dm.matrix.shape == (0, 0)


@pytest.mark.parametrize("kwargs", [{"ksize": 0}, {"ksize": -1}, {"scaled": 0}, {"scaled": -5}])
def test_rejects_nonpositive_parameters(tmp_path, kwargs):
    fasta = tmp_path / "a.fna"
    fasta.write_text(">a\nACGT\n", encoding="utf-8")
    with pytest.raises(ValueError):
        compute_sourmash_distances({"a": fasta}, **kwargs)


# --------------------------------------------------------------------------- #
# real sourmash
# --------------------------------------------------------------------------- #

def _write_genomes(tmp_path: Path) -> dict[str, Path]:
    """Two near-identical genomes plus one unrelated, ~20 kb each."""
    rng = random.Random(1234)
    base = "".join(rng.choice("ACGT") for _ in range(20_000))
    # `near` differs from `base` in its last 2% only.
    near = base[:19_600] + "".join(rng.choice("ACGT") for _ in range(400))
    far = "".join(rng.choice("ACGT") for _ in range(20_000))

    paths = {}
    for name, seq in (("base", base), ("near", near), ("far", far)):
        p = tmp_path / f"{name}.fna"
        p.write_text(f">{name}_contig1\n{seq}\n", encoding="utf-8")
        paths[name] = p
    return paths


@needs_sourmash
def test_returns_wellformed_distance_matrix(tmp_path):
    dm = compute_sourmash_distances(_write_genomes(tmp_path), ksize=21, scaled=100)
    assert isinstance(dm, DistanceMatrix)
    assert sorted(dm.strain_ids) == ["base", "far", "near"]
    assert dm.matrix.shape == (3, 3)
    assert np.allclose(np.diag(dm.matrix), 0.0), "diagonal must be zero"
    assert np.allclose(dm.matrix, dm.matrix.T), "matrix must be symmetric"
    assert np.isfinite(dm.matrix).all(), "no NaN/inf may reach the matrix"
    assert (dm.matrix >= 0.0).all() and (dm.matrix <= 1.0).all()


@needs_sourmash
def test_near_identical_pair_is_closer_than_the_unrelated_one(tmp_path):
    dm = compute_sourmash_distances(_write_genomes(tmp_path), ksize=21, scaled=100)
    idx = {s: i for i, s in enumerate(dm.strain_ids)}
    d_near = dm.matrix[idx["base"], idx["near"]]
    d_far = dm.matrix[idx["base"], idx["far"]]
    assert d_near < d_far, f"near={d_near} should be < far={d_far}"


@needs_sourmash
def test_seam_feeds_greedy_representative_clustering(tmp_path):
    """THE seam test: the producer's output must drop straight into the clusterer.

    This is the assertion that proves compatibility with the existing clonality
    path -- the whole reason the producer returns `phylogeny.DistanceMatrix` rather
    than a new type of its own.
    """
    genomes = _write_genomes(tmp_path)
    dm = compute_sourmash_distances(genomes, ksize=21, scaled=100)
    clusters = greedy_representative_clusters_from_matrix(dm, threshold=0.005)
    assert isinstance(clusters, dict)
    assert set(clusters) == set(genomes), "every input id must be assigned"
    assert all(isinstance(v, int) for v in clusters.values())
    # A permissive threshold must put everything in one cluster; a strict one must not.
    wide = greedy_representative_clusters_from_matrix(dm, threshold=1.0)
    assert len(set(wide.values())) == 1


@needs_sourmash
def test_sketch_too_small_is_raised_not_coerced(tmp_path):
    """A 20 kb genome at scaled=1000 is exactly sourmash's unreliable region.

    The producer must refuse rather than let `.ani is None` become a NaN distance.
    """
    genomes = _write_genomes(tmp_path)
    with pytest.raises(SketchTooSmall):
        compute_sourmash_distances(genomes, ksize=21, scaled=1000)


@needs_sourmash
def test_empty_fasta_raises_empty_sketch(tmp_path):
    p = tmp_path / "empty.fna"
    p.write_text(">only_a_header\n", encoding="utf-8")
    with pytest.raises(EmptySketch):
        compute_sourmash_distances({"empty": p}, ksize=21, scaled=100)


@needs_sourmash
def test_non_acgt_characters_do_not_break_the_sketch(tmp_path):
    """Real assemblies carry N runs and IUPAC codes; `force=True` must cover them.

    OPERATIONAL CAVEAT, measured while writing this test: ambiguity codes do not
    merely survive, they THIN the sketch sharply, because every k-mer overlapping a
    non-ACGT base is dropped. A 20.5 kb sequence that is ~half IUPAC/N sketched to
    89 hashes where a clean 20 kb sequence gave 195 -- enough to tip sourmash into
    `size_is_inaccurate` and trip `SketchTooSmall` at the same `scaled`. So a
    heavily-N draft assembly needs a DENSER sketch than its length alone suggests.
    `scaled=10` here is what makes the fixture exercise the distance path rather
    than the refusal path.
    """
    rng = random.Random(7)
    seq = "".join(rng.choice("ACGT") for _ in range(10_000))
    seq += "N" * 500
    seq += "".join(rng.choice("ACGTRYKMSWN") for _ in range(10_000))
    a = tmp_path / "ambig.fna"
    a.write_text(f">c1\n{seq}\n", encoding="utf-8")
    b = tmp_path / "clean.fna"
    b.write_text(">c1\n" + "".join(rng.choice("ACGT") for _ in range(20_000)) + "\n",
                 encoding="utf-8")
    dm = compute_sourmash_distances({"ambig": a, "clean": b}, ksize=21, scaled=10)
    assert np.isfinite(dm.matrix).all()
    assert dm.matrix.shape == (2, 2)


@needs_sourmash
def test_gzipped_fasta_is_read(tmp_path):
    import gzip

    rng = random.Random(99)
    seq = "".join(rng.choice("ACGT") for _ in range(20_000))
    plain = tmp_path / "p.fna"
    plain.write_text(f">c\n{seq}\n", encoding="utf-8")
    gz = tmp_path / "g.fna.gz"
    gz.write_bytes(gzip.compress(f">c\n{seq}\n".encode()))

    dm = compute_sourmash_distances({"plain": plain, "gz": gz}, ksize=21, scaled=100)
    idx = {s: i for i, s in enumerate(dm.strain_ids)}
    # Same sequence through both paths -> distance 0.
    assert dm.matrix[idx["plain"], idx["gz"]] == pytest.approx(0.0, abs=1e-9)


@needs_sourmash
def test_one_minus_ani_matches_sourmash_own_dist(tmp_path):
    """Pins that `ani_to_distance` is sourmash's own quantity, not a reinterpretation."""
    from sourmash import MinHash

    genomes = _write_genomes(tmp_path)
    sketches = {}
    for sid, path in genomes.items():
        mh = MinHash(n=0, ksize=21, scaled=100)
        mh.add_sequence(path.read_text(encoding="utf-8").split("\n", 1)[1].replace("\n", ""),
                        force=True)
        sketches[sid] = mh
    res = sketches["base"].jaccard_ani(sketches["near"])
    assert ani_to_distance(res.ani) == pytest.approx(res.dist, abs=1e-12)


# --------------------------------------------------------------------------- #
# the remediation advice must name a real extra
# --------------------------------------------------------------------------- #

def test_the_remediation_advice_names_an_extra_that_actually_exists():
    """`SourmashUnavailable` tells the reader to run `uv sync --extra lineage`.

    An advertised install command is a promise exactly like an advertised CLI
    invocation, and this repo has a recorded history of those dangling -- `--help`
    examples no parser accepts, a `--genome-fasta` flag that never existed. The
    advice is the FIRST thing a blocked user acts on, so a renamed or dropped extra
    would send them to a command that does nothing.

    Reads the module source rather than triggering the exception, so the guard holds
    whether or not sourmash happens to be installed in the running interpreter.
    """
    import tomllib

    src = (REPO_ROOT / "dna_decode" / "eval" / "sourmash_distance.py").read_text(
        encoding="utf-8")
    extras = tomllib.loads(
        (REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8")
    )["project"]["optional-dependencies"]

    assert "uv sync --extra lineage" in src, "the advice string moved or changed"
    assert "lineage" in extras, (
        "the advice names `--extra lineage` but pyproject has no such extra"
    )
    assert any(dep.replace(" ", "").startswith("sourmash") for dep in extras["lineage"])

    # The floor quoted in the message must be the floor actually declared, or the
    # pip fallback installs a version the module was never probed against (the
    # verified surface in its docstring is sourmash 4.9.4).
    assert "sourmash>=4.9" in src
    assert any("sourmash>=4.9" in dep.replace(" ", "") for dep in extras["lineage"]), (
        f"message pins sourmash>=4.9; pyproject declares {extras['lineage']}"
    )
