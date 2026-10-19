"""Genome distance via sourmash FracMinHash, as a drop-in `DistanceMatrix` producer.

WHY THIS EXISTS
---------------
`eval/phylogeny.compute_mash_distances` is the repo's genome-distance producer and it
needs the Mash binary. There is no native Windows Mash build (the v2.3 release ships
only `mash-Linux64` / `mash-OSX64`, and bioconda publishes no `win-64`), so on this
host Mash has exactly two routes -- Docker or WSL2 -- and BOTH need the same
virtualization layer, which is gated in BIOS. `sourmash` ships a `win_amd64` wheel,
so it is the one genome-distance route that runs natively here.

FracMinHash IS NOT MASH'S BOTTOM-SKETCH MINHASH
-----------------------------------------------
Mash draws a fixed-size bottom sketch (`-s`); sourmash keeps every hash below
`max_hash = 2**64 / scaled`, so sketch size scales with genome length. The resulting
ANI estimators are therefore NOT the same estimator, and the fixed `(0.001, 0.005)`
thresholds this repo's clonality work uses were calibrated against MASH distances.
**Threshold semantics under FracMinHash are ASSUMED, not measured, until the
sourmash-vs-Mash crosscheck runs and is scored against its frozen bar.** Do not read
a sourmash partition as interchangeable with a Mash partition before then.

VERIFIED SURFACE (probed 2026-10-19, sourmash 4.9.4, py3.11 win_amd64)
----------------------------------------------------------------------
Recorded because the planning pass marked this `[unverified]` and one plan
assumption turned out to be wrong:

  * `MinHash(n, ksize, *, scaled=0, ...)` -- `n` is POSITIONAL and required;
    `n=0` with `scaled=N` selects FracMinHash.
  * `MinHash.add_sequence(seq, force=False)` -- `force=True` is required for real
    assemblies, which carry non-ACGT characters (N runs, IUPAC ambiguity codes).
  * `MinHash.jaccard_ani(other)` returns a `jaccardANIResult` carrying
    `dist`, `size_is_inaccurate`, and `ani` **as a property**.
  * `jaccardANIResult.ani` IS `None` WHENEVER `size_is_inaccurate` is True. The
    plan assumed `.ani` was always a float. Measured: a 50 kb sequence at
    `scaled=1000` sketches to 49 hashes and returns `ani=None`; ~207 hashes is
    already enough to return a float. A real 5 Mb bacterial genome at
    `scaled=1000` sketches to ~4,900 hashes and is comfortably reliable.
  * `1 - ani` reproduces `.dist` exactly (checked: ani 0.9995408641127369 against
    dist 0.00045913588726309573), so `ani_to_distance` is not a reinterpretation
    of sourmash's own number -- it is the same quantity, named.

The `None` case is why `SketchTooSmall` exists. Letting it through would put
`None`/`NaN` into the matrix, and a NaN distance silently clusters as "far from
everything" -- a broken measurement wearing the costume of a biological result.
"""
from __future__ import annotations

import gzip
import itertools
from pathlib import Path

import numpy as np

from dna_decode.eval.phylogeny import DistanceMatrix

__all__ = [
    "SourmashUnavailable",
    "SketchTooSmall",
    "EmptySketch",
    "ani_to_distance",
    "compute_sourmash_distances",
]


class SourmashUnavailable(RuntimeError):
    """`sourmash` is not importable.

    Named, and named narrowly, so a caller can distinguish "the optional extra is
    not installed" from "the measurement ran and failed". Mirrors
    `forward/prosst_scorer.StructureMethodUnavailable`.
    """


class SketchTooSmall(RuntimeError):
    """sourmash reported `size_is_inaccurate`, so its own ANI estimate is `None`.

    Raised rather than coerced. The remedy is a smaller `scaled` (a denser sketch),
    not a fabricated distance.
    """


class EmptySketch(RuntimeError):
    """A genome produced zero hashes -- an empty, unreadable, or too-short FASTA."""


def ani_to_distance(ani: float) -> float:
    """Convert an ANI in [0, 1] to a distance in [0, 1] as `1 - ani`.

    One named, tested helper rather than an inline expression, because the UNITS are
    the whole point: the fixed `(0.001, 0.005)` clustering thresholds are distances,
    and silently feeding a similarity where a distance is expected inverts every
    comparison while still producing a plausible-looking matrix.
    """
    if ani is None or not isinstance(ani, (int, float)) or isinstance(ani, bool):
        raise ValueError(f"ani must be a real number in [0, 1], got {ani!r}")
    if not (0.0 <= float(ani) <= 1.0):
        raise ValueError(f"ani must be in [0, 1], got {ani!r}")
    return 1.0 - float(ani)


def _read_sequences(path: Path) -> list[str]:
    """Sequences from a FASTA, transparently handling `.gz`.

    Delegates parsing to the repo's existing `pathotype.detect.parse_fasta` rather
    than adding a fourth FASTA reader; the only thing added here is decompression,
    which that function does not do.
    """
    from dna_decode.pathotype.detect import parse_fasta

    path = Path(path)
    if path.suffix == ".gz":
        text = gzip.decompress(path.read_bytes()).decode("utf-8", errors="replace")
        out, name, buf = [], None, []
        for line in text.splitlines():
            if line.startswith(">"):
                if name is not None:
                    out.append("".join(buf))
                name, buf = line, []
            else:
                buf.append(line.strip())
        if name is not None:
            out.append("".join(buf))
        return [s for s in out if s]
    return [seq for _, seq in parse_fasta(path) if seq]


def _minhash_cls():
    try:
        from sourmash import MinHash
    except ImportError as exc:  # pragma: no cover - exercised via monkeypatch
        raise SourmashUnavailable(
            "sourmash is not installed. It is an optional extra: "
            "`uv sync --extra lineage` (or `pip install 'sourmash>=4.9'`)."
        ) from exc
    return MinHash


def compute_sourmash_distances(
    strain_genomes: dict[str, Path],
    *,
    ksize: int = 21,
    scaled: int = 1000,
) -> DistanceMatrix:
    """Pairwise genome distances for `{strain_id: fasta_path}` as a `DistanceMatrix`.

    The return type is `eval.phylogeny.DistanceMatrix` imported by reference, so the
    dataclass keeps exactly one definition and the result feeds
    `clonality.greedy_representative_clusters_from_matrix` unchanged.

    Raises `SourmashUnavailable`, `EmptySketch`, or `SketchTooSmall` -- never a bare
    `Exception`, and never a matrix containing `NaN`.
    """
    if ksize <= 0:
        raise ValueError(f"ksize must be positive, got {ksize}")
    if scaled <= 0:
        raise ValueError(f"scaled must be positive, got {scaled}")

    strain_ids = list(strain_genomes)
    n = len(strain_ids)
    if n == 0:
        return DistanceMatrix(strain_ids=[], matrix=np.zeros((0, 0), dtype=float))

    minhash_cls = _minhash_cls()

    sketches = {}
    for sid in strain_ids:
        sketch = minhash_cls(n=0, ksize=ksize, scaled=scaled)
        for seq in _read_sequences(Path(strain_genomes[sid])):
            # force=True: real assemblies carry N runs and IUPAC ambiguity codes.
            sketch.add_sequence(seq, force=True)
        if not len(sketch.hashes):
            raise EmptySketch(
                f"{sid!r} ({strain_genomes[sid]}) produced zero hashes at "
                f"ksize={ksize}, scaled={scaled}"
            )
        sketches[sid] = sketch

    matrix = np.zeros((n, n), dtype=float)
    for i, j in itertools.combinations(range(n), 2):
        a, b = strain_ids[i], strain_ids[j]
        result = sketches[a].jaccard_ani(sketches[b])
        if getattr(result, "size_is_inaccurate", False) or result.ani is None:
            raise SketchTooSmall(
                f"sourmash reports its ANI estimate unreliable for {a!r} vs {b!r} "
                f"at ksize={ksize}, scaled={scaled} "
                f"({len(sketches[a].hashes)} and {len(sketches[b].hashes)} hashes). "
                "Lower `scaled` for a denser sketch; do not substitute a distance."
            )
        dist = ani_to_distance(result.ani)
        matrix[i, j] = matrix[j, i] = dist

    return DistanceMatrix(strain_ids=strain_ids, matrix=matrix)
