"""Locate peak tiles in the genome BY SEQUENCE, and refuse the ambiguous. Pure, offline, no torch.

**Why sequence search rather than the coordinates the substrate ships.** The substrate's reference
is `U00096.2`; the cached annotation is `GCF_000005845.2` / `NC_000913.3`. Measured across every
uniquely-located tile, the stated `start` column agrees with the genome at **1193/43061 = 2.77%**,
a coincidence rate. Located by sequence instead, the top stated-vs-located offsets are +1978 (387),
+1977 (374), +1976 (294), +777 (286), +776 (85), 0 (83) -- offset FAMILIES, i.e. cumulative indel
drift between assembly versions, not a constant shift that could be corrected. So reading the
stated columns would yield a frame that is fully self-consistent and entirely wrong, which is the
worst available failure: it produces numbers rather than an error.

This module therefore never reads a coordinate column, and `tests/test_tile_locate.py` asserts that
STRUCTURALLY with `ast` rather than by text search -- a substring scan cannot tell an executable
reference from the prose above, and this docstring is full of the very words it would match. (The
`Tile` dataclass does not even expose `start`/`end`, so the module would have to go back to the raw
file to get them; the structural assertion is what keeps that true as the file changes.)

**The sweep direction is the whole performance story.** Indexing the ~43k query sequences into a
dict and sliding both strands of the 4,641,652 bp genome once is measured at **5 s**; 43k
independent `str.find` calls over 4.6 Mb is not viable. Reverse-strand hits are normalised back to
FORWARD genome coordinates, so every returned position indexes the forward genome.

**Ambiguity is excluded and COUNTED, never resolved by taking the first hit.** A sequence occurring
more than once -- including a sequence that occurs once on each strand -- has no single frame, and
silently keeping hit one would assign it a position that is wrong half the time. `pos0` and
`strand` are `None` for such a sequence, so arithmetic on them raises instead of quietly
succeeding. Measured on the real cohort: 43,386 of 43,417 distinct sequences located (0.9993), of
which 43,061 are UNIQUE (0.9918), 325 ambiguous, 31 not found.

**The floor does not crash the run, and that is deliberate.** `LOCATE_FLOOR = 0.95` sits below the
measured 0.9918 with real headroom rather than being fitted to it. A report below the floor sets
`frame_claim_permitted` to False and carries the `LOCATE_BELOW_FLOOR` verdict; it does NOT raise,
because the experiment's grader owns an `INDETERMINATE_LOCATE_RATE` outcome and a module that
crashed first would make that branch dead code -- the exact defect the grader's reachability proof
exists to prevent. `require_floor()` is available for a caller that genuinely wants to stop.

**Provenance.** `scripts/frame_window_ordering_probe.py` carried this sweep first and measured the
figures above; this module is the extraction its own docstring names as the fix for the duplication.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

#: Below this UNIQUE-locate rate a frame claim is not permitted. Measured 0.9918.
LOCATE_FLOOR = 0.95

_COMPLEMENT = str.maketrans("ACGT", "TGCA")


class TileLocateError(ValueError):
    """The inputs cannot be swept."""


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


@dataclass(frozen=True)
class LocateResult:
    """Where one sequence sits in the forward genome, or why it has no single position.

    `pos0` and `strand` are populated ONLY when `n_occurrences == 1`. For an ambiguous sequence
    they are `None` so that using them raises a `TypeError` rather than silently computing a frame
    from an arbitrary hit.
    """

    n_occurrences: int
    pos0: int | None = None
    strand: str | None = None

    @property
    def is_unique(self) -> bool:
        return self.n_occurrences == 1


@dataclass(frozen=True)
class LocateReport:
    """Accounting over a sweep. Every distinct input sequence lands in exactly one bucket."""

    n_distinct: int
    n_unique: int
    n_ambiguous: int
    n_not_found: int
    tile_length: int
    genome_length: int

    def __post_init__(self) -> None:
        total = self.n_unique + self.n_ambiguous + self.n_not_found
        if total != self.n_distinct:
            raise TileLocateError(
                "buckets sum to %d but %d distinct sequences were swept; a reduced denominator "
                "reads as a better locate rate rather than as an error"
                % (total, self.n_distinct)
            )

    @property
    def n_located(self) -> int:
        return self.n_unique + self.n_ambiguous

    @property
    def unique_rate(self) -> float:
        return self.n_unique / float(self.n_distinct) if self.n_distinct else 0.0

    @property
    def located_rate(self) -> float:
        return self.n_located / float(self.n_distinct) if self.n_distinct else 0.0

    @property
    def meets_floor(self) -> bool:
        return self.unique_rate >= LOCATE_FLOOR

    @property
    def frame_claim_permitted(self) -> bool:
        return self.meets_floor

    @property
    def verdict(self) -> str:
        return "LOCATE_OK" if self.meets_floor else "LOCATE_BELOW_FLOOR"

    def require_floor(self) -> None:
        """Raise if the unique-locate rate is below the floor. For a caller that wants to stop."""
        if not self.meets_floor:
            raise TileLocateError(
                "unique-locate rate %.4f < floor %.2f: the frame is unknown for too much of the "
                "cohort, so the scored subset is no longer the cohort the arm measured"
                % (self.unique_rate, LOCATE_FLOOR)
            )

    def as_dict(self) -> dict:
        return {
            "n_distinct": self.n_distinct,
            "n_located": self.n_located,
            "n_unique": self.n_unique,
            "n_ambiguous": self.n_ambiguous,
            "n_not_found": self.n_not_found,
            "unique_rate": round(self.unique_rate, 6),
            "located_rate": round(self.located_rate, 6),
            "tile_length": self.tile_length,
            "genome_length": self.genome_length,
            "locate_floor": LOCATE_FLOOR,
            "verdict": self.verdict,
            "frame_claim_permitted": self.frame_claim_permitted,
        }


def read_fasta_single(path: Path | str) -> tuple[str, str]:
    """Return `(contig_name, sequence)` for a single-record FASTA, uppercased.

    REFUSES a multi-record file: the window replica and this sweep both assume one replicon, and
    concatenating two contigs would place tiles at positions that exist in neither.
    """
    name = None
    chunks: list[str] = []
    n_records = 0
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            n_records += 1
            if n_records == 1:
                name = line[1:].split()[0]
        else:
            chunks.append(line.strip())
    if n_records != 1:
        raise TileLocateError(
            "%s holds %d FASTA records; this sweep assumes ONE replicon and concatenating contigs "
            "would place tiles at positions that exist in neither" % (path, n_records)
        )
    return str(name), "".join(chunks).upper()


def locate_sequences(seqs, genome: str) -> dict[str, LocateResult]:
    """Sweep both strands of `genome` once, returning a result per DISTINCT input sequence.

    All inputs must share one length: the sweep is a fixed-width window, and a ragged cohort is
    REFUSED rather than swept at one width, because the mismatched sequences would silently come
    back not-found and read as a locate-rate problem.
    """
    distinct = {s.upper() for s in seqs}
    if not distinct:
        raise TileLocateError("locate_sequences got zero sequences")
    lengths = {len(s) for s in distinct}
    if len(lengths) != 1:
        raise TileLocateError(
            "ragged input (%d distinct lengths, %d-%d); a single-width sweep would report the "
            "other widths as not-found, which reads as a locate-rate problem rather than as a "
            "caller error" % (len(lengths), min(lengths), max(lengths))
        )
    L = lengths.pop()
    if L <= 0:
        raise TileLocateError("tile length must be positive")
    G = genome.upper()
    if len(G) < L:
        raise TileLocateError(
            "genome is %d bp, shorter than the %d bp tiles" % (len(G), L)
        )

    hits: dict[str, list[tuple[int, str]]] = {}
    for strand, strand_seq in (("+", G), ("-", reverse_complement(G))):
        n = len(strand_seq)
        for p in range(0, n - L + 1):
            key = strand_seq[p:p + L]
            if key in distinct:
                # normalise a reverse-strand hit back to FORWARD genome coordinates
                fwd = p if strand == "+" else n - p - L
                hits.setdefault(key, []).append((fwd, strand))

    out: dict[str, LocateResult] = {}
    for s in distinct:
        found = hits.get(s, [])
        if len(found) == 1:
            pos0, strand = found[0]
            out[s] = LocateResult(1, pos0, strand)
        else:
            out[s] = LocateResult(len(found))
    return out


def locate_report(located: dict[str, LocateResult], *, tile_length: int,
                  genome_length: int) -> LocateReport:
    """Bucket a sweep's results. The buckets are asserted to sum to the input count."""
    n_unique = sum(1 for r in located.values() if r.n_occurrences == 1)
    n_amb = sum(1 for r in located.values() if r.n_occurrences > 1)
    n_nf = sum(1 for r in located.values() if r.n_occurrences == 0)
    return LocateReport(
        n_distinct=len(located),
        n_unique=n_unique,
        n_ambiguous=n_amb,
        n_not_found=n_nf,
        tile_length=tile_length,
        genome_length=genome_length,
    )


def unique_positions(located: dict[str, LocateResult]) -> dict[str, tuple[int, str]]:
    """Just the sequences with exactly one hit, as `{seq: (pos0, strand)}`."""
    return {
        s: (r.pos0, r.strand)
        for s, r in located.items()
        if r.n_occurrences == 1 and r.pos0 is not None and r.strand is not None
    }


# ------------------------------------------------------------------------------------- caching


CACHE_SCHEMA = "tile-locate-cache-v1"


def save_locate_cache(path: Path | str, located: dict[str, LocateResult], *,
                      report: LocateReport) -> Path:
    """Persist a sweep so the experiment does not re-sweep. Gitignored data, not a wiki artifact.

    Stores every sequence including the ambiguous and not-found ones, because their COUNTS are part
    of the experiment's accounting and a cache holding only the usable rows could not reconcile.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "schema": CACHE_SCHEMA,
        "report": report.as_dict(),
        "located": {
            s: [r.n_occurrences, r.pos0, r.strand] for s, r in sorted(located.items())
        },
    }
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


def load_locate_cache(path: Path | str) -> dict[str, LocateResult]:
    """Read a sweep cache back. Refuses an unexpected schema rather than mis-parsing it."""
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    if payload.get("schema") != CACHE_SCHEMA:
        raise TileLocateError(
            "locate cache schema is %r, expected %r" % (payload.get("schema"), CACHE_SCHEMA)
        )
    out: dict[str, LocateResult] = {}
    for s, row in payload["located"].items():
        n, pos0, strand = row
        out[s] = LocateResult(int(n), pos0 if pos0 is None else int(pos0), strand)
    return out
