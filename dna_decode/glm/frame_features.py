"""Positional frame features on SEQUENCE-LOCATED tile coordinates. Pure, offline, no torch.

These are the features handed to the cheap baseline in the additive half of the exposure
investigation: if corpus-window POSITION is what the frozen PromoterAtlas encoder's margin is made
of, then telling a GC baseline where the corpus windows are should absorb that margin.

**Every coordinate here comes from `tile_locate`, never from the substrate's columns.** The stated
column agrees with the cached genome at a coincidence rate, so a frame built on it would be fully
self-consistent and entirely wrong. A tile that is not UNIQUELY located gets no feature row at all
-- it is excluded and counted, never imputed -- because imputing a frame for an ambiguous tile
invents the very quantity under test.

**Two granularities ship, because neither is obviously right and the choice must not be made after
seeing which answers better.** Measured on the real cohort: `in_window_any` fires on 28,622 of
43,061 uniquely-located tiles (66.47%) and `in_window_full` on 8,239 (19.13%). The frozen bar names
which one is PRIMARY in advance and the other is reported beside it.

**`in_window_any`'s threshold is NAMED, not left implicit** (`IN_WINDOW_ANY_MIN_OVERLAP_BASES = 1`,
i.e. strictly more than zero overlapping bases). The pre-execution review asked for this because a
small epsilon would move the count. The related worry -- that the flag is dominated by trivial
one-base overlaps -- was MEASURED AND REFUTED: only 8.2% of in-window tiles overlap by 0.2 or less.
`overlap_fraction` ships alongside as the graded form, so a consumer need not accept the binary.

**`in_window_full` means fully covered by the UNION of windows, which is the measured definition.**
The technical plan's prose says "entirely inside one", but the figure it pins (8,239) comes from
`scripts/frame_plan_feasibility_probe.py`, which computes `c == TILE_LEN` against a single union
coverage mask. Union coverage and inside-a-single-window are not the same predicate -- two
overlapping windows can jointly cover a tile that neither contains -- so this module implements the
one the committed number was measured with, and says so rather than letting the two drift under one
name.

**The sign convention is the part most likely to be got wrong silently, so it is explicit.**
`signed_distance_to_next_cds_start` is POSITIVE when the tile lies UPSTREAM of the gene start IN
THE GENE'S OWN ORIENTATION, on both strands. For a plus-strand gene that means lower genome
coordinates (`start0 - midpoint`); for a minus-strand gene it means HIGHER genome coordinates
(`midpoint - end0`), because a minus-strand gene reads the other way along the genome and its
corpus window sits at `[end0, end0 + 200)`. So a tile sitting in a corpus window scores positive on
either strand, and a sign error would show up as in-window tiles scoring negative.

The nearest same-strand gene start is used, signed -- not "the next start ahead" -- because an
ahead-only definition is undefined for a tile with no gene ahead of it and would need a sentinel,
and a sentinel in a feature matrix becomes a magic number the ridge happily fits.
"""

from __future__ import annotations

import bisect
from dataclasses import dataclass, field

from dna_decode.glm.corpus_window import window_coverage

#: Row order of the emitted matrix. Fixed, because the bar names a primary by NAME.
FEATURE_NAMES: tuple[str, ...] = (
    "in_window_any",
    "in_window_full",
    "overlap_fraction",
    "signed_distance_to_next_cds_start",
    "inside_cds",
)

#: `in_window_any` fires at strictly MORE than zero overlapping bases. Named because an epsilon
#: here moves the count, and the review asked for the threshold to be stated rather than implied.
IN_WINDOW_ANY_MIN_OVERLAP_BASES = 1


class FrameFeatureError(ValueError):
    """The inputs cannot produce a feature table."""


@dataclass(frozen=True)
class FrameFeatureTable:
    """Features for the uniquely-located sequences, plus the accounting for everything else.

    `seqs` is the row order of `matrix`. Keyed by SEQUENCE rather than by tile, because 689 of the
    44,106 scored rows share a sequence with another row; a consumer maps tile -> sequence -> row.
    """

    names: tuple[str, ...]
    seqs: list[str]
    matrix: list[list[float]]
    excluded: dict[str, str] = field(default_factory=dict)
    n_input: int = 0

    def __post_init__(self) -> None:
        if len(self.seqs) != len(self.matrix):
            raise FrameFeatureError(
                "%d sequences but %d feature rows" % (len(self.seqs), len(self.matrix))
            )
        total = len(self.seqs) + len(self.excluded)
        if total != self.n_input:
            raise FrameFeatureError(
                "emitted %d + excluded %d = %d, but %d sequences were supplied. An earlier "
                "abstention arm in this repo turned true positives into an apparent bias through "
                "a quietly reduced denominator, and the defect was caught by exactly this "
                "arithmetic." % (len(self.seqs), len(self.excluded), total, self.n_input)
            )
        for row in self.matrix:
            if len(row) != len(self.names):
                raise FrameFeatureError(
                    "a row has %d values for %d feature names" % (len(row), len(self.names))
                )

    def column(self, name: str) -> list[float]:
        """One feature across every emitted row."""
        try:
            j = self.names.index(name)
        except ValueError:
            raise FrameFeatureError(
                "no feature %r; have %r" % (name, self.names)
            ) from None
        return [row[j] for row in self.matrix]

    def is_constant(self, name: str) -> bool:
        """A constant feature cannot absorb anything, so the bar checks this before grading."""
        col = self.column(name)
        return not col or len(set(col)) == 1

    @property
    def n_excluded_ambiguous(self) -> int:
        return sum(1 for v in self.excluded.values() if v == "ambiguous")

    @property
    def n_excluded_not_found(self) -> int:
        return sum(1 for v in self.excluded.values() if v == "not_found")

    def accounting(self) -> dict:
        return {
            "n_input_distinct": self.n_input,
            "n_emitted": len(self.seqs),
            "n_excluded": len(self.excluded),
            "n_excluded_ambiguous": self.n_excluded_ambiguous,
            "n_excluded_not_found": self.n_excluded_not_found,
            "feature_names": list(self.names),
            "in_window_any_min_overlap_bases": IN_WINDOW_ANY_MIN_OVERLAP_BASES,
            "in_window_full_definition": "fully covered by the UNION of corpus windows",
        }


def _orientation_starts(cds_rows) -> tuple[list[int], list[int]]:
    """Sorted gene-start positions in each strand's own orientation.

    A plus-strand gene starts at its low coordinate; a minus-strand gene starts at its high one.
    """
    plus = sorted(int(r[0]) for r in cds_rows if r[2] == "+")
    minus = sorted(int(r[1]) for r in cds_rows if r[2] == "-")
    return plus, minus


def _signed_distance(starts: list[int], midpoint: int, strand: str) -> float:
    """Signed distance to the NEAREST same-strand gene start. Positive = tile is upstream.

    Uses bisect rather than a scan: 43k tiles against 4,340 genes would otherwise be ~188M
    comparisons.
    """
    if not starts:
        raise FrameFeatureError(
            "no CDS rows on strand %r, so no distance is defined; a sentinel here would become a "
            "magic number the ridge fits" % strand
        )
    i = bisect.bisect_left(starts, midpoint)
    best = None
    for j in (i - 1, i):
        if 0 <= j < len(starts):
            s = starts[j]
            d = float(s - midpoint) if strand == "+" else float(midpoint - s)
            if best is None or abs(d) < abs(best):
                best = d
    return float(best)


def frame_features(located, mask: bytes, cds_rows, *, tile_length: int) -> FrameFeatureTable:
    """Build the frame feature table for every UNIQUELY located sequence.

    `located` is `{seq: LocateResult}` from `tile_locate.locate_sequences`; `mask` is the corpus
    window coverage from `corpus_window.window_coverage`, whose length defines the genome length.

    The plan's signature also listed the window LIST. It is not taken, deliberately: the only
    predicate that would use it is "inside a single window", and the committed `in_window_full`
    figure was measured against the union mask. Accepting both inputs would invite two different
    definitions to drift under one name.
    """
    if tile_length <= 0:
        raise FrameFeatureError("tile_length must be positive, got %r" % (tile_length,))
    genome_length = len(mask)
    if genome_length <= 0:
        raise FrameFeatureError("the coverage mask is empty")

    cds_mask = window_coverage(cds_rows, genome_length)
    plus_starts, minus_starts = _orientation_starts(cds_rows)

    seqs: list[str] = []
    matrix: list[list[float]] = []
    excluded: dict[str, str] = {}

    for seq in sorted(located):
        res = located[seq]
        if res.n_occurrences != 1 or res.pos0 is None or res.strand is None:
            excluded[seq] = "ambiguous" if res.n_occurrences > 1 else "not_found"
            continue

        pos0, strand = int(res.pos0), str(res.strand)
        if pos0 < 0 or pos0 + tile_length > genome_length:
            raise FrameFeatureError(
                "located tile at %d..%d falls outside the %d-base mask; a truncated denominator "
                "would read as a lower overlap rather than as an error"
                % (pos0, pos0 + tile_length, genome_length)
            )

        covered = sum(mask[pos0:pos0 + tile_length])
        frac = covered / float(tile_length)
        midpoint = pos0 + tile_length // 2
        starts = plus_starts if strand == "+" else minus_starts

        matrix.append([
            1.0 if covered >= IN_WINDOW_ANY_MIN_OVERLAP_BASES else 0.0,
            1.0 if covered == tile_length else 0.0,
            frac,
            _signed_distance(starts, midpoint, strand),
            1.0 if cds_mask[midpoint] else 0.0,
        ])
        seqs.append(seq)

    return FrameFeatureTable(
        names=FEATURE_NAMES,
        seqs=seqs,
        matrix=matrix,
        excluded=excluded,
        n_input=len(located),
    )
