"""Replica of the PromoterAtlas corpus-extraction window. Pure, offline, no torch.

**What this is for.** The frozen PromoterAtlas encoder was pretrained on a corpus built by
`extract_promoter_regions(gbs, min_intergenic_length=5, regulatory_region_length=200)` --
`src/promoter_atlas/utils/genomics.py:6` in the MIT-licensed clone -- which takes the 200 bp
immediately upstream of each annotated gene start, strand-aware, gated on a 5 bp minimum
intergenic distance. To ask whether the arm's 0.5076 is explained by corpus-window POSITION, the
rule has to be reconstructed rather than guessed at.

**Two upstream quirks are reproduced DELIBERATELY, and a future tidy-up must fail loudly.**
Verified against the upstream source, not paraphrased from a memo:

1. On both strands the `continue` guard sits BEFORE the running `last_end` assignment, so a gene
   skipped for sitting too close to the replicon end does NOT update `last_end`. The next gene is
   therefore gated against an older neighbour than the one immediately preceding it.
2. The reverse branch tests `last_end_reverse - end >= min_intergenic_length` while
   `last_end_reverse` holds the previous gene's **start** (`last_end_reverse = start` closes the
   loop body). On the forward branch the same variable holds the previous gene's **end**.

Both look like defects. The question this answers is what the corpus CONTAINED, not what it should
have contained -- a corrected replica would measure a rule nobody trained on.
`tests/test_corpus_window.py` pins each quirk with a case whose window set CHANGES if it is
corrected.

**This does NOT claim exact reproduction of the published corpus, and the reason is measured.**
The upstream function threads `last_end` through features in FILE order, and the cached GFF3's
coding-sequence rows are not ascending by start on either strand, so "replicate faithfully" and
"sort first" are different programs: file order gives 3,490 windows over 624,812 bases, sorted
gives 3,496 over 625,255, a symmetric difference of 26 windows
(`wiki/frame_plan_window_ordering_2026-11-01.json`). Compounding it, the authors read GenBank
rather than this GFF3, so "file order" here is a different file's order and faithful reproduction
is not well-defined from this input at all. `order` is therefore an EXPLICIT parameter defaulting
to `"sorted"`, and every figure pinned off this module is a sorted-order figure. The claim it
supports is scoped to *the corpus rule as specified at 200 bp with a 5 bp intergenic gate*, never
*as the authors ran it*.

**Provenance.** `scripts/frame_window_ordering_probe.py` carried this replica first and measured
the figures above; its own docstring records that if the two copies drift "the fix on record is to
extract the replica into the package rather than to re-merge the scripts". This module is that
extraction, so the drift hazard is discharged rather than restated.
"""

from __future__ import annotations

from typing import NamedTuple

#: The upstream defaults. Not tunable knobs -- the corpus was built at these values.
REGULATORY_REGION_LENGTH = 200
MIN_INTERGENIC_LENGTH = 5

#: Accepted values for `order`. `"sorted"` is the default and the one every pinned figure uses.
ORDERS = ("sorted", "file")


class CdsRow(NamedTuple):
    """One coding-sequence feature in 0-BASED HALF-OPEN coordinates.

    The conversion from GFF3's 1-based-inclusive columns happens in exactly one place --
    `cds_rows_from_annotation_table` -- and is asserted by test rather than assumed, because an
    off-by-one here shifts every window by a base and still produces a fully self-consistent
    result.
    """

    start0: int
    end0: int          # exclusive
    strand: str        # "+" or "-"
    locus_tag: str


class Window(NamedTuple):
    """One corpus window in 0-based half-open coordinates, with the gene it was taken for."""

    start0: int
    end0: int          # exclusive
    strand: str
    locus_tag: str


class CorpusWindowError(ValueError):
    """The inputs cannot be used to build a window set."""


def cds_rows_from_annotation_table(tab) -> list[CdsRow]:
    """Adapt `dna_decode.data.annotations.parse_gff3` output to `CdsRow`, in FILE order.

    Filters to `type == "CDS"` carrying a NON-EMPTY `locus_tag`, because the upstream extractor
    filters on `"locus_tag" in feature.qualifiers`. Row order is preserved exactly as the parser
    returned it, which is what makes `order="file"` meaningful.

    **This is the single 1-based-inclusive to 0-based-half-open boundary.** GFF3 `start` is
    1-based inclusive and `end` is inclusive; Biopython's `location.start` is 0-based and `.end`
    is exclusive, which is what the upstream function consumed. So `start0 = start - 1` and
    `end0 = end`.
    """
    sub = tab[(tab["type"] == "CDS") & (tab["locus_tag"] != "")]
    out: list[CdsRow] = []
    for r in sub.itertuples():
        strand = str(r.strand)
        if strand not in ("+", "-"):
            continue
        out.append(CdsRow(int(r.start) - 1, int(r.end), strand, str(r.locus_tag)))
    return out


def corpus_windows(
    cds_rows,
    genome_length: int,
    *,
    order: str = "sorted",
    regulatory_region_length: int = REGULATORY_REGION_LENGTH,
    min_intergenic_length: int = MIN_INTERGENIC_LENGTH,
) -> list[Window]:
    """The ordered corpus windows for one replicon.

    `cds_rows` is a sequence of `CdsRow` in 0-based half-open coordinates. `order` decides how
    `last_end` is threaded and is an INPUT to the result, not a presentation detail: see the
    module docstring. Returns forward windows then reverse windows, each in the order emitted.
    """
    if order not in ORDERS:
        raise CorpusWindowError(
            "order=%r is not one of %r; the threading order changes which windows exist at all, "
            "so it may not be guessed" % (order, ORDERS)
        )
    if genome_length <= 0:
        raise CorpusWindowError("genome_length must be positive, got %r" % (genome_length,))

    rows = [CdsRow(*r) for r in cds_rows]
    fwd = [r for r in rows if r.strand == "+"]
    rev = [r for r in rows if r.strand == "-"]
    if order == "sorted":
        # Sort on (start0, end0) ONLY -- a stable sort then preserves file order among ties, and
        # including locus_tag in the key would make the result depend on gene naming.
        fwd = sorted(fwd, key=lambda r: (r.start0, r.end0))
        rev = sorted(rev, key=lambda r: (r.start0, r.end0))

    windows: list[Window] = []

    # ---- forward strand: the window is the RL bases immediately upstream of the gene start ----
    last_end = 0
    for r in fwd:
        if r.start0 < regulatory_region_length:
            continue                      # QUIRK 1: precedes the last_end update, deliberately
        if r.start0 - last_end >= min_intergenic_length:
            windows.append(
                Window(r.start0 - regulatory_region_length, r.start0, "+", r.locus_tag)
            )
        last_end = r.end0

    # ---- reverse strand: iterated end-of-genome-first; the window sits DOWNSTREAM in genome
    # coordinates, which is upstream relative to the gene's own direction ----------------------
    last_end_rev = genome_length
    for r in reversed(rev):
        if genome_length - r.end0 < regulatory_region_length:
            continue                      # QUIRK 1 again, same shape
        if last_end_rev - r.end0 >= min_intergenic_length:
            windows.append(
                Window(r.end0, r.end0 + regulatory_region_length, "-", r.locus_tag)
            )
        last_end_rev = r.start0            # QUIRK 2: the previous gene's START, not its end

    return windows


def window_coverage(windows, genome_length: int) -> bytes:
    """Per-base mask: 1 where any corpus window covers the base. Windows may overlap."""
    if genome_length <= 0:
        raise CorpusWindowError("genome_length must be positive, got %r" % (genome_length,))
    cov = bytearray(genome_length)
    for w in windows:
        a = max(int(w[0]), 0)
        b = min(int(w[1]), genome_length)
        if b > a:
            cov[a:b] = bytes([1]) * (b - a)
    return bytes(cov)


def overlap_fraction(mask: bytes, pos0: int, length: int) -> float:
    """What fraction of a `length`-bp tile at `pos0` sits inside a corpus window.

    REFUSES a tile that would run off the end of the mask rather than silently measuring a short
    window, because a truncated denominator reads as a lower overlap instead of as an error.
    """
    if length <= 0:
        raise CorpusWindowError("length must be positive, got %r" % (length,))
    if pos0 < 0 or pos0 + length > len(mask):
        raise CorpusWindowError(
            "tile at %d..%d falls outside the %d-base mask" % (pos0, pos0 + length, len(mask))
        )
    return sum(mask[pos0:pos0 + length]) / float(length)
