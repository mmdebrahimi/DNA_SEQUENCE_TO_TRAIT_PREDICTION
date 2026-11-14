"""Circular shift of a whole CDS annotation -- the structurally-matched null for the frame arm.

WHY THIS EXISTS, and why it is not the control that was originally proposed.

The frame arm's required matched-informativeness control failed its own premise: the non-positional
comparators (dinucleotide 0.2895, k-mer4 0.2811) are WEAKER standalone than the frame (0.3170), so
they were never a matched control and `attributable_absorption` is only weakly controlled. The
proposed replacement was a SHUFFLED-WINDOW control -- identical feature construction with the
window set relocated along the genome.

**That replacement was measured invalid before it was built**
(`wiki/promoteratlas_frame_column_split_2026-11-05.json`). `frame_features` receives window
information ONLY as the union `mask`, so relocating windows changes three of its five columns and
leaves two -- `signed_distance_to_next_cds_start` and `inside_cds` -- byte-identical. Measured over
the same 10 committed splits, those two CDS columns ALONE score 0.3183 standalone against the full
five-column table's 0.3170 (paired: the CDS pair is at least the full table on 6 of 10 seeds), and
handed to the GC baseline they absorb +0.0705 against the frame's +0.0687. A shuffled-window
control would therefore absorb slightly MORE than the treatment it controls, driving
`attributable_absorption` negative and firing an indeterminate verdict that is an artifact of
shared columns rather than a fact about position.

**A circular shift of the whole annotation is the construction that actually matches.** Every
`CdsRow` is translated by one offset modulo the genome length, and then BOTH halves of the feature
table are regenerated from the shifted rows -- the windows via `corpus_windows` and the CDS columns
via `frame_features`'s own `cds_rows` argument. Nothing is left invariant.

What it is matched on, exactly, and by construction rather than by hoping:

  * gene COUNT -- a translation is a bijection on rows
  * every gene LENGTH -- translation preserves `end0 - start0` per row, so the length multiset is
    identical, not merely similar
  * strand composition -- the strand field is untouched
  * the intergenic spacing multiset -- preserved, rotated
  * total CDS-covered bases, and therefore coding density

What it DESTROYS is the only thing it should: the phase between the annotation and the measured
tiles. It answers "is the frame table's absorption specific to the TRUE annotation positions, or
would any annotation of this structure absorb as much?"

**Why a circular shift and not independent per-window relocation.** E. coli is gene-dense. The
committed window set is 3,496 windows x 200 bp = 625,255 bases, about 13.5 percent of the genome,
and the genes themselves cover most of the rest. Relocating each window independently and
uniformly would land most of them inside coding sequence, so `in_window_any` would begin tracking
"inside a gene body" instead of "upstream of a gene start" -- a DIFFERENT confound substituted for
the one being removed, not a null. A rigid rotation cannot do that: the density profile travels
with the annotation.

**The genome is circular, so the rotation has no edge.** It does, however, have an ORIGIN
bookkeeping problem, and this module refuses rather than papers over it. A gene `[s, e)` shifted by
`delta` straddles the coordinate origin exactly when `(L - delta) mod L` falls strictly inside it,
because the shifted interval then crosses a multiple of `L`. `CdsRow` is a half-open interval with
`start0 < end0` and carries no wrap flag, so a straddling row cannot be represented.

Three ways to handle that, and the choice is load-bearing:

  (a) DROP straddling rows -- loses 1 or 2 genes and breaks the exact gene-count match, which is
      the property the whole construction rests on.
  (b) SPLIT them in two -- raises the gene count and invents two short genes that are not in any
      real annotation.
  (c) CHOOSE an offset at which nothing straddles -- exact on every matched quantity, with no edge
      handling at all.

This module does (c). `valid_shift_offsets` is precisely the set of offsets whose origin does not
land STRICTLY INSIDE a gene, read off a CDS-interior mask built with the same `window_coverage`
primitive `frame_features` uses, so the rule stays checkable against that primitive. The strictness
matters: an offset that aligns a gene's START to the origin is valid, and an earlier version of
this module rejected one such offset per gene. `circular_shift_cds` RAISES on an
offset that would straddle; it never silently drops or splits a row.
"""

from __future__ import annotations

import random

from dna_decode.glm.corpus_window import CdsRow, window_coverage


class AnnotationShiftError(ValueError):
    """The annotation cannot be shifted as requested."""


def circular_shift_cds(cds_rows, genome_length: int, delta: int) -> list[CdsRow]:
    """Translate every CDS row by `delta` bases modulo `genome_length`.

    Row ORDER is preserved, so a caller passing `order="file"` to `corpus_windows` still gets a
    meaningful file order. Strand and `locus_tag` are untouched -- the locus tag travels with its
    gene so a shifted window remains attributable to a named row in diagnostics.

    RAISES `AnnotationShiftError` if any row would straddle the coordinate origin. Use
    `valid_shift_offsets` (or `sample_shift_offset`) to pick an offset for which it cannot happen;
    dropping or splitting a straddling row would break the exact gene-count and gene-length match
    that is this control's entire basis.
    """
    if genome_length <= 0:
        raise AnnotationShiftError("genome_length must be positive, got %r" % (genome_length,))
    rows = [CdsRow(*r) for r in cds_rows]
    if not rows:
        raise AnnotationShiftError("no CDS rows to shift")

    d = int(delta) % genome_length
    out: list[CdsRow] = []
    for r in rows:
        if not (0 <= r.start0 < r.end0 <= genome_length):
            raise AnnotationShiftError(
                "row %r is not a half-open interval inside a %d-base genome; shifting it would "
                "produce coordinates that look valid and are not"
                % ((r.start0, r.end0, r.strand, r.locus_tag), genome_length)
            )
        s = (r.start0 + d) % genome_length
        e = r.end0 + d
        if e > genome_length:
            e -= genome_length
        if s >= e:
            raise AnnotationShiftError(
                "offset %d makes %s (%d..%d) straddle the origin. Pick an offset from "
                "valid_shift_offsets(): dropping the row would break the gene-COUNT match and "
                "splitting it would invent two genes that are in no real annotation."
                % (d, r.locus_tag, r.start0, r.end0)
            )
        out.append(CdsRow(s, e, r.strand, r.locus_tag))
    return out


def valid_shift_offsets(cds_rows, genome_length: int) -> list[int]:
    """Every offset at which NO row straddles the coordinate origin.

    A row `[s, e)` straddles under `delta` exactly when `(L - delta) mod L` lies STRICTLY inside
    it -- in the open interval `(s, e)`, not the half-open `[s, e)`. The boundary case is real and
    was a bug here until an independent interval-arithmetic derivation caught it: when
    `L - delta == s`, the shifted row begins exactly at the coordinate origin, which is a perfectly
    ordinary interval and straddles nothing. Testing plain CDS coverage rejects exactly one valid
    offset per gene -- the one that aligns that gene's start to the origin.

    So the mask built here covers each row's INTERIOR, `[s + 1, e)`, rather than the row itself.
    That keeps the derivation on the same `window_coverage` primitive `frame_features` uses for
    `inside_cds` (so the two cannot drift apart) while being exact rather than merely conservative.
    A length-1 row contributes an empty interior, which is correct: no integer offset can land
    strictly inside it.

    On a gene-dense genome this set is small as a FRACTION and large in absolute terms: E. coli is
    roughly 86 percent coding, which still leaves several hundred thousand usable offsets.
    """
    if genome_length <= 0:
        raise AnnotationShiftError("genome_length must be positive, got %r" % (genome_length,))
    interior = window_coverage([(r[0] + 1, r[1]) for r in cds_rows], genome_length)
    # delta is valid iff position (L - delta) % L is not strictly inside a row.
    return [d for d in range(genome_length) if not interior[(genome_length - d) % genome_length]]


def sample_shift_offset(cds_rows, genome_length: int, *, seed: int,
                        min_abs_shift: int = 50_000) -> int:
    """One reproducible valid offset, drawn uniformly from the offsets that strand nothing.

    `min_abs_shift` keeps the draw away from the identity and its neighbourhood measured
    CIRCULARLY -- an offset of 20 bases, or of `genome_length - 20`, is a near-copy of the true
    annotation and would make the control agree with the treatment for a trivial reason. It is a
    floor on a genuine rotation, not a tuning knob.
    """
    valid = valid_shift_offsets(cds_rows, genome_length)
    far = [d for d in valid if min(d, genome_length - d) >= min_abs_shift]
    if not far:
        raise AnnotationShiftError(
            "no offset is both non-straddling and at least %d bases from the identity on a "
            "%d-base genome (%d non-straddling offsets exist); the annotation is too dense for "
            "this control" % (min_abs_shift, genome_length, len(valid))
        )
    return int(random.Random(seed).choice(far))


def shift_match_report(cds_rows, shifted_rows, genome_length: int) -> dict:
    """What the shift preserved, MEASURED rather than asserted.

    The claim "matched by construction" is only worth making if it is checked, so this returns the
    quantities a reader would otherwise have to take on faith. Every `*_match` field must be True
    for the control to be the one this module documents.
    """
    a = [CdsRow(*r) for r in cds_rows]
    b = [CdsRow(*r) for r in shifted_rows]
    la = sorted(r.end0 - r.start0 for r in a)
    lb = sorted(r.end0 - r.start0 for r in b)
    ma = window_coverage([(r.start0, r.end0) for r in a], genome_length)
    mb = window_coverage([(r.start0, r.end0) for r in b], genome_length)
    sa = sorted(r.strand for r in a)
    sb = sorted(r.strand for r in b)
    return {
        "n_rows": len(a),
        "n_rows_shifted": len(b),
        "gene_count_match": len(a) == len(b),
        "gene_length_multiset_match": la == lb,
        "strand_composition_match": sa == sb,
        "cds_covered_bases": int(sum(ma)),
        "cds_covered_bases_shifted": int(sum(mb)),
        "cds_covered_bases_match": int(sum(ma)) == int(sum(mb)),
        "coding_density": round(sum(ma) / float(genome_length), 4),
        "note": (
            "A translation is a bijection on rows and preserves each row's length, so gene count, "
            "the gene-length multiset, strand composition and total covered bases are exact "
            "matches, not approximations. Window COUNT is reported separately by the caller and "
            "may differ by a few windows: corpus_windows skips genes closer to the linear "
            "coordinate boundary than one regulatory-region length, and a rotation moves which "
            "genes those are."
        ),
    }
