"""The non-tile control cohorts the peak-tile rung has discarded on every arm to date.

`genomewide.load_peak_tiles` parses all 46,713 rows of
`GSE144621_peak_tile_expression_formatted_std.txt.gz`, and every consumer then filters
`category == 'tile'`. That filter throws away **1,406 control sequences that are all exactly 150 bp**
-- 936 `random` and 470 `neg_control` -- i.e. the arm's own primary geometry, unchanged. They are the
substrate's own answer to "does this encoder separate real regulatory sequence from matched
non-regulatory sequence", which is the question the PromoterAtlas arm's load-bearing honest limit
asks and could not answer.

**This module does NOT open the substrate. It CALLS the owner and filters the result.** That is a
correction: the first version parsed the file itself with a hardcoded tab split, which made it a
second reader of a substrate that already has one, and `tests/test_substrate_loader_ownership.py`
registers `glm/genomewide.py` as that owner. The trap the guard names is real -- the fragment files
in the same series are SPACE-delimited while the peak-tile file is TAB-delimited, so a second loader
that does not inherit `sniff_delimiter` can parse one of them into zero rows, which reads as "no
signal" rather than "no data". Routing through the owner also means this module inherits its
non-ACGT and unparseable-value filters for free rather than repeating them.

**Measured, which is why the correction costs nothing:** the owner already returns the control rows.
`load_peak_tiles` yields all four categories, and filtering its output for 150 bp gives 936 `random`
and 470 `neg_control` -- byte-identical to what the second parser produced.

**`pos_control` is refused, not padded.** All 106 of its rows are 48-75 bp, so ZERO survive at 150 bp.
Padding was measured on this exact encoder to contaminate the pooled embedding (max |delta| 1.232,
min cosine 0.282, against a between-tile spread of about 0.17) because attention lets the all-zero
columns reach the real positions, which masked pooling cannot undo. So the strongest positive control
is simply unavailable in the native-150 cell, and `load_control_tiles` raises rather than coercing it.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

#: The categories that survive at 150 bp. Not a default to be widened casually: `pos_control` is
#: geometry-incompatible (see the module docstring) and `tile` is the thing being controlled FOR.
CONTROL_CATEGORIES: tuple[str, ...] = ("random", "neg_control")

#: Measured 2026-11-01 THROUGH THE OWNER, after its own ACGT and value filters plus the 150 bp
#: length filter applied here. A loader that silently returns a different count is a defect, so the
#: real-data test asserts these exactly.
EXPECTED_AT_150: dict[str, int] = {"random": 936, "neg_control": 470}


class ControlCohortError(RuntimeError):
    """A requested control category cannot be served at the requested geometry."""


@dataclass(frozen=True)
class ControlSeq:
    """One non-tile control sequence, carrying the same fields the arm scores tiles on.

    Deliberately NOT the owner's `Tile` type. `Tile` carries a `peak` field, which is the SPLIT UNIT
    for tiles -- tiles inside one peak overlap heavily, so the arm holds out whole peaks -- and a
    control sequence has no peak. Returning `Tile` would hand callers a field that is meaningless
    here and invite someone to split controls on it.
    """

    seq: str
    expression: float
    category: str
    rep1: float = 0.0
    rep2: float = 0.0


def load_control_tiles(
    cache_dir: Path | str,
    *,
    categories: tuple[str, ...] = CONTROL_CATEGORIES,
    length: int = 150,
    allow_download: bool = False,
) -> list[ControlSeq]:
    """Load the non-tile control sequences at exactly `length` bp, via the substrate's owner.

    `cache_dir` is a **DIRECTORY** and is required. Two reasons, and the narrowing is deliberate.
    The owner resolves its own path by joining a cache directory with the file name, so a file path
    is no longer meaningful here. And `genomewide.DEFAULT_CACHE` is `D:/dna_decode_cache/mpra`, a
    volume that is not attached on this host, so a defaulted call would either attempt a download or
    fail with a confusing device error. Every call site passes the path explicitly.

    `allow_download` exists only to make the offline intent explicit; this loader never fetches.

    RAISES `ControlCohortError` when a requested category yields zero rows at `length`. That is the
    load-bearing refusal: an empty list reads downstream as "the control found nothing", which is a
    statement about the encoder, when the truth would be "this category does not exist at this
    geometry" -- a statement about the substrate.
    """
    from dna_decode.glm.genomewide import GenomeWideDataError, load_peak_tiles

    if allow_download:
        raise ControlCohortError(
            "load_control_tiles never downloads; the substrate is expected to be cached already"
        )
    d = Path(cache_dir)
    if d.exists() and not d.is_dir():
        raise ControlCohortError(
            f"{d} is a file; this function takes the cache DIRECTORY because the owner "
            f"(genomewide.load_peak_tiles) resolves the substrate's filename itself"
        )

    try:
        rows, _stats = load_peak_tiles(d, allow_download=False)
    except GenomeWideDataError as exc:
        # Translate the owner's error into this module's public failure type so existing callers
        # keep a stable contract. The owner raises this for an unreachable volume AND for a missing
        # file with downloading disabled, which are the two cases a caller here can hit.
        raise ControlCohortError(str(exc)) from exc

    want = set(categories)
    out: list[ControlSeq] = []
    per_cat: dict[str, int] = {c: 0 for c in want}
    for t in rows:
        if t.category not in want or len(t.seq) != length:
            continue
        per_cat[t.category] += 1
        out.append(
            ControlSeq(
                seq=t.seq,
                expression=t.expression,
                category=t.category,
                rep1=t.rep1,
                rep2=t.rep2,
            )
        )

    empty = sorted(c for c, n in per_cat.items() if n == 0)
    if empty:
        raise ControlCohortError(
            f"categories {empty} yield ZERO rows at {length} bp. This is a SUBSTRATE fact, not an "
            f"encoder result: `pos_control` rows are 48-75 bp, so none survives at 150. Padding "
            f"them is refused -- it was measured to contaminate this encoder's embedding (max "
            f"|delta| 1.232). Request a geometry-compatible category, or accept the control is "
            f"unavailable."
        )
    return out


def control_counts(controls: list[ControlSeq]) -> dict[str, int]:
    """Per-category counts, for the reconcile assertion and the artifact."""
    out: dict[str, int] = {}
    for c in controls:
        out[c.category] = out.get(c.category, 0) + 1
    return out
