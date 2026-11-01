"""The non-tile control cohorts the peak-tile rung has discarded on every arm to date.

`genomewide.load_peak_tiles` parses all 46,713 rows of
`GSE144621_peak_tile_expression_formatted_std.txt.gz`, and every consumer then filters
`category == 'tile'`. That filter throws away **1,406 control sequences that are all exactly 150 bp**
-- 936 `random` and 470 `neg_control` -- i.e. the arm's own primary geometry, unchanged. They are the
substrate's own answer to "does this encoder separate real regulatory sequence from matched
non-regulatory sequence", which is the question the PromoterAtlas arm's load-bearing honest limit
asks and could not answer.

**This is a NEW module rather than an edit to `genomewide.py`, and that is not a style choice.**
`scripts/kaggle/build_promoteratlas_tile_arm_kernel.py` embeds `genomewide.py` as base64+gzip and
`scripts/promoteratlas_arm_dryrun.py` asserts the embedded payload is BYTE-IDENTICAL to the repo
copy. Adding a loader there would break that assert and silently fork project code -- the same
reasoning the arm kernel records for putting its GEO retry in the kernel instead.

**`pos_control` is refused, not padded.** All 106 of its rows are 48-75 bp, so ZERO survive at 150 bp.
Padding was measured on this exact encoder to contaminate the pooled embedding (max |delta| 1.232,
min cosine 0.282, against a between-tile spread of about 0.17) because attention lets the all-zero
columns reach the real positions, which masked pooling cannot undo. So the strongest positive control
is simply unavailable in the native-150 cell, and `load_control_tiles` raises rather than coercing it.
That narrows what the exposure experiment can show and is recorded as a limit, never worked around.
"""

from __future__ import annotations

import gzip
from dataclasses import dataclass
from pathlib import Path

#: The categories that survive at 150 bp. Not a default to be widened casually: `pos_control` is
#: geometry-incompatible (see the module docstring) and `tile` is the thing being controlled FOR.
CONTROL_CATEGORIES: tuple[str, ...] = ("random", "neg_control")

#: Measured 2026-11-01 from the cached substrate, after the same `ACGT`-only and length filters
#: `load_peak_tiles` applies. A loader that silently returns a different count is a defect, so the
#: real-data test asserts these exactly.
EXPECTED_AT_150: dict[str, int] = {"random": 936, "neg_control": 470}

PEAK_TILE_FILE = "GSE144621_peak_tile_expression_formatted_std.txt.gz"


class ControlCohortError(RuntimeError):
    """A requested control category cannot be served at the requested geometry."""


@dataclass(frozen=True)
class ControlSeq:
    """One non-tile control sequence, carrying the same fields the arm scores tiles on."""

    seq: str
    expression: float
    category: str
    rep1: float = 0.0
    rep2: float = 0.0


def _rows(path: Path):
    with gzip.open(path, "rt", errors="replace") as fh:
        hdr = next(fh).rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(hdr)}
        need = ("variant", "expn_med", "category", "RNA_exp_sum_1", "RNA_exp_sum_2")
        missing = [c for c in need if c not in ix]
        if missing:
            raise ControlCohortError(
                f"{path.name} lacks required columns {missing}; the delimiter sniff or the file is wrong"
            )
        for line in fh:
            q = line.rstrip("\n").split("\t")
            if len(q) == len(hdr):
                yield q, ix


def load_control_tiles(
    cache_dir: Path | str,
    *,
    categories: tuple[str, ...] = CONTROL_CATEGORIES,
    length: int = 150,
    allow_download: bool = False,
) -> list[ControlSeq]:
    """Load the non-tile control sequences at exactly `length` bp.

    `cache_dir` is REQUIRED and has no default on purpose. `genomewide.DEFAULT_CACHE` is
    `D:/dna_decode_cache/mpra`, a volume that is not attached on this host, so a defaulted call would
    either attempt a download or fail with a confusing device error. Every call site passes the path.

    `allow_download` exists only to make the offline intent explicit; this loader never fetches.

    RAISES `ControlCohortError` when a requested category yields zero rows at `length`. That is the
    load-bearing refusal: an empty list reads downstream as "the control found nothing", which is a
    statement about the encoder, when the truth would be "this category does not exist at this
    geometry" -- a statement about the substrate.
    """
    if allow_download:
        raise ControlCohortError(
            "load_control_tiles never downloads; the substrate is expected to be cached already"
        )
    d = Path(cache_dir)
    path = d / PEAK_TILE_FILE if d.is_dir() else d
    if not path.is_file():
        raise ControlCohortError(f"substrate not cached at {path}")

    want = set(categories)
    out: list[ControlSeq] = []
    per_cat: dict[str, int] = {c: 0 for c in want}
    for q, ix in _rows(path):
        cat = q[ix["category"]]
        if cat not in want:
            continue
        seq = q[ix["variant"]].upper()
        if len(seq) != length or (set(seq) - set("ACGT")):
            continue
        try:
            e = float(q[ix["expn_med"]])
            r1 = float(q[ix["RNA_exp_sum_1"]])
            r2 = float(q[ix["RNA_exp_sum_2"]])
        except ValueError:
            continue
        per_cat[cat] += 1
        out.append(ControlSeq(seq=seq, expression=e, category=cat, rep1=r1, rep2=r2))

    empty = sorted(c for c, n in per_cat.items() if n == 0)
    if empty:
        raise ControlCohortError(
            f"categories {empty} yield ZERO rows at {length} bp. This is a SUBSTRATE fact, not an "
            f"encoder result: `pos_control` rows are 48-75 bp, so none survives at 150. Padding them "
            f"is refused -- it was measured to contaminate this encoder's embedding (max |delta| "
            f"1.232). Request a geometry-compatible category, or accept the control is unavailable."
        )
    return out


def control_counts(controls: list[ControlSeq]) -> dict[str, int]:
    """Per-category counts, for the reconcile assertion and the artifact."""
    out: dict[str, int] = {}
    for c in controls:
        out[c.category] = out.get(c.category, 0) + 1
    return out
