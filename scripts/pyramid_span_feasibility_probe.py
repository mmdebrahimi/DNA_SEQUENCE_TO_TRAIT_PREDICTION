#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Measure the maximum HONEST pyramid span each cached rung permits, and compare it to 4 decades.

WHY THIS EXISTS -- it is the `multiresolution-pyramid-4-decades-UNRUN` arm's OWN EXECUTED EVIDENCE.

That arm sat `planned` on the ceiling criterion, so `ceiling_verdict` returned `blocked:untried-arm`
forever. It cannot be retired by a sibling's result or on a redundancy reason: R6 permits a kill only
on the arm's own executed evidence, and R7 requires a scope narrow enough that the evidence supports
it. A prose note saying "499 bp is less than 10,000 bp" is an ARGUMENT, not evidence.

So this RUNS. It loads every sequence-bearing rung through its canonical loader, measures the real
length distribution, and computes the span each rung can actually support. The arm's verdict then
rests on a measurement of the substrate it would have to run on -- which is the one kind of executed
evidence an unrunnable arm can produce.

WHAT "HONEST SPAN" MEANS, and why the floor is not a free parameter.
`dna_decode/glm/encoder.py` already settled this, and its reasoning is IMPORTED rather than restated:

    WHY NOT (1, 10, 100, 150), which would score 2.18 decades: a width-1 kernel cannot see a motif at
    all, so buying span by shrinking the finest branch to 1 bp inflates the decade metric without
    testing the hypothesis the metric stands for.

So the finest branch is read from `PYRAMID_WIDTHS[0]`, NOT written here as a literal. If someone
changes the encoder's finest scale, this probe follows it, and a 1 bp floor cannot be smuggled in to
make an arm look reachable.

WHAT WAS SEARCHED (R19: a substrate claim is a claim about the search).
Every rung with a canonical loader in `dna_decode.glm.genomewide`: fragments LB, fragments M9, peak
tiles. The probe REPORTS the set it covered and the rungs it could not reach, so a reader can tell
"measured infeasible" from "did not look" -- the distinction the 2026-10-26 memo got wrong.

Read-only. Touches no ledger. Exit 0 on a completed measurement whatever the verdict, 2 on a fault --
an infeasible verdict is a RESULT, not an error, and must not be confused with a broken run.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from dataclasses import dataclass, field
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from dna_decode.glm.encoder import PYRAMID_WIDTHS  # noqa: E402
from dna_decode.glm.genomewide import (  # noqa: E402
    DEFAULT_CACHE,
    GenomeWideDataError,
    load_fragments,
    load_peak_tiles,
)

#: The arm under test is named for FOUR decades. Not a tunable -- it is the arm's own declaration.
REQUIRED_DECADES = 4.0

#: The finest branch that can still see a motif, READ from the encoder rather than restated.
FINEST_HONEST_WIDTH = PYRAMID_WIDTHS[0]

ARM = "multiresolution-pyramid-4-decades-UNRUN"


class ProbeFault(RuntimeError):
    """A measurement could not be completed. Distinct from a completed INFEASIBLE verdict."""


def honest_span_decades(max_len: int, finest: int = FINEST_HONEST_WIDTH) -> float:
    """log10(max_len / finest). The span a pyramid can cover with every branch a real scale."""
    if max_len <= 0 or finest <= 0:
        raise ProbeFault(f"span undefined for max_len={max_len}, finest={finest}")
    if max_len < finest:
        raise ProbeFault(
            f"max_len {max_len} is below the finest honest width {finest}: no pyramid fits at all"
        )
    return math.log10(max_len / finest)


def required_bp(decades: float = REQUIRED_DECADES, finest: int = FINEST_HONEST_WIDTH) -> int:
    """How long a sequence must be for `decades` of span above `finest`."""
    return int(round(finest * (10.0 ** decades)))


@dataclass
class RungSpan:
    """One rung's measured length distribution and the span it permits."""

    rung: str
    n_sequences: int
    length_min: int
    length_median: int
    length_max: int
    n_distinct_lengths: int
    honest_span_decades: float
    feasible_for_4_decades: bool
    shortfall_factor: float | None
    source_file: str

    def to_dict(self) -> dict:
        d = {k: getattr(self, k) for k in (
            "rung", "n_sequences", "length_min", "length_median", "length_max",
            "n_distinct_lengths", "source_file", "feasible_for_4_decades", "shortfall_factor")}
        d["honest_span_decades"] = round(self.honest_span_decades, 4)
        return d


def _rung_from_stats(rung: str, n: int, st) -> RungSpan:
    span = honest_span_decades(st.length_max)
    need = required_bp()
    feasible = st.length_max >= need
    return RungSpan(
        rung=rung,
        n_sequences=n,
        length_min=st.length_min,
        length_median=st.length_median,
        length_max=st.length_max,
        n_distinct_lengths=st.n_distinct_lengths,
        honest_span_decades=span,
        feasible_for_4_decades=feasible,
        shortfall_factor=None if feasible else round(need / st.length_max, 2),
        source_file=Path(st.path).name,
    )


@dataclass
class ProbeReport:
    rungs: list = field(default_factory=list)
    unreachable: list = field(default_factory=list)

    @property
    def best(self) -> RungSpan | None:
        return max(self.rungs, key=lambda r: r.honest_span_decades) if self.rungs else None

    def verdict(self) -> str:
        """FEASIBLE / INFEASIBLE / INDETERMINATE.

        INDETERMINATE when NO rung could be measured -- that is the "did not look" case and must not
        read as infeasible. This is the distinction the prose version of this claim collapsed.
        """
        if not self.rungs:
            return "INDETERMINATE_NO_RUNG_MEASURED"
        if any(r.feasible_for_4_decades for r in self.rungs):
            return "FEASIBLE"
        return "INFEASIBLE_SUBSTRATE_GATED"

    def to_dict(self) -> dict:
        best = self.best
        return {
            "schema": "pyramid-span-feasibility-v1",
            "arm_under_test": ARM,
            "required_decades": REQUIRED_DECADES,
            "finest_honest_width_bp": FINEST_HONEST_WIDTH,
            "finest_width_source": "dna_decode.glm.encoder.PYRAMID_WIDTHS[0] (imported, not restated)",
            "required_sequence_length_bp": required_bp(),
            "verdict": self.verdict(),
            "rungs_measured": [r.to_dict() for r in self.rungs],
            "rungs_unreachable": self.unreachable,
            "best_achievable": None if best is None else {
                "rung": best.rung,
                "honest_span_decades": round(best.honest_span_decades, 4),
                "max_sequence_length_bp": best.length_max,
                "n_sequences": best.n_sequences,
            },
            "what_was_searched": (
                "every rung with a canonical loader in dna_decode.glm.genomewide: fragments LB, "
                "fragments M9, peak tiles. Reported so 'measured infeasible' is distinguishable "
                "from 'did not look'."
            ),
            "scope": (
                "Bounds the SUBSTRATE CACHED ON THIS HOST only. It is NOT a claim that no such "
                "dataset exists in the world, and NOT a claim about the arm's predictive merit -- "
                "the arm was never run, and this says it cannot be instantiated here, not that it "
                "would have failed."
            ),
            "revival_condition": (
                "any rung reaching >= %d bp becomes available -- a long-fragment or whole-locus "
                "expression assay, or a different label type (e.g. per-gene RNA-seq, where upstream "
                "context length is unbounded by the assay)." % required_bp()
            ),
        }


def probe(cache_dir: Path | str = DEFAULT_CACHE, *, limit: int | None = None,
          allow_download: bool = False) -> ProbeReport:
    """Measure every reachable rung. A rung that cannot load is RECORDED, never silently skipped."""
    rep = ProbeReport()

    for which, label in (("LB", "fragments_LB"), ("M9", "fragments_M9")):
        try:
            frags, st = load_fragments(cache_dir, which=which, limit=limit,
                                       allow_download=allow_download)
            rep.rungs.append(_rung_from_stats(label, len(frags), st))
        except (GenomeWideDataError, OSError) as exc:
            rep.unreachable.append({"rung": label, "reason": f"{type(exc).__name__}: {exc}"})

    # NOTE: load_peak_tiles takes NO `limit` kwarg (unlike load_fragments) -- verified by reading it,
    # not assumed from the sibling's signature. Passing one raises TypeError. The tile file is small
    # (46,713 rows) so it always loads whole, which makes its max length exact even in a smoke run.
    try:
        tiles, st = load_peak_tiles(cache_dir, allow_download=allow_download)
        rep.rungs.append(_rung_from_stats("peak_tiles", len(tiles), st))
    except (GenomeWideDataError, OSError) as exc:
        rep.unreachable.append({"rung": "peak_tiles", "reason": f"{type(exc).__name__}: {exc}"})

    return rep


def render(rep: ProbeReport, limit: int | None = None) -> str:
    """Render the report. `limit` is REQUIRED to be passed through.

    The first version of this took no `limit`, so a capped smoke printed
    `VERDICT: INFEASIBLE_SUBSTRATE_GATED` to stdout while the artifact correctly recorded
    `INDETERMINATE_ROW_LIMITED` -- two surfaces of one run disagreeing, with the louder one wrong.
    Measured on the real files: at `--limit 2000` the fragment max reads 324 bp against a true 499,
    so a capped run genuinely cannot ground this verdict.
    """
    out = [
        "pyramid span feasibility -- arm: %s" % ARM,
        "  required: %.2f decades  =>  %d bp at a %d bp finest branch"
        % (REQUIRED_DECADES, required_bp(), FINEST_HONEST_WIDTH),
        "",
        "  %-16s %9s %7s %7s %9s  %s"
        % ("rung", "n", "max bp", "span", "short by", "verdict"),
    ]
    for r in sorted(rep.rungs, key=lambda x: -x.honest_span_decades):
        out.append("  %-16s %9d %7d %7.4f %9s  %s" % (
            r.rung, r.n_sequences, r.length_max, r.honest_span_decades,
            "-" if r.shortfall_factor is None else "%.1fx" % r.shortfall_factor,
            "FEASIBLE" if r.feasible_for_4_decades else "infeasible",
        ))
    for u in rep.unreachable:
        out.append("  %-16s %s" % (u["rung"], "UNREACHABLE: " + u["reason"][:70]))
    best = rep.best
    out.append("")
    if limit is not None:
        out.append("  VERDICT: INDETERMINATE_ROW_LIMITED (--limit %d)" % limit)
        out.append("  a row cap can only UNDER-report max length, so it cannot ground a verdict;")
        out.append("  the uncapped verdict this run would otherwise print is %s" % rep.verdict())
    else:
        out.append("  VERDICT: %s" % rep.verdict())
    if best is not None:
        out.append("  best achievable here: %.4f decades on %s (max %d bp, n=%d)"
                   % (best.honest_span_decades, best.rung, best.length_max, best.n_sequences))
    return "\n".join(out)


def main(argv: list | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--cache-dir", default=str(DEFAULT_CACHE))
    ap.add_argument("--limit", type=int, default=None,
                    help="cap rows per rung (SMOKE ONLY -- a capped run can UNDER-report max length, "
                         "so it must never back a verdict; the artifact records the cap)")
    ap.add_argument("--allow-download", action="store_true")
    ap.add_argument("--out", default=None, help="write the JSON artifact here")
    a = ap.parse_args(argv)

    try:
        rep = probe(a.cache_dir, limit=a.limit, allow_download=a.allow_download)
    except ProbeFault as exc:
        print("PROBE FAULT: %s" % exc, file=sys.stderr)
        return 2

    print(render(rep, limit=a.limit))

    if not rep.rungs:
        print("\nno rung could be measured; refusing to write a verdict artifact", file=sys.stderr)
        return 2

    if a.out:
        payload = rep.to_dict()
        if a.limit is not None:
            payload["row_limit_applied"] = a.limit
            payload["verdict"] = "INDETERMINATE_ROW_LIMITED"
            payload["limit_note"] = (
                "a row cap can only UNDER-report the maximum length, so a capped run cannot "
                "establish infeasibility; re-run uncapped for a verdict."
            )
        Path(a.out).write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        print("\nwrote %s" % a.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
