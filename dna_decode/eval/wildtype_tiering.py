"""Wild-type / non-wild-type classification against an ECOFF.

WHY THIS IS A SEPARATE FUNCTION AND NOT A REUSE OF `mic_tiers.classify_tier`.
`classify_tier` answers a CLINICAL question -- will treatment succeed -- and its vocabulary reflects
that: R/S tiers with 4x safety margins, an intermediate BORDERLINE zone, and a CLSI-vs-EUCAST
disagreement check. An ECOFF answers a different question -- does this isolate carry an acquired
mechanism -- and its vocabulary is BINARY: wild-type or not, with no intermediate zone and no safety
margin.

Passing an ECOFF into `classify_tier` as `clsi_s` would typecheck, run, and return plausible tier
strings for every isolate. It would also be meaningless, and nothing downstream could tell. That is
the same failure shape as the log2-unit trap in the Oxford deposit, where applying clinical
breakpoints to dilution indices returned zero resistant isolates on all three drugs and read as a
clean, powered, all-susceptible cohort rather than as an error.

WT/NWT IS NOT S/R. A non-wild-type isolate is one whose MIC sits above the wild-type distribution;
that is evidence of a mechanism, NOT a clinical resistance call, and it must never be rendered to a
user as one. Callers that report these values are responsible for the distinction.

The median-of-MICs convention matches `mic_tiers.classify_tier` deliberately, so that an ECOFF arm
and a clinical arm scored on the same isolate differ ONLY in the anchor.

NOT FROZEN, and imports nothing from the frozen surface.
"""
from __future__ import annotations

from math import isnan
from statistics import median
from typing import Iterable, Optional

WT = "WT"            # at or below the ECOFF -- within the wild-type distribution
NWT = "NWT"          # above the ECOFF -- non-wild-type, i.e. mechanism present
NO_MIC = "NO_MIC"    # nothing numeric to classify


def classify_wildtype(mics: Iterable[Optional[float]], ecoff: float) -> str:
    """Classify an isolate as `WT`, `NWT` or `NO_MIC` against a single ECOFF.

    AT the ECOFF is wild-type: EUCAST defines the ECOFF as the HIGHEST MIC of the wild-type
    distribution, so the boundary value belongs to the wild type and only values strictly above it are
    non-wild-type. An off-by-one here silently reclassifies an entire dilution.
    """
    if ecoff is None or (isinstance(ecoff, float) and isnan(ecoff)) or ecoff <= 0:
        raise ValueError(f"ECOFF must be a positive number in mg/L; got {ecoff!r}")
    valid = [float(m) for m in (mics or [])
             if m is not None and not (isinstance(m, float) and isnan(m))]
    if not valid:
        return NO_MIC
    return NWT if median(valid) > ecoff else WT


def wildtype_counts(per_isolate: dict[str, Iterable[Optional[float]]], ecoff: float) -> dict[str, int]:
    """Convenience roll-up: {WT: n, NWT: n, NO_MIC: n} over a {isolate: mics} map."""
    out = {WT: 0, NWT: 0, NO_MIC: 0}
    for mics in per_isolate.values():
        out[classify_wildtype(mics, ecoff)] += 1
    return out


def ecoff_is_resolvable_on_grid(ecoff: float, measured_mics: Iterable[float]) -> bool:
    """False when the ECOFF sits at or below the lowest MIC the panel can report.

    In that regime EVERY isolate is non-wild-type and the anchor discriminates nothing -- the
    `DEGENERATE_ECOFF_BELOW_PANEL` outcome. This is a live risk rather than a theoretical one: the
    Oxford ciprofloxacin panel bottoms out at 0.125 mg/L.
    """
    vals = [float(m) for m in measured_mics
            if m is not None and not (isinstance(m, float) and isnan(m))]
    if not vals:
        return False
    return ecoff >= min(vals)
