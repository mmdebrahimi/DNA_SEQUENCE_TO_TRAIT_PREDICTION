"""External-cohort MIC -> tiered R/S labels.

The frozen decoder is scored against MIC-derived binary labels. Naive
`MIC >= breakpoint -> R` thresholding is WRONG: it silently forces
intermediate / borderline / censored MICs into R/S and corrupts sens/spec
(the brainstorm C3 finding). Instead we route every isolate's MIC through
`mic_tiers.classify_tier` (the project's label-quality authority) and score
only DECISIVE tiers, reporting the excluded buckets so the dropped fraction
is visible BEFORE any sens/spec is computed:

  STRICT  pass = HIGH_R / HIGH_S            (4x safety margin; primary metric)
  RELAXED pass = + DECISIVE_R / DECISIVE_S  (clearly R/S, no margin; secondary)
  EXCLUDED     = BORDERLINE / AMBIGUOUS / CONFLICT / NO_MIC

Drug aliases (CIP/CRO/CN/…) are normalized to the canonical mic_tiers names
before any breakpoint lookup or filename. Drugs outside the cipro/cef/gent
pilot are rejected.
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

from dna_decode.data.mic_tiers import breakpoints_for, classify_tier

_LOWER_OPS = (">=", ">")   # lower bound -> true value is AT OR ABOVE
_UPPER_OPS = ("<=", "<")   # upper bound -> true value is AT OR BELOW

# Breakpoint provenance (matches dna_decode/data/mic_tiers.py DRUG_BREAKPOINTS).
BREAKPOINT_VERSION = "CLSI 2024 + EUCAST 14.0 (E. coli)"

# Pilot drugs only. Aliases -> canonical mic_tiers drug name.
CANONICAL_DRUG: dict[str, str] = {
    "ciprofloxacin": "ciprofloxacin", "cipro": "ciprofloxacin", "cip": "ciprofloxacin",
    "ceftriaxone": "ceftriaxone", "cro": "ceftriaxone", "cef": "ceftriaxone", "axo": "ceftriaxone",
    "gentamicin": "gentamicin", "gen": "gentamicin", "gent": "gentamicin", "cn": "gentamicin", "gm": "gentamicin",
}

PILOT_DRUGS = frozenset({"ciprofloxacin", "ceftriaxone", "gentamicin"})

STRICT_TIERS = {"HIGH_R": "R", "HIGH_S": "S", "CENSORED_HIGH_R": "R", "CENSORED_HIGH_S": "S"}
RELAXED_EXTRA = {"DECISIVE_R": "R", "DECISIVE_S": "S"}
EXCLUDED_TIERS = ("BORDERLINE", "AMBIGUOUS", "CONFLICT", "NO_MIC", "CENSORED_EXCLUDED")


def canonical_drug(name: str) -> str | None:
    """Normalize a cohort drug alias to a canonical pilot drug name, or None to skip."""
    if not name:
        return None
    return CANONICAL_DRUG.get(name.strip().lower())


def parse_mic_token(token) -> float | None:
    """Parse a MIC cell tolerant of censored operators + units.

    Strips a leading comparison operator (`>`, `>=`, `<`, `<=`, `=`) and a
    trailing unit (mg/L, ug/mL, µg/mL). The censored value is returned as its
    NUMERIC BOUND (a documented approximation; callers keep the raw token).
    Returns None for empty / NA / unparseable cells.
    """
    if token is None:
        return None
    if isinstance(token, (int, float)):
        return float(token)
    s = str(token).strip()
    if not s or s.upper() in ("NA", "N/A", "NULL", "NONE", "."):
        return None
    for op in (">=", "<=", ">", "<", "="):
        if s.startswith(op):
            s = s[len(op):].strip()
            break
    # strip a trailing unit
    for unit in ("mg/l", "ug/ml", "µg/ml", "mg/L", "ug/mL"):
        if s.lower().endswith(unit):
            s = s[: -len(unit)].strip()
            break
    try:
        return float(s)
    except ValueError:
        return None


@dataclass(frozen=True)
class MicValue:
    """A parsed MIC cell that PRESERVES the censoring operator (vs parse_mic_token,
    which discards it). `operator` is one of '', '=', '>', '>=', '<', '<='."""
    value: float
    operator: str
    raw: str


def parse_mic_value(token) -> MicValue | None:
    """Parse a MIC cell into a MicValue (numeric bound + operator + raw token).

    Same numeric tolerance as parse_mic_token (units, NA) but RETAINS the operator
    so the censor DIRECTION can be modeled by tier_for_isolate. None if unparseable.
    """
    if token is None:
        return None
    if isinstance(token, (int, float)):
        return MicValue(float(token), "=", str(token))
    raw = str(token).strip()
    s = raw
    if not s or s.upper() in ("NA", "N/A", "NULL", "NONE", "."):
        return None
    op = ""
    for cand in (">=", "<=", ">", "<", "="):
        if s.startswith(cand):
            op = cand
            s = s[len(cand):].strip()
            break
    for unit in ("mg/l", "ug/ml", "µg/ml", "mg/L", "ug/mL"):
        if s.lower().endswith(unit):
            s = s[: -len(unit)].strip()
            break
    try:
        return MicValue(float(s), op or "=", raw)
    except ValueError:
        return None


def _censored_tier(values: list[MicValue], bps) -> str:
    """Operator-aware tier for an all-censored isolate (no plain `=` value).

    A lower bound (`>`/`>=`) can only support R, and only if its bound itself lands
    HIGH_R (true value is >= bound, so at least that resistant). An upper bound
    (`<`/`<=`) can only support S, and only if its bound lands HIGH_S. Anything else
    (mid-range bound, contradictory directions, or DECISIVE-only) is interval-censored
    -> CENSORED_EXCLUDED. NEVER calls R from an upper bound or S from a lower bound.
    """
    lowers = [v.value for v in values if v.operator in _LOWER_OPS]
    uppers = [v.value for v in values if v.operator in _UPPER_OPS]
    if lowers and not uppers:
        # most-resistant lower bound; classify the bound itself
        if classify_tier([max(lowers)], set(), bps) == "HIGH_R":
            return "CENSORED_HIGH_R"
        return "CENSORED_EXCLUDED"
    if uppers and not lowers:
        if classify_tier([min(uppers)], set(), bps) == "HIGH_S":
            return "CENSORED_HIGH_S"
        return "CENSORED_EXCLUDED"
    return "CENSORED_EXCLUDED"   # both directions or neither -> ambiguous


def tier_for_isolate(mic_tokens, distinct_calls, drug: str) -> str:
    """Classify one isolate's MICs for `drug` into a mic_tiers tier string.

    Plain (`=`) values go through classify_tier (the existing path). When an isolate
    has ONLY censored values, censor direction is modeled by _censored_tier so an
    upper bound can't be falsely called R (and a lower bound can't be falsely S).
    """
    bps = breakpoints_for(drug)
    calls = {str(c).strip().upper() for c in (distinct_calls or set()) if str(c).strip()}
    # contradictory categorical calls dominate (matches classify_tier's CONFLICT).
    if (calls & {"R", "RESISTANT"}) and (calls & {"S", "SUSCEPTIBLE"}):
        return "CONFLICT"
    parsed = [p for p in (parse_mic_value(t) for t in mic_tokens) if p is not None]
    plain = [p.value for p in parsed if p.operator in ("", "=")]
    if plain:
        return classify_tier(plain, calls, bps)        # existing numeric path
    if not parsed:
        return classify_tier([], calls, bps)           # NO_MIC
    return _censored_tier(parsed, bps)                 # all-censored, operator-aware


def build_drug_labels(isolate_to_mics: dict[str, list], drug: str,
                      isolate_to_calls: dict[str, set] | None = None) -> dict:
    """Tier-label every isolate for `drug`. Returns strict/relaxed label maps + buckets.

    `isolate_to_mics[acc]` is a list of raw MIC tokens (one isolate may have
    multiple AST rows). `isolate_to_calls[acc]` is an optional set of categorical
    R/S calls (feeds classify_tier's CONFLICT detection).
    """
    canon = canonical_drug(drug)
    if canon is None or canon not in PILOT_DRUGS:
        raise ValueError(f"drug {drug!r} is not a pilot drug ({sorted(PILOT_DRUGS)})")
    isolate_to_calls = isolate_to_calls or {}
    strict: dict[str, str] = {}
    relaxed: dict[str, str] = {}
    buckets: dict[str, int] = {}
    for acc, tokens in isolate_to_mics.items():
        tier = tier_for_isolate(tokens, isolate_to_calls.get(acc, set()), canon)
        buckets[tier] = buckets.get(tier, 0) + 1
        if tier in STRICT_TIERS:
            strict[acc] = STRICT_TIERS[tier]
            relaxed[acc] = STRICT_TIERS[tier]
        elif tier in RELAXED_EXTRA:
            relaxed[acc] = RELAXED_EXTRA[tier]
        # EXCLUDED tiers contribute only to buckets
    return {
        "drug": canon,
        "breakpoint_version": BREAKPOINT_VERSION,
        "strict": strict,
        "relaxed": relaxed,
        "buckets": buckets,
        "n_total": len(isolate_to_mics),
        "n_strict": len(strict),
        "n_relaxed": len(relaxed),
    }


def _write_selected(path: Path, labels: dict[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(f"{a}\t{rs}\n" for a, rs in sorted(labels.items())), encoding="utf-8")


def write_labels(out_dir: str | Path, result: dict) -> dict:
    """Write selected_strict.tsv + selected_relaxed.tsv + buckets_<drug>.json. Returns paths."""
    out_dir = Path(out_dir)
    drug = result["drug"]
    strict_p = out_dir / "selected_strict.tsv"
    relaxed_p = out_dir / "selected_relaxed.tsv"
    buckets_p = out_dir / f"buckets_{drug}.json"
    _write_selected(strict_p, result["strict"])
    _write_selected(relaxed_p, result["relaxed"])
    buckets_p.parent.mkdir(parents=True, exist_ok=True)
    buckets_p.write_text(json.dumps(
        {"drug": drug, "breakpoint_version": result["breakpoint_version"],
         "buckets": result["buckets"], "n_total": result["n_total"],
         "n_strict": result["n_strict"], "n_relaxed": result["n_relaxed"]},
        indent=2), encoding="utf-8")
    return {"strict": str(strict_p), "relaxed": str(relaxed_p), "buckets": str(buckets_p)}
