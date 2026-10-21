"""POST-RUN validity gates for a sampling-sweep grid.

Deliberately separate from `knob_contract.py`. That module's job is PRECONDITIONS — refuse a
run before it spends 22 minutes. These are the gates applied AFTER the points come back, when
the question is which of them may enter a verdict. Mixing the two would make neither name true.

WHY (2026-10-20/21). The sampling-sweep plan's Step 4 says a point whose knob is unverified is
"excluded from the verdict rather than averaged in." Read literally that excludes the T=1.0
ANCHOR -- the comparator the sharpened points are measured against -- because verification is
gated on a non-default knob having been REQUESTED, and at T=1.0 none was. The real in-flight
anchor artifact shows the consequence exactly:

    "knob_contract": {"requested_non_default": false,
                      "knob_wired": false,
                      "knob_verified_end_to_end": false}

A default point is not UNVERIFIED. It is NOTHING-TO-VERIFY. Collapsing those two into one
boolean is the defect, and it is the same three-state distinction the doubt layer already
makes (not-applicable / not-assessable / assessed) -- a distinction this repo has lost before.

The plan's own Verification Signal gate 2 was already right ("at every NON-DEFAULT grid
point"); the rule simply lived in three places and disagreed with itself. So it lives here
once, and the prose cites it.
"""
from __future__ import annotations


class GridGateError(ValueError):
    """A grid that cannot be assembled into a verdict."""


# Returned as REASON STRINGS, not booleans. A bool is what collapsed "nothing to verify" into
# "unverified"; a caller that has to name which case it is in cannot make that mistake silently.
ELIGIBLE_DEFAULT = "eligible:default-sampling"
ELIGIBLE_VERIFIED = "eligible:verified"
EXCLUDED_UNVERIFIED = "excluded:requested-but-unverified"
EXCLUDED_KNOB_DEAD = "excluded:knob-dead"

ELIGIBLE = (ELIGIBLE_DEFAULT, ELIGIBLE_VERIFIED)


def verdict_eligibility(knob_record: dict) -> str:
    """Whether a grid point's sampling may enter the verdict, and why.

    Ordered; first match wins:

      1. no non-default knob was REQUESTED -> `eligible:default-sampling`
         Nothing to verify. This is the anchor, and it is the comparator -- excluding it would
         leave the sharpened points with nothing to be sharper than.
      2. requested but `knob_wired` false -> `excluded:knob-dead`
         Checked BEFORE the verification result, because a knob that installed no warper is
         inert whatever a verification says; this is the 2026-10-20 regression's own signature.
      3. requested, wired, verified     -> `eligible:verified`
      4. requested, wired, NOT verified -> `excluded:requested-but-unverified`

    A missing key is treated as falsy rather than raising: an older artifact predating the
    knob-contract stamp records default sampling, which is case 1 and is honestly eligible.
    """
    if not isinstance(knob_record, dict):
        raise GridGateError(f"knob_record must be a dict, got {type(knob_record).__name__}")

    requested = bool(knob_record.get("requested_non_default"))
    if not requested:
        return ELIGIBLE_DEFAULT
    if not bool(knob_record.get("knob_wired")):
        return EXCLUDED_KNOB_DEAD
    if bool(knob_record.get("knob_verified_end_to_end")):
        return ELIGIBLE_VERIFIED
    return EXCLUDED_UNVERIFIED


def null_reproduced(point_null, anchor_null, recorded_null=None, *, tol: float = 0.0) -> dict:
    """The cross-point invariant: temperature must not move the Markov null.

    Temperature touches only the GLM candidate, so every grid point's null must equal the
    anchor's. If it moves, that point is not a single-variable comparison and its GLM number is
    uninterpretable -- a VALIDITY failure of that point, never something to average in.

    HONEST SCOPE. The sweep plan already specifies a null cross-check, against the RECORDED
    0.5990 (`kmer`) / 0.6679 (`both`). Two things were missing and both are covered here: it was
    prose rather than a coded gate, and it compared only against the recorded baseline -- never
    anchor-vs-grid-point, which is the comparison that detects drift WITHIN one sweep.

    `tol=0.0` by default because the null is empirically deterministic: it reproduced to four
    decimals across two independent runs on different devices. A caller that wants slack must
    ask for it, so the strict case cannot be lost by omission.
    """
    out = {"point_null": point_null, "anchor_null": anchor_null,
           "recorded_null": recorded_null, "tol": tol}

    if point_null is None or anchor_null is None:
        out.update(ok=False, reason="null_missing")
        return out

    d_anchor = abs(float(point_null) - float(anchor_null))
    out["delta_vs_anchor"] = d_anchor
    if d_anchor > tol:
        out.update(ok=False, reason="null_moved_vs_anchor")
        return out

    if recorded_null is not None:
        d_rec = abs(float(point_null) - float(recorded_null))
        out["delta_vs_recorded"] = d_rec
        if d_rec > tol:
            out.update(ok=False, reason="null_moved_vs_recorded")
            return out

    out.update(ok=True, reason="null_reproduced")
    return out


def assemble_grid(points, *, anchor_null=None, recorded_null=None, tol: float = 0.0) -> dict:
    """Partition grid points into eligible / excluded, every exclusion carrying its reason.

    `points` is an iterable of dicts with at least `label` and `knob_contract`, optionally
    `null`. Nothing is ever dropped silently -- a point that leaves the verdict leaves a record
    saying which gate removed it, because an unexplained absence reads as "not run".

    The anchor's null is taken from the first point whose eligibility is `eligible:default-
    sampling` when `anchor_null` is not supplied, since that point IS the default-sampling run.
    """
    pts = list(points)
    eligible, excluded = [], []

    if anchor_null is None:
        for p in pts:
            if verdict_eligibility(p.get("knob_contract") or {}) == ELIGIBLE_DEFAULT:
                anchor_null = p.get("null")
                break

    for p in pts:
        label = p.get("label")
        why = verdict_eligibility(p.get("knob_contract") or {})
        if why not in ELIGIBLE:
            excluded.append({"label": label, "reason": why})
            continue
        if "null" in p or anchor_null is not None:
            nr = null_reproduced(p.get("null"), anchor_null, recorded_null, tol=tol)
            if not nr["ok"]:
                excluded.append({"label": label, "reason": f"excluded:{nr['reason']}",
                                 "null_check": nr})
                continue
            eligible.append({"label": label, "reason": why, "null_check": nr})
        else:
            eligible.append({"label": label, "reason": why})

    return {"n_points": len(pts), "eligible": eligible, "excluded": excluded,
            "anchor_null": anchor_null, "recorded_null": recorded_null, "tol": tol}
