"""Is a decay-curve bin a DISTANCE bin, or is it just one held-out GROUP wearing a distance label?

    uv run python scripts/decay_curve_confound.py --self-check
    uv run python scripts/decay_curve_confound.py --axis Minus10 --seed 0

**WHY THIS EXISTS.** `scripts/oracle_decay_curve.py` reported two scored bins declining 0.8004 ->
0.5926 and its frozen monotonicity function duly returned `declines_monotonically`. Decomposing those
bins by their held-out element value showed bin 0-1 was 100.0% the SINGLE-substitution group and bin
2-3 was 98.6% the DOUBLE -- so the two bins were two GROUPS, not two distances, and the decline could
not be attributed to distance at all.

That decomposition originally lived in a throwaway script. The verdict
`confounded:distance_collinear_with_mutation_order` therefore rested on code nobody could re-run,
which is exactly the shape this repo refuses elsewhere: a claim whose evidence cannot be re-derived
from the artifact. This module is that decomposition, committed, with the arithmetic in a PURE
function so it is testable with no substrate, no checkpoint and no torch.

**WHAT IT DOES AND DOES NOT DECIDE.** It measures COLLINEARITY between the distance axis and a
grouping variable, and it reports whether the confound is BREAKABLE -- i.e. whether any single group
spans two distances at n >= the floor, which is the only way to get a within-group distance contrast.
It does NOT decide that distance has no effect; a confounded design is uninformative about the
question, not evidence for the null.

**THE FLOOR IS AN INPUT, NEVER TUNED HERE.** It is passed in from the caller's pre-registration.
Lowering it until a within-group cell qualifies is the move this module exists to make visible, so
`separability` is reported against the floor it was given and the floor is echoed into the result.
"""
from __future__ import annotations

import argparse
import json
from collections import Counter
from datetime import date
from pathlib import Path

VERDICTS = (
    "confounded:collinear",            # no group spans 2 distances at the floor
    "separable:within-group-possible",  # at least one group does
    "indeterminate:single-group",       # only one group present -- nothing to be collinear with
    "indeterminate:no-scored-bins",     # nothing cleared the floor
)


def decompose(
    distances,
    labels,
    *,
    bins,
    min_bin_n: int,
) -> dict:
    """PURE. Cross-tab distance against a grouping label and report whether the confound breaks.

    `distances` and `labels` are parallel sequences. `bins` is the caller's pre-registered
    ((lo,hi),...) and `min_bin_n` its floor. Nothing here reads a file or a model.
    """
    distances = [int(d) for d in distances]
    labels = [str(x) for x in labels]
    if len(distances) != len(labels):
        raise ValueError(
            "distances and labels must be parallel: got %d vs %d -- a misaligned zip would "
            "silently attribute one variant's distance to another's group"
            % (len(distances), len(labels))
        )

    groups = sorted(set(labels))
    seen_d = sorted(set(distances))
    pair = Counter(zip(labels, distances))

    cross_tab = {g: {"d=%d" % d: pair.get((g, d), 0) for d in seen_d} for g in groups}

    # Composition of every pre-registered bin, and which of them cleared the floor.
    bin_rows = []
    for k, (lo, hi) in enumerate(bins):
        idx = [i for i, d in enumerate(distances) if lo <= d <= hi]
        n = len(idx)
        comp = Counter(labels[i] for i in idx)
        row = {
            "bin": "%d-%d" % (lo, hi),
            "bin_index": k,
            "n": n,
            "scored": n >= min_bin_n,
            "composition": {
                g: {"n": c, "pct_of_bin": round(100.0 * c / n, 1)} for g, c in sorted(comp.items())
            } if n else {},
        }
        # The headline number: how dominated is this bin by its single commonest group?
        row["dominant_group"] = comp.most_common(1)[0][0] if comp else None
        row["dominant_pct"] = round(100.0 * comp.most_common(1)[0][1] / n, 1) if n else None
        bin_rows.append(row)

    # Can the confound be broken? Only if some group spans >= 2 distances at the floor.
    within = {}
    for g in groups:
        per_d = Counter(d for d, l in zip(distances, labels) if l == g)
        at_floor = sorted(d for d, c in per_d.items() if c >= min_bin_n)
        within[g] = {
            "cells": {"d=%d" % d: c for d, c in sorted(per_d.items())},
            "distances_at_or_above_floor": at_floor,
            "separable": len(at_floor) >= 2,
        }

    scored = [r for r in bin_rows if r["scored"]]
    any_separable = any(w["separable"] for w in within.values())

    if not scored:
        verdict = "indeterminate:no-scored-bins"
    elif len(groups) < 2:
        verdict = "indeterminate:single-group"
    elif any_separable:
        verdict = "separable:within-group-possible"
    else:
        verdict = "confounded:collinear"

    assert verdict in VERDICTS, "verdict escaped the closed enum"

    # Max dominance across scored bins -- one number saying how bad the collinearity is.
    worst = max((r["dominant_pct"] for r in scored if r["dominant_pct"] is not None), default=None)

    return {
        "verdict": verdict,
        "verdict_enum_closed": list(VERDICTS),
        "min_bin_n_supplied_by_caller": min_bin_n,
        "bins_supplied_by_caller": [list(b) for b in bins],
        "n_total": len(distances),
        "n_groups": len(groups),
        "groups": groups,
        "distances_observed": seen_d,
        "cross_tab": cross_tab,
        "bins": bin_rows,
        "n_scored_bins": len(scored),
        "max_dominant_pct_across_scored_bins": worst,
        "within_group_separation": within,
        "reading": _reading(verdict, scored, worst),
        "what_this_does_not_establish": (
            "A confounded design is UNINFORMATIVE about the distance question -- it is not evidence "
            "that distance has no effect. Nor does a separable verdict mean the within-group contrast "
            "was actually computed; it means one is possible at the supplied floor."),
    }


def _reading(verdict: str, scored: list, worst) -> str:
    if verdict == "indeterminate:no-scored-bins":
        return "No bin cleared the floor, so there is nothing to decompose."
    if verdict == "indeterminate:single-group":
        return ("Only one group is present, so the distance axis cannot be collinear with a grouping "
                "variable here. This is NOT a clean bill of health for other confounds.")
    if verdict == "separable:within-group-possible":
        return ("At least one group spans two distances at the floor, so a within-group distance "
                "contrast is computable and would separate distance from group.")
    return ("No group spans two distances at the floor, so no within-group distance contrast exists. "
            "The scored bins are group contrasts wearing distance labels"
            + (" (worst bin is %.1f%% one group)." % worst if worst is not None else "."))


def _self_check() -> int:
    """Exercise the pure function on synthetic data. No substrate, no torch."""
    ok = True
    BINS = ((0, 1), (2, 3), (4, 7))

    # Perfectly collinear: group A only at d=1, group B only at d=2.
    d = [1] * 40 + [2] * 40
    g = ["A"] * 40 + ["B"] * 40
    r = decompose(d, g, bins=BINS, min_bin_n=30)
    print("  perfectly collinear      -> %s" % r["verdict"])
    ok &= r["verdict"] == "confounded:collinear"
    ok &= r["max_dominant_pct_across_scored_bins"] == 100.0
    ok &= not any(w["separable"] for w in r["within_group_separation"].values())

    # Separable: group A spans d=1 AND d=2, both above the floor.
    d = [1] * 40 + [2] * 40
    g = ["A"] * 80
    r = decompose(d, g, bins=BINS, min_bin_n=30)
    print("  one group, two distances -> %s" % r["verdict"])
    ok &= r["verdict"] == "indeterminate:single-group"

    d = [1] * 40 + [2] * 40 + [1] * 40
    g = ["A"] * 80 + ["B"] * 40
    r = decompose(d, g, bins=BINS, min_bin_n=30)
    print("  A spans 2 distances      -> %s" % r["verdict"])
    ok &= r["verdict"] == "separable:within-group-possible"
    ok &= r["within_group_separation"]["A"]["separable"] is True
    ok &= r["within_group_separation"]["B"]["separable"] is False

    # Floor is respected: the same data with a higher floor stops being separable.
    r = decompose(d, g, bins=BINS, min_bin_n=50)
    print("  same data, floor 50      -> %s" % r["verdict"])
    ok &= r["verdict"] == "confounded:collinear"

    # Nothing clears the floor.
    r = decompose([1] * 5, ["A"] * 5, bins=BINS, min_bin_n=30)
    print("  below floor everywhere   -> %s" % r["verdict"])
    ok &= r["verdict"] == "indeterminate:no-scored-bins"

    # Misaligned input must raise, not silently zip short.
    try:
        decompose([1, 2, 3], ["A", "B"], bins=BINS, min_bin_n=1)
        print("  misaligned input did NOT raise -- BAD"); ok = False
    except ValueError:
        print("  misaligned input raises, as intended")

    print("SELF-CHECK %s" % ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--self-check", action="store_true", help="exercise the pure function and exit")
    ap.add_argument("--cache-dir", default="data/mpra_gse108535")
    ap.add_argument("--axis", default="Minus10")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    if args.self_check:
        return _self_check()

    # Imported here so --self-check needs none of this.
    from scripts.oracle_decay_curve import (
        BINS,
        MIN_BIN_N,
        encode,
        hamming_min_to_set,
        leave_element_out,
        load_pairs,
    )

    pairs, _rep = load_pairs(args.cache_dir, allow_download=False)
    train, test, held = leave_element_out(pairs, args.axis, seed=args.seed)
    print("%s: held %s -> train %d / test %d" % (args.axis, held, len(train), len(test)))

    dists = hamming_min_to_set(encode([p.seq for p in test]), encode([p.seq for p in train]))
    labels = [p.elements[args.axis] for p in test]

    res = decompose(dists.tolist(), labels, bins=BINS, min_bin_n=MIN_BIN_N)

    art = {
        "schema": "decay-curve-confound-v1",
        "date": str(date.today()),
        "question": ("are the decay-curve's scored bins DISTANCE bins, or single held-out groups "
                     "wearing distance labels?"),
        "substrate": {
            "dataset": "GEO GSE108535 (Urtecho/Kosuri sigma70 variants)",
            "axis": args.axis, "seed": args.seed, "held_out_values": held,
            "n_train": len(train), "n_test": len(test),
        },
        "floor_and_bins_are_INPUTS": (
            "Both come from scripts/oracle_decay_curve.py, which fixed them before any data was read. "
            "They are echoed into the result so a later run cannot quietly use a different floor."),
        "result": res,
    }
    out = Path(args.out) if args.out else Path("wiki/decay_curve_confound_%s.json" % date.today())
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(art, indent=2), encoding="utf-8")

    print("VERDICT: %s" % res["verdict"])
    print("  %s" % res["reading"])
    for r in res["bins"]:
        if r["n"]:
            print("  bin %-6s n=%-5d scored=%-5s dominant=%s at %.1f%%" % (
                r["bin"], r["n"], r["scored"], r["dominant_group"], r["dominant_pct"]))
    for g, w in res["within_group_separation"].items():
        print("  within %-18s cells=%s at_floor=%s separable=%s" % (
            g, w["cells"], w["distances_at_or_above_floor"], w["separable"]))
    print("wrote %s" % out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
