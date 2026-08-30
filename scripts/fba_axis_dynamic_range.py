"""Does an axis's own DYNAMIC RANGE predict how flat the model's knockout ratios are on it?

The within-gene switch measurement produced a flatness ordering across three axes -- 4 media 61.2%,
25 carbon sources 68.2%, 13 nitrogen sources 75.5% -- and the nitrogen run was fired as a PRE-REGISTERED
prediction that it would be flattest, because six of its conditions were known to give identical wildtype
growth. The prediction held. But dynamic range had only ever been measured on nitrogen, so "flatness
tracks dynamic range" was an ordering with one measured predictor, not a relationship.

This measures the predictor on ALL THREE axes with one yardstick, and it is cheap: WILDTYPE growth per
condition only -- 4 + 25 + 13 = 42 LP solves, no deletion panel.

Two summaries per axis, because they answer different questions:
  * distinct_fraction -- how many conditions the model can tell apart AT ALL (exact ties)
  * cv               -- how much the growths actually spread (a tiny spread still counts as "distinct")

n=3 axes cannot establish a relationship. What it can do is say whether the pre-registered direction
survives a common yardstick, or whether nitrogen was flattest for some other reason.

Run: uv run python scripts/fba_axis_dynamic_range.py
"""
from __future__ import annotations

import json
import statistics
import sys
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

# Flatness measured by scripts/fba_within_gene_ranking.py on each axis, same metric, same day.
OBSERVED_FLATNESS = {"media4": 0.6119, "carbon": 0.6820, "nitrogen": 0.7548}
TOL = 1e-9


def _summarise(name: str, growths: list[float]) -> dict:
    distinct = {round(g, 9) for g in growths}
    mean = statistics.fmean(growths) if growths else 0.0
    sd = statistics.pstdev(growths) if len(growths) > 1 else 0.0
    return {"axis": name, "n_conditions": len(growths),
            "n_distinct_growths": len(distinct),
            "distinct_fraction": round(len(distinct) / len(growths), 4) if growths else None,
            "cv": round(sd / mean, 4) if mean > TOL else None,
            "min": round(min(growths), 5) if growths else None,
            "max": round(max(growths), 5) if growths else None,
            "flatness_of_knockout_ratios": OBSERVED_FLATNESS.get(name)}


def main() -> int:
    from dna_decode.fba.conditional_essentiality import CONDITIONS, apply_condition
    from dna_decode.fba.fitness_browser import apply_carbon_condition, carbon_conditions, open_db
    from dna_decode.fba.model import load_model, wildtype_growth
    from dna_decode.fba.nitrogen import apply_nitrogen_condition, nitrogen_conditions

    model = load_model()
    rows = []

    g = []
    for c in sorted(CONDITIONS):
        with model:
            apply_condition(model, c)
            g.append(wildtype_growth(model))
    rows.append(_summarise("media4", g))

    conn = open_db(None)
    for name, get_conds, apply_fn, kw in (
            ("carbon", carbon_conditions, apply_carbon_condition, "all_carbon"),
            ("nitrogen", nitrogen_conditions, apply_nitrogen_condition, "all_nitrogen")):
        conds = get_conds(conn, model)
        all_ex = tuple(conds.values())
        g = []
        for c in sorted(conds):
            with model:
                apply_fn(model, conds[c], **{kw: all_ex})
                g.append(wildtype_growth(model))
        rows.append(_summarise(name, g))

    rows.sort(key=lambda r: r["flatness_of_knockout_ratios"] or 0)
    print(f"{'axis':10} {'n':>3} {'distinct':>9} {'frac':>7} {'cv':>8}   flatness")
    for r in rows:
        print(f"{r['axis']:10} {r['n_conditions']:>3} {r['n_distinct_growths']:>9} "
              f"{r['distinct_fraction']:>7} {str(r['cv']):>8}   {r['flatness_of_knockout_ratios']}")

    fl = [r["flatness_of_knockout_ratios"] for r in rows]
    frac = [r["distinct_fraction"] for r in rows]
    cvs = [r["cv"] for r in rows]
    mono_frac = all(a > b for a, b in zip(frac, frac[1:]))     # flatness up => distinct_fraction down
    mono_cv = all(a > b for a, b in zip(cvs, cvs[1:])) if all(c is not None for c in cvs) else None

    out = {"schema": "fba-axis-dynamic-range-v1", "generated": date.today().isoformat(),
           "axes": rows,
           "flatness_ascending": fl,
           "distinct_fraction_monotonically_decreasing": mono_frac,
           "cv_monotonically_decreasing": mono_cv,
           "verdict": ("dynamic range tracks flatness on both summaries" if (mono_frac and mono_cv)
                       else "dynamic range tracks flatness on one summary only" if (mono_frac or mono_cv)
                       else "dynamic range does NOT track flatness"),
           "honest_limit": ("n=3 axes. This cannot establish a relationship; it says only whether the "
                            "pre-registered direction survives a common yardstick. The three axes also "
                            "differ in substrate, label source and gene set, any of which could drive "
                            "flatness instead.")}
    (ROOT / "wiki" / f"fba_axis_dynamic_range_{out['generated']}.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8")
    print(f"\n  distinct_fraction falls as flatness rises : {mono_frac}")
    print(f"  cv falls as flatness rises               : {mono_cv}")
    print(f"  VERDICT: {out['verdict']}")
    print(f"  LIMIT: n=3 axes, and they differ in substrate/labels/gene set too.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
