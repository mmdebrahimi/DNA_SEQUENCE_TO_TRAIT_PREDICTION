"""Which eta^2 governs the pooling gain -- the LABEL's, the PREDICTOR's, or neither separably?

WHY THIS EXISTS. `wiki/forward_epistasis_h2_confirmed_2026-08-25.md` claimed H2 CONFIRMED: that eta^2(k),
the share of FITNESS variance sitting between mutation orders, governs the additive score's pooling gain.
It correlated the gain against eta^2 computed on FITNESS ONLY.

But a pooling gain needs BOTH the label and the predictor to separate by group. If the additive score had
identical distributions across orders, pooling groups that differ in fitness MEAN would add no rank
correlation at all -- the gain would be zero however large eta^2(fitness) grew. So eta^2(fitness) is
NECESSARY, not obviously SUFFICIENT, and the score-side term was never measured.

This measures it. Per order-subset it computes eta^2 on fitness AND eta^2 on the additive score, their
collinearity, and the partial Spearman of gain~eta^2_fitness controlling for eta^2_score -- plus a
LEAVE-ONE-ORDER-OUT sweep, because the effective sample size here is the NUMBER OF ORDERS (5), not the
number of subsets (26): the subsets are nested and share most of their variants.

WHAT IT FOUND (2026-08-25). The two eta^2 are near-collinear (rho +0.98 / +0.999 / +0.60), so their
contributions are NOT separably identifiable on this data. On GFP the fitness-side partial collapses to
~+0.02 once the score side is controlled; on ParD the score side is a perfect rank predictor while fitness
manages +0.600; only HIS7 leaves fitness-side signal standing. The supportable claim is therefore the JOINT
one -- the gain is governed by aligned between-order separation in label AND predictor -- and NOT "eta^2 of
fitness governs it". Swapping in "eta^2 of the score governs it" would just be the next overclaim.

IT ALSO GROUNDS THE MECHANISM. The additive score is a SUM of k per-mutation log-ratios, so its mean scales
with k BY CONSTRUCTION -> its between-order separation is structurally nonzero. The joint score (one forward
pass on the multi-mutant) has no such construction. That is why the sweep saw the additive score's pooling
gain swing 4x across proteins while the joint score's stayed near-constant -- the asymmetry the correction
memo flagged but could not explain. `mean_score_by_order` reports the scaling directly.

LIMITATION, STATED. The joint score's own eta^2 CANNOT be computed here: only aggregate per-order Spearman
values were persisted by the 2026-07-27 sweep, not per-variant joint scores. The mechanism above is
therefore grounded on the additive side and INFERRED on the joint side.

ZERO NEW COMPUTE -- cached ESM2 per-position log-prob matrices + cached assay CSVs. No GPU, no network.

Run: uv run python scripts/epistasis_eta2_identifiability.py
"""
from __future__ import annotations

import itertools
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

MIN_SUBSETS = 4          # below this a within-protein correlation is not worth reporting
OUT = ROOT / "wiki" / "forward_epistasis_eta2_identifiability_2026-08-25.json"


def eta2_of(groups: list[list[float]]) -> float:
    """Between-group share of total variance. PURE. Re-exported shape of the H2 script's `eta2`."""
    from epistasis_pooling_h2_test import eta2
    return eta2(groups)


def score_groups(by_order: dict, orders) -> list[list[float]]:
    """The additive SCORE values, grouped by mutation order. PURE."""
    return [[s for _, s in by_order[k]] for k in orders]


def fitness_groups(by_order: dict, orders) -> list[list[float]]:
    """The FITNESS values, grouped by mutation order. PURE."""
    return [[f for f, _ in by_order[k]] for k in orders]


def partial_spearman(r_xy: float, r_zy: float, r_xz: float) -> float:
    """Partial rank correlation of x~y controlling z, from the three pairwise rhos. PURE.

    Returns nan when the controlling variable explains essentially all of either margin -- which is the
    HONEST output under near-collinearity, not a failure. Do not paper over it with a fallback value.
    """
    den = ((1.0 - r_zy ** 2) * (1.0 - r_xz ** 2)) ** 0.5
    if den <= 1e-12:
        return float("nan")
    return (r_xy - r_zy * r_xz) / den


def subset_rows(by_order: dict) -> list[dict]:
    """One row per subset of mutation orders, carrying BOTH eta^2 values and the pooling gain. PURE."""
    from epistasis_pooling_h2_test import subset_stats
    orders = sorted(by_order)
    rows = []
    for size in range(2, len(orders) + 1):
        for combo in itertools.combinations(orders, size):
            st = subset_stats(by_order, combo)
            if st is None:
                continue
            st["eta2_fitness"] = st.pop("eta2")
            st["eta2_score"] = round(eta2_of(score_groups(by_order, combo)), 4)
            rows.append(st)
    return rows


def correlate(rows: list[dict]) -> dict:
    """rho(gain, eta2_fitness), rho(gain, eta2_score), their collinearity, and both partials. PURE."""
    from scipy.stats import spearmanr
    g = [r["pooling_gain"] for r in rows]
    ef = [r["eta2_fitness"] for r in rows]
    es = [r["eta2_score"] for r in rows]
    r_fit = float(spearmanr(ef, g).statistic)
    r_score = float(spearmanr(es, g).statistic)
    r_collin = float(spearmanr(ef, es).statistic)
    return {
        "n_subsets": len(rows),
        "rho_gain_eta2_fitness": round(r_fit, 4),
        "rho_gain_eta2_score": round(r_score, 4),
        "rho_eta2_fitness_eta2_score": round(r_collin, 4),
        "partial_fitness_given_score": round(partial_spearman(r_fit, r_score, r_collin), 4),
        "partial_score_given_fitness": round(partial_spearman(r_score, r_fit, r_collin), 4),
    }


def leave_one_order_out(by_order: dict) -> list[dict]:
    """Re-correlate using only subsets that EXCLUDE each order in turn. PURE.

    This is the honest resampling unit. The 26 subsets are nested and share most of their variants, so a
    p-value over them assumes an independence they do not have; the effective sample size is the number of
    ORDERS. Dropping one order at a time is order-level sensitivity, not variant-level significance.
    """
    orders = sorted(by_order)
    out = []
    for drop in orders:
        kept = {k: v for k, v in by_order.items() if k != drop}
        if len(kept) < 2:
            continue
        rows = subset_rows(kept)
        if len(rows) < MIN_SUBSETS:
            out.append({"dropped_order": drop, "status": "too_few_subsets", "n_subsets": len(rows)})
            continue
        c = correlate(rows)
        c["dropped_order"] = drop
        c["status"] = "ok"
        out.append(c)
    return out


def order_slates(sizes: dict) -> dict:
    """Alternative 5-order slates for a protein that carries more than 5 orders. PURE.

    `epistasis_pooling_h2_test.MAX_ORDERS=5` keeps the MOST POPULOUS orders. That began as a fix for a
    runtime blow-up (GFP carries orders 2..12 -> 2036 subsets) and was THEN justified as "also the
    statistically right call". A post-hoc justification of a performance fix deserves a falsifier: it is
    selection on SAMPLE SIZE, and per-order Spearman precision scales with n.
    """
    by_pop = sorted(sizes, key=lambda k: -sizes[k])
    by_k = sorted(sizes)
    return {
        "most_populous_5_SHIPPED": sorted(by_pop[:5]),
        "least_populous_5": sorted(by_pop[-5:]),
        "every_other_k": by_k[::2][:5],
    }


def slate_sensitivity(by_all: dict) -> dict | None:
    """Re-correlate on alternative order slates -- does popularity-selection drive the result? PURE.

    MEASURED on GFP (2026-08-25): the relationship stays POSITIVE on every slate but is markedly weaker on
    the all-sparse tail -- most-populous [2-6] n=45623 rho +0.973, least-populous [8-12] n=2451 rho +0.490,
    every-other-k [2,4,6,8,10] n=28125 rho +0.963. `every_other_k` is the discriminating slate: it INCLUDES
    sparse orders yet holds, so the weakening concentrates where EVERY order is sparse. Two readings remain
    entangled -- genuinely weaker at high k, vs noisier at low n -- and this design cannot separate them.
    Note leave-one-order-out never reaches this regime: it drops ONE order from the populous slate.
    """
    sizes = {k: len(v) for k, v in by_all.items()}
    if len(sizes) <= 5:
        return None
    out = {"order_sizes": {str(k): v for k, v in sorted(sizes.items())}, "slates": {}}
    for label, orders in order_slates(sizes).items():
        sub = {k: by_all[k] for k in orders if k in by_all}
        if len(sub) < 3:
            continue
        rows = subset_rows(sub)
        if len(rows) < MIN_SUBSETS:
            continue
        c = correlate(rows)
        c["orders"] = sorted(sub)
        c["n_variants"] = sum(len(v) for v in sub.values())
        out["slates"][label] = c
    return out


def mean_score_by_order(by_order: dict) -> dict:
    """Mean additive score per mutation order -- the STRUCTURAL claim, checkable directly. PURE.

    An additive score is a sum of k per-mutation log-ratios (each typically negative), so its mean should
    fall roughly linearly in k. That construction is why the additive score separates by order at all, and
    why the joint score need not.
    """
    return {str(k): round(sum(s for _, s in v) / len(v), 4) for k, v in sorted(by_order.items())}


def main() -> int:
    from epistasis_pooling_h2_test import PROTEINS, load_protein
    report = {"_schema": "epistasis-eta2-identifiability-v1", "proteins": {},
              "honest_scope": (
                  "The two eta^2 are near-collinear, so label-side and predictor-side contributions are "
                  "NOT separably identifiable on this data. The supportable claim is the JOINT condition "
                  "(aligned between-order separation in BOTH label and predictor), not that fitness eta^2 "
                  "governs the gain -- and not that score eta^2 does either. Effective sample size is the "
                  "number of ORDERS, so leave-one-order-out is reported instead of a subset p-value. The "
                  "joint score's own eta^2 is NOT computable here: the 2026-07-27 sweep persisted only "
                  "aggregate per-order rho, not per-variant joint scores.")}

    print("WHICH eta^2 GOVERNS THE POOLING GAIN? (label-side, predictor-side, or not separable)\n")
    hdr = (f"{'protein':26s} {'sub':>4s} {'rho(g,FIT)':>11s} {'rho(g,SCORE)':>13s} "
           f"{'collin':>8s} {'part FIT|sc':>12s} {'part SCORE|fit':>15s}")
    print(hdr)

    import epistasis_pooling_h2_test as H2
    for name in PROTEINS:
        by_order = load_protein(name)
        # Load EVERY order (not just the most populous 5) purely to falsify the popularity selection.
        saved, H2.MAX_ORDERS = H2.MAX_ORDERS, 99
        by_all = load_protein(name)
        H2.MAX_ORDERS = saved
        if not by_order:
            print(f"{name[:26]:26s} (assay or ESM cache absent -- skipped)")
            report["proteins"][name] = {"status": "cache_absent"}
            continue
        rows = subset_rows(by_order)
        if len(rows) < MIN_SUBSETS:
            print(f"{name[:26]:26s} only {len(rows)} subsets -- not reported")
            report["proteins"][name] = {"status": "too_few_subsets", "n_subsets": len(rows)}
            continue
        c = correlate(rows)
        loo = leave_one_order_out(by_order)
        slates = slate_sensitivity(by_all) if by_all else None
        report["proteins"][name] = {
            "status": "ok", "orders": sorted(by_order), **c,
            "mean_additive_score_by_order": mean_score_by_order(by_order),
            "leave_one_order_out": loo, "slate_sensitivity": slates, "subsets": rows,
        }
        print(f"{name[:26]:26s} {c['n_subsets']:4d} {c['rho_gain_eta2_fitness']:11.3f} "
              f"{c['rho_gain_eta2_score']:13.3f} {c['rho_eta2_fitness_eta2_score']:8.3f} "
              f"{c['partial_fitness_given_score']:12.3f} {c['partial_score_given_fitness']:15.3f}")
        ok = [r for r in loo if r.get("status") == "ok"]
        if ok:
            lo = min(r["rho_gain_eta2_fitness"] for r in ok)
            hi = max(r["rho_gain_eta2_fitness"] for r in ok)
            los = min(r["rho_gain_eta2_score"] for r in ok)
            his = max(r["rho_gain_eta2_score"] for r in ok)
            print(f"{'':26s}      leave-one-order-out: rho(g,FIT) {lo:+.3f}..{hi:+.3f} | "
                  f"rho(g,SCORE) {los:+.3f}..{his:+.3f}   ({len(ok)} drops)")
        print(f"{'':26s}      mean additive score by order: "
              f"{report['proteins'][name]['mean_additive_score_by_order']}")
        if slates and slates["slates"]:
            print(f"{'':26s}      SLATE SENSITIVITY (does popularity-selection drive it?):")
            for lbl, s in slates["slates"].items():
                print(f"{'':30s}{lbl:26s} k={str(s['orders']):20s} n={s['n_variants']:6d} "
                      f"rho(g,FIT)={s['rho_gain_eta2_fitness']:+.3f} "
                      f"rho(g,SCORE)={s['rho_gain_eta2_score']:+.3f}")

    OUT.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"\n-> {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
