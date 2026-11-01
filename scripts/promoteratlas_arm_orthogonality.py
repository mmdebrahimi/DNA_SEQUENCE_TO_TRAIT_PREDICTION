"""Join the PromoterAtlas arm's per-tile predictions against the incumbent channels.

THIS IS THE MEASUREMENT THAT SETTLES THE ARM'S OWN MERIT QUESTION, and it matters more than the
arm's headline number. The standing argument against this arm is an INFERENCE recorded in
wiki/promoteratlas_circularity_2026-10-31.md: PromoterAtlas would be a THIRD sequence-learned
channel beside the conv already in the stack, where the next channel should carry a different PRIOR
rather than a different neural net. That memo says the measurement which would settle it "needs
per-tile vectors the arm has not produced". The arm produces them, so the inference becomes
testable here.

What it computes, per seed, on the incumbents' own held-out tiles:
  1. Spearman between the atlas predictions and each incumbent channel (conv / gc / max_tx /
     decomposition). PREDICTION on record before running: atlas-vs-conv HIGH (both sequence-learned),
     atlas-vs-max_tx LOWER (different prior).
  2. The rank-average of max_tx + atlas, against the committed two-channel hybrid at 0.4254. If
     atlas correlates highly with conv AND adds nothing over the hybrid, the inferred merit-negative
     becomes MEASURED rather than argued.
  3. A paired comparison against conv on identical test sets, with a distribution-free interval.

Usage:  uv run python -I scripts/promoteratlas_arm_orthogonality.py <arm_artifact.json>
"""
from __future__ import annotations

import gzip
import json
import statistics as st
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
TABLES = REPO / "data" / "promcalc" / "tile_predictions"
#: The committed two-channel hybrid this arm would have to improve on to earn a place in the stack.
HYBRID_2CH = 0.4254
STACK = 0.4404


def spearman(a, b):
    from scipy.stats import spearmanr
    r = spearmanr(a, b).statistic
    return float(r) if r == r else 0.0


def rank_average(*vecs):
    """Mean of within-vector ranks. Matches the rank-average the hybrid/stack arms used: it is
    scale-free, so a channel on a different unit cannot dominate by magnitude alone."""
    n = len(vecs[0])
    out = [0.0] * n
    for v in vecs:
        order = sorted(range(n), key=lambda i: v[i])
        rk = [0.0] * n
        for pos, i in enumerate(order):
            rk[i] = pos
        for i in range(n):
            out[i] += rk[i] / len(vecs)
    return out


def median_ci95(vals):
    """Exact distribution-free order-statistic interval for the MEDIAN, reporting its ACHIEVED
    coverage. The same discipline the fragment arm uses -- and deliberately NOT a normal
    approximation around the mean, which is a different estimand and has been misread here before.
    """
    from math import comb
    v = sorted(vals)
    n = len(v)
    best = None
    for k in range(n // 2 + 1):
        cov = 1.0 - 2.0 * sum(comb(n, i) for i in range(k)) / (2 ** n)
        if cov >= 0.95:
            best = (k, cov)
    if best is None:
        return None, None, None
    k, cov = best
    return v[k], v[n - 1 - k], round(cov, 4)


def main() -> int:
    if len(sys.argv) < 2:
        print(__doc__)
        return 2
    art = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
    if art.get("smoke"):
        print("REFUSED: %s is a SMOKE artifact. A subsampled pool did not reproduce the "
              "incumbents' split, so no join is valid." % sys.argv[1])
        return 3
    if art.get("verdict", "").startswith("INVALID"):
        print("REFUSED: arm verdict is %s -- nothing may be joined." % art["verdict"])
        return 3

    pertile = art["per_tile_predictions"]
    rows, orth, ranks = {}, {}, {}
    print("%-6s %-9s %-9s %-9s %-9s %-9s %s"
          % ("seed", "atlas", "conv", "gc", "max_tx", "decomp", "atlas-vs: conv / max_tx"))
    print("-" * 104)
    for s in sorted(pertile, key=int):
        with gzip.open(TABLES / ("seed-%02d.json.gz" % int(s)), "rt", encoding="utf-8") as fh:
            tab = json.load(fh)
        idx_arm = pertile[s]["tile_indices"]
        if list(idx_arm) != list(tab["tile_indices"]):
            print("REFUSED: seed %s tile_indices differ from the committed table. The arm was not "
                  "scored on the incumbents' test set." % s)
            return 3
        atlas = pertile[s]["pred"]
        meas = tab["measured"]
        pr = tab["predictions"]
        o = {c: round(spearman(atlas, pr[c]), 4) for c in ("conv", "gc", "max_tx", "decomposition")}
        orth[s] = o
        ra = rank_average(pr["max_tx"], atlas)
        ranks[s] = round(spearman(ra, meas), 4)
        rows[s] = {"atlas": art["arms"][art["geometry_primary"]]["per_seed"][s],
                   **{c: round(tab["spearman_recomputed"][c], 4)
                      for c in ("conv", "gc", "max_tx", "decomposition")},
                   "rank_avg_maxtx_atlas": ranks[s]}
        print("%-6s %-9.4f %-9.4f %-9.4f %-9.4f %-9.4f %.4f / %.4f"
              % (s, rows[s]["atlas"], rows[s]["conv"], rows[s]["gc"], rows[s]["max_tx"],
                 rows[s]["decomposition"], o["conv"], o["max_tx"]))

    seeds = sorted(rows, key=int)
    med = lambda key: round(st.median(rows[s][key] for s in seeds), 4)
    paired = [rows[s]["atlas"] - rows[s]["conv"] for s in seeds]
    lo, hi, cov = median_ci95(paired)
    wins = sum(1 for d in paired if d > 0)

    out = {
        "schema": "promoteratlas-orthogonality-v1",
        "source_artifact": sys.argv[1],
        "n_seeds": len(seeds),
        "medians": {k: med(k) for k in
                    ("atlas", "conv", "gc", "max_tx", "decomposition", "rank_avg_maxtx_atlas")},
        "orthogonality_per_seed": orth,
        "orthogonality_median": {c: round(st.median(orth[s][c] for s in seeds), 4)
                                 for c in ("conv", "gc", "max_tx", "decomposition")},
        "paired_vs_conv": {"per_seed_delta": [round(d, 4) for d in paired],
                           "median": round(st.median(paired), 4),
                           "median_ci95": [round(lo, 4), round(hi, 4)] if lo is not None else None,
                           "achieved_coverage": cov, "seeds_won": wins, "n": len(paired)},
        "stack_comparison": {"rank_avg_maxtx_plus_atlas": med("rank_avg_maxtx_atlas"),
                             "committed_2ch_hybrid_maxtx_plus_conv": HYBRID_2CH,
                             "committed_3ch_stack": STACK,
                             "delta_vs_2ch_hybrid": round(med("rank_avg_maxtx_atlas") - HYBRID_2CH, 4)},
    }

    om = out["orthogonality_median"]
    d2 = out["stack_comparison"]["delta_vs_2ch_hybrid"]
    if om["conv"] > om["max_tx"] and d2 <= 0:
        out["merit_inference"] = (
            "MEASURED, AND IT CONFIRMS THE STANDING INFERENCE. The atlas predictions correlate more "
            "with the conv channel (%.4f) than with the thermodynamic one (%.4f), and rank-averaging "
            "atlas with max_tx does NOT beat the committed max_tx+conv hybrid (%.4f vs %.4f, delta "
            "%+.4f). So PromoterAtlas behaves as a third SEQUENCE-LEARNED channel and the next "
            "channel should carry a different PRIOR, exactly as inferred. The inference is no longer "
            "an inference." % (om["conv"], om["max_tx"], out["medians"]["rank_avg_maxtx_atlas"],
                               HYBRID_2CH, d2))
    elif d2 > 0:
        out["merit_inference"] = (
            "MEASURED, AND IT REFUTES THE STANDING INFERENCE. Rank-averaging atlas with max_tx beats "
            "the committed max_tx+conv hybrid by %+.4f, so this encoder carries something the conv "
            "channel does not. The claim that the next channel must carry a different PRIOR does not "
            "survive." % d2)
    else:
        out["merit_inference"] = (
            "MIXED: atlas correlates more with the thermodynamic channel (%.4f) than with conv "
            "(%.4f), which the inference did not predict, yet it still does not beat the 2-channel "
            "hybrid (delta %+.4f). Reported as mixed rather than forced either way."
            % (om["max_tx"], om["conv"], d2))

    dest = REPO / "wiki" / "promoteratlas_orthogonality_2026-11-01.json"
    dest.write_text(json.dumps(out, indent=1) + "\n", encoding="utf-8")

    print("\n--- medians ---")
    for k, v in out["medians"].items():
        print("   %-26s %.4f" % (k, v))
    print("\n--- orthogonality (median Spearman between atlas and each channel) ---")
    for c, v in om.items():
        print("   atlas vs %-16s %.4f" % (c, v))
    p = out["paired_vs_conv"]
    print("\n--- paired vs conv, identical test sets ---")
    print("   median %+.4f  ci95 %s (achieved %.4f)  seeds won %d/%d"
          % (p["median"], p["median_ci95"], p["achieved_coverage"] or 0, p["seeds_won"], p["n"]))
    print("\n--- merit inference ---")
    print("   " + out["merit_inference"].replace(". ", ".\n   "))
    print("\nwrote %s" % dest)
    return 0


if __name__ == "__main__":
    sys.exit(main())
