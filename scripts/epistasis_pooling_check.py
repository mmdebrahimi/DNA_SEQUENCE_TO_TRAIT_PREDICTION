"""Is the epistasis sweep's per-protein Delta a MUTATION-ORDER POOLING ARTIFACT?

The 2026-07-27 sweep reported per-protein joint-vs-additive Deltas POOLED across mutation orders, and its
headline anomaly was ParD at Delta = -0.283 (vs ~+-0.005 everywhere else). Pooling across a variable that
BOTH scores track is exactly the shape that inflates a metric -- the same trap as clonality inflation on
the AMR side, with mutation order in place of the clone.

This recomputes each protein's Delta WITHIN order (n-weighted over per-order Spearman) and contrasts it
with the pooled figure. Pure re-analysis of the committed sweep JSON -- no GPU, no network, no re-scoring.

The assay CSVs (cached, gitignored) additionally allow two explanation tests; both are reported honestly,
and the first is FALSIFIED by its own numbers:
  H1  additive is a SUM of k terms, so it should harvest cross-order signal where fitness declines with k
      -> measured by Spearman(k, fitness). FALSIFIED: GFP has the steeper slope and the smaller gain.
  H2  pooling inflation should track BETWEEN-order variance as a fraction of total -> eta^2(k).
      CONSISTENT but n=3 (a monotone ordering of 3 points arises ~17% of the time by chance).

Run: uv run python scripts/epistasis_pooling_check.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

SWEEP = ROOT / "wiki" / "forward_epistasis_sweep_2026-07-27.json"
ASSAYS = Path("D:/dna_decode_cache/epistasis")


def within_order(per_order: dict, metric: str) -> float:
    """n-weighted mean of the per-order Spearman for `metric`. PURE."""
    num = sum(v[metric] * v["n"] for v in per_order.values())
    den = sum(v["n"] for v in per_order.values())
    return num / den if den else float("nan")


def pooling_report(results: list[dict]) -> list[dict]:
    """Per protein: within-order vs pooled rho for both scorers, and both Deltas. PURE."""
    out = []
    for r in results:
        po = r.get("per_order") or {}
        if len(po) < 2:            # a single order cannot show a pooling effect
            continue
        wa, wj = within_order(po, "additive"), within_order(po, "joint")
        out.append({
            "protein": r["dms"], "seqlen": r.get("seqlen"), "orders": sorted(int(k) for k in po),
            "within_additive": round(wa, 4), "pooled_additive": round(r["additive"], 4),
            "within_joint": round(wj, 4), "pooled_joint": round(r["joint"], 4),
            "additive_pooling_gain": round(r["additive"] - wa, 4),
            "joint_pooling_gain": round(r["joint"] - wj, 4),
            "within_delta": round(wj - wa, 4), "pooled_delta": round(r["joint"] - r["additive"], 4),
        })
    return out


def order_structure(csv_path: Path) -> dict | None:
    """Spearman(k, fitness) [H1] + eta^2(k) [H2] for one assay. None when the CSV is absent."""
    if not csv_path.exists():
        return None
    import pandas as pd
    from scipy.stats import spearmanr
    df = pd.read_csv(csv_path)
    df["k"] = df["mutant"].str.count(":") + 1
    s = df[df.k.between(2, 6)]
    if s.empty or s.k.nunique() < 2:
        return None
    rho, _ = spearmanr(s["k"], s["DMS_score"])
    grand = s["DMS_score"].mean()
    ss_b = sum(len(v) * (v["DMS_score"].mean() - grand) ** 2 for _, v in s.groupby("k"))
    ss_t = float(((s["DMS_score"] - grand) ** 2).sum())
    return {"spearman_k_fitness": round(float(rho), 4),
            "eta2_k": round(ss_b / ss_t, 4) if ss_t else None,
            "n": int(len(s))}


def main() -> int:
    rows = pooling_report(json.loads(SWEEP.read_text(encoding="utf-8"))["results"])
    print("WITHIN-ORDER vs POOLED (the pooled figure is what the 2026-07-27 sweep published)\n")
    print(f"{'protein':30s} {'within Δ':>9s} {'pooled Δ':>9s} {'add gain':>9s} {'joint gain':>11s}")
    for r in rows:
        print(f"{r['protein'][:30]:30s} {r['within_delta']:+9.3f} {r['pooled_delta']:+9.3f} "
              f"{r['additive_pooling_gain']:+9.3f} {r['joint_pooling_gain']:+11.3f}")

    print("\nEXPLANATION TESTS (assay CSVs; skipped when the cache is absent)\n")
    print(f"{'protein':30s} {'rho(k,fit) H1':>14s} {'eta2(k) H2':>11s} {'add gain':>9s}")
    struct = {}
    for r in rows:
        st = order_structure(ASSAYS / f"{r['protein']}.csv")
        if not st:
            print(f"{r['protein'][:30]:30s} {'(assay not cached)':>14s}")
            continue
        struct[r["protein"]] = st
        print(f"{r['protein'][:30]:30s} {st['spearman_k_fitness']:14.3f} {st['eta2_k']:11.3f} "
              f"{r['additive_pooling_gain']:+9.3f}")

    if len(struct) >= 3:
        print("\n  H1 (fitness-vs-k slope): FALSIFIED where the steepest slope does NOT carry the largest "
              "gain.")
        print(f"  H2 (eta^2): consistent if monotone -- but n={len(struct)}; a monotone ordering of 3 "
              "points arises ~17% of the time by chance. NOT a mechanism.")

    out = ROOT / "wiki" / "forward_epistasis_pooling_correction_2026-08-25.json"
    out.write_text(json.dumps({
        "_schema": "epistasis-pooling-correction-v1",
        "source_sweep": SWEEP.name,
        "honest_scope": ("within-order Delta is the honest number when a sweep spans mutation orders; the "
                         "pooled Delta flatters any scorer that scales with k. H1 falsified, H2 "
                         "underpowered (n=3), density confounded (n=1 protein above 3% k/L)."),
        "per_protein": rows, "order_structure": struct}, indent=2), encoding="utf-8")
    print(f"\n-> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
