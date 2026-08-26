"""Does BETWEEN-ORDER variance (eta^2) predict the additive score's pooling gain? -- the H2 test.

CONTEXT. `wiki/forward_epistasis_pooling_correction_2026-08-25.md` established that the epistasis sweep's
headline ParD anomaly was a mutation-order POOLING artifact, and that the inflated half is the ADDITIVE
score's pooling gain (+0.284 on ParD vs +0.069/+0.082 elsewhere). Two explanations were offered:
  H1  additive is a SUM of k terms, so it harvests cross-order signal where fitness declines with k.
      FALSIFIED -- GFP has the steeper fitness-vs-k slope and the SMALLER gain.
  H2  the governing quantity is how much fitness variance sits BETWEEN orders -- eta^2(k).
      CONSISTENT but on n=3 PROTEINS, where a monotone ordering arises ~17% of the time by chance.

WHY THIS DESIGN. Adding proteins is not possible from cache: of the 8 cached epistasis assays only 3 carry
3+ mutation orders (GFP k2-6, HIS7 k2-6, ParD k2-4). But n=3 was never the real problem -- the problem was
that eta^2 was CONFOUNDED WITH PROTEIN IDENTITY (each protein contributed exactly one point, so "high
eta^2" and "is ParD" were the same statement). That is the same confound that killed H3 (density).

So this tests H2 **WITHIN protein**: for each protein and each SUBSET of its mutation orders, compute
eta^2(k) and the additive pooling gain on that subset. Protein identity is then held FIXED while eta^2
varies, which is the only way to see whether eta^2 is doing the work. 26 subsets for GFP, 26 for HIS7,
4 for ParD.

ZERO NEW COMPUTE. The additive score needs only the cached ESM2 per-position log-prob matrices (the JOINT
score, which would need a forward pass per mutant, is not used -- the pooling gain is a property of the
ADDITIVE score alone). Fitness comes from the cached assay CSVs. No GPU, no network.

It also uses the FULL assays rather than the sweep's 300-per-order subsample, so each estimate rests on
thousands of variants instead of hundreds.

Run: uv run python scripts/epistasis_pooling_h2_test.py
"""
from __future__ import annotations

import itertools
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

ASSAYS = Path("D:/dna_decode_cache/epistasis")
ESM = Path("D:/dna_decode_cache/esm")
PROTEINS = ["GFP_AEQVI_Sarkisyan_2016", "HIS7_YEAST_Pokusaeva_2019", "F7YBW8_MESOW_Aakre_2015"]
PER_ORDER_CAP = 20_000     # far above the sweep's 300; keeps the string parsing tractable
MIN_PER_ORDER = 30         # an order with fewer variants than this is dropped from a subset
MAX_ORDERS = 5             # see below -- bounds the subset count at 2^5-6 = 26 per protein

# WHY MAX_ORDERS EXISTS. GFP's assay carries mutation orders 2..12, not the 2..6 the sweep sampled, so
# "every subset of orders" is 2^11-12 = 2036 subsets -- a combinatorial blow-up that ran >25 min without
# finishing. Restricting to the MOST POPULOUS orders is also the statistically right call: the sparse
# high-k tails give noisy per-order Spearmans, and a noisy within_rho is exactly what would manufacture
# a spurious pooling gain.

_MUT = re.compile(r"^([A-Z])(\d+)([A-Z])$")


def additive_scores(mutants, table: dict) -> list[float | None]:
    """Sum of per-mutation ESM2 log-ratios. PURE. None when any token is unparseable/out of range."""
    out = []
    for m in mutants:
        tot, ok = 0.0, True
        for tok in str(m).split(":"):
            g = _MUT.match(tok.strip())
            if not g:
                ok = False
                break
            wt, pos, mut = g.group(1), g.group(2), g.group(3)
            col = table.get(pos)
            if not col or wt not in col or mut not in col:
                ok = False
                break
            tot += col[mut] - col[wt]
        out.append(tot if ok else None)
    return out


def eta2(groups: list) -> float:
    """Between-group share of total variance. PURE."""
    vals = [v for g in groups for v in g]
    if len(vals) < 2:
        return float("nan")
    grand = sum(vals) / len(vals)
    ss_b = sum(len(g) * (sum(g) / len(g) - grand) ** 2 for g in groups if g)
    ss_t = sum((v - grand) ** 2 for v in vals)
    return ss_b / ss_t if ss_t else float("nan")


def subset_stats(by_order: dict, orders: tuple) -> dict | None:
    """eta^2(k) + pooled/within additive rho for one subset of mutation orders. PURE."""
    from scipy.stats import spearmanr
    groups = [by_order[k] for k in orders]
    fit_groups = [[f for f, _ in g] for g in groups]
    pooled_f = [f for g in groups for f, _ in g]
    pooled_s = [s for g in groups for _, s in g]
    if len(pooled_f) < 3:
        return None
    pooled_rho = float(spearmanr(pooled_s, pooled_f).statistic)
    num = den = 0.0
    for g in groups:
        if len(g) < 3:
            return None
        r = float(spearmanr([s for _, s in g], [f for f, _ in g]).statistic)
        num += r * len(g)
        den += len(g)
    within = num / den
    return {"orders": list(orders), "n": int(den), "eta2": round(eta2(fit_groups), 4),
            "pooled_rho": round(pooled_rho, 4), "within_rho": round(within, 4),
            "pooling_gain": round(pooled_rho - within, 4)}


def load_protein(name: str):
    """-> {order: [(fitness, additive_score), ...]} using cached assay + cached ESM2 matrix."""
    import pandas as pd
    csv, esm = ASSAYS / f"{name}.csv", ESM / f"esm2_t33_650M_UR50D__{name}.json"
    if not (csv.exists() and esm.exists()):
        return None
    table = json.loads(esm.read_text(encoding="utf-8"))
    df = pd.read_csv(csv)
    df["k"] = df["mutant"].astype(str).str.count(":") + 1
    by_order = {}
    for k, grp in df[df.k >= 2].groupby("k"):
        if len(grp) > PER_ORDER_CAP:
            grp = grp.sample(PER_ORDER_CAP, random_state=0)
        sc = additive_scores(grp["mutant"], table)
        pairs = [(f, s) for f, s in zip(grp["DMS_score"], sc) if s is not None]
        if len(pairs) >= MIN_PER_ORDER:
            by_order[int(k)] = pairs
    if len(by_order) > MAX_ORDERS:      # keep the most populous orders, then restore ascending order
        keep = sorted(sorted(by_order, key=lambda k: -len(by_order[k]))[:MAX_ORDERS])
        by_order = {k: by_order[k] for k in keep}
    return by_order or None


def controls(by_order: dict, seed: int = 0) -> dict:
    """Two INTERVENTIONS that isolate where the pooling gain comes from. PURE (given by_order).

    A -- permute ORDER LABELS across all variants. Keeps every marginal distribution, destroys the
         fitness<->order association. H2 predicts the gain collapses to ~0.
    B -- permute FITNESS WITHIN each order. Destroys the within-order score<->fitness signal, keeps the
         per-order means. H2 predicts the gain SURVIVES, because it was never made of within-order signal.

    Together these turn a correlation into an intervention: if the gain is entirely between-order
    structure, A must kill it and B must not.
    """
    import random
    orders = tuple(sorted(by_order))
    rnd = random.Random(seed)

    allv = [v for k in orders for v in by_order[k]]
    shuf = allv[:]
    rnd.shuffle(shuf)
    a, i = {}, 0
    for k in orders:                       # re-deal the same group SIZES from the shuffled pool
        n = len(by_order[k])
        a[k] = shuf[i:i + n]
        i += n

    b = {}
    for k in orders:
        fs = [f for f, _ in by_order[k]]
        ss = [s for _, s in by_order[k]]
        rnd.shuffle(fs)
        b[k] = list(zip(fs, ss))

    return {"real": subset_stats(by_order, orders),
            "ctl_a_order_labels_shuffled": subset_stats(a, orders),
            "ctl_b_fitness_shuffled_within_order": subset_stats(b, orders)}


def main() -> int:
    from scipy.stats import spearmanr
    report = {"_schema": "epistasis-pooling-h2-v1", "per_order_cap": PER_ORDER_CAP, "proteins": {}}
    print = lambda *a, **k: (__import__("builtins").print(*a, **{**k, "flush": True}))
    print("H2 WITHIN-PROTEIN: does eta^2(k) predict the ADDITIVE score's pooling gain?")
    print("(protein identity held FIXED -- the confound that made the n=3 cross-protein read useless)\n")

    for name in PROTEINS:
        by_order = load_protein(name)
        if not by_order:
            print(f"{name[:34]:34s} (assay or ESM cache absent -- skipped)")
            continue
        orders = sorted(by_order)
        rows = []
        for size in range(2, len(orders) + 1):
            for combo in itertools.combinations(orders, size):
                st = subset_stats(by_order, combo)
                if st:
                    rows.append(st)
        if len(rows) < 4:
            print(f"{name[:34]:34s} only {len(rows)} subsets -- too few to test")
            continue
        rho = float(spearmanr([r["eta2"] for r in rows], [r["pooling_gain"] for r in rows]).statistic)
        p = float(spearmanr([r["eta2"] for r in rows], [r["pooling_gain"] for r in rows]).pvalue)
        ctl = controls(by_order)
        report["proteins"][name] = {"n_subsets": len(rows), "orders": orders,
                                    "spearman_eta2_vs_gain": round(rho, 4), "p": round(p, 5),
                                    "controls": ctl, "subsets": rows}
        print(f"{name[:34]:34s} orders={orders} subsets={len(rows):3d}  "
              f"spearman(eta2, pooling_gain) = {rho:+.3f}  (p={p:.4g})")
        for lbl, key in (("real", "real"), ("ctlA order-labels shuffled", "ctl_a_order_labels_shuffled"),
                         ("ctlB fitness shuffled in-order", "ctl_b_fitness_shuffled_within_order")):
            c = ctl[key]
            print(f"{'':36s}{lbl:32s} eta2={c['eta2']:.3f} within={c['within_rho']:+.3f} "
                  f"gain={c['pooling_gain']:+.3f}")

    out = ROOT / "wiki" / "forward_epistasis_h2_within_protein_2026-08-25.json"
    out.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"\n-> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
