"""Does the E-Flux null's mechanism hold at scale, or was it a selection artifact?

WHY THIS EXISTS
---------------
Commit 5b61e1e killed its own first explanation of the E-Flux null and replaced it with a better one:

    "The growth-ratio distribution is BIMODAL -- genes sit at ~1.0 or ~0.0 with near-nothing between.
     E-Flux moves magnitudes by up to 0.43 but almost never across the 0.01 line, because that line
     sits in an empty region."

That mechanism now carries real weight: it produced the conclusion "expression constraints can only
help if the readout changes to something graded -- a metric redesign, not more constraint work",
which redirects the whole workstream.

But it was measured on `conditionally_essential_genes(recs)[:40]` -- 40 genes SELECTED FOR BEING
CONDITIONALLY ESSENTIAL, out of 1,516 in the model. Conditionally essential genes are by construction
the ones that flip between lethal and viable across conditions, i.e. exactly the population most
likely to sit at the extremes. So the observed bimodality may be a property of the SAMPLE rather than
of the readout, and the conclusion drawn from it would not follow.

WHAT THIS MEASURES
------------------
Two DISTINCT claims that the commit message runs together:

  (A) GLOBAL BIMODALITY -- are growth ratios concentrated at ~0 and ~1 with an empty middle?
  (B) LOCAL EMPTINESS AT THE THRESHOLD -- is the 0.01 line specifically in a sparse region?

(B) is the one that actually carries the mechanism. A distribution can be strongly bimodal and STILL
have genes packed just above 0.01 (in the lower mode's shoulder), and those are precisely the genes a
flux rescale could flip. Only (B) licenses "the classification is robust to large flux changes".

Baseline arm only: the question is about the READOUT's shape, which does not require the E-Flux arm.

Run: .venv/Scripts/python.exe scripts/fba_bimodality_at_scale.py
"""
import sys

sys.path.insert(0, '.')
sys.path.insert(0, 'scripts')

import numpy as np
from cobra.flux_analysis import single_gene_deletion

SOLVER_NOISE = []   # every raw (pre-clamp) ratio, so the clamp can be audited

from dna_decode.fba.conditional_essentiality import conditionally_essential_genes
from dna_decode.fba.fitness_browser import (ESSENTIAL_FITNESS, apply_carbon_condition, carbon_conditions,
                                            load_records, open_db)
from dna_decode.fba.model import load_model, wildtype_growth
from fba_eflux_bridge import FRAC, build_condition_expression

# Window around the essentiality threshold. A gene inside it is one a flux rescale could plausibly
# flip between lethal and viable -- the population the mechanism claims is empty.
NEAR_LO, NEAR_HI = 0.001, 0.10
MID_LO, MID_HI = 0.05, 0.95          # the "empty middle" of a bimodal distribution


def ratios_for_condition(model, cond, allex, genes):
    with model:
        apply_carbon_condition(model, cond, all_carbon=allex)
        wt = wildtype_growth(model)
        res = single_gene_deletion(model, gene_list=genes, processes=1)
        out = {}
        for _, r in res.iterrows():
            g = next(iter(r["ids"]))
            growth = r["growth"]
            ratio = 0.0 if growth != growth else (growth / wt if wt > 0 else 0.0)
            # The LP returns tiny NEGATIVE growth for lethal knockouts (observed min ~-1.1e-12).
            # A negative growth rate is not physical - it is solver noise - so clamp at zero. The
            # magnitude is asserted below so a real negative could never be silently absorbed.
            SOLVER_NOISE.append(ratio)
            out[g] = max(0.0, ratio)
        return out


def main():
    model = load_model()
    conn = open_db()
    conds_all = carbon_conditions(conn, model)
    expr, _ = build_condition_expression(conds_all)
    keys = sorted(expr)
    conds = {k: conds_all[k] for k in keys}
    recs = load_records(conn, conds, gene_filter={g.id for g in model.genes}, threshold=ESSENTIAL_FITNESS)

    ce40 = [r.gene_id for r in conditionally_essential_genes(recs)][:40]
    all_genes = [g for g in model.genes]
    allex = tuple(conds_all.values())

    print("BIMODALITY AT SCALE - is the E-Flux null's mechanism real or a selection artifact?")
    print("=" * 92)
    print(f"model genes {len(all_genes)} | conditions {len(conds)} | essentiality threshold FRAC={FRAC}")
    print(f"the mechanism was measured on {len(ce40)} CONDITIONALLY-ESSENTIAL genes; this runs ALL genes\n")

    pooled_all, pooled_ce = [], []
    print(f"  {'condition':>22} {'n':>6} {'frac<thresh':>12} {'in mid band':>12} {'NEAR thresh':>12}")
    for name in keys:
        r = ratios_for_condition(model, conds[name], allex, all_genes)
        v = np.array([r[g.id] for g in all_genes])
        pooled_all.append(v)
        ce = np.array([r[g] for g in ce40 if g in r])
        pooled_ce.append(ce)
        mid = ((v > MID_LO) & (v < MID_HI)).mean()
        near = ((v > NEAR_LO) & (v < NEAR_HI)).sum()
        print(f"  {name[:22]:>22} {len(v):6d} {(v < FRAC).mean():11.1%} {mid:11.1%} {near:9d}   ")

    A = np.concatenate(pooled_all)
    C = np.concatenate(pooled_ce)

    print(f"\n[1] GLOBAL SHAPE - pooled over all conditions (n={len(A)} gene x condition cells)")
    edges = [0.0, 0.001, 0.01, 0.05, 0.2, 0.5, 0.8, 0.95, 0.999, 1.0001]
    hist, _ = np.histogram(A, bins=edges)
    for i in range(len(hist)):
        print(f"  [{edges[i]:<6.3f},{edges[i+1]:>6.3f})  {hist[i]:7d}  {hist[i]/len(A):7.2%}")
    mid_all = ((A > MID_LO) & (A < MID_HI)).mean()
    mid_ce = ((C > MID_LO) & (C < MID_HI)).mean()
    print(f"\n  intermediate band ({MID_LO}-{MID_HI}):  ALL genes {mid_all:.2%}   |   the CE-40 sample {mid_ce:.2%}")

    print(f"\n[2] THE LOAD-BEARING CLAIM - is the {FRAC} threshold in an EMPTY region?")
    near_all = ((A > NEAR_LO) & (A < NEAR_HI)).sum()
    print(f"  genes within ({NEAR_LO}, {NEAR_HI}) of the threshold: {near_all} of {len(A)} "
          f"= {near_all/len(A):.3%}")
    print(f"  ALL genes exactly 0.0 : {(A <= 0.0).mean():.1%}")
    print(f"  ALL genes >= 0.95     : {(A >= 0.95).mean():.1%}")

    print("\n[3] VERDICT")
    artifact = mid_all > 3 * max(mid_ce, 1e-9) and mid_all > 0.05
    if near_all / len(A) < 0.01:
        print(f"  The threshold region IS sparse ({near_all/len(A):.3%} of cells). The mechanism's")
        print("  load-bearing claim SURVIVES at scale: a flux rescale has almost no genes to flip.")
    else:
        print(f"  The threshold region is NOT empty ({near_all/len(A):.3%} of cells sit beside it).")
        print("  The mechanism's load-bearing claim does NOT hold at scale.")
    if artifact:
        print(f"  AND the intermediate band is {mid_all/max(mid_ce,1e-9):.1f}x denser over ALL genes than")
        print("  over the CE-40 sample -> the observed bimodality was partly a SELECTION ARTIFACT.")
    else:
        print("  Bimodality is NOT merely a selection artifact: the full gene set shows it too.")
    return A, C, near_all, mid_all, mid_ce


def _verify(A, C, near_all, mid_all, mid_ce):
    assert len(A) > 10000, f"expected >10k gene x condition cells, got {len(A)}"
    assert len(C) > 0, "the CE-40 comparison sample must be non-empty"
    raw = np.array(SOLVER_NOISE)
    worst_neg = raw.min()
    assert worst_neg > -1e-6, (
        f"clamped a NON-trivial negative growth ratio ({worst_neg:.3g}) - that would be a real "
        f"result being hidden, not solver noise")
    assert (A >= 0.0).all() and (A <= 1.5).all(), "clamped ratios outside a sane range"
    # the sample the mechanism was measured on must be a small minority of the model
    assert len(set(range(len(C)))) < len(A) / 10, "CE sample should be far smaller than the full set"
    print(f"\n[verify-in-batch] PASS: scanned {len(A)} gene x condition cells over the FULL "
          f"{len(A)//11}-gene model (the mechanism was derived from 40 SELECTED genes); the only "
          f"out-of-range values were solver noise (worst raw negative {worst_neg:.2g}, clamped to 0, "
          f"magnitude asserted < 1e-6 so a real negative cannot hide); and the threshold-window count "
          f"({near_all}) plus both intermediate-band fractions (ALL {mid_all:.2%} vs CE-40 {mid_ce:.2%}) "
          f"come from the SAME executed scan, so the artifact check compares like with like.")


if __name__ == "__main__":
    _verify(*main())
