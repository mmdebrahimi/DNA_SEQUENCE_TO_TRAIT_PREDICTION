"""Is the +0.040 (E over B) real, or an artifact of ONE arbitrary position split?

Two checks, both required before the headline is quotable:

1. STRIDE SWEEP. The main run holds out every 3rd sorted position. That is one arbitrary partition. If the
   ordering of the five variants flips across strides 2/3/4/5, nothing here is a finding.

2. POSITION-CLUSTERED BOOTSTRAP. n=496 variants but only 88 INDEPENDENT positions -- ~19 variants share a
   position, so a variant-level bootstrap would overstate precision by treating clustered observations as
   independent. Resample POSITIONS with replacement instead (the same discipline the lineage layer applies
   to clonal isolates).

The paired comparison is E-minus-B on the SAME held-out variants, never a difference of separately-computed
medians -- differencing two independently-computed numbers is a documented trap in this repo.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from pear_genotype_alphabet import (  # noqa: E402
    CODON, assert_sign_convention, blosum62_fn, esm2_marginals, load_ref, onehot,
    parse_variants, predict, ridge, spearman,
)

OUT = ROOT / "wiki" / "pear_genotype_alphabet_robustness_2026-09-29.json"


def run_split(mis, esm, bl, stride):
    positions = sorted({r["pos"] for r in mis})
    test_pos = set(positions[::stride])
    tr = [r for r in mis if r["pos"] not in test_pos]
    te = [r for r in mis if r["pos"] in test_pos]
    y_tr = np.array([r["fitness"] for r in tr])
    y_te = np.array([r["fitness"] for r in te])

    sA = np.array([float(bl(r["wt"], r["mut"])) for r in te])
    sB = np.array([esm[str(r["pos"])][r["mut"]] - esm[str(r["pos"])][r["wt"]] for r in te])
    pair, permut = {}, {}
    for r in tr:
        pair.setdefault((r["wt"], r["mut"]), []).append(r["fitness"])
        permut.setdefault(r["mut"], []).append(r["fitness"])
    pm = {k: float(np.mean(v)) for k, v in pair.items()}
    mm = {k: float(np.mean(v)) for k, v in permut.items()}
    g = float(np.mean(y_tr))
    sC = np.array([pm.get((r["wt"], r["mut"]), mm.get(r["mut"], g)) for r in te])
    sD = predict(ridge(onehot(tr, esm, False), y_tr), onehot(te, esm, False))
    sE = predict(ridge(onehot(tr, esm, True), y_tr), onehot(te, esm, True))
    return {
        "stride": stride, "n_test_positions": len(test_pos), "n_test": len(te),
        "A_blosum62": spearman(sA, y_te), "B_esm2": spearman(sB, y_te),
        "C_pairmeans": spearman(sC, y_te), "D_ridge_onehot": spearman(sD, y_te),
        "E_ridge_plus_esm2": spearman(sE, y_te),
    }, (te, y_te, sB, sE)


def clustered_bootstrap(te, y, sB, sE, n=2000, seed=0):
    """Resample POSITIONS with replacement -- variants within a position are not independent."""
    rng = np.random.default_rng(seed)
    by_pos: dict[int, list[int]] = {}
    for i, r in enumerate(te):
        by_pos.setdefault(r["pos"], []).append(i)
    pos = list(by_pos)
    deltas = []
    for _ in range(n):
        draw = rng.choice(len(pos), size=len(pos), replace=True)
        idx = [i for d in draw for i in by_pos[pos[d]]]
        if len({y[i] for i in idx}) < 3:
            continue
        deltas.append(spearman(sE[idx], y[idx]) - spearman(sB[idx], y[idx]))
    d = np.array(deltas)
    return {
        "n_draws": int(len(d)), "mean_delta": float(d.mean()),
        "ci95": [float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))],
        "frac_positive": float((d > 0).mean()),
    }


def main() -> int:
    cds = load_ref()
    rows = parse_variants(cds)
    assert_sign_convention(rows)
    mis = [r for r in rows if r["kind"] == "missense"]
    protein = "".join(CODON[cds[i:i + 3]] for i in range(0, len(cds) - len(cds) % 3, 3)).rstrip("*")
    esm = esm2_marginals(protein)
    bl = blosum62_fn()

    sweep, keep = [], None
    for stride in (2, 3, 4, 5):
        res, payload = run_split(mis, esm, bl, stride)
        sweep.append(res)
        if stride == 3:
            keep = payload
        order = " > ".join(k.split("_")[0] for k, _ in
                           sorted(((k, v) for k, v in res.items() if len(k) > 2 and k[1] == "_"),
                                  key=lambda kv: -kv[1]))
        print(f"stride {stride}: pos_test {res['n_test_positions']:>3} n {res['n_test']:>4} | "
              f"A {res['A_blosum62']:+.3f} B {res['B_esm2']:+.3f} C {res['C_pairmeans']:+.3f} "
              f"D {res['D_ridge_onehot']:+.3f} E {res['E_ridge_plus_esm2']:+.3f}   [{order}]")

    # stability of the two claims across strides
    b1 = [max(s["C_pairmeans"], s["D_ridge_onehot"]) > s["A_blosum62"] for s in sweep]
    b3 = [s["E_ridge_plus_esm2"] > s["B_esm2"] for s in sweep]
    print(f"\nB1 (learned alphabet > BLOSUM62) across strides: {b1}  -> {'STABLE' if len(set(b1))==1 else 'UNSTABLE'}")
    print(f"B3 (ESM2+head > ESM2 alone)      across strides: {b3}  -> {'STABLE' if len(set(b3))==1 else 'UNSTABLE'}")

    te, y, sB, sE = keep
    boot = clustered_bootstrap(te, y, sB, sE)
    print(f"\nposition-clustered bootstrap of (E - B), stride 3, {boot['n_draws']} draws:")
    print(f"  mean delta {boot['mean_delta']:+.4f}   95% CI [{boot['ci95'][0]:+.4f}, {boot['ci95'][1]:+.4f}]"
          f"   P(delta>0) = {boot['frac_positive']:.3f}")
    excl = boot["ci95"][0] > 0
    print(f"  CI excludes zero: {excl}")

    OUT.write_text(json.dumps({
        "schema": "pear-genotype-alphabet-robustness-v1", "analysis_date": "2026-09-29",
        "stride_sweep": sweep,
        "B1_stable_across_strides": len(set(b1)) == 1, "B1_values": b1,
        "B3_stable_across_strides": len(set(b3)) == 1, "B3_values": b3,
        "clustered_bootstrap_E_minus_B": boot,
        "bootstrap_unit": "POSITION (not variant) -- ~19 variants share a position and are not independent",
        "honest_limits": [
            "Strides are nested partitions of the SAME 264 positions, so the four rows are not four "
            "independent experiments -- they test partition-sensitivity, not replication.",
            "ESM2 sees no measured fitness; the supervised variants see the training positions. That "
            "asymmetry favours the supervised side, so E>B is SUGGESTIVE while C/D<A is DECISIVE.",
            "One protein, one assay, one lab.",
        ],
    }, indent=1), encoding="utf-8")
    print(f"\nwrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
