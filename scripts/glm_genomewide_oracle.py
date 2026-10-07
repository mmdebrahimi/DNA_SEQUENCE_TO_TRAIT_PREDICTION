"""Does the oracle's rho 0.59 survive OFF the designed grid? And does growth condition carry signal?

    uv run python scripts/glm_genomewide_oracle.py [--limit 20000] [--no-download]

Two questions on one substrate (GSE144621, 321k sheared *E. coli* genomic fragments in LB and M9):

**Q1 -- generalisation.** The oracle's 0.59 is on a fully-crossed DESIGNED grid, with the standing caveat
that such a grid is not a natural promoter distribution. These are real genomic fragments.

**Q2 -- condition.** The bacterial analogue of cell-state conditioning: the SAME fragments are measured in
two media, so whether condition carries signal beyond noise is directly measurable rather than assumed.

**The comparison is CEILING-RELATIVE, and that is the point.** This assay ships two replicates, so its own
reproducibility is computable and bounds any sequence model. Measured ceiling ~0.79. A raw
0.59-versus-X comparison against the designed grid would be unfair to X because the grid's measurements are
far cleaner -- and GSE108535 offers NO derivable ceiling (its two published tables correlate at exactly
1.0, one monotone transform of the other, not replicates). So the grid number cannot be ceiling-normalised
and the two are **not** directly comparable. Said plainly rather than papered over.
"""
from __future__ import annotations

import argparse
import json
from datetime import date
from pathlib import Path

from dna_decode.glm.genomewide import (
    DEFAULT_CACHE,
    condition_effect,
    gc_features,
    kmer_features,
    leave_peak_out,
    load_fragments,
    load_peak_tiles,
    noise_ceiling,
    position_blocked_split,
    positional_onehot,
    spearman,
)

#: Frozen before any number was looked at.
PREREGISTERED = {
    "q1_question": "does a learned sequence feature set predict expression on REAL GENOMIC fragments, "
                   "under a position-blocked split, and how much of the DERIVED noise ceiling does it reach?",
    "q1_bars": {
        "beats_crude_baseline_by": 0.05,
        "fraction_of_ceiling_for_GENERALISES": 0.50,
    },
    "q1_verdicts": {
        "GENERALISES_OFF_THE_DESIGNED_GRID": "learned beats GC by >= 0.05 AND reaches >= 50% of the ceiling",
        "WEAK_BUT_REAL": "learned beats GC by >= 0.05 but reaches < 50% of the ceiling",
        "DOES_NOT_GENERALISE": "learned does not beat GC by 0.05",
    },
    "q2_question": "does growth condition carry signal beyond measurement noise?",
    "q2_rule": "the comparison MUST be noise-matched (single-replicate vs single-replicate both within and "
               "across media). The unmatched average-vs-average form INVERTS the answer and is reported "
               "alongside only to show that it does.",
    "thresholds_are": "ASSERTED; the noise ceiling that normalises them is DERIVED from the assay's own "
                      "replicates.",
    "designed_grid_comparator": {
        "gc": 0.1414, "kmer4": 0.4137, "kmer6": 0.5872, "onehot": 0.5914,
        "caveat": "leave-ELEMENT-out on a designed grid; NOT ceiling-normalisable (no replicates), so not "
                  "directly comparable to a ceiling-relative number here",
    },
}

FEATURES = {"gc": gc_features, "kmer4": lambda s: kmer_features(s, 4),
            "kmer6": lambda s: kmer_features(s, 6), "onehot": positional_onehot}


def _fit(train, test, feature: str) -> dict:
    import numpy as np
    from sklearn.linear_model import Ridge
    from sklearn.preprocessing import StandardScaler
    fn = FEATURES[feature]
    Xtr = np.asarray(fn([f.seq for f in train]), dtype=np.float32)
    Xte = np.asarray(fn([f.seq for f in test]), dtype=np.float32)
    ytr = np.log1p(np.asarray([f.expression for f in train], dtype=np.float64))
    yte = [f.expression for f in test]
    sc = StandardScaler().fit(Xtr)
    m = Ridge(alpha=1.0).fit(sc.transform(Xtr), ytr)
    pred = m.predict(sc.transform(Xte))
    return {"feature": feature, "n_features": int(Xtr.shape[1]),
            "n_train": len(train), "n_test": len(test),
            "spearman": round(spearman(list(pred), yte), 4)}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cache-dir", default=str(DEFAULT_CACHE))
    ap.add_argument("--limit", type=int, default=20000,
                    help="subsample size for the MODEL arms (memory: 6-mers are 4096 features)")
    ap.add_argument("--n-blocks", type=int, default=10)
    ap.add_argument("--held-blocks", type=int, default=2)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--no-download", action="store_true")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    dl = not args.no_download
    lb, st_lb = load_fragments(args.cache_dir, which="LB", allow_download=dl)
    m9, st_m9 = load_fragments(args.cache_dir, which="M9", allow_download=dl)
    print(f"LB {st_lb.n_parsed} fragments ({st_lb.delimiter}-delimited, {st_lb.n_header_columns} cols, "
          f"len {st_lb.length_min}-{st_lb.length_max} median {st_lb.length_median})")
    print(f"M9 {st_m9.n_parsed} fragments ({st_m9.delimiter}-delimited)")

    # --- the ceiling, on the FULL data (free, and it reframes everything below) -----------------------
    ceil_lb = noise_ceiling(lb)
    ceil_m9 = noise_ceiling(m9)
    print(f"\nNOISE CEILING  LB: single-rep {ceil_lb['single_replicate_agreement']} -> "
          f"model ceiling {ceil_lb['ceiling_for_a_sequence_model']}")
    print(f"NOISE CEILING  M9: single-rep {ceil_m9['single_replicate_agreement']} -> "
          f"model ceiling {ceil_m9['ceiling_for_a_sequence_model']}")

    # --- Q2: condition, on the FULL shared set -------------------------------------------------------
    cond = condition_effect(lb, m9)
    print(f"\nCONDITION (n={cond['n_shared']}): within {cond['within_condition_single_vs_single']['mean']} "
          f"vs cross {cond['cross_condition_single_vs_single']['mean']} "
          f"(gap {cond['noise_matched_gap']:+.4f}) -> {cond['verdict']}")
    print(f"  unmatched ave-vs-ave would read {cond['unmatched_ave_vs_ave']} "
          f"(ceiling if irrelevant {cond['ceiling_if_condition_were_irrelevant']}); "
          f"disattenuated true cross-condition {cond['disattenuated_true_cross_condition']}")

    # --- Q1: the model, on a position-blocked split --------------------------------------------------
    import random as _r
    sub = lb if len(lb) <= args.limit else _r.Random(args.seed).sample(lb, args.limit)
    train, test, held = position_blocked_split(sub, n_blocks=args.n_blocks,
                                               held_blocks=args.held_blocks, seed=args.seed)
    print(f"\nposition-blocked split: {len(train)} train / {len(test)} test, held blocks {held} "
          f"of {args.n_blocks} (subsampled to {len(sub)} of {len(lb)})")
    arms = [_fit(train, test, f) for f in ("gc", "kmer4", "kmer6")]
    for a in arms:
        frac = a["spearman"] / ceil_lb["ceiling_for_a_sequence_model"]
        a["fraction_of_ceiling"] = round(frac, 4)
        print(f"  {a['feature']:7s} ({a['n_features']:5d} feats) rho={a['spearman']:+.4f}  "
              f"= {frac:.1%} of the {ceil_lb['ceiling_for_a_sequence_model']} ceiling")

    # --- Q1b: the LIKE-FOR-LIKE arm -- peak tiles selected from called peaks --------------------------
    # TWO filters, both measured, both load-bearing:
    #  (a) DESIGNED CONTROLS OUT. 106 pos_control + 470 neg_control + 936 random are chosen to be extreme,
    #      so leaving them in inflates any correlation with spread the model never had to earn.
    #  (b) FIXED 150 bp ONLY. Peak tiles are 48-150 (NOT a fixed 150, which I assumed and the module's own
    #      guard refused), but 97.4% are exactly 150 -- so this costs 2.6% and is what makes the positional
    #      feature valid, and therefore what makes this arm like-for-like with the designed grid.
    tiles_all, st_t = load_peak_tiles(args.cache_dir, allow_download=dl)
    n_ctrl = sum(1 for t in tiles_all if t.category != "tile")
    n_ragged = sum(1 for t in tiles_all if t.category == "tile" and len(t.seq) != 150)
    tiles = [t for t in tiles_all if t.category == "tile" and len(t.seq) == 150]
    ttrain, ttest, n_held = leave_peak_out(tiles, seed=args.seed)
    print(f"\npeak TILES: {st_t.n_parsed} loaded ({st_t.delimiter}-delimited, len {st_t.length_min}-"
          f"{st_t.length_max}); dropped {n_ctrl} designed controls + {n_ragged} non-150bp -> "
          f"{len(tiles)} kept")
    print(f"  leave-peak-out: {len(ttrain)} train / {len(ttest)} test ({n_held} peaks held)")
    tile_arms = []
    for f in ("gc", "kmer4", "kmer6", "onehot"):
        a = _fit(ttrain, ttest, f)
        tile_arms.append(a)
        print(f"  {a['feature']:7s} ({a['n_features']:5d} feats) rho={a['spearman']:+.4f}")

    verdict, detail = _verdict(arms, ceil_lb["ceiling_for_a_sequence_model"])
    # EACH ARM USES ITS OWN CEILING. The frag and tile subsets have very different reproducibility
    # (single-replicate 0.458 vs 0.795), so normalising tiles by the frag ceiling FLATTERS them -- it read
    # 34.9% of 0.7928 when the honest figure is 29.4% of 0.9410.
    ceil_t = noise_ceiling(tiles)
    print(f"  TILE noise ceiling: single-rep {ceil_t['single_replicate_agreement']} -> "
          f"model ceiling {ceil_t['ceiling_for_a_sequence_model']} "
          f"(NOT the frag ceiling {ceil_lb['ceiling_for_a_sequence_model']})")
    for a in tile_arms:
        a["fraction_of_ceiling"] = round(a["spearman"] / ceil_t["ceiling_for_a_sequence_model"], 4)
    tverdict, tdetail = _verdict(tile_arms, ceil_t["ceiling_for_a_sequence_model"])
    print(f"\nQ1b (peak tiles, like-for-like) VERDICT: {tverdict}\n  {tdetail}")
    artifact = {
        "schema": "glm-genomewide-oracle-v1",
        "date": str(date.today()),
        "substrate": "GEO GSE144621, sheared E. coli MG1655 genomic fragments, LB + M9",
        "load_stats": {"LB": st_lb.as_dict(), "M9": st_m9.as_dict()},
        "noise_ceiling": {"LB": ceil_lb, "M9": ceil_m9},
        "condition_effect": cond,
        "preregistered": PREREGISTERED,
        "split": {"kind": "position_blocked", "n_blocks": args.n_blocks, "held_blocks": held,
                  "n_train": len(train), "n_test": len(test), "subsample": len(sub),
                  "why": "91.5% of consecutive fragments OVERLAP, so a random split puts near-duplicates "
                         "on both sides; measured, not assumed"},
        "model_arms_random_fragments": arms,
        "q1_verdict": verdict,
        "q1_detail": detail,
        "peak_tiles": {"load_stats": st_t.as_dict(), "n_peaks_held": n_held,
                       "noise_ceiling": ceil_t,
                       "n_designed_controls_dropped": n_ctrl, "n_non_150bp_dropped": n_ragged,
                       "n_kept": len(tiles),
                       "n_train": len(ttrain), "n_test": len(ttest),
                       "split": "leave_peak_out -- tiles within a peak overlap heavily",
                       "arms": tile_arms, "verdict": tverdict, "detail": tdetail,
                       "why_this_is_the_like_for_like_arm":
                           "fixed 150 bp like the designed grid (so positional-onehot applies) AND "
                           "selected from called peaks, so it asks the grid's question -- how strong is "
                           "this promoter -- rather than the frag library's rare-event detection question"},
        "honest_limits": [
            "The model arms run on a SUBSAMPLE (6-mers are 4096 features); the ceiling and the condition "
            "effect use the FULL data.",
            "Fragments are 48-475 bp with 280 distinct lengths, so the designed grid's winning "
            "positional-onehot feature is NOT applicable without an anchor, and no anchor is obviously "
            "right. Only length-invariant k-mer features are used -- the comparison to the grid's 0.5914 "
            "one-hot number is therefore NOT like-for-like.",
            "The designed-grid numbers cannot be ceiling-normalised (GSE108535 has no replicates), so "
            "'fraction of ceiling' exists only on this substrate. Do not compare the two fractions.",
            "A fragment straddling a held-block boundary is assigned by midpoint, so a little sequence can "
            "cross the seam; with 10 blocks that is 2 seams over ~4.64 Mb.",
            "Expression is RNA_exp_ave, which is ALREADY DNA-normalised (measured: spearman with DNA_ave "
            "= -0.039). It is NOT divided by DNA_ave again.",
            "97% of fragments carry a single integrated barcode, which is why the replicate-derived "
            "ceiling is low (~0.79) and why it must be quoted with every model number.",
        ],
    }
    out = Path(args.out) if args.out else Path(f"wiki/glm_genomewide_oracle_{date.today()}.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(artifact, indent=2), encoding="utf-8")
    print(f"\nQ1 VERDICT: {verdict}\n  {detail}\nwrote {out}")
    return 0


def _verdict(arms: list[dict], ceiling: float) -> tuple[str, str]:
    gc = next(a for a in arms if a["feature"] == "gc")["spearman"]
    best = max((a for a in arms if a["feature"] != "gc"), key=lambda a: a["spearman"])
    margin = best["spearman"] - gc
    frac = best["spearman"] / ceiling if ceiling else 0.0
    d = (f"best learned = {best['feature']} rho {best['spearman']:+.4f} ({frac:.1%} of the {ceiling} "
         f"ceiling); GC {gc:+.4f}; margin {margin:+.4f}")
    if margin < PREREGISTERED["q1_bars"]["beats_crude_baseline_by"]:
        return "DOES_NOT_GENERALISE", d
    if frac < PREREGISTERED["q1_bars"]["fraction_of_ceiling_for_GENERALISES"]:
        return "WEAK_BUT_REAL", d
    return "GENERALISES_OFF_THE_DESIGNED_GRID", d


if __name__ == "__main__":
    raise SystemExit(main())
