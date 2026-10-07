"""Train + evaluate the bacterial expression oracle on real measured MPRA data. Writes an artifact.

    uv run python scripts/glm_expression_oracle.py

Scores every feature set under BOTH a random split and leave-one-element-out splits. Both are reported
side by side on purpose: this repo has measured three separate times that random splits inflate results
~3x, so showing only the random number would be the exact failure the comparison exists to expose.

Substrate: GEO GSE108535 (Urtecho et al. 2019), 10,898 sequence->expression pairs, verified by download.
No GPU, no network after the first run (~227 KB + 8.5 MB cached on D:).
"""
from __future__ import annotations

import argparse
import json
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
WIKI = ROOT / "wiki"


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cache-dir", default="D:/dna_decode_cache/mpra")
    p.add_argument("--features", default="gc,kmer4,kmer6,onehot")
    p.add_argument("--alpha", type=float, default=1.0)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--out", default=None)
    a = p.parse_args(argv)

    from dna_decode.glm.expression import (
        ELEMENT_COLUMNS,
        fit_predict_ridge,
        leave_element_out,
        load_pairs,
        random_split,
    )

    pairs, rep = load_pairs(a.cache_dir)
    print("=" * 78)
    print(f"substrate: GSE108535  ->  {rep.n_joined} sequence->expression pairs "
          f"(join rate {rep.join_rate:.3f})")
    print(f"  a naive join without the hyphen/underscore fix would give {rep.naive_join} "
          f"({100 * rep.naive_join / rep.n_expression_rows:.1f}%)")
    vals = [x.expression for x in pairs]
    print(f"  expression: min {min(vals):.4f}  median {sorted(vals)[len(vals)//2]:.4f}  "
          f"max {max(vals):.4f}  ({max(vals)/min(vals):.0f}x)")
    for c in ELEMENT_COLUMNS:
        print(f"  {c:12} {len({x.elements.get(c,'') for x in pairs}):4d} distinct values")
    print("=" * 78)

    feats = [f.strip() for f in a.features.split(",")]
    results = []

    # ---- random split: reported for comparison, NOT the headline ----------------------------------
    tr, te = random_split(pairs, seed=a.seed)
    print(f"\n[RANDOM split]  train {len(tr)} / test {len(te)}   <- INFLATED, for comparison only")
    for f in feats:
        r = fit_predict_ridge(tr, te, feature=f, alpha=a.alpha)
        r["split"] = "random"
        r["split_axis"] = None
        results.append(r)
        print(f"   {f:8} n_feat {r['n_features']:6d}   spearman {r['spearman']:+.4f}")

    # ---- leave-one-element-out: the honest generalisation test ------------------------------------
    for col in ELEMENT_COLUMNS:
        try:
            tr, te, held = leave_element_out(pairs, col, seed=a.seed)
        except ValueError as e:
            print(f"\n[LEAVE-OUT {col}] SKIPPED: {e}")
            continue
        if len(te) < 50 or len(tr) < 200:
            print(f"\n[LEAVE-OUT {col}] SKIPPED: train {len(tr)} / test {len(te)} too small")
            continue
        print(f"\n[LEAVE-OUT {col}]  train {len(tr)} / test {len(te)}   held out {len(held)} values")
        for f in feats:
            r = fit_predict_ridge(tr, te, feature=f, alpha=a.alpha)
            r["split"] = "leave_element_out"
            r["split_axis"] = col
            r["held_out_values"] = held
            results.append(r)
            print(f"   {f:8} n_feat {r['n_features']:6d}   spearman {r['spearman']:+.4f}")

    # ---- the comparison that matters -------------------------------------------------------------
    print("\n" + "=" * 78)
    print("INFLATION: random vs held-out element, per feature")
    print("=" * 78)
    for f in feats:
        rnd = next((r["spearman"] for r in results if r["feature"] == f and r["split"] == "random"), None)
        lo = [r["spearman"] for r in results if r["feature"] == f and r["split"] == "leave_element_out"]
        if rnd is None or not lo:
            continue
        mean_lo = sum(lo) / len(lo)
        print(f"  {f:8} random {rnd:+.4f}   held-out mean {mean_lo:+.4f}   "
              f"inflation {rnd - mean_lo:+.4f}"
              + ("   <- random OVERSTATES" if rnd > mean_lo else ""))

    today = date.today().isoformat()
    art = {
        "record": "glm-expression-oracle-v1",
        "date": today,
        "substrate": {
            "accession": "GSE108535",
            "citation": ("Urtecho G, Tripp AD, Insigne KD, Kim H, Kosuri S. Biochemistry "
                         "2019;58(11):1539-1551. DOI 10.1021/acs.biochem.7b01069"),
            "assay": "genomically-integrated MPRA, RNA/DNA barcode ratio, fully-crossed sigma70 grid",
            **rep.as_dict(),
        },
        "alpha": a.alpha,
        "results": results,
        "honest_limits": [
            "ONE organism (E. coli), ONE assay, ONE designed grid. A designed grid is not a natural "
            "promoter distribution, so this bounds grid generalisation, not genome-wide performance.",
            "Random-split numbers are reported ONLY for comparison. This repo has measured ~3x inflation "
            "from random splits three separate times; the leave-element-out numbers are the headline.",
            "Ridge on fixed features is deliberately the TRIVIAL HEAD half of AlphaGenome's design "
            "principle. No learned representation is trained here, so these are a FLOOR, not a ceiling.",
            "Spearman is rank-based and so invariant to the log transform; the FIT is not, which is why "
            "log_target defaults on (expression spans 615x).",
            "No confidence intervals on the Spearman differences.",
        ],
    }
    out = Path(a.out) if a.out else WIKI / f"glm_expression_oracle_{today}.json"
    out.write_text(json.dumps(art, indent=2), encoding="utf-8")
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
