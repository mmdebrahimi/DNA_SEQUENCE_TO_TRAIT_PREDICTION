"""The LAST planned arm on the tile-chain ceiling criterion: change the SUBSTRATE, not the method.

Graded by `wiki/eukaryotic_arm_acceptance_bar.json`, frozen before any number here existed. The
question is whether the bacterial rung's ceiling (best arm 0.5076 against a 0.941 replicate
ceiling, after ten method arms) was the METHOD's limit or the SUBSTRATE's. So this runs the SAME
learned-encoder-vs-dumb-baseline contrast on a constructed human lentiMPRA panel and compares the
CONTRAST, because a raw Spearman cannot cross a substrate boundary -- two panels have different
noise floors, assays and dynamic range.

Bacterial anchor, read from the committed arms ledger rather than retyped:
    conv-3width-13k-params-150bp 0.3234 - gc-single-feature 0.3112 = +0.0122

**Three traps this script exists to not fall into.**

1. `activity` is `log2(RNA/DNA)` and GOES NEGATIVE, while both incumbent fitters apply `np.log1p`
   to the target -- undefined at <= -1, and `spearman` maps NaN to 0.0, so the failure would
   surface as a clean-looking null rather than as an error. Every record is therefore shifted by a
   single constant to `expression = activity - min + 1.0`, computed over the WHOLE cell once (never
   per fold, which would make folds incomparable). Spearman is invariant to monotone transforms of
   the target, so the EVALUATION is unaffected; only the fitting scale moves, and that is named in
   the artifact rather than left implicit.

2. `tile_peak_key` REFUSES a record with no `.peak`, which is correct -- a silent fallback would
   produce a RANDOM inner split on a substrate that needs a blocked one. The leakage unit here is
   the published fold, so `group_key` is passed explicitly.

3. The published incumbents name no metric, so the frozen bar refuses the head-to-head. They are
   carried through verbatim with `metric_unstated: true` and inform no predicate.

The verdict is computed by a pure function of the measurements plus the frozen anchor, so a
reader can re-derive it from the artifact without re-running anything.
"""

from __future__ import annotations

import argparse
import json
import time
from dataclasses import dataclass
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
BAR_PATH = REPO / "wiki/eukaryotic_arm_acceptance_bar.json"

N_BOOT = 10_000
BOOT_SEED = 606


@dataclass(frozen=True)
class Row:
    """Adapter onto the shape both incumbent fitters expect (`.seq`, `.expression`)."""
    name: str
    seq: str
    expression: float
    fold: int


def shift_to_positive(records) -> tuple[list, float]:
    """One constant over the WHOLE cell, so folds stay comparable. See trap 1 in the docstring."""
    lo = min(r.activity for r in records)
    off = 1.0 - lo
    rows = [Row(name=r.name, seq=r.seq, expression=r.activity + off, fold=r.fold)
            for r in records]
    assert min(x.expression for x in rows) > 0.0, "shift failed; log1p would produce NaN"
    return rows, off


def paired_bootstrap(gains: list[float], n_boot: int = N_BOOT, seed: int = BOOT_SEED) -> dict:
    """Bootstrap the FOLD-LEVEL gains, which is the pairing the frozen bar names."""
    import numpy as np

    g = np.asarray(gains, dtype=float)
    rng = np.random.default_rng(seed)
    draws = rng.choice(g, size=(n_boot, len(g)), replace=True).mean(axis=1)
    return {
        "mean": float(g.mean()),
        "ci_lo": float(np.percentile(draws, 2.5)),
        "ci_hi": float(np.percentile(draws, 97.5)),
        "n_folds": len(g),
        "n_boot": n_boot,
        "seed": seed,
    }


def classify(ci_lo: float, ci_hi: float, anchor: float) -> tuple[str, str]:
    """Pure verdict. Ordered, first match wins; mirrors the frozen bar's predicates exactly."""
    if ci_lo > anchor:
        return ("SUBSTRATE_LIMITED",
                "the CI for the human gain lies entirely ABOVE the bacterial gain %.4f" % anchor)
    if ci_hi < anchor:
        return ("METHOD_LIMITED",
                "the CI lies entirely BELOW the bacterial gain %.4f" % anchor)
    return ("SUBSTRATE_NEUTRAL",
            "the CI contains the bacterial gain %.4f, so no difference is detectable" % anchor)


def run_fold(rows: list, fold: int, *, device: str, seed: int, widths) -> dict:
    """GC-only ridge and a retrained conv encoder on IDENTICAL train/test rows."""
    from dna_decode.glm.encoder import fit_encoder
    from dna_decode.glm.expression import fit_predict_ridge

    train = [r for r in rows if r.fold != fold]
    test = [r for r in rows if r.fold == fold]
    if not train or not test:
        raise RuntimeError("fold %r leaves an empty side (train=%d test=%d)"
                           % (fold, len(train), len(test)))

    gc = fit_predict_ridge(train, test, feature="gc")
    conv = fit_encoder(train, test, widths=widths, seed=seed, device=device,
                       group_key=lambda r: str(r.fold))
    conv.pop("predictions", None)  # one float per test row; never serialised

    return {
        "fold": fold,
        "n_train": len(train),
        "n_test": len(test),
        "gc_spearman": round(float(gc["spearman"]), 4),
        "conv_spearman": round(float(conv["spearman"]), 4),
        "gain": round(float(conv["spearman"]) - float(gc["spearman"]), 4),
        "conv_protocol": {k: v for k, v in conv.items()
                          if k in ("epochs_run", "n_train", "n_val", "widths", "channels",
                                   "device", "early_stopped")},
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cell", default="K562", choices=("K562", "HepG2"))
    ap.add_argument("--folds", type=int, default=0,
                    help="cap the number of folds (0 = all). A capped run is stamped PARTIAL and "
                         "its verdict is withheld -- the bar names a 10-fold pairing.")
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args(argv)

    bar = json.loads(BAR_PATH.read_text(encoding="utf-8"))
    anchor = bar["anchors"]["conv_minus_gc"]

    from dna_decode.glm.encoder import validate_widths
    from dna_decode.glm.lentimpra import load_cell

    t0 = time.time()
    records, stats = load_cell(args.cell)
    rows, offset = shift_to_positive(records)
    length = len(rows[0].seq)
    widths = validate_widths((6, 20, 75), length)

    folds = sorted({r.fold for r in rows})
    if args.folds:
        folds = folds[:args.folds]
    partial = bool(args.folds) and len(folds) < len({r.fold for r in rows})

    print("cell=%s n=%d len=%d folds=%s join_rate=%.4f widths=%s offset=%+.4f"
          % (args.cell, len(rows), length, folds, stats.join_rate, widths, offset), flush=True)

    per_fold = []
    for f in folds:
        r = run_fold(rows, f, device=args.device, seed=args.seed, widths=widths)
        per_fold.append(r)
        print("  fold %-3s gc %.4f  conv %.4f  gain %+.4f   (%.0fs)"
              % (f, r["gc_spearman"], r["conv_spearman"], r["gain"], time.time() - t0), flush=True)

    boot = paired_bootstrap([r["gain"] for r in per_fold])
    verdict, detail = classify(boot["ci_lo"], boot["ci_hi"], anchor)

    rec = {
        "schema": "eukaryotic-arm-result-v1",
        "date": "2026-11-01",
        "bar": "wiki/eukaryotic_arm_acceptance_bar.json",
        "criterion": bar["criterion"],
        "arm": bar["arm"],
        "cell": args.cell,
        "substrate": {
            "n_rows": len(rows),
            "sequence_length": length,
            "join_rate": round(stats.join_rate, 4),
            "n_designed": stats.n_designed,
            "n_measured": stats.n_measured,
            "n_joined": stats.n_joined,
            "sequence_column": stats.sequence_column,
            "n_folds_published": stats.n_folds,
        },
        "target_reparameterization": {
            "applied": "expression = activity + %.6f" % offset,
            "why": ("activity is log2(RNA/DNA) and goes negative; both incumbent fitters apply "
                    "np.log1p to the target, which is undefined at <= -1, and spearman maps NaN to "
                    "0.0 -- so without this the failure would read as a clean null"),
            "computed_over": "the whole cell, once (never per fold, which would make folds "
                             "incomparable)",
            "evaluation_unaffected": ("Spearman is invariant to monotone transforms of the target; "
                                      "only the fitting scale moves"),
        },
        "per_fold": per_fold,
        "medians": {
            "gc": round(sorted(r["gc_spearman"] for r in per_fold)[len(per_fold) // 2], 4),
            "conv": round(sorted(r["conv_spearman"] for r in per_fold)[len(per_fold) // 2], 4),
        },
        "gain_bootstrap": {k: (round(v, 4) if isinstance(v, float) else v)
                           for k, v in boot.items()},
        "bacterial_anchor": bar["anchors"],
        "verdict": "PARTIAL_WITHHELD" if partial else verdict,
        "verdict_detail": ("capped at %d of %d folds; the bar names a 10-fold pairing, so the "
                           "verdict is withheld rather than reported on a different pairing"
                           % (len(folds), stats.n_folds)) if partial else detail,
        "provisional_verdict_if_complete": verdict if partial else None,
        "refused": bar["refused_comparisons"],
        "honest_limits": bar["honest_limits_fixed_in_advance"],
        "wall_seconds": round(time.time() - t0, 1),
    }

    out = REPO / ("wiki/eukaryotic_arm_result_%s_2026-11-01.json" % args.cell)
    out.write_text(json.dumps(rec, indent=2, sort_keys=True), encoding="utf-8")

    print()
    print("gc median   %.4f" % rec["medians"]["gc"])
    print("conv median %.4f" % rec["medians"]["conv"])
    print("gain        mean %+.4f  CI [%+.4f, %+.4f]  over %d folds"
          % (boot["mean"], boot["ci_lo"], boot["ci_hi"], boot["n_folds"]))
    print("anchor      %+.4f (bacterial conv - gc)" % anchor)
    print("VERDICT     %s" % rec["verdict"])
    print("            %s" % rec["verdict_detail"])
    print("-> %s" % out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
