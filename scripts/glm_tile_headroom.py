"""Is there POSITIONAL headroom on genomic promoter tiles, or is 0.3000 a composition ceiling?

    uv run python scripts/glm_tile_headroom.py [--seeds 10] [--no-download]

**This is a STOP GATE and it can close G-A in an afternoon.** G-A proposes a learned convolutional encoder,
and a conv encoder is fundamentally a *composition-plus-position* model. On real genomic peak tiles a single
GC number already scores 0.3000 and BEATS every fixed learned feature (one-hot 0.2691, 4-mer 0.2769, 6-mer
0.2022). So before building an encoder, measure the thing the encoder's whole advantage rests on: **does
adding POSITION to composition buy anything at all here?**

If it does not, the encoder's main lever is dead before it is built, and the honest move is to say so rather
than spend the build and discover it.

**FOUR ARMS, chosen so the composition/position contrast is clean:**

| arm | features | what it isolates |
|---|---|---|
| `gc` | 1 | the incumbent bar, and the reconcile anchor |
| `gc_dinuc` | 17 | composition only, richer (GC + 16 dinucleotide frequencies) |
| `gc_posgc` | 11 | position only, minimal (GC + GC-per-bin) |
| `dinuc_posgc` | 26 | composition AND position together |

**The pre-registered question:** `best_position_aware - best_composition_only`. If position adds less than
`MIN_POSITION_GAIN`, hypothesis GA-H3 is supported and `recommend_user_ratification_before_cnn` is set — a
MACHINE-READABLE flag, so the gate is not a judgement buried in prose.

**RECONCILE GATE (decision D1's spirit, applied to an input rather than an output).** Before any arm is
interpreted, the `gc` arm on seed 0 must reproduce the committed **0.3000 ± 0.01**. If it does not, this
script is reading a different split or different data than the published number, and it REFUSES rather than
reporting — a divergence here would make every comparison below meaningless while looking fine.

Multi-seed per decision D2: compute buys rigour, not search. Verdict on the median, dispersion reported.
"""
from __future__ import annotations

import argparse
import json
import statistics
from datetime import date
from pathlib import Path

from dna_decode.glm.genomewide import (
    DEFAULT_CACHE,
    gc_features,
    kmer_features,
    leave_peak_out,
    load_peak_tiles,
    noise_ceiling,
    spearman,
)

#: Frozen BEFORE the run.
PREREGISTERED = {
    "question": "does adding POSITION to composition improve prediction of genomic promoter-tile "
                "expression, i.e. does a conv encoder's core advantage exist on this substrate at all?",
    "min_position_gain": 0.02,
    "meaning": "if (best position-aware arm) - (best composition-only arm) < 0.02, position buys nothing "
               "measurable here, GA-H3 is supported, and building a conv encoder is a predicted-negative",
    "reconcile_anchor": {"arm": "gc", "seed": 0, "expected": 0.3000, "tolerance": 0.01,
                         "source": "wiki/glm_genomewide_oracle_2026-10-07.json"},
    "composition_only_arms": ["gc", "gc_dinuc"],
    "position_aware_arms": ["gc_posgc", "dinuc_posgc"],
    "verdicts": {
        "POSITION_HAS_HEADROOM": "position adds >= 0.02 -- the encoder's lever is real, proceed to build it",
        "NO_POSITIONAL_HEADROOM": "position adds < 0.02 -- GA-H3 supported, encoder is a predicted-negative, "
                                  "recommend user ratification before spending the build",
        "INDETERMINATE_RECONCILE_FAILED": "the gc arm did not reproduce the committed 0.3000 on seed 0, so "
                                          "this script is not reading the same data as the published number "
                                          "and NOTHING is interpretable",
    },
    "thresholds_are": "ASSERTED; the 0.9410 ceiling that contextualises them is DERIVED from the tiles' own "
                      "replicates",
}

N_BINS = 10


def positional_gc(seqs: list[str], *, bins: int = N_BINS) -> list[list[float]]:
    """GC content per positional bin. The MINIMAL positional feature.

    Deliberately minimal: if even this adds nothing over plain GC, a 600-feature one-hot or a conv stack is
    very unlikely to find positional signal that this cannot see. It is the cheapest possible test of the
    encoder's premise.
    """
    out = []
    for s in seqs:
        n = len(s)
        v = []
        for b in range(bins):
            lo, hi = (b * n) // bins, ((b + 1) * n) // bins
            seg = s[lo:hi]
            v.append(((seg.count("G") + seg.count("C")) / len(seg)) if seg else 0.0)
        out.append(v)
    return out


def _concat(*blocks: list[list[float]]) -> list[list[float]]:
    return [sum(parts, []) for parts in zip(*blocks)]


ARMS = {
    "gc": lambda s: gc_features(s),
    "gc_dinuc": lambda s: _concat(gc_features(s), kmer_features(s, 2)),
    "gc_posgc": lambda s: _concat(gc_features(s), positional_gc(s)),
    "dinuc_posgc": lambda s: _concat(kmer_features(s, 2), positional_gc(s)),
}


def _fit(train, test, arm: str) -> float:
    import numpy as np
    from sklearn.linear_model import Ridge
    from sklearn.preprocessing import StandardScaler
    fn = ARMS[arm]
    Xtr = np.asarray(fn([t.seq for t in train]), dtype=np.float64)
    Xte = np.asarray(fn([t.seq for t in test]), dtype=np.float64)
    ytr = np.log1p(np.asarray([t.expression for t in train], dtype=np.float64))
    yte = [t.expression for t in test]
    sc = StandardScaler().fit(Xtr)
    m = Ridge(alpha=1.0).fit(sc.transform(Xtr), ytr)
    return spearman(list(m.predict(sc.transform(Xte))), yte)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cache-dir", default=str(DEFAULT_CACHE))
    ap.add_argument("--seeds", type=int, default=10)
    ap.add_argument("--no-download", action="store_true")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    tiles_all, st = load_peak_tiles(args.cache_dir, allow_download=not args.no_download)
    # the SAME two filters the baseline used, or the reconcile gate cannot possibly pass
    tiles = [t for t in tiles_all if t.category == "tile" and len(t.seq) == 150]
    ceil = noise_ceiling(tiles)["ceiling_for_a_sequence_model"]
    print(f"tiles: {len(tiles)} kept of {st.n_parsed} | ceiling {ceil}")

    seeds = list(range(args.seeds))
    per_arm: dict[str, dict] = {}
    for arm in ARMS:
        vals = []
        for s in seeds:
            tr, te, _ = leave_peak_out(tiles, seed=s)
            vals.append(_fit(tr, te, arm))
        med = statistics.median(vals)
        per_arm[arm] = {
            "n_features": len(ARMS[arm]([tiles[0].seq])[0]),
            "median": round(med, 4),
            "stdev": round(statistics.stdev(vals), 4) if len(vals) > 1 else 0.0,
            "min": round(min(vals), 4), "max": round(max(vals), 4),
            "fraction_of_ceiling": round(med / ceil, 4) if ceil else None,
            "per_seed": [round(v, 4) for v in vals],
        }
        print(f"  {arm:14s} ({per_arm[arm]['n_features']:3d} feats) median {med:+.4f} "
              f"sd {per_arm[arm]['stdev']:.4f}  = {per_arm[arm]['fraction_of_ceiling']:.1%} of ceiling")

    verdict, detail, ratify = _verdict(per_arm, ceil)
    artifact = {
        "schema": "glm-tile-headroom-v1",
        "date": str(date.today()),
        "substrate": "GEO GSE144621 peak tiles, category=tile, len=150, leave-peak-out",
        "n_tiles": len(tiles),
        "load_stats": st.as_dict(),
        "noise_ceiling": ceil,
        "seeds": seeds,
        "preregistered": PREREGISTERED,
        "arms": per_arm,
        "verdict": verdict,
        "verdict_detail": detail,
        "recommend_user_ratification_before_cnn": ratify,
        "honest_limits": [
            "Positional GC in 10 bins is the MINIMAL positional feature. A conv encoder could in principle "
            "find positional structure this cannot see -- a motif at a precise offset rather than a bulk "
            "shift in one bin. So a null here makes the encoder a predicted-negative, it does not PROVE "
            "the encoder must fail.",
            "Ridge is linear in these features; an interaction between composition and position is not "
            "representable here and is one thing a conv encoder could add.",
            "The 2/3-inactive composition of the tile set applies to every arm equally, so it cannot "
            "explain a DIFFERENCE between arms -- but it does bound all of them.",
            "min_position_gain = 0.02 is ASSERTED, not derived.",
        ],
    }
    out = Path(args.out) if args.out else Path(f"wiki/glm_tile_headroom_{date.today()}.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(artifact, indent=2), encoding="utf-8")
    print(f"\nVERDICT: {verdict}\n  {detail}")
    print(f"recommend_user_ratification_before_cnn: {ratify}")
    print(f"wrote {out}")
    return 0


def _verdict(per_arm: dict, ceil: float) -> tuple[str, str, bool]:
    """RECONCILE FIRST. A script reading different data than the published number cannot grade anything."""
    anc = PREREGISTERED["reconcile_anchor"]
    gc_seed0 = per_arm[anc["arm"]]["per_seed"][anc["seed"]]
    if abs(gc_seed0 - anc["expected"]) > anc["tolerance"]:
        return ("INDETERMINATE_RECONCILE_FAILED",
                f"gc arm on seed {anc['seed']} gave {gc_seed0:.4f}, expected {anc['expected']} "
                f"+/-{anc['tolerance']} from {anc['source']}. This script is not reading the same split or "
                f"data as the published number; nothing below is interpretable.", False)
    comp = max(PREREGISTERED["composition_only_arms"], key=lambda a: per_arm[a]["median"])
    pos = max(PREREGISTERED["position_aware_arms"], key=lambda a: per_arm[a]["median"])
    gain = per_arm[pos]["median"] - per_arm[comp]["median"]
    d = (f"best composition-only = {comp} {per_arm[comp]['median']:+.4f}; best position-aware = {pos} "
         f"{per_arm[pos]['median']:+.4f}; POSITION GAIN {gain:+.4f} against a bar of "
         f"{PREREGISTERED['min_position_gain']}; reconcile OK (gc seed0 {gc_seed0:.4f})")
    if gain >= PREREGISTERED["min_position_gain"]:
        return "POSITION_HAS_HEADROOM", d, False
    return "NO_POSITIONAL_HEADROOM", d, True


if __name__ == "__main__":
    raise SystemExit(main())
