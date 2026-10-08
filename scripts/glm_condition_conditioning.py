"""Does ONE model conditioned on growth medium beat TWO independent per-medium models?

    uv run python scripts/glm_condition_conditioning.py [--seeds 10] [--no-download]

G-C. The bacterial analogue of AlphaGenome's cell-state conditioning, and the one item from that critique
that was both new and already evidenced: on 297,868 coordinate-keyed fragments shared between LB and M9 the
noise-matched cross-condition agreement is **0.4062** against within-condition **0.4676** (gap **-0.0614**),
and the disattenuated true cross-condition correlation is **0.7887** — so roughly 38% of the real variance
is condition-specific. Something is there to capture. The question is whether a shared model captures it.

**THE REVIEW FINDING THAT RESHAPED THIS SCRIPT.** The first design made the primary arm
`[sequence ‖ medium indicator]`. With a linear head that learns only a GLOBAL per-medium offset — it cannot
make sequence effects differ by medium. Meanwhile the two-model comparator fits *separate weight vectors per
medium*, so it is **strictly more expressive on exactly the axis conditioning is supposed to exploit**. The
design was rigged against itself: with the data-matched arm removing the one remaining advantage (the shared
model sees ~2x the rows), it could essentially only return `GAIN_IS_DATA_VOLUME` or `NO_GAIN`.

So the PRIMARY arm now carries `sequence × medium` INTERACTIONS, and the indicator arm is retained
deliberately as a diagnostic — to show that an indicator alone cannot win.

**FOUR ARMS:**

| arm | features | role |
|---|---|---|
| `two_model` | per-medium weights | the comparator to beat |
| `interaction` | `[seq ‖ medium ‖ seq×medium]` | **PRIMARY** — can represent condition-specific response |
| `indicator` | `[seq ‖ medium]` | diagnostic; expected to fail, and that failure is the point |
| `shuffled` | interaction, medium label permuted | **validity** — must show no gain |
| `data_matched` | interaction, row-count matched | **validity** — isolates conditioning from data volume |

**WHY A 3-MER SEQUENCE REPRESENTATION (64 features).** Not a compromise — it is what the substrate
indicated. Step 3 measured on this project's own tiles that ADDING features makes prediction WORSE (1
feature beat 11, 17 and 26). A 4-mer interaction model over 595k rows would also need ~1.2 GB for the design
matrix alone, which would force a subsample and violate decision D2 (use the full data). 3-mers keep the
full 297,868 shared fragments at ~300 MB.

**TWO METRICS, both reported, and the distinction matters.** `pooled` Spearman over all (fragment, medium)
test pairs is where a global offset can show value; `per_medium_mean` is the within-medium average, where
only condition-specific SEQUENCE response can show value. A gain that appears only in `pooled` is an offset
effect, not the conditioning this family is about.
"""
from __future__ import annotations

import argparse
import json
import statistics
from datetime import date
from pathlib import Path

from dna_decode.glm.genomewide import (
    DEFAULT_CACHE,
    condition_effect,
    kmer_features,
    load_fragments,
    noise_ceiling,
    position_blocked_split,
    spearman,
)

#: Frozen BEFORE any run.
PREREGISTERED = {
    "question": "does one model taking (sequence, medium) beat two independent per-medium models on "
                "held-out position-blocked data?",
    "primary_arm": "interaction",
    "comparator": "two_model",
    "primary_metric": "per_medium_mean",
    "secondary_metric": "pooled",
    "min_gain": 0.02,
    "validity": {
        "shuffled_max_gain": 0.02,
        "data_matched_must_also_clear": True,
    },
    "verdicts": {
        "CONDITIONING_HELPS": "interaction beats two_model by >= 0.02 AND shuffled shows no gain AND "
                              "data_matched also clears -- conditioning, not data volume",
        "GAIN_IS_DATA_VOLUME": "interaction beats two_model but data_matched does NOT clear -- the gain is "
                               "attributable to seeing ~2x the rows",
        "NO_GAIN": "interaction does not beat two_model by the margin",
        "INDETERMINATE_NULL_NOT_CLEAN": "the medium-shuffled arm itself gains >= 0.02, so the pipeline "
                                        "fabricates signal and NOTHING is graded",
    },
    "min_shared_fragments": 100_000,
    "sequence_features": "kmer3 (64) -- chosen because Step 3 measured that MORE features hurt on this "
                         "project's substrate, and because 4-mers would force a subsample (decision D2 says "
                         "use the full data)",
    "thresholds_are": "ASSERTED; the per-medium ceilings that contextualise them are DERIVED",
}

K = 3
MEDIA = ("LB", "M9")


def _seq_features(seqs: list[str]) -> list[list[float]]:
    return kmer_features(seqs, K)


def build_design(frags, medium_flags: list[float], *, arm: str):
    """Feature matrix per arm. `medium_flags` is 0.0 for LB and 1.0 for M9, aligned with `frags`.

    `interaction` is the load-bearing one: appending `seq * medium` lets the model give a sequence feature a
    DIFFERENT coefficient per medium, which an indicator alone cannot do.
    """
    seqs = [f.seq for f in frags]
    S = _seq_features(seqs)
    if arm == "sequence_only":
        return S
    if arm == "indicator":
        return [s + [m] for s, m in zip(S, medium_flags)]
    if arm == "interaction":
        return [s + [m] + [x * m for x in s] for s, m in zip(S, medium_flags)]
    raise ValueError(f"unknown arm {arm!r}")


def fit_ridge(Xtr, ytr, Xte):
    import numpy as np
    from sklearn.linear_model import Ridge
    from sklearn.preprocessing import StandardScaler
    Xtr = np.asarray(Xtr, dtype=np.float32)
    Xte = np.asarray(Xte, dtype=np.float32)
    y = np.log1p(np.asarray(ytr, dtype=np.float64))
    sc = StandardScaler().fit(Xtr)
    m = Ridge(alpha=1.0).fit(sc.transform(Xtr), y)
    return [float(v) for v in m.predict(sc.transform(Xte))]


def load_shared(cache_dir: str, *, allow_download: bool):
    """LB and M9 intersected on the COORDINATE key, so a repeated genomic element cannot collapse.

    Measured: 340 LB and 302 M9 sequences occur at more than one locus, and every one maps to a different
    coordinate triple. A sequence-only key would block-assign such a fragment by whichever locus survived,
    letting the two loci land on opposite sides of a position-blocked split.
    """
    lb, st_lb = load_fragments(cache_dir, which="LB", allow_download=allow_download)
    m9, st_m9 = load_fragments(cache_dir, which="M9", allow_download=allow_download)
    key = lambda f: (f.start, f.end, f.strand, f.seq)
    ka = {key(f): f for f in lb}
    kb = {key(f): f for f in m9}
    shared = sorted(set(ka) & set(kb))
    if len(shared) < PREREGISTERED["min_shared_fragments"]:
        raise SystemExit(f"REFUSE: only {len(shared)} shared fragments, floor is "
                         f"{PREREGISTERED['min_shared_fragments']}. A degraded intersection must not be "
                         f"scored.")
    return [ka[k] for k in shared], [kb[k] for k in shared], st_lb, st_m9


def per_medium_and_pooled(pred_by_medium: dict, truth_by_medium: dict) -> dict:
    """Both metrics. A gain visible only in `pooled` is a global-offset effect, not conditioning."""
    per = {m: spearman(pred_by_medium[m], truth_by_medium[m]) for m in MEDIA}
    pooled_p = [v for m in MEDIA for v in pred_by_medium[m]]
    pooled_t = [v for m in MEDIA for v in truth_by_medium[m]]
    return {
        "per_medium": {m: round(per[m], 4) for m in MEDIA},
        "per_medium_mean": round(sum(per.values()) / len(per), 4),
        "pooled": round(spearman(pooled_p, pooled_t), 4),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cache-dir", default=str(DEFAULT_CACHE))
    ap.add_argument("--seeds", type=int, default=10)
    ap.add_argument("--n-blocks", type=int, default=10)
    ap.add_argument("--held-blocks", type=int, default=2)
    ap.add_argument("--no-download", action="store_true")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    A, B = {}, {}
    lb, m9, st_lb, st_m9 = load_shared(args.cache_dir, allow_download=not args.no_download)
    frag = {"LB": lb, "M9": m9}
    print(f"shared fragments (coordinate key): {len(lb)}")

    # PER-MEDIUM ceilings, never a shared normaliser -- LB 0.7928 vs M9 0.8014 differ
    ceilings = {m: noise_ceiling(frag[m]) for m in MEDIA}
    for m in MEDIA:
        print(f"  {m} ceiling {ceilings[m]['ceiling_for_a_sequence_model']} "
              f"(single-replicate {ceilings[m]['single_replicate_agreement']})")
    cond = condition_effect(lb, m9)
    print(f"  condition effect: gap {cond['noise_matched_gap']:+.4f}, "
          f"disattenuated {cond['disattenuated_true_cross_condition']}, verdict {cond['verdict']}")

    seeds = list(range(args.seeds))
    two_model_runs = []
    for s in seeds:
        # ONE split per seed, derived from LB coordinates and applied to BOTH media by index, so every arm
        # sees the identical partition and arms differ only in the model.
        tr_lb, te_lb, held = position_blocked_split(frag["LB"], n_blocks=args.n_blocks,
                                                    held_blocks=args.held_blocks, seed=s)
        tr_idx = {(f.start, f.end, f.strand) for f in tr_lb}
        idx_tr = [i for i, f in enumerate(frag["LB"]) if (f.start, f.end, f.strand) in tr_idx]
        idx_te = [i for i, f in enumerate(frag["LB"]) if (f.start, f.end, f.strand) not in tr_idx]

        pred, truth = {}, {}
        for m in MEDIA:
            fr = frag[m]
            tr = [fr[i] for i in idx_tr]
            te = [fr[i] for i in idx_te]
            X_tr = build_design(tr, [0.0] * len(tr), arm="sequence_only")
            X_te = build_design(te, [0.0] * len(te), arm="sequence_only")
            pred[m] = fit_ridge(X_tr, [f.expression for f in tr], X_te)
            truth[m] = [f.expression for f in te]
        r = per_medium_and_pooled(pred, truth)
        r["held_blocks"] = held
        r["n_train"], r["n_test"] = len(idx_tr), len(idx_te)
        two_model_runs.append(r)
        print(f"  seed {s}: two_model per-medium mean {r['per_medium_mean']:+.4f} "
              f"pooled {r['pooled']:+.4f} (held {held}, {len(idx_te)} test)")

    baseline = {
        "arm": "two_model",
        "per_medium_mean_median": round(statistics.median(r["per_medium_mean"] for r in two_model_runs), 4),
        "pooled_median": round(statistics.median(r["pooled"] for r in two_model_runs), 4),
        "per_seed": two_model_runs,
    }
    artifact = {
        "schema": "glm-condition-conditioning-v1",
        "date": str(date.today()),
        "stage": "step4-baselines-only",
        "substrate": "GEO GSE144621 LB + M9 sheared fragments, coordinate-keyed intersection",
        "n_shared": len(lb),
        "load_stats": {"LB": st_lb.as_dict(), "M9": st_m9.as_dict()},
        "per_medium_ceilings": ceilings,
        "condition_effect": cond,
        "preregistered": PREREGISTERED,
        "seeds": seeds,
        "two_model_baseline": baseline,
    }
    out = Path(args.out) if args.out else Path(f"wiki/glm_condition_conditioning_{date.today()}.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(artifact, indent=2), encoding="utf-8")
    print(f"\ntwo_model baseline: per-medium mean {baseline['per_medium_mean_median']:+.4f} | "
          f"pooled {baseline['pooled_median']:+.4f}")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
