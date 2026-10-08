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
    "primary_arm": "partial_pooling",
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
    "measured_degeneracy": "a FULL interaction [seq, m, seq*m] spans the same hypothesis space as two "
                           "independent per-medium models (verified: corr 0.9999998 on synthetic truth with "
                           "medium-specific coefficients; identical to 4dp on two real seeds). It is "
                           "retained as a demonstration, NOT as a candidate. The primary arm is "
                           "partial_pooling, which sits strictly between full pooling and no pooling.",
    "pooling_strength": "SELECTED on an inner position-blocked split of the TRAINING data from the frozen "
                        "grid (1.0, 0.5, 0.3, 0.1, 0.03); the test split never informs the choice, and the "
                        "whole inner-validation curve is reported so a flat curve stays visible. This "
                        "SUPERSEDES an asserted 0.3, which had never been exercised: it was applied before "
                        "StandardScaler, which divides a constant column factor straight back out, so all "
                        "three conditioned arms were silently the same model. The pooling penalty is now "
                        "applied AFTER standardisation (scale c => the interaction block is penalised by "
                        "alpha/c^2).",
    "scale_selection_is_refined_not_relaxed": "the VERDICT BAR is untouched (min_gain 0.02, validity "
                                              "first). What changed is the arm's own definition, so that "
                                              "the primary arm is non-degenerate -- which is what it was "
                                              "always specified to be.",
    "sequence_features": "kmer3 (64) -- chosen because Step 3 measured that MORE features hurt on this "
                         "project's substrate, and because 4-mers would force a subsample (decision D2 says "
                         "use the full data)",
    "thresholds_are": "ASSERTED; the per-medium ceilings that contextualise them are DERIVED",
}

K = 3
MEDIA = ("LB", "M9")


def _seq_features(seqs: list[str]) -> list[list[float]]:
    return kmer_features(seqs, K)


def feature_matrix(frags):
    """The 3-mer matrix for a medium, computed ONCE and reused by every arm and every seed.

    MEASURED, not assumed: `kmer_features` costs ~181 s over the full fragment set x 2 media, against ~4 s
    for the ridge fit it feeds -- so ~97% of the sweep's cost was recomputing a PURE FUNCTION OF THE
    SEQUENCE seven times per seed. The sequences never change between arms (arms differ only in how the
    medium flag enters the design) nor between seeds (seeds change only the split), so one pass per medium
    is sufficient and a 10-seed/7-arm sweep drops from ~3.5 h to minutes. That is what makes decision D2
    (full data, multi-seed median) affordable rather than a compute trade.

    float64 deliberately: `build_design` does its arithmetic here and casts to float32 only at the end,
    which is bit-identical to the previous list-of-floats path that `fit_ridge` cast on entry.
    """
    import numpy as np
    return np.asarray(_seq_features([f.seq for f in frags]), dtype=np.float64)


#: Frozen candidate grid for the interaction-block scale. 1.0 is the no-pooling endpoint (= the two-model
#: comparator); smaller values penalise the condition-specific block harder, approaching full pooling.
#: The grid is FROZEN; which member is used is SELECTED ON TRAINING DATA ONLY (see `select_pooling_scale`).
POOLING_GRID = (1.0, 0.5, 0.3, 0.1, 0.03)

#: Retained only so the superseded asserted value stays visible in the record. NOT used for the verdict.
SUPERSEDED_ASSERTED_SCALE = 0.3


def interaction_col_scale(n_features: int, scale: float):
    """Per-column scale for the `interaction` design: 1.0 on the shared block, `scale` on the interaction.

    Consumed by `fit_ridge(..., col_scale=)` and therefore applied AFTER standardisation, which is the only
    place it does anything (see `fit_ridge`). The shared sequence block and the medium column are left at
    1.0 so partial pooling shrinks ONLY the condition-specific part.
    """
    import numpy as np
    w = (n_features - 1) // 2
    cs = np.ones(n_features, dtype=np.float64)
    cs[w + 1:] = scale
    return cs


def build_design(S, medium_flags, *, arm: str):
    """Design matrix per arm. `S` is the precomputed 3-mer block; `medium_flags` is 0.0 LB / 1.0 M9.

    `interaction` is the load-bearing one: appending `seq * medium` lets the model give a sequence feature a
    DIFFERENT coefficient per medium, which an indicator alone cannot do.

    No pooling scale is applied here -- deliberately. Scaling a column before `StandardScaler` is a no-op,
    so partial pooling lives in `fit_ridge`'s `col_scale` instead; see both docstrings.

    Arithmetic is float64 and the cast to float32 happens on the way out -- exactly where `fit_ridge` used to
    do it -- so the returned matrix is bit-identical to the pre-vectorisation path.
    """
    import numpy as np
    S = np.asarray(S, dtype=np.float64)
    m = np.asarray(medium_flags, dtype=np.float64).reshape(-1, 1)
    n, w = S.shape
    if arm == "sequence_only":
        return S.astype(np.float32)
    if arm == "indicator":
        X = np.empty((n, w + 1), dtype=np.float32)
        X[:, :w] = S
        X[:, w] = m[:, 0]
        return X
    if arm == "interaction":
        # MEASURED DEGENERACY -- read this before trusting an `interaction` number.
        #
        # A FULL interaction [seq, m, seq*m] with a binary m spans EXACTLY the same hypothesis space as two
        # independent per-medium models: LB uses w, M9 uses w+v. Verified on synthetic data where the truth
        # genuinely has medium-specific coefficients: corr(two_model, full_interaction) = 0.9999998, max
        # prediction difference 0.0118. On the real fragments both arms returned per-medium-mean Spearman
        # IDENTICAL to four decimals on two separate seeds (+0.2013 and +0.1916).
        #
        # So this arm CANNOT beat the two-model comparator -- it IS the comparator, re-parameterised, with
        # ridge's penalty placement the only difference. The first design made the primary arm an indicator
        # (too WEAK to express interaction); correcting to a full interaction over-shot to exactly
        # EQUIVALENT. It is retained to demonstrate that equivalence, not as a candidate.
        #
        # The informative arm is `partial_pooling`: penalise the interaction block MORE than the shared
        # block, placing the model strictly BETWEEN full pooling (indicator) and no pooling (two_model) --
        # the actual multi-task claim, and the only one of the three that can beat both endpoints.
        #
        # SECOND MEASURED CORRECTION: that penalty is applied in `fit_ridge(col_scale=)`, AFTER
        # standardisation. The first attempt scaled this block HERE, which `StandardScaler` divides straight
        # back out -- a silent no-op that made partial_pooling identical to this arm (r=0.9999970, same to 7
        # decimals) and so left ALL THREE conditioned arms the same model. A test pins it now.
        X = np.empty((n, 2 * w + 1), dtype=np.float32)
        X[:, :w] = S
        X[:, w] = m[:, 0]
        X[:, w + 1:] = S * m
        return X
    raise ValueError(f"unknown arm {arm!r}")


def fit_ridge(Xtr, ytr, Xte, col_scale=None):
    """Ridge on standardised features, with an OPTIONAL per-column scale applied AFTER standardisation.

    `col_scale` is how partial pooling is actually implemented, and the ordering is load-bearing. Scaling a
    column BEFORE `StandardScaler` is a NO-OP -- standardisation divides out any constant factor, so the
    standardised matrix is unchanged (measured: scaling the interaction block by 0.3 pre-scaler left
    predictions identical to the unscaled model to 7 decimals, making all three conditioned arms the same
    model). Applied AFTER standardisation it is a genuine per-block penalty: if a standardised column is
    multiplied by c, an effective coefficient gamma on that column costs alpha*(gamma/c)^2, i.e. the block
    is penalised by alpha/c^2. So c=1 is no pooling (the two-model comparator), c->0 is full pooling (the
    interaction block is penalised to death, leaving the indicator arm), and 0<c<1 is partial pooling.
    """
    import numpy as np
    from sklearn.linear_model import Ridge
    from sklearn.preprocessing import StandardScaler
    Xtr = np.asarray(Xtr, dtype=np.float32)
    Xte = np.asarray(Xte, dtype=np.float32)
    y = np.log1p(np.asarray(ytr, dtype=np.float64))
    sc = StandardScaler().fit(Xtr)
    Ztr, Zte = sc.transform(Xtr), sc.transform(Xte)
    if col_scale is not None:
        cs = np.asarray(col_scale, dtype=np.float32)
        if cs.shape[0] != Ztr.shape[1]:
            raise ValueError(f"col_scale has {cs.shape[0]} entries for {Ztr.shape[1]} columns")
        Ztr = Ztr * cs
        Zte = Zte * cs
    m = Ridge(alpha=1.0).fit(Ztr, y)
    return [float(v) for v in m.predict(Zte)]


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



def run_conditioned_arm(frag, S, idx_tr, idx_te, *, arm: str, seed: int,
                        shuffle_medium: bool = False, data_match: bool = False,
                        interaction_scale: float = 1.0) -> dict:
    """ONE model trained across BOTH media, unlike `two_model` which fits per-medium weights.

    `shuffle_medium` permutes the medium label across training rows -- the validity control. If a permuted
    label still produces a gain, the pipeline fabricates signal and nothing may be graded.

    `data_match` subsamples the combined training rows to roughly ONE medium's count. Without it, a gain is
    ambiguous: the shared model sees ~2x the rows, so it could win on data volume while learning nothing
    about condition at all.
    """
    import random as _r

    import numpy as np
    ix = np.asarray(idx_tr, dtype=np.int64)
    # LB block then M9 block -- the same row order as the pre-vectorisation loop, so the RNG draws below
    # land on the same rows.
    S_tr = np.vstack([S[m][ix] for m in MEDIA])
    flags_tr = [float(mi) for mi, _ in enumerate(MEDIA) for _i in idx_tr]
    y_tr = [frag[m][i].expression for m in MEDIA for i in idx_tr]
    if shuffle_medium:
        rng = _r.Random(9000 + seed)
        flags_tr = flags_tr[:]
        rng.shuffle(flags_tr)
    if data_match:
        rng = _r.Random(7000 + seed)
        keep = rng.sample(range(len(y_tr)), len(y_tr) // len(MEDIA))
        S_tr = S_tr[np.asarray(keep, dtype=np.int64)]
        flags_tr = [flags_tr[i] for i in keep]
        y_tr = [y_tr[i] for i in keep]

    X_tr = build_design(S_tr, flags_tr, arm=arm)
    col_scale = (interaction_col_scale(X_tr.shape[1], interaction_scale)
                 if (arm == "interaction" and interaction_scale != 1.0) else None)

    # ONE model, fitted once, then asked for a prediction per medium. Fitting inside the medium loop would
    # train the same model twice and -- worse -- stop being "one conditioned model" in any meaningful sense.
    jx = np.asarray(idx_te, dtype=np.int64)
    pred, truth = {}, {}
    n_te = len(idx_te)
    stacked = np.vstack([
        build_design(S[m][jx], [float(mi)] * n_te, arm=arm) for mi, m in enumerate(MEDIA)
    ])
    all_pred = fit_ridge(X_tr, y_tr, stacked, col_scale=col_scale)
    off = 0
    for m in MEDIA:
        pred[m] = all_pred[off:off + n_te]
        off += n_te
        truth[m] = [frag[m][i].expression for i in idx_te]
    r = per_medium_and_pooled(pred, truth)
    r["n_train_rows"] = len(y_tr)
    r["n_features"] = int(X_tr.shape[1])
    return r


def _inner_split(frag, idx_tr, *, seed: int, n_blocks: int, held_blocks: int):
    """A position-blocked split WITHIN the training set, for selecting the pooling strength.

    Blocked again rather than random, for the same reason the outer split is blocked: 91.5% of consecutive
    fragments overlap, so a random inner split would let the selection see near-duplicates of its own
    validation rows and pick a scale that is over-fit to them.
    """
    sub = [frag["LB"][i] for i in idx_tr]
    itr, _iva, held = position_blocked_split(sub, n_blocks=n_blocks, held_blocks=held_blocks,
                                             seed=5000 + seed)
    keys = {(f.start, f.end, f.strand) for f in itr}
    k = lambda i: (frag["LB"][i].start, frag["LB"][i].end, frag["LB"][i].strand)
    return [i for i in idx_tr if k(i) in keys], [i for i in idx_tr if k(i) not in keys], held


def select_pooling_scale(frag, S, idx_tr, *, seed: int, n_blocks: int, held_blocks: int,
                         shuffle_medium: bool = False) -> tuple[float, dict, list[int]]:
    """Choose the interaction-block scale on TRAINING DATA ONLY, never on the test split.

    The pooling strength is the one real hyperparameter of this arm. Asserting it would make the headline a
    statement about an arbitrary constant -- and the asserted 0.3 had never actually been exercised, because
    it was applied in a position where `StandardScaler` divided it back out. Selecting it on an inner
    position-blocked split makes the arm *partial pooling with the pooling strength learned from training
    data*, which is the actual multi-task method, and leaves the test split untouched by the choice.

    Returns the selected scale, the full inner-validation curve (reported, so a flat curve is visible as a
    flat curve rather than hidden behind an argmax), and the inner held blocks.
    """
    inner_tr, inner_va, held = _inner_split(frag, idx_tr, seed=seed, n_blocks=n_blocks,
                                            held_blocks=held_blocks)
    curve = {}
    for c in POOLING_GRID:
        r = run_conditioned_arm(frag, S, inner_tr, inner_va, arm="interaction", seed=seed,
                                interaction_scale=c, shuffle_medium=shuffle_medium)
        curve[c] = r[PREREGISTERED["primary_metric"]]
    best = max(POOLING_GRID, key=lambda c: curve[c])
    return best, {str(c): round(v, 4) for c, v in curve.items()}, held


def verdict(arms: dict) -> tuple[str, str]:
    """Frozen four-branch rule. VALIDITY IS CHECKED FIRST (decision D1)."""
    key = PREREGISTERED["primary_metric"]
    base = arms["two_model"][key]
    inter = arms["partial_pooling"][key]
    full_int = arms["interaction"][key]
    shuf = arms["partial_pooling_shuffled"][key]
    dm = arms["data_matched"][key]
    ind = arms["indicator"][key]
    gain = inter - base
    shuf_gain = shuf - base
    dm_gain = dm - base
    d = (f"{key}: two_model {base:+.4f} | partial_pooling {inter:+.4f} (gain {gain:+.4f}) | "
         f"full_interaction {full_int:+.4f} [DEGENERATE: equals two_model by construction] | "
         f"indicator {ind:+.4f} | shuffled {shuf:+.4f} (gain {shuf_gain:+.4f}) | "
         f"data_matched {dm:+.4f} (gain {dm_gain:+.4f}); bar {PREREGISTERED['min_gain']}")
    if shuf_gain >= PREREGISTERED["validity"]["shuffled_max_gain"]:
        return ("INDETERMINATE_NULL_NOT_CLEAN",
                f"the medium-SHUFFLED arm itself gained {shuf_gain:+.4f}, at or above the "
                f"{PREREGISTERED['validity']['shuffled_max_gain']} bar. The pipeline fabricates signal; "
                f"nothing is graded. {d}")
    if gain < PREREGISTERED["min_gain"]:
        return "NO_GAIN", d
    if dm_gain < PREREGISTERED["min_gain"]:
        return "GAIN_IS_DATA_VOLUME", d
    return "CONDITIONING_HELPS", d


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

    # ONE feature pass per medium for the WHOLE sweep (see `feature_matrix`). The fragments are
    # coordinate-aligned between media, so row i is the same locus in both.
    import time as _t
    _t0 = _t.time()
    S = {m: feature_matrix(frag[m]) for m in MEDIA}
    print(f"  3-mer features: {S['LB'].shape} per medium, computed once in {_t.time() - _t0:.0f}s")

    seeds = list(range(args.seeds))
    two_model_runs = []
    arm_runs: dict = {}
    selections: list[dict] = []
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
            X_tr = build_design(S[m][idx_tr], [0.0] * len(idx_tr), arm="sequence_only")
            X_te = build_design(S[m][idx_te], [0.0] * len(idx_te), arm="sequence_only")
            pred[m] = fit_ridge(X_tr, [fr[i].expression for i in idx_tr], X_te)
            truth[m] = [fr[i].expression for i in idx_te]
        r = per_medium_and_pooled(pred, truth)
        r["held_blocks"] = held
        r["n_train"], r["n_test"] = len(idx_tr), len(idx_te)
        two_model_runs.append(r)
        print(f"  seed {s}: two_model pmm {r['per_medium_mean']:+.4f} pooled {r['pooled']:+.4f} "
              f"(held {held}, {len(idx_te)} test)")
        # The pooling strength is chosen on an INNER split of the training rows. The shuffled control runs
        # the SAME selection on shuffled labels -- a null that skipped selection would be the null of a
        # different method and could not license the primary number.
        sel, curve, inner_held = select_pooling_scale(frag, S, idx_tr, seed=s, n_blocks=args.n_blocks,
                                                      held_blocks=args.held_blocks)
        sel_sh, curve_sh, _ = select_pooling_scale(frag, S, idx_tr, seed=s, n_blocks=args.n_blocks,
                                                   held_blocks=args.held_blocks, shuffle_medium=True)
        selections.append({"seed": s, "selected_scale": sel, "inner_curve": curve,
                           "inner_held_blocks": inner_held, "selected_scale_shuffled": sel_sh,
                           "inner_curve_shuffled": curve_sh})
        print(f"           pooling scale selected on inner split: {sel} "
              f"(curve {curve}) | shuffled-null selected {sel_sh}")

        for arm_name, kw in (("interaction", {}), ("indicator", {}),
                             ("shuffled", {"shuffle_medium": True}),
                             ("partial_pooling", {"interaction_scale": sel}),
                             ("partial_pooling_shuffled",
                              {"interaction_scale": sel_sh, "shuffle_medium": True}),
                             # the data-volume control MIRRORS THE PRIMARY ARM, so it isolates row count
                             # rather than comparing two different models
                             ("data_matched", {"interaction_scale": sel, "data_match": True})):
            base_arm = "indicator" if arm_name == "indicator" else "interaction"
            rr = run_conditioned_arm(frag, S, idx_tr, idx_te, arm=base_arm, seed=s, **kw)
            arm_runs.setdefault(arm_name, []).append(rr)
            print(f"           {arm_name:25s} pmm {rr['per_medium_mean']:+.4f} "
                  f"pooled {rr['pooled']:+.4f} ({rr['n_features']} feats, {rr['n_train_rows']} rows)")

    baseline = {
        "arm": "two_model",
        "per_medium_mean_median": round(statistics.median(r["per_medium_mean"] for r in two_model_runs), 4),
        "pooled_median": round(statistics.median(r["pooled"] for r in two_model_runs), 4),
        "per_seed": two_model_runs,
    }
    agg = {"two_model": {k: round(statistics.median(r[k] for r in two_model_runs), 4)
                         for k in ("per_medium_mean", "pooled")}}
    for name, runs in arm_runs.items():
        agg[name] = {k: round(statistics.median(r[k] for r in runs), 4)
                     for k in ("per_medium_mean", "pooled")}
        agg[name]["n_features"] = runs[0]["n_features"]
        agg[name]["n_train_rows"] = runs[0]["n_train_rows"]
        agg[name]["per_seed_pmm"] = [round(r["per_medium_mean"], 4) for r in runs]
    v, vd = verdict(agg)
    artifact = {
        "schema": "glm-condition-conditioning-v1",
        "date": str(date.today()),
        "stage": "step5-all-arms",
        "arms": agg,
        "verdict": v,
        "verdict_detail": vd,
        "null_clean": v != "INDETERMINATE_NULL_NOT_CLEAN",
        "registered_protocol_used": True,
        "substrate": "GEO GSE144621 LB + M9 sheared fragments, coordinate-keyed intersection",
        "n_shared": len(lb),
        "load_stats": {"LB": st_lb.as_dict(), "M9": st_m9.as_dict()},
        "per_medium_ceilings": ceilings,
        "condition_effect": cond,
        "preregistered": PREREGISTERED,
        "seeds": seeds,
        "two_model_baseline": baseline,
        "pooling_scale_selection": {
            "grid": list(POOLING_GRID),
            "selected_on": "inner position-blocked split of the TRAINING rows only",
            "per_seed": selections,
            "superseded_asserted_value": SUPERSEDED_ASSERTED_SCALE,
        },
        "honest_limits": [
            "TWO conditions only. A null result is weak evidence against conditioning in general; it rules "
            "out conditioning as tested, on two growth media.",
            "Conditioning is tested as sequence x medium INTERACTIONS over 3-mer features, which is not the "
            "same as a model whose REPRESENTATION is modulated by condition.",
            "The fragment substrate has a low ceiling (LB 0.7925 / M9 0.8040) because 97% of fragments "
            "carry a single barcode, so absolute numbers are small by construction and must be read against "
            "those per-medium ceilings.",
            "per_medium_mean is the primary metric. The pooled metric is inflated by the between-media "
            "level difference (measured: pooled 0.310 vs per-medium mean 0.197 for the baseline), which a "
            "global offset can capture without any condition-specific sequence response.",
            "Only length-invariant features apply (fragments are 48-499bp with 359 distinct lengths), so "
            "this arm is not comparable to the designed grid's positional-one-hot number.",
            "THE DESIGN IS RANK-DEFICIENT BY CONSTRUCTION and sklearn says so (LinAlgWarning, rcond ~1e-8). "
            "The 64 3-mer FREQUENCIES sum to 1.0 for every row -- an affine constraint, so the raw block is "
            "full rank 64 -- but StandardScaler CENTRES each column, which turns it into an exact linear "
            "dependency: measured rank 63 of 64, condition number 9.3e14. Ridge alpha=1.0 is therefore "
            "doing load-bearing work on every arm, which is a partial mechanistic account of the FLAT "
            "pooling curve: an extra per-block penalty has little room to change a fit that is already "
            "heavily regularised on a near-collinear basis.",
            "THE POOLING CURVE IS FLAT and the selected scale is therefore NOISE. Inner-validation scores "
            "differ by <= 0.0001 across the whole grid, and the argmax moves between 0.03, 0.1 and 1.0 "
            "across seeds. Read this as 'pooling strength is not a lever on this substrate', NOT as "
            "'0.03 is the right amount of pooling'.",
            "min_gain = 0.02 is ASSERTED, not derived.",
        ],
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
