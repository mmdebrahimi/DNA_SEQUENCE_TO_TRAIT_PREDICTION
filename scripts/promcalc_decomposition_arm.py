"""Stage C: does the FULL mechanistic decomposition beat the single scalar reduced from it?

The 2026-10-12 measurement reduced ~73 candidate promoters x 10 numbers per tile to ONE scalar
(`max_tx`) and scored median Spearman 0.3968 -- beating every learned arm and the only arm on this
substrate whose paired CI against GC excludes zero. This asks whether the ~730 discarded numbers carry
anything the single strongest-promoter rate does not.

**What makes this the right next question rather than more capacity.** The learned arms here were fed
RAW SEQUENCE and lost: 13,217 params scored 0.3234 and 27,649 params scored 0.3167, both with paired
CIs against a one-feature GC baseline that include zero. The literature's first-ranked lever is
training-set size and it is EXHAUSTED on this rung (46,713 rows is the entire file). So the live lever
is the REPRESENTATION, not the capacity or the row count -- and the 2026-10-12 feature-permutation
control established that the signal lives in the sequence-to-mechanism pairing specifically. This arm
gives a learner the mechanism instead of asking it to rediscover sigma70 physics from 44k rows.

## Pre-registered, and the threshold is DERIVED

Frozen before the fit. The bar is read off the substrate's own measured effect sizes, because a
pre-registered number that was never analysed is just an assertion with a timestamp:

* Gains that turned out to be NOISE here: +0.0069 (pyramid) and +0.0125 (encoder-multi) over GC, both
  CI-including-zero over 10 seeds.
* The one gain that was REAL: +0.0686 (biophysical over GC), CI excluding zero.

So `+0.02` sits above the band that swallowed both learned arms and well below the single real effect.
It is the smallest gain this substrate has shown itself able to resolve, not a round number.

| verdict | condition |
|---|---|
| `DECOMPOSITION_ADDS` | paired median gain over `max_tx` >= +0.02 AND >= 8/10 seeds win AND control clean |
| `WEAK` | gain > 0 but below either the margin or the seed count |
| `NO_GAIN` | paired median gain <= 0 |
| `INDETERMINATE_RECONCILE_FAILED` | `gc` seed 0 is not 0.3000 +- 0.01, or `max_tx` is not 0.3968 +- 0.02 |
| `INDETERMINATE_CONTROL_NOT_CLEAN` | the feature-permutation control itself scores |

**The control is FEATURE permutation, never shuffled labels.** Shuffling training labels randomises a
fitted coefficient's sign but not its ranking power: the prediction stays a monotone function of a
feature that genuinely correlates with the target, so the control returns +-0.4 by construction and can
never collapse. That error was made and caught on this exact substrate on 2026-10-12. Permuting the
tile->feature assignment destroys the sequence-to-mechanism pairing, which is the thing under test.

**Both baselines are reproduced INSIDE this pipeline, not quoted.** `gc` is the gate's own `ridge_gc`
imported, and `max_tx` is refit here through the identical StandardScaler -> Ridge path. Comparing a
rich-feature fit against a number carried over from another script would confound the feature set with
the pipeline; the two reconcile gates are what make the comparison mean anything.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import statistics
import sys
from collections import defaultdict
from datetime import date
from pathlib import Path

DEFAULT_DIR = Path("D:/dna_decode_cache/promcalc")

#: Frozen BEFORE the fit. See the module docstring for how +0.02 is derived.
PREREGISTERED = {
    "question": "does the full per-site mechanistic decomposition beat the single max_tx scalar "
                "reduced from it, on the same tiles under the same registered protocol?",
    "primary_metric": "spearman",
    "comparator": "max_tx",
    "secondary_comparator": "gc",
    "margin": 0.02,
    "min_seeds_winning": 8,
    "registered_seeds": list(range(10)),
    "control": "feature_permutation",
    "control_max_abs": 0.05,
    "reconcile": {
        "gc_seed0": {"expected": 0.3000, "tolerance": 0.01,
                     "source": "glm_encoder_gate.PREREGISTERED.reconcile_anchor"},
        "max_tx_median": {"expected": 0.3968, "tolerance": 0.02,
                          "source": "wiki/promcalc_vs_learned_2026-10-12.json"},
    },
    "margin_is": "DERIVED from this substrate's measured effect sizes: above the +0.0069/+0.0125 "
                 "CI-including-zero band of the learned arms, below the +0.0686 of the one real gain",
}

#: The nine component energies, in the order stage B writes them.
DG_FIELDS = ("dG_total", "dG_10", "dG_35", "dG_ext10", "dG_spacer", "dG_disc", "dG_ITR",
             "dG_UP", "dG_bind")


# --------------------------------------------------------------------------------------------------
# features (pure)
# --------------------------------------------------------------------------------------------------
def tile_features(sites: list[dict]) -> list[float]:
    """Fixed-length mechanism features for ONE tile from its per-site decomposition.

    Three families, each answering a different question the single scalar cannot:
      1. the WINNING promoter's own energy breakdown -- how its strength is distributed across
         -10 / -35 / spacer / discriminator / UP, which `max_tx` collapses to a single rate;
      2. distributional summaries over all candidates -- whether many medium promoters differ from
         one strong one (`sum_tx` tested this crudely and did not separate from `max`);
      3. the winner's POSITION and spacer length -- geometry the rate does not carry.

    Order is fixed and documented by `feature_names()`; a permutation control depends on the vector
    being positionally meaningful.
    """
    if not sites:
        # A tile with no candidate promoter is real (not an error). Zeros plus an explicit flag, so
        # the model can separate "no promoter" from "a promoter with zero energy".
        return [0.0] * (len(DG_FIELDS) + 4) + [0.0] * (len(DG_FIELDS) * 4 + 4) + [0.0, 1.0]

    tx = [s["Tx_rate"] for s in sites]
    best = sites[max(range(len(sites)), key=lambda i: tx[i])]

    f: list[float] = []
    # (1) the winning site's decomposition + geometry
    f += [best[k] for k in DG_FIELDS]
    f += [math.log1p(max(best["Tx_rate"], 0.0)), float(best["spacer_len"]),
          float(best["hex35_start"]), float(best["hex10_start"])]

    # (2) distribution over every candidate
    for k in DG_FIELDS:
        v = [s[k] for s in sites]
        f += [min(v), max(v), statistics.fmean(v),
              statistics.pstdev(v) if len(v) > 1 else 0.0]
    txp = [max(t, 0.0) for t in tx]
    srt = sorted(txp, reverse=True)
    f += [math.log1p(sum(txp)), math.log1p(srt[0]),
          math.log1p(srt[1]) if len(srt) > 1 else 0.0,
          math.log1p(srt[2]) if len(srt) > 2 else 0.0]

    # (3) how many candidates there were at all
    f += [float(len(sites)), 0.0]
    return f


def feature_names() -> list[str]:
    n = [f"best_{k}" for k in DG_FIELDS]
    n += ["best_log1p_tx", "best_spacer_len", "best_hex35_start", "best_hex10_start"]
    for k in DG_FIELDS:
        n += [f"{k}_min", f"{k}_max", f"{k}_mean", f"{k}_std"]
    n += ["log1p_sum_tx", "log1p_tx1", "log1p_tx2", "log1p_tx3"]
    n += ["n_sites", "no_sites_flag"]
    return n


def load_decomposition(in_dir: Path | str, *, smoke: bool = False) -> tuple[dict[int, list[dict]], list[str]]:
    """tile_id -> per-site dicts, read from the main file AND every resume shard.

    **Reads a GLOB, not one filename, and that is load-bearing.** Stage B writes
    `decomposition.part-NN.csv.gz` when resuming an interrupted sweep, so a loader hardcoded to
    `decomposition.csv.gz` would silently score whichever fraction the first run happened to finish
    (88.4% on the real interruption) and report a perfectly clean verdict over a partial substrate.
    That is the project's recorded failure shape: a consumer that parsed only part of its input looks
    exactly like one whose input was complete.

    **Truncation is tolerated; MISSING COVERAGE is not — they are different questions.** An earlier
    version of this function let the `EOFError` propagate, reasoning that a short read means scoring
    an unknown fraction. That was wrong in a way only running it showed: a killed sweep leaves
    `decomposition.csv.gz` PERMANENTLY truncated, so after a successful resume brought coverage to
    100% the loader would still refuse forever. The tail of a gzip stream says nothing about whether
    the substrate is complete. So a truncated member is read as far as it goes and counted, and the
    caller's **coverage check against the tile list** is the real gate — it answers the question that
    actually matters (are any tiles missing?) rather than a proxy for it (did every file close?).
    """
    in_dir = Path(in_dir)
    if smoke:
        files = sorted(in_dir.glob("decomposition_smoke*.csv.gz"))
    else:
        # EXPLICIT, not `decomposition*` - that glob also matches decomposition_smoke.csv.gz and
        # merged the 20-tile smoke run into the real substrate.
        files = [p for p in [in_dir / "decomposition.csv.gz",
                             *sorted(in_dir.glob("decomposition.part-*.csv.gz"))] if p.exists()]
    if not files:
        raise SystemExit(f"no decomposition shards in {in_dir} - run stage B first")
    need = (*DG_FIELDS, "Tx_rate", "spacer_len", "hex35_start", "hex10_start", "tile_id")
    out: dict[int, list[dict]] = defaultdict(list)
    partial = 0
    truncated: list[str] = []
    for p in files:
        try:
            with gzip.open(p, "rt", encoding="utf-8") as fh:
                for r in csv.DictReader(fh):
                    # A killed run ends MID-RECORD, and csv.DictReader pads the stub with None rather
                    # than raising -- so an unguarded int() crashes on the tail of every interrupted
                    # sweep. Skip a stub, but COUNT it: gzip can only lose the TAIL, so more than one
                    # per file means something other than truncation is wrong and must not be absorbed.
                    if any(r.get(k) is None for k in need):
                        partial += 1
                        continue
                    out[int(r["tile_id"])].append({
                        **{k: float(r[k]) for k in DG_FIELDS},
                        "Tx_rate": float(r["Tx_rate"]),
                        "spacer_len": int(r["spacer_len"]),
                        "hex35_start": int(r["hex35_start"]),
                        "hex10_start": int(r["hex10_start"]),
                    })
        except (EOFError, OSError):
            truncated.append(p.name)  # interrupted writer; the coverage gate decides, not this
    if truncated:
        print(f"  note: truncated tail on {truncated} (interrupted run) - coverage gate decides",
              flush=True)
    # Malformed rows are COUNTED and skipped, never used to infer completeness.
    # An earlier version refused above a threshold of one-per-file, on the premise that "truncation
    # can only damage one row". MEASURED FALSE: a cut deflate stream replays its sliding window, and
    # the real interruption produced 2,199 copies of a single mangled row (a float written as
    # `999999999993236`) all belonging to the one tile being written. The threshold was arbitrary and
    # the premise was wrong, so completeness is decided ONLY by the caller's coverage check against
    # the tile list - which is the question that matters and does not depend on this count.
    if partial:
        print(f"  note: skipped {partial} malformed row(s) (interrupted-writer garbage)", flush=True)
    return dict(out), [p.name for p in files]


# --------------------------------------------------------------------------------------------------
# fitting
# --------------------------------------------------------------------------------------------------
def _fit_ridge(Xtr, ytr, Xte):
    import numpy as np
    from sklearn.linear_model import Ridge
    from sklearn.preprocessing import StandardScaler
    sc = StandardScaler().fit(np.asarray(Xtr, dtype=np.float64))
    m = Ridge(alpha=1.0).fit(sc.transform(np.asarray(Xtr, dtype=np.float64)),
                             np.asarray(ytr, dtype=np.float64))
    return [float(v) for v in m.predict(sc.transform(np.asarray(Xte, dtype=np.float64)))]


def _agg(vals: list[float]) -> dict:
    return {"median": round(statistics.median(vals), 4),
            "min": round(min(vals), 4), "max": round(max(vals), 4),
            "per_seed": [round(v, 4) for v in vals]}


def paired(a: list[float], b: list[float]) -> dict:
    """Per-seed paired margin a-b with a normal-approx CI95; the CI is the point, not the mean."""
    d = [x - y for x, y in zip(a, b)]
    mean = statistics.fmean(d)
    if len(d) > 1:
        se = statistics.stdev(d) / math.sqrt(len(d))
        lo, hi = mean - 1.96 * se, mean + 1.96 * se
    else:
        lo = hi = mean
    return {"mean": round(mean, 4), "median": round(statistics.median(d), 4),
            "ci95": [round(lo, 4), round(hi, 4)], "includes_zero": bool(lo <= 0.0 <= hi),
            "wins": sum(1 for v in d if v > 0), "n": len(d)}


def verdict(arms: dict, pre: dict = PREREGISTERED) -> tuple[str, str]:
    rec = pre["reconcile"]
    gc0 = arms["gc"]["per_seed"][0]
    if abs(gc0 - rec["gc_seed0"]["expected"]) > rec["gc_seed0"]["tolerance"]:
        return ("INDETERMINATE_RECONCILE_FAILED",
                f"gc seed 0 = {gc0:.4f}, not {rec['gc_seed0']['expected']} "
                f"+-{rec['gc_seed0']['tolerance']}; this script is not reading the same data as the "
                f"published numbers and nothing here is interpretable")
    mt = arms["max_tx"]["median"]
    if abs(mt - rec["max_tx_median"]["expected"]) > rec["max_tx_median"]["tolerance"]:
        return ("INDETERMINATE_RECONCILE_FAILED",
                f"max_tx median = {mt:.4f}, not {rec['max_tx_median']['expected']} "
                f"+-{rec['max_tx_median']['tolerance']}; the scalar did not reproduce, so a "
                f"decomposition-vs-scalar comparison would confound features with pipeline")
    ctrl = arms.get("control_feature_permuted", {}).get("median")
    if ctrl is not None and abs(ctrl) >= pre["control_max_abs"]:
        return ("INDETERMINATE_CONTROL_NOT_CLEAN",
                f"the feature-permutation control scores {ctrl:+.4f} (bar {pre['control_max_abs']}); "
                f"the pipeline fabricates signal and nothing is graded")
    p = arms["_paired_vs_max_tx"]
    if p["median"] <= 0:
        return ("NO_GAIN", f"decomposition {p['median']:+.4f} vs the max_tx scalar at the median "
                           f"({p['wins']}/{p['n']} seeds) -- the ~730 discarded numbers add nothing "
                           f"the strongest-promoter rate does not already carry")
    if p["median"] >= pre["margin"] and p["wins"] >= pre["min_seeds_winning"]:
        return ("DECOMPOSITION_ADDS",
                f"decomposition beats the scalar by {p['median']:+.4f} at the median, "
                f"{p['wins']}/{p['n']} seeds, CI95 {p['ci95']}")
    return ("WEAK", f"decomposition {p['median']:+.4f} over the scalar, {p['wins']}/{p['n']} seeds -- "
                    f"positive but under the {pre['margin']} margin or the "
                    f"{pre['min_seeds_winning']}-seed bar; not enough to carry the representation")


def _require_decomposition_dir(in_dir: Path | str, *, flag: str = "--in-dir") -> Path:
    """`in_dir`, or a refusal naming the flag and where the data actually is.

    Mirrors the shape `scripts/glm_tile_predictions.py::require_external_dir` already ships, and for
    the same reason: `allow_discovered` stays at its default False, because adopting a same-named root
    automatically would convert a loud missing-ROOT into a quiet missing-DATASET -- a directory of the
    right NAME is not evidence of the right CONTENTS. This diagnoses and refuses; the operator passes
    the path.

    Unlike stage A's output dir, this one must ALREADY exist and hold shards, so the resolver's
    `configured`-means-is_dir semantics are exactly right here.

    Without this, an absent volume reached `load_decomposition` and surfaced as
    "no decomposition shards in <path> - run stage B first", which points the operator at stage B when
    the real wall is the missing volume. The D: default is unchanged.
    """
    from dna_decode.data.external_root import describe, resolve_external_root

    path, status = resolve_external_root(in_dir)
    if status == "configured":
        return path
    raise SystemExit(
        f"{flag} {in_dir} is not reachable [{status}].\n\n{describe(in_dir)}\n\n"
        f"Pass {flag} explicitly. The root is NOT substituted automatically. This is a different wall "
        f"from 'stage B has not run yet' -- the directory itself is absent."
    )


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--in-dir", default=str(DEFAULT_DIR))
    ap.add_argument("--mpra-cache", default=None,
                    help="GSE144621 cache directory. Defaults to genomewide.DEFAULT_CACHE, a "
                         "hardcoded D: path. Pass this when that volume is absent or carries the "
                         "cache at a different layout; the refusal names every candidate it found "
                         "and what each one holds, so the value to pass comes from the error rather "
                         "than from a guess. Do not name a drive letter from memory -- letters move.")
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--seeds", type=int, default=10)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import numpy as np
    from glm_encoder_gate import ridge_gc
    from dna_decode.glm.genomewide import (
        DEFAULT_CACHE,
        GenomeWideDataError,
        leave_peak_out,
        load_peak_tiles,
        noise_ceiling,
        spearman,
    )

    # Until 2026-10-22 this read `load_peak_tiles()` bare, so it took the hardcoded `DEFAULT_CACHE`
    # unconditionally and died with an unhandled WinError on any host without that volume. The
    # sibling consumer (`glm_tile_predictions.py`) was parameterised on 2026-10-21 and this one was
    # MISSED, because the AST guard that was supposed to catch exactly this read a single file while
    # being named `..._in_the_tile_chain`. The guard now DERIVES its file set from the callers, so a
    # third consumer cannot be added without being covered.
    #
    # The loader RAISES (typed) and this script PRESENTS -- same layering as the sibling. A traceback
    # here reads as "this tool is broken" rather than "the volume is not attached".
    mpra_cache = a.mpra_cache or str(DEFAULT_CACHE)
    try:
        tiles_all, st = load_peak_tiles(cache_dir=mpra_cache)
    except GenomeWideDataError as exc:
        raise SystemExit(f"--mpra-cache {mpra_cache}\n\n{exc}") from None
    tiles = [t for t in tiles_all if t.category == "tile" and len(t.seq) == 150]
    ceil = noise_ceiling(tiles)["ceiling_for_a_sequence_model"]

    dec, shards = load_decomposition(_require_decomposition_dir(a.in_dir), smoke=a.smoke)
    cov = len(dec) / len(tiles)
    print(f"decomposition: {len(dec)} tiles | {sum(len(v) for v in dec.values())} site rows | "
          f"coverage {cov:.4f} | shards {shards}", flush=True)

    # REFUSE a partial substrate rather than score it. A missing tile is not an empty tile: it enters
    # the feature matrix as the all-zero `no_sites_flag` row, which is a real observation ("this tile
    # has no candidate promoter") and would be silently wrong for a tile that was simply never swept.
    # 0.999 rather than 1.0 leaves room for a genuinely site-free tile, and the message says which.
    if not a.smoke and cov < 0.999:
        raise SystemExit(
            f"coverage {cov:.4f} over {len(tiles)} tiles - the sweep is INCOMPLETE "
            f"({len(tiles) - len(dec)} tiles missing). Scoring now would feed every missing tile in as "
            f"a no-promoter observation and report a clean verdict over a partial substrate. "
            f"Resume stage B: D:/tools/promcalc_env/Scripts/python.exe scripts/promcalc_decompose.py "
            f"--workers 7 --resume"
        )

    # Feature matrix aligned to the FILTERED tile list by position (stage A wrote that index).
    feats = [tile_features(dec.get(i, [])) for i in range(len(tiles))]
    scalar = [[math.log1p(max((max(s["Tx_rate"] for s in dec[i]) if dec.get(i) else 0.0), 0.0))]
              for i in range(len(tiles))]
    nF = len(feature_names())
    assert all(len(f) == nF for f in feats), "ragged feature matrix"

    _pos = {id(t): i for i, t in enumerate(tiles)}
    rng = np.random.default_rng(0)
    runs: dict[str, list[float]] = defaultdict(list)
    for s in range(a.seeds):
        tr, te, held = leave_peak_out(tiles, seed=s)
        # `Tile` is a dataclass with no id field, so position in the filtered list IS the join key --
        # the same key stage A exported. VERIFIED: `leave_peak_out` returns the very objects it was
        # given (not copies), the split is exhaustive and the halves are disjoint. Asserted rather
        # than assumed, because a splitter that started copying would otherwise mis-join every
        # feature row to the wrong tile and still produce a plausible number.
        itr = [_pos[id(t)] for t in tr]
        ite = [_pos[id(t)] for t in te]
        assert not set(itr) & set(ite), f"seed {s}: train/test index overlap"
        ytr = [math.log1p(t.expression) for t in tr]
        yte = [t.expression for t in te]

        g, _ = ridge_gc(tr, te)
        runs["gc"].append(g)

        runs["max_tx"].append(spearman(
            _fit_ridge([scalar[i] for i in itr], ytr, [scalar[i] for i in ite]), yte))

        runs["decomposition"].append(spearman(
            _fit_ridge([feats[i] for i in itr], ytr, [feats[i] for i in ite]), yte))

        # control: destroy the tile->feature pairing, keep the feature DISTRIBUTION intact
        perm = rng.permutation(len(itr))
        runs["control_feature_permuted"].append(spearman(
            _fit_ridge([feats[itr[j]] for j in perm], ytr, [feats[i] for i in ite]), yte))

        print(f"  seed {s}: gc {g:+.4f} | max_tx {runs['max_tx'][-1]:+.4f} | "
              f"decomp {runs['decomposition'][-1]:+.4f} | "
              f"ctrl {runs['control_feature_permuted'][-1]:+.4f} "
              f"(held {held} peaks, {len(te)} test)", flush=True)

    arms = {k: _agg(v) for k, v in runs.items()}
    arms["_paired_vs_max_tx"] = paired(runs["decomposition"], runs["max_tx"])
    arms["_paired_vs_gc"] = paired(runs["decomposition"], runs["gc"])
    arms["_paired_max_tx_vs_gc"] = paired(runs["max_tx"], runs["gc"])
    v, why = verdict(arms)

    out = {
        "schema": "promcalc-decomposition-arm-v1",
        "measured_on": str(date.today()),
        "substrate": {"name": "GSE144621 peak tiles", "n": len(tiles),
                      "length_bp": 150, "replicate_ceiling": ceil,
                      "protocol": "leave_peak_out, registered seeds, ~11,000 test tiles/seed"},
        "n_features": nF, "feature_names": feature_names(),
        "arms": arms, "verdict": v, "verdict_detail": why,
        "preregistered": PREREGISTERED,
        "smoke": a.smoke,
    }
    p = Path(a.out) if a.out else Path(f"wiki/promcalc_decomposition_arm_{date.today()}.json")
    if a.smoke:
        p = p.with_name(p.stem + "_smoke.json")
    p.write_text(json.dumps(out, indent=1), encoding="utf-8")

    print()
    print(f"VERDICT: {v}")
    print(f"  {why}")
    print(f"  ceiling {ceil} | gc {arms['gc']['median']} | max_tx {arms['max_tx']['median']} | "
          f"decomposition {arms['decomposition']['median']}")
    print(f"  wrote {p}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
