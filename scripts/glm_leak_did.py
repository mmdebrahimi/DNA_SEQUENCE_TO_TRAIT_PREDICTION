"""Is the leak real, or is it just too small to see at 120 windows? The DiD, at any n.

WHY THIS EXISTS. `wiki/glm_representative_repretrain_graded_2026-11-25.json` records the one
quantity that separates a TRAINING LEAK from ordinary CROSS-ORGANISM cost, and records that it could
not be separated:

    DiD        = mean paired delta on `ortholog_bearing` minus the same on `withheld_orthology_free`
    delta      = committed_CE - assessed_CE per window, so POSITIVE = the assessed model is better
    predicted  = -0.0054 nats   (removing the draft removes the near-ortholog from the bearing half
                                 only, while cross-organism cost acts on both halves and cancels)
    measured   = -0.001897 nats, CI [-0.006801, +0.002813] -- INCLUDES ZERO
    noise floor=  0.00281 nats

The interval is ~0.0096 wide and the predicted effect is 0.0054, so at 120 windows the measurement
cannot distinguish the predicted leak from nothing. That is a POWER limit, not a null result, and it
was dismissed once as "more eval windows is trivial" -- which was right about headroom on the
cross-entropy itself and WRONG about this DiD.

THE HEADROOM IS REAL AND WAS MEASURED, NOT ASSUMED. `ortholog_bearing` has 2,830 distinct windows
available and `withheld_orthology_free` has 519; both were scored at 120. The clean half is the
binding constraint at 519, which is 4.3x the current n, so a balanced re-score should shrink the
interval by about sqrt(4.3) ~ 2.1x -- enough that the predicted -0.0054 falls outside an interval
centred near the measured -0.0019.

WHAT THIS SCRIPT DOES, AND WHAT IT DOES NOT. It does NOT score anything: scoring is
`scripts/glm_beyond_kmer.py --eval-set <set> --n-windows <n> --out <path>`, which already persists
per-window rows. This reads two scored artifacts per set -- one for the committed checkpoint, one
for an assessed checkpoint -- and computes the DiD with an interval. So it is torch-free and the
scoring cost is paid once and reused.

THE GEOMETRY GATE IS THE BASIS OF TRUST, AND IT IS STRICTER THAN IT LOOKS. `select_eval_windows`
draws its sample with `random.Random(selection_seed).sample(pool, n_windows)`, and Python's `sample`
has NO prefix property: the 120-window draw is not the first 120 of the 500-window draw. So an
artifact scored at a different `n_windows` describes DIFFERENT WINDOWS, and pairing across them is
meaningless. Every artifact in one DiD must therefore agree on `n_windows`, `positions`, `seed` and
`selection_seed`, and the script REFUSES (exit 2) otherwise rather than silently pairing window i of
one draw against window i of another.

THE BOOTSTRAP RESAMPLES THE TWO SETS INDEPENDENTLY, which is what the design calls for: the halves
are DISJOINT sets of windows, so there is no pairing between them to preserve. The within-set delta
IS paired (same window, two models) and that pairing is kept.

HONEST LIMITS.
* The DiD identifies the leak only under the stated assumption that cross-organism cost acts equally
  on both halves and cancels. That assumption is the design's, not this script's, and it is not
  testable from these numbers.
* A resolved DiD bounds the leak on THIS corpus and THIS pair of checkpoints. It is not a general
  statement about leakage.
* Narrowing an interval does not make a point estimate correct; an interval that excludes the
  predicted size says the prediction was wrong, not that the true value is the new centre.
"""
from __future__ import annotations

import argparse
import json
import random
import statistics
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

BEARING = "ortholog_bearing"
CLEAN = "withheld_orthology_free"

#: From the graded artifact. Quoted, not re-derived, so this script cannot drift from the claim it
#: is testing.
PREDICTED_DID_NATS = -0.0054
PRIOR_DID_NATS = -0.001897
PRIOR_CI = [-0.006801, 0.002813]
NOISE_FLOOR_NATS = 0.00281
GRADED = "wiki/glm_representative_repretrain_graded_2026-11-25.json"

#: FROZEN BEFORE ANY NUMBER AT THE NEW n EXISTS. Four outcomes, and the one that matters most is
#: the third: at a larger n the interval can RULE OUT the predicted leak size without resolving a
#: nonzero effect, which is a real finding and not a null.
PREREGISTERED = {
    "quantity": ("DiD = mean(delta on ortholog_bearing) - mean(delta on withheld_orthology_free), "
                 "delta = committed_CE - assessed_CE per window (positive = assessed better)"),
    "predicted_nats": PREDICTED_DID_NATS,
    "prior_measurement": {"n_windows_per_set": 120, "did_nats": PRIOR_DID_NATS, "ci95": PRIOR_CI,
                          "noise_floor_nats": NOISE_FLOOR_NATS, "source": GRADED},
    "outcomes": {
        "LEAK_RESOLVED_AT_PREDICTED_SIZE":
            "CI excludes zero AND contains the predicted -0.0054",
        "LEAK_RESOLVED_BUT_SMALLER_THAN_PREDICTED":
            "CI excludes zero AND lies entirely between zero and the predicted size -- a real leak, "
            "smaller than the design expected",
        "LEAK_EXCLUDED_AT_PREDICTED_SIZE":
            "CI contains zero AND excludes the predicted -0.0054 -- the predicted leak size is "
            "ruled out, which the 120-window run could not do",
        "STILL_UNRESOLVED":
            "CI contains BOTH zero and the predicted size -- more windows did not settle it",
    },
    "resolution_rule": "a percentile bootstrap 95% CI, matching the committed scorer's convention",
    "why_this_is_not_a_relabelled_null": (
        "the 120-window interval contained both zero and the predicted size, so it distinguished "
        "nothing. An interval that contains zero while EXCLUDING -0.0054 is a different and "
        "stronger statement, and it is reachable only at larger n."),
}


class DidError(RuntimeError):
    pass


GEOMETRY_KEYS = ("n_windows", "positions", "seed", "selection_seed")


def load_artifact(path: Path) -> dict:
    d = json.loads(path.read_text(encoding="utf-8"))
    for k in ("per_window", "run_args"):
        if not d.get(k):
            raise DidError(f"{path.name} lacks {k}")
    return d


def deltas(committed: dict, assessed: dict, *, arm: str) -> list[float]:
    """committed_CE - assessed_CE per window. Positive = the assessed model is better.

    Both artifacts must describe the SAME windows; `check_geometry` is what establishes that, and
    it must be called first. Length equality alone is not sufficient -- two different 120-window
    draws are both length 120.
    """
    c, s = committed["per_window"], assessed["per_window"]
    if len(c) != len(s):
        raise DidError(f"window counts differ: {len(c)} vs {len(s)}")
    return [cw[arm] - sw[arm] for cw, sw in zip(c, s)]


def check_geometry(arts: dict) -> dict:
    """Every artifact must agree on the draw, and each must sit on its intended eval set.

    Refuses rather than pairing across draws. See the module docstring on `sample`'s lack of a
    prefix property -- this is the check that makes the whole comparison meaningful.
    """
    geoms = {}
    for label, (path, d) in arts.items():
        ra = d["run_args"]
        geoms[label] = {k: ra.get(k) for k in GEOMETRY_KEYS}
        geoms[label]["eval_set"] = ra.get("eval_set")
        geoms[label]["n_per_window_rows"] = len(d["per_window"])
    draws = {tuple(g[k] for k in GEOMETRY_KEYS) for g in geoms.values()}
    if len(draws) != 1:
        raise DidError(f"artifacts disagree on the draw geometry: {geoms}")
    for label, g in geoms.items():
        want = BEARING if "bearing" in label else CLEAN
        if g["eval_set"] != want:
            raise DidError(f"{label} is on eval_set {g['eval_set']!r}, expected {want!r}")
        if g["n_per_window_rows"] != g["n_windows"]:
            raise DidError(f"{label}: {g['n_per_window_rows']} rows vs n_windows {g['n_windows']}")
    return geoms


def did_with_ci(d_bearing: list[float], d_clean: list[float], *, seed: int, boot: int) -> dict:
    """The DiD and a percentile interval, resampling the two DISJOINT sets independently."""
    if not d_bearing or not d_clean:
        raise DidError("both sets need at least one window")
    point = statistics.fmean(d_bearing) - statistics.fmean(d_clean)
    r = random.Random(seed)
    draws = []
    for _ in range(boot):
        b = statistics.fmean([d_bearing[r.randrange(len(d_bearing))] for _ in d_bearing])
        c = statistics.fmean([d_clean[r.randrange(len(d_clean))] for _ in d_clean])
        draws.append(b - c)
    draws.sort()
    lo, hi = draws[int(0.025 * len(draws))], draws[int(0.975 * len(draws))]
    return {
        "did_nats": round(point, 6),
        "ci95": [round(lo, 6), round(hi, 6)],
        "ci_width": round(hi - lo, 6),
        "excludes_zero": bool(lo > 0 or hi < 0),
        "excludes_predicted": bool(lo > PREDICTED_DID_NATS or hi < PREDICTED_DID_NATS),
        "mean_delta_bearing": round(statistics.fmean(d_bearing), 6),
        "mean_delta_clean": round(statistics.fmean(d_clean), 6),
        "n_bearing": len(d_bearing),
        "n_clean": len(d_clean),
        "boot": boot,
        "seed": seed,
    }


def verdict_for(stat: dict) -> dict:
    """The frozen bar, applied by code. Ordered; first match wins."""
    zero, pred = stat["excludes_zero"], stat["excludes_predicted"]
    if zero and not pred:
        v = "LEAK_RESOLVED_AT_PREDICTED_SIZE"
    elif zero and pred:
        v = "LEAK_RESOLVED_BUT_SMALLER_THAN_PREDICTED"
    elif not zero and pred:
        v = "LEAK_EXCLUDED_AT_PREDICTED_SIZE"
    else:
        v = "STILL_UNRESOLVED"
    return {"VERDICT": v, "reading": PREREGISTERED["outcomes"][v],
            "excludes_zero": zero, "excludes_predicted_size": pred}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--committed-bearing", type=Path, required=True)
    ap.add_argument("--committed-clean", type=Path, required=True)
    ap.add_argument("--assessed-bearing", type=Path, required=True)
    ap.add_argument("--assessed-clean", type=Path, required=True)
    ap.add_argument("--label", default="assessed", help="name for the assessed checkpoint")
    ap.add_argument("--arm", default="left_only", choices=("left_only", "bidirectional"))
    ap.add_argument("--boot", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()

    arts = {}
    for label, p in (("committed_bearing", a.committed_bearing),
                     ("committed_clean", a.committed_clean),
                     ("assessed_bearing", a.assessed_bearing),
                     ("assessed_clean", a.assessed_clean)):
        if not p.exists():
            print(f"REFUSING: missing {p}", file=sys.stderr)
            return 2
        try:
            arts[label] = (p, load_artifact(p))
        except DidError as exc:
            print(f"REFUSING: {exc}", file=sys.stderr)
            return 2

    try:
        geoms = check_geometry(arts)
        d_bearing = deltas(arts["committed_bearing"][1], arts["assessed_bearing"][1], arm=a.arm)
        d_clean = deltas(arts["committed_clean"][1], arts["assessed_clean"][1], arm=a.arm)
    except DidError as exc:
        print(f"REFUSING: {exc}", file=sys.stderr)
        return 2

    stat = did_with_ci(d_bearing, d_clean, seed=a.seed, boot=a.boot)
    verdict = verdict_for(stat)
    n = stat["n_bearing"]
    prior_width = PRIOR_CI[1] - PRIOR_CI[0]

    out = {
        "record": "glm-leak-did-v1",
        "arm": a.arm,
        "assessed_label": a.label,
        "PREREGISTERED": PREREGISTERED,
        "geometry": geoms,
        "inputs": {k: str(v[0]) for k, v in arts.items()},
        "checkpoint_steps": {k: v[1].get("checkpoint_step") for k, v in arts.items()},
        "result": stat,
        "VERDICT": verdict,
        "comparison_to_120_windows": {
            "prior_did_nats": PRIOR_DID_NATS,
            "prior_ci_width": round(prior_width, 6),
            "this_ci_width": stat["ci_width"],
            "width_ratio": round(stat["ci_width"] / prior_width, 3) if prior_width else None,
            "n_per_set_now": n,
            "note": ("a narrower interval is expected from n alone; the VERDICT, not the width, is "
                     "what changed or did not"),
        },
        "honest_limits": [
            "The DiD identifies the leak only if cross-organism cost acts equally on both halves "
            "and cancels. That is the design's assumption and is not testable from these numbers.",
            "Scoped to this corpus and this checkpoint pair; not a general claim about leakage.",
            "An interval excluding the predicted size says the prediction was wrong, NOT that the "
            "true value is the new centre.",
            "The two halves are disjoint, so the DiD is a difference of independent means; the "
            "within-set delta is paired per window and that pairing is preserved.",
        ],
    }
    a.out.write_text(json.dumps(out, indent=1) + "\n", encoding="utf-8")

    print(f"arm {a.arm} | n per set: bearing {stat['n_bearing']}, clean {stat['n_clean']}")
    print(f"  mean delta bearing {stat['mean_delta_bearing']:+.6f}  "
          f"clean {stat['mean_delta_clean']:+.6f}")
    print(f"  DiD {stat['did_nats']:+.6f}  CI {stat['ci95']}  width {stat['ci_width']:.6f} "
          f"(was {prior_width:.6f} at n=120)")
    print(f"  excludes zero: {stat['excludes_zero']} | "
          f"excludes predicted {PREDICTED_DID_NATS}: {stat['excludes_predicted']}")
    print(f"VERDICT {verdict['VERDICT']}")
    print(f"wrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
