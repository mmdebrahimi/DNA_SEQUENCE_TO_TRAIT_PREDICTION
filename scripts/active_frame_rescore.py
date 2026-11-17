"""Re-score the committed PromoterAtlas arms on the active frame. No training, no `D:` volume.

WHAT THIS ANSWERS
-----------------
The published bacterial conv-over-GC margin is +0.0122 on the POOLED frame, and a 2026-11-16
scratch measurement found that restricting to the 14,913 active tiles flips the learned-versus-GC
comparison sign unanimously -- but it did so with a 4-mer ridge proxy and never touched the
STRONGEST arm. The strongest arm's per-tile predictions are already committed, so its active-frame
margin can be read with no model trained and no detached volume: the realistic upper bound on
this substrate.

WHAT IT REFUSES TO DO
---------------------
Arm identity is assigned by RECONCILE, never by filename. `data/promoteratlas/` holds two
prediction files with different shapes, and one of them carries a self-declared `rung` label.
Neither is trusted: each file's pooled Spearman is recomputed from its own stored predictions
and matched against the committed per-seed series in `wiki/asgiven_arm_result_2026-11-02.json`.
If neither series matches, the script exits non-zero with
`INDETERMINATE_ARM_IDENTITY_UNRESOLVED` and reports nothing.

Every active-frame figure carries `not_a_ledger_candidate`. The reason is intrinsic rather than
bookkeeping: `active` is a threshold ON THE LABEL, so an active-frame score is range-restricted
on the outcome and is not a candidate for "best achievable prediction" of the population the
ledger's `ceiling:` criterion names. What IS readable is the differential between two arms under
the SAME restriction, which is why the margin is the headline and the absolutes are not.
"""

from __future__ import annotations

import argparse
import json
import statistics as st
import sys
import time
from datetime import date
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from dna_decode.glm.active_frame import (  # noqa: E402
    ActiveFrameError,
    reproducibility_curve,
    restrict_predictions,
    threshold_audit,
)
from dna_decode.glm.genomewide import leave_peak_out, load_peak_tiles, spearman  # noqa: E402
from scripts.glm_encoder_gate import ridge_gc  # noqa: E402

#: Candidate prediction files. Order is irrelevant: identity comes from the reconcile.
CANDIDATE_FILES = (
    "data/promoteratlas/asgiven_arm_pertile_2026-11-02.json",
    "data/promoteratlas/per_tile_predictions_2026-11-01.json",
)

COMMITTED = "wiki/asgiven_arm_result_2026-11-02.json"

#: Tight enough to reject a 0.01 perturbation, which a test pins so the gate is not vacuous.
RECONCILE_TOLERANCE = 0.001

SEEDS = tuple(range(10))
SPLIT_FRAC = 0.25
N_BOOT = 1000
BOOT_SEED = 707

EXIT_OK = 0
EXIT_IDENTITY = 3
EXIT_SPLIT = 4
EXIT_DATA = 5


class RescoreError(RuntimeError):
    """A refusal that must not be reported as a result."""


# --------------------------------------------------------------------------------------
# loading
# --------------------------------------------------------------------------------------

def normalise_prediction_file(obj: dict) -> dict[int, dict]:
    """Both committed shapes to one: ``{seed: {"tile_indices": [...], "pred": [...]}}``.

    One file wraps seeds under ``per_seed`` and carries a self-declared ``rung``; the other puts
    seeds at the top level. The wrapper's label is deliberately DISCARDED here -- identity is
    settled by reconcile, and reading the label would reintroduce the filename trust this script
    exists to avoid.
    """
    block = obj.get("per_seed") if isinstance(obj.get("per_seed"), dict) else obj
    out: dict[int, dict] = {}
    for k, v in block.items():
        if not isinstance(v, dict) or "pred" not in v or "tile_indices" not in v:
            continue
        try:
            seed = int(k)
        except (TypeError, ValueError):
            continue
        out[seed] = {"tile_indices": list(v["tile_indices"]), "pred": list(v["pred"])}
    if not out:
        raise RescoreError(
            "no per-seed {tile_indices, pred} blocks found under either shape; an empty "
            "normalisation would score nothing and report success"
        )
    return out


def load_committed(path: str | Path) -> dict:
    """The two committed per-seed series plus the split gate they were recorded with."""
    d = json.loads(Path(path).read_text(encoding="utf-8"))
    rung = d["rung"]["per_seed"]
    asgiven = {int(k): float(v["spearman"]) for k, v in rung.items()}
    atlas = {int(k): float(v) for k, v in d["comparator"]["atlas_per_seed"].items()}
    split = {
        int(k): {"n_test": int(v["n_test"]), "n_train": int(v["n_train"]),
                 "held_out_peaks": int(v["held_out_peaks"])}
        for k, v in d["gates"]["split_equality"]["per_seed"].items()
    }
    return {
        "series": {
            "asgiven-position-resolved-promoteratlas-150bp": asgiven,
            "promoteratlas-base-embeddings-frozen-150bp": atlas,
        },
        "medians": {
            "asgiven-position-resolved-promoteratlas-150bp": round(st.median(asgiven.values()), 4),
            "promoteratlas-base-embeddings-frozen-150bp": round(st.median(atlas.values()), 4),
        },
        "split": split,
        "comparator_provenance": {
            "recomputed_in_run": d["comparator"].get("recomputed_in_run"),
            "source": d["comparator"].get("source"),
            "note": "the atlas series was NOT recomputed in its own run; this script recomputes "
                    "it from the stored predictions, which is a stronger check than the one the "
                    "committed artifact performed on itself",
        },
    }


# --------------------------------------------------------------------------------------
# the blocking gates
# --------------------------------------------------------------------------------------

def recompute_pooled(norm: dict[int, dict], tiles) -> dict[int, float]:
    out = {}
    for seed, blk in sorted(norm.items()):
        y, p = restrict_predictions(blk["tile_indices"], blk["pred"], tiles, lambda _t: True)
        out[seed] = round(float(spearman(p, y)), 4)
    return out


def resolve_arm_identity(recomputed: dict[int, float], committed: dict,
                         tolerance: float = RECONCILE_TOLERANCE) -> dict:
    """Assign the arm by whichever committed series the recompute matches on EVERY seed."""
    report = {}
    for name, series in committed["series"].items():
        shared = sorted(set(recomputed) & set(series))
        deltas = {s: round(recomputed[s] - series[s], 4) for s in shared}
        worst = max((abs(v) for v in deltas.values()), default=float("inf"))
        report[name] = {
            "n_seeds_compared": len(shared),
            "max_abs_delta": round(worst, 4),
            "per_seed_delta": deltas,
            "matches": bool(shared) and worst <= tolerance,
        }
    hits = [n for n, r in report.items() if r["matches"]]
    return {
        "tolerance": tolerance,
        "candidates": report,
        "resolved_to": hits[0] if len(hits) == 1 else None,
        "n_matching_series": len(hits),
    }


def check_split(norm: dict[int, dict], tiles, committed: dict) -> dict:
    """Two checks: the committed ``n_test``, and the EXACT held-out index set.

    The plan declares the ``n_test`` cross-check. The exact-set check is strictly stronger and
    in the same spirit: two different splits can agree on size, and a prediction vector aligned
    to the wrong rows scores perfectly happily.
    """
    idx_of = {id(t): i for i, t in enumerate(tiles)}
    rows = {}
    ok = True
    for seed, blk in sorted(norm.items()):
        stored = list(blk["tile_indices"])
        want = committed["split"].get(seed, {}).get("n_test")
        _train, test, held = leave_peak_out(tiles, frac=SPLIT_FRAC, seed=seed)
        reproduced = sorted(idx_of[id(t)] for t in test)
        size_ok = (want is None) or (len(stored) == want)
        set_ok = sorted(stored) == reproduced
        ok &= bool(size_ok and set_ok)
        rows[seed] = {
            "n_stored": len(stored),
            "n_test_committed": want,
            "size_matches_committed": bool(size_ok),
            "n_reproduced_by_splitter": len(reproduced),
            "held_out_peaks_reproduced": int(held),
            "exact_index_set_matches_splitter": bool(set_ok),
        }
    return {"passed": bool(ok), "per_seed": rows,
            "note": "n_test is the plan's declared gate; the exact index-set check is a "
                    "strengthening, since equal sizes do not imply the same rows"}


# --------------------------------------------------------------------------------------
# the active-frame rescore
# --------------------------------------------------------------------------------------

def _paired_bootstrap(margins: list[float], n_boot: int = N_BOOT, seed: int = BOOT_SEED) -> dict:
    rng = np.random.default_rng(seed)
    arr = np.asarray(margins, dtype=float)
    draws = arr[rng.integers(0, arr.size, size=(n_boot, arr.size))].mean(axis=1)
    return {
        "mean": round(float(arr.mean()), 4),
        "ci_lo": round(float(np.percentile(draws, 2.5)), 4),
        "ci_hi": round(float(np.percentile(draws, 97.5)), 4),
        "n_boot": n_boot,
        "excludes_zero": bool(np.percentile(draws, 2.5) > 0 or np.percentile(draws, 97.5) < 0),
    }


def _binned_margin(bin_by, arm_pred, gc_pred, y, n_bins: int) -> dict:
    """The arm-versus-GC margin inside equal-sized rank bins of ``bin_by``."""
    v = np.asarray(bin_by, dtype=float)
    a = np.asarray(arm_pred, dtype=float)
    g = np.asarray(gc_pred, dtype=float)
    t = np.asarray(y, dtype=float)
    order = np.argsort(v, kind="stable")
    rows = []
    for i, idx in enumerate(np.array_split(order, n_bins)):
        ra = spearman(a[idx].tolist(), t[idx].tolist())
        rg = spearman(g[idx].tolist(), t[idx].tolist())
        rows.append({"bin": i, "n": int(idx.size),
                     "arm": round(float(ra), 4), "gc": round(float(rg), 4),
                     "margin": round(float(ra - rg), 4)})
    return {"n_bins": n_bins, "bins": rows,
            "n_total": int(sum(r["n"] for r in rows))}


def rescore(tiles, norm: dict[int, dict], *, n_bins: int = 10) -> dict:
    """Pooled and active-frame scores per seed, plus the GC comparator on the same rows."""
    active = [t.active == "active" for t in tiles]
    per_seed = {}
    for seed, blk in sorted(norm.items()):
        train, test, _held = leave_peak_out(tiles, frac=SPLIT_FRAC, seed=seed)
        gc_r, gc_pred = ridge_gc(train, test)
        idx_of = {id(t): i for i, t in enumerate(tiles)}
        test_idx = [idx_of[id(t)] for t in test]
        gc_by_tile = dict(zip(test_idx, gc_pred))

        stored = list(blk["tile_indices"])
        arm_by_tile = dict(zip(stored, blk["pred"]))

        rows = [i for i in stored if i in gc_by_tile]
        y_all = [float(tiles[i].expression) for i in rows]
        arm_all = [float(arm_by_tile[i]) for i in rows]
        gc_all = [float(gc_by_tile[i]) for i in rows]
        act = [i for i in rows if active[i]]
        y_act = [float(tiles[i].expression) for i in act]
        arm_act = [float(arm_by_tile[i]) for i in act]
        gc_act = [float(gc_by_tile[i]) for i in act]

        per_seed[seed] = {
            "n_test": len(rows),
            "n_test_active": len(act),
            "pooled_arm": round(float(spearman(arm_all, y_all)), 4),
            "pooled_gc": round(float(spearman(gc_all, y_all)), 4),
            "pooled_gc_from_ridge_gc_return": round(float(gc_r), 4),
            "active_arm": round(float(spearman(arm_act, y_act)), 4),
            "active_gc": round(float(spearman(gc_act, y_act)), 4),
            "not_a_ledger_candidate": True,
            "why_not_a_ledger_candidate": (
                "the active frame is a threshold ON THE LABEL, so these are range-restricted "
                "on the outcome and are not candidates for best achievable prediction of the "
                "population the ceiling: criterion names"
            ),
        }
        per_seed[seed]["pooled_margin"] = round(
            per_seed[seed]["pooled_arm"] - per_seed[seed]["pooled_gc"], 4)
        per_seed[seed]["active_margin"] = round(
            per_seed[seed]["active_arm"] - per_seed[seed]["active_gc"], 4)
        if seed == 0:
            per_seed[seed]["binned_margin_by_label"] = dict(
                _binned_margin(y_all, arm_all, gc_all, y_all, n_bins),
                binning_variable="expn_med, the SCORED label",
                conditions_on_the_outcome=True,
                reads_as="binning on the scored label conditions on the outcome; the MARGIN is "
                         "still readable because both arms face the identical restriction, but "
                         "neither absolute is comparable to a pooled number",
            )
            per_seed[seed]["binned_margin_by_group1"] = dict(
                _binned_margin([float(tiles[i].rep1) for i in rows], arm_all, gc_all, y_all,
                               n_bins),
                binning_variable="RNA_exp_sum_1, group 1's published summary",
                conditions_on_the_outcome=False,
                reads_as="PREFERRED: bins are assigned by a measurement rather than by the "
                         "scored label, so the bin assignment does not condition on the outcome",
            )

    pooled_m = [per_seed[s]["pooled_margin"] for s in sorted(per_seed)]
    active_m = [per_seed[s]["active_margin"] for s in sorted(per_seed)]
    return {
        "per_seed": per_seed,
        "pooled": {
            "arm_median": round(st.median(per_seed[s]["pooled_arm"] for s in per_seed), 4),
            "gc_median": round(st.median(per_seed[s]["pooled_gc"] for s in per_seed), 4),
            "margin_median": round(st.median(pooled_m), 4),
            "n_seeds_arm_wins": sum(1 for m in pooled_m if m > 0),
            "bootstrap": _paired_bootstrap(pooled_m),
            "conditions_on_the_outcome": False,
        },
        "active": {
            "arm_median": round(st.median(per_seed[s]["active_arm"] for s in per_seed), 4),
            "gc_median": round(st.median(per_seed[s]["active_gc"] for s in per_seed), 4),
            "margin_median": round(st.median(active_m), 4),
            "n_seeds_arm_wins": sum(1 for m in active_m if m > 0),
            "bootstrap": _paired_bootstrap(active_m),
            "conditions_on_the_outcome": True,
            "not_a_ledger_candidate": True,
        },
    }


# --------------------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------------------

def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--cache", default="data/mpra",
                    help="tile cache dir; the D:/ default in genomewide is detached here")
    ap.add_argument("--committed", default=COMMITTED)
    ap.add_argument("--files", nargs="*", default=list(CANDIDATE_FILES))
    ap.add_argument("--n-bins", type=int, default=10)
    ap.add_argument("--out", default=None)
    args = ap.parse_args(argv)

    t0 = time.time()
    try:
        loaded, _st = load_peak_tiles(args.cache, allow_download=False)
    except Exception as e:  # noqa: BLE001
        print(f"INDETERMINATE_SUBSTRATE_UNREADABLE: {e}")
        return EXIT_DATA
    tiles = [t for t in loaded if t.category == "tile" and len(t.seq) == 150]
    print(f"loaded {len(tiles)} length-150 tiles from {args.cache}")

    audit = threshold_audit(tiles)
    print(f"active flag: {audit['n_active']} active / {audit['n_inactive']} inactive, "
          f"implied cut {audit['implied_cut']}, overlap {audit['overlap_either_side']}")

    committed = load_committed(args.committed)
    print(f"committed medians: {committed['medians']}")

    arms = {}
    for path in args.files:
        p = Path(path)
        if not p.exists():
            print(f"INDETERMINATE_PREDICTIONS_ABSENT: {p}")
            return EXIT_DATA
        norm = normalise_prediction_file(json.loads(p.read_text(encoding="utf-8")))
        pooled = recompute_pooled(norm, tiles)
        ident = resolve_arm_identity(pooled, committed)
        print(f"\n{p.name}: recomputed pooled median "
              f"{round(st.median(pooled.values()), 4)} over {len(pooled)} seeds")
        for name, r in ident["candidates"].items():
            print(f"   vs {name:<48s} max|delta| {r['max_abs_delta']:.4f}  "
                  f"{'MATCH' if r['matches'] else 'no'}")
        if ident["resolved_to"] is None:
            print(f"INDETERMINATE_ARM_IDENTITY_UNRESOLVED for {p.name}: "
                  f"{ident['n_matching_series']} series matched within {ident['tolerance']}")
            return EXIT_IDENTITY
        split = check_split(norm, tiles, committed)
        if not split["passed"]:
            bad = [s for s, r in split["per_seed"].items()
                   if not (r["size_matches_committed"] and r["exact_index_set_matches_splitter"])]
            print(f"INDETERMINATE_SPLIT_MISMATCH for {p.name}: seeds {bad}")
            return EXIT_SPLIT
        print(f"   resolved to {ident['resolved_to']}; split gate passed on all "
              f"{len(split['per_seed'])} seeds")
        arms[ident["resolved_to"]] = {
            "file": str(p), "identity": ident, "split_gate": split,
            "recomputed_pooled_per_seed": pooled,
            "scores": rescore(tiles, norm, n_bins=args.n_bins),
        }

    out = {
        "schema": "active-frame-rescore-v1",
        "date": str(date.today()),
        "substrate": {"cache": args.cache, "n_tiles": len(tiles)},
        "threshold_audit": audit,
        "committed_source": args.committed,
        "comparator_provenance": committed["comparator_provenance"],
        "committed_medians": committed["medians"],
        "arms": arms,
        "no_model_was_trained": True,
        "reads_as": (
            "the MARGIN between an arm and GC under the SAME restriction is the readable "
            "quantity. Active-frame absolutes are range-restricted on the outcome and carry "
            "not_a_ledger_candidate; do not compare them to pooled numbers or to the 0.941 "
            "ceiling, whose own frame is what the parent plan is investigating"
        ),
        "wall_seconds": round(time.time() - t0, 1),
    }
    dest = Path(args.out) if args.out else Path(f"wiki/active_frame_rescore_{date.today()}.json")
    dest.write_text(json.dumps(out, indent=2), encoding="utf-8")

    print("\n" + "=" * 74)
    for name, a in arms.items():
        s = a["scores"]
        print(f"{name}")
        print(f"   pooled  arm {s['pooled']['arm_median']:+.4f}  gc {s['pooled']['gc_median']:+.4f}"
              f"  margin {s['pooled']['margin_median']:+.4f}"
              f"  wins {s['pooled']['n_seeds_arm_wins']}/10  CI {s['pooled']['bootstrap']['ci_lo']:+.4f}"
              f",{s['pooled']['bootstrap']['ci_hi']:+.4f}")
        print(f"   active  arm {s['active']['arm_median']:+.4f}  gc {s['active']['gc_median']:+.4f}"
              f"  margin {s['active']['margin_median']:+.4f}"
              f"  wins {s['active']['n_seeds_arm_wins']}/10  CI {s['active']['bootstrap']['ci_lo']:+.4f}"
              f",{s['active']['bootstrap']['ci_hi']:+.4f}")
    print("=" * 74)
    print(f"wrote {dest}")
    return EXIT_OK


if __name__ == "__main__":
    raise SystemExit(main())
