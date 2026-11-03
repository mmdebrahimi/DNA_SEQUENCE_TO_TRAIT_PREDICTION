"""Score the shipped-expression-heads arm against its frozen bar.

The arm applies five trained prediction heads ZERO-SHOT to the GSE144621 peak-tile rung and asks
whether any of them beats the frozen encoder's median Spearman 0.5076 by the inherited 0.02 margin.

**Everything decidable was decided in the bar, and this runner only executes it.** The windowing
route (the faithful per-window recompute), the primary reduction (mean over 65 windows), the
ineligible secondary (max), the correction (alpha/5 over five heads), the scoring definition (direct
Spearman, sign as given) and the null (fixed predictor, permuted TEST labels, no refit) are all read
from `wiki/expression_heads_acceptance_bar.json`. Nothing here chooses any of them.

**Why this STREAMS instead of computing the latents once and reusing them.** The obvious shape --
one forward pass into an array, then five heads over it -- would hold
`44,106 x 128 x 150 x 4 bytes = 3.39 GB` of float32 latent at once. This host has around 6 GB free
and has already had a background job reaped for system memory pressure, so that shape is not
available. Instead each batch computes its latent ONCE, immediately applies all five heads to it,
reduces to scalars, and discards the latent. What persists is
`44,106 x 5 heads x 2 reductions` doubles, which is under 4 MB. The "one forward pass per tile"
property that justifies sharing the encoder across heads is preserved exactly; only the lifetime of
the latent changes.

**The true route is the expensive half and that is a bar decision, not an accident.** Per tile per
head it recomputes the convolutional trunk 65 times, which the feasibility probe measured at 15.5x
to 23.3x the cheap route and projected at 4.64 h over all five heads. The cheap route was rejected
because one head's rank agreement with the faithful route (0.9888) fell below the margin-derived
tolerance (0.999).

**A capped run is stamped PARTIAL and the loader refuses it**, so a smoke run can never be read as
the arm's result.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from dna_decode.glm import atlas_heads as AH
from dna_decode.glm.atlas_embed import DEFAULT_PA_REPO
from scripts.freeze_expression_heads_bar import (
    ARM,
    CHECK_ORDER,
    NULL_SIGNAL_RATIO,
    OUTCOMES,
    PRIMARY_REDUCTION,
    ROUTE_TRUE,
    SECONDARY_REDUCTION,
    classify,
)
from scripts.positional_arm_probe import free_ram_bytes, load_arm_cohort
from scripts.positional_arm_score import (
    ScoreError,
    leave_peak_out_indices,
    median_of,
    paired_bootstrap,
)

REPO = Path(__file__).resolve().parents[1]
BAR_PATH = REPO / "wiki/expression_heads_acceptance_bar.json"
OUT = REPO / "wiki/expression_heads_result_2026-11-02.json"
#: per-tile scalars go under data/, never wiki/ -- a large decimal dump in wiki/ has previously
#: polluted the contract-number audit's decoy pool and broken two of its guards
BLOCK_DIR = REPO / "data/promoteratlas/expression_heads_block"
BLOCK_NAME = "head_scalars.npy"
SIDECAR_NAME = "head_scalars_sidecar.json"

REDUCTIONS = (PRIMARY_REDUCTION, SECONDARY_REDUCTION)
BATCH = 64
N_BOOT = 10_000
NULL_SEEDS = 10

#: Gate -> the outcome a failure returns, derived from the frozen bar's own check order.
GATE_SEQUENCE = tuple((g, o) for g, o in CHECK_ORDER if isinstance(o, str))


# ---------------------------------------------------------------------------------------------------
# streaming extraction -- resumable, PARTIAL-stamped
# ---------------------------------------------------------------------------------------------------
def column_layout(head_names: list[str]) -> dict:
    """`(head, reduction) -> column index`, plus `_total`. Pure and order-stable."""
    layout: dict = {}
    k = 0
    for name in sorted(head_names):
        for red in REDUCTIONS:
            layout["%s::%s" % (name, red)] = k
            k += 1
    layout["_total"] = k
    return layout


def _sidecar_path(block_dir: Path) -> Path:
    return Path(block_dir) / SIDECAR_NAME


def read_sidecar(block_dir: Path) -> dict | None:
    p = _sidecar_path(block_dir)
    if not p.is_file():
        return None
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return None


def load_scalars(block_dir: Path, *, allow_partial: bool = False):
    """Load the per-tile head scalars, REFUSING a PARTIAL block unless explicitly allowed.

    A capped run cannot be read as the arm's result, which is the whole reason the stamp exists.
    """
    side = read_sidecar(block_dir)
    if side is None:
        raise ScoreError("no sidecar at %s; run the extraction first" % _sidecar_path(block_dir))
    if side.get("status") != "COMPLETE" and not allow_partial:
        raise ScoreError(
            "block status is %r, not COMPLETE (%d of %d rows). Refusing: a capped run is not the "
            "arm's result." % (side.get("status"), side.get("n_done", -1), side.get("n_rows", -1))
        )
    arr = np.load(Path(block_dir) / BLOCK_NAME, mmap_mode="r")
    if arr.shape != (side["n_rows"], side["layout"]["_total"]):
        raise ScoreError("block is %r, sidecar says %r"
                         % (arr.shape, (side["n_rows"], side["layout"]["_total"])))
    return arr, side


def extract(*, tiles, heads, block_dir: Path, limit: int | None = None,
            batch: int = BATCH, encode=None, verbose: bool = True, heartbeat=None) -> dict:
    """One streaming pass: per batch, latent once, all heads, reduce, discard the latent."""
    import torch

    block_dir = Path(block_dir)
    block_dir.mkdir(parents=True, exist_ok=True)
    names = [h.name for h in heads]
    layout = column_layout(names)
    n_rows = len(tiles)
    n_target = n_rows if limit is None else min(limit, n_rows)

    path = block_dir / BLOCK_NAME
    side = read_sidecar(block_dir)
    resume_from = 0
    if side and side.get("n_rows") == n_rows and side.get("layout") == layout and path.is_file():
        resume_from = int(side.get("n_done", 0))
        arr = np.lib.format.open_memmap(path, mode="r+")
    else:
        arr = np.lib.format.open_memmap(path, mode="w+", dtype="float64",
                                        shape=(n_rows, layout["_total"]))

    if encode is None:
        from scripts.positional_arm_extract import make_atlas_encoder

        encode = make_atlas_encoder()

    t0 = time.time()
    done = resume_from
    while done < n_target:
        hi = min(done + batch, n_target)
        seqs = [tiles[i].seq for i in range(done, hi)]
        z = torch.from_numpy(np.asarray(encode(seqs), dtype=np.float32))
        if not torch.isfinite(z).all():
            raise ScoreError("batch starting at %d produced non-finite latents" % done)
        for h in heads:
            per_window = AH.head_forward_per_window(z, h.state_dict)   # the FAITHFUL route
            for red in REDUCTIONS:
                col = layout["%s::%s" % (h.name, red)]
                arr[done:hi, col] = AH.reduce_windows(per_window, red).detach().numpy()
        del z
        done = hi
        arr.flush()
        _sidecar_path(block_dir).write_text(json.dumps({
            "schema": "expression-heads-scalars-v1",
            "status": "COMPLETE" if done >= n_rows else "PARTIAL",
            "n_rows": n_rows, "n_done": done, "layout": layout,
            "heads": sorted(names), "reductions": list(REDUCTIONS),
            "route": ROUTE_TRUE, "batch": batch,
            "pa_repo": str(DEFAULT_PA_REPO),
            "wall_seconds": round(time.time() - t0, 1),
            "why_streamed": ("the full-cohort latent would be 3.39 GB of float32; each batch "
                             "computes its latent once, applies every head, reduces, and discards"),
        }, indent=1, sort_keys=True), encoding="utf-8")
        # The heartbeat is the surface an UNATTENDED reader checks. Without it a job running for
        # hours classifies NOT_STARTED, which is exactly the conflation the state set exists to
        # prevent. Written per batch; the Heartbeat object enforces its own cadence.
        if heartbeat is not None:
            heartbeat.beat(done, force=(done >= n_target or done == resume_from + batch),
                           resumed_from=resume_from, n_heads=len(heads), route=ROUTE_TRUE)
        if verbose and (done % (batch * 25) == 0 or done >= n_target):
            el = time.time() - t0
            rate = (done - resume_from) / el if el > 0 else 0.0
            left = (n_target - done) / rate if rate > 0 else float("nan")
            print("  %6d/%6d tiles  %.1f tiles/s  eta %.1f min"
                  % (done, n_target, rate, left / 60.0), flush=True)

    return read_sidecar(block_dir) or {}


# ---------------------------------------------------------------------------------------------------
# gates
# ---------------------------------------------------------------------------------------------------
def gate_split_equality(ctx: dict) -> dict:
    """The committed per-seed test-set sizes must reproduce exactly."""
    want, got, bad = ctx["committed_n_test"], {}, []
    for s in ctx["seeds"]:
        _tr, te, _ = ctx["splits"][s]
        got[str(s)] = len(te)
        if int(want[str(s)]) != len(te):
            bad.append((s, int(want[str(s)]), len(te)))
    return {"passed": not bad, "n_test_per_seed": got,
            "detail": ("all %d seeds reproduce the committed test-set sizes" % len(ctx["seeds"]))
            if not bad else "seeds %r differ from the committed sizes" % bad}


def gate_geometry(ctx: dict) -> dict:
    """Every head's fc1 is (64, 11008) and the window count is 65. Re-checked at SCORE time."""
    bad = [h.name for h in ctx["heads"]
           if tuple(h.state_dict["fc1.weight"].shape) != AH.FC1_SHAPE]
    nw = AH.n_windows(ctx["length"])
    ok = (not bad) and nw == ctx["bar"]["windowing"]["n_windows"]
    return {"passed": bool(ok), "n_windows": nw, "heads_with_wrong_fc1": bad,
            "detail": ("%d heads at fc1 %r, %d windows" % (len(ctx["heads"]), AH.FC1_SHAPE, nw))
            if ok else "wrong geometry: heads %r, windows %d" % (bad, nw)}


def gate_head_validity(ctx: dict) -> dict:
    """Finite and NON-CONSTANT per head. A collapsed head would score as a clean null."""
    rows, bad = {}, []
    for h in ctx["heads"]:
        col = ctx["layout"]["%s::%s" % (h.name, PRIMARY_REDUCTION)]
        v = np.asarray(ctx["arr"][:, col], dtype=float)
        finite = bool(np.isfinite(v).all())
        sd = float(np.std(v))
        rows[h.name] = {"finite": finite, "std": sd, "n_distinct": int(len(np.unique(v)))}
        if not finite or sd == 0.0:
            bad.append(h.name)
    return {"passed": not bad, "per_head": rows,
            "detail": "all %d heads finite and non-constant" % len(ctx["heads"]) if not bad
            else "head(s) %r are non-finite or constant" % bad}


def gate_permutation_null(ctx: dict) -> dict:
    """FIXED predictor, PERMUTED test labels, NO refit -- per the amended bar.

    Both conditions must hold: the median must be near zero, AND the magnitude must be small
    relative to the arm's own real signal. The second is what catches a null that cannot destroy
    the signal it exists to destroy, which is what the inherited shuffled-TRAINING-label null did
    on a one-feature predictor.
    """
    from dna_decode.glm.genomewide import spearman

    margin = ctx["bar"]["inherited_not_chosen"]["margin"]["value"]
    best_name, best_med = ctx["best_head"], ctx["best_median"]
    col = ctx["layout"]["%s::%s" % (best_name, PRIMARY_REDUCTION)]
    pred = np.asarray(ctx["arr"][:, col], dtype=float)

    nulls = []
    for s in range(NULL_SEEDS):
        tr, te, _ = ctx["splits"][ctx["seeds"][s % len(ctx["seeds"])]]
        y = np.asarray(ctx["ys"](te), dtype=float)
        rng = np.random.default_rng(90_000 + s)
        yp = y.copy()
        rng.shuffle(yp)
        nulls.append(round(float(spearman(list(pred[te]), list(yp))), 4))

    med = float(np.median(nulls))
    mx = float(np.max(np.abs(nulls)))
    ratio = (mx / abs(best_med)) if best_med else float("inf")
    ok_med = abs(med) < margin
    ok_mag = ratio < NULL_SIGNAL_RATIO
    n_test = len(ctx["splits"][ctx["seeds"][0]][1])
    return {
        "passed": bool(ok_med and ok_mag),
        "form": "fixed predictor, permuted TEST labels, NO refit",
        "head": best_name, "per_seed": nulls,
        "median": round(med, 4), "max_abs": round(mx, 4),
        "mean_abs": round(float(np.mean(np.abs(nulls))), 4),
        "real_median_spearman": round(float(best_med), 4),
        "max_abs_over_real": round(ratio, 4),
        "condition_median_passed": bool(ok_med),
        "condition_magnitude_passed": bool(ok_mag),
        "sampling_scale_1_over_sqrt_n_test": round(1.0 / np.sqrt(max(1, n_test)), 4),
        "detail": ("median %+.4f (bar |median|<%.2f) and max|null| %.4f = %.3fx the real %.4f "
                   "(bar <%.2fx)" % (med, margin, mx, ratio, best_med, NULL_SIGNAL_RATIO))
        if (ok_med and ok_mag) else
        ("null FAILED: median %+.4f, max|null| %.4f = %.3fx the real %.4f"
         % (med, mx, ratio, best_med)),
    }


GATE_FUNCS = {
    "split_equality": gate_split_equality,
    "geometry": gate_geometry,
    "head_validity": gate_head_validity,
    "permutation_null": gate_permutation_null,
}


# ---------------------------------------------------------------------------------------------------
# scoring
# ---------------------------------------------------------------------------------------------------
def score_head(ctx: dict, head_name: str, reduction: str) -> dict:
    """DIRECT Spearman per seed, sign AS GIVEN. No fitting of any kind."""
    from dna_decode.glm.genomewide import spearman

    col = ctx["layout"]["%s::%s" % (head_name, reduction)]
    pred = np.asarray(ctx["arr"][:, col], dtype=float)
    per_seed, gains, wins = {}, [], 0
    for s in ctx["seeds"]:
        _tr, te, _ = ctx["splits"][s]
        rho = round(float(spearman(list(pred[te]), list(ctx["ys"](te)))), 4)
        want = float(ctx["committed_atlas"][str(s)])
        per_seed[str(s)] = {"spearman": rho, "atlas_committed": want,
                            "gain": round(rho - want, 4)}
        gains.append(rho - want)
        wins += int(rho > want)
    boot = paired_bootstrap(gains, percentiles=ctx["pct"], n_boot=ctx["n_boot"])
    return {"head": head_name, "reduction": reduction,
            "median_spearman": median_of(v["spearman"] for v in per_seed.values()),
            "per_seed": per_seed, "seed_wins_over_atlas": wins,
            "n_seeds": len(ctx["seeds"]), "gain_bootstrap": boot,
            "eligible_as_the_arms_result": reduction == PRIMARY_REDUCTION}


def verdict_from_artifact(rec: dict, bar: dict) -> tuple[str, str]:
    """Re-derive the verdict from the emitted numbers plus the frozen bar, and nothing else."""
    gates = rec.get("gates", {})
    margin = bar["inherited_not_chosen"]["margin"]["value"]
    for gate, outcome in GATE_SEQUENCE:
        g = gates.get(gate)
        if g is None:
            return (outcome, "gate %r was never evaluated" % gate)
        if not g.get("passed"):
            return (outcome, g.get("detail", "gate %r failed" % gate))
    best = rec.get("best_primary")
    if not best:
        return ("HEADS_TIES", "no primary rung was scored")
    b = best["gain_bootstrap"]
    return classify(b["ci_lo"], b["ci_hi"], margin=margin)


def run(*, bar: dict, block_dir: Path, limit: int | None = None, cohort=None,
        verbose: bool = True, job: str | None = None) -> dict:
    t0 = time.time()
    inh = bar["inherited_not_chosen"]
    tiles, st = cohort if cohort is not None else load_arm_cohort()
    length = len(tiles[0].seq)

    base_sd = AH.load_base_state_dict()
    heads = AH.load_heads(Path(DEFAULT_PA_REPO) / "trained_weights/expression_prediction",
                          base_state_dict=base_sd)
    if verbose:
        print("EXTRACT  %d tiles x %d heads, route=%s" % (len(tiles), len(heads), ROUTE_TRUE),
              flush=True)
    hb = None
    if job:
        from scripts.detached_launch import Heartbeat

        n_target = len(tiles) if limit is None else min(limit, len(tiles))
        hb = Heartbeat(job, units_total=n_target, label="tiles extracted (%d heads, %s route)"
                                                        % (len(heads), ROUTE_TRUE))
    extract(tiles=tiles, heads=heads, block_dir=block_dir, limit=limit, verbose=verbose,
            heartbeat=hb)
    arr, side = load_scalars(block_dir, allow_partial=False)

    seeds = list(range(inh["split"]["n_seeds"]))
    frac = inh["split"]["frac"]
    splits = {s: leave_peak_out_indices(tiles, frac=frac, seed=s) for s in seeds}
    ctx = {
        "tiles": tiles, "heads": heads, "arr": arr, "layout": side["layout"],
        "seeds": seeds, "splits": splits, "bar": bar, "length": length,
        "committed_atlas": inh["comparator"]["atlas_per_seed"],
        "committed_n_test": inh["split"]["n_test_per_seed"],
        "pct": bar["multiplicity"]["ci_percentiles"], "n_boot": N_BOOT,
        "ys": lambda idx: [float(tiles[i].expression) for i in idx],
    }

    rec: dict = {
        "schema": "expression-heads-result-v1", "date": "2026-11-02",
        "bar": "wiki/expression_heads_acceptance_bar.json", "arm": ARM,
        "criterion": bar["criterion"], "route": ROUTE_TRUE,
        "scoring": bar["scoring"], "block": side,
        "substrate": {"n_tiles": len(tiles), "bp": length, "path": st.path,
                      "n_peaks": len({t.peak for t in tiles})},
        "free_ram_bytes_at_start": free_ram_bytes(), "gates": {},
    }

    # primaries first, because the null gate needs the arm's real signal to scale against
    primaries = {h.name: score_head(ctx, h.name, PRIMARY_REDUCTION) for h in heads}
    best_name = max(primaries, key=lambda k: primaries[k]["median_spearman"])
    ctx["best_head"] = best_name
    ctx["best_median"] = primaries[best_name]["median_spearman"]

    for gate_name, outcome in GATE_SEQUENCE:
        g = GATE_FUNCS[gate_name](ctx)
        rec["gates"][gate_name] = g
        if not g["passed"]:
            rec["verdict"], rec["verdict_detail"] = outcome, g["detail"]
            rec["wall_seconds"] = round(time.time() - t0, 1)
            return rec

    rec["primaries"] = primaries
    rec["secondaries"] = {h.name: score_head(ctx, h.name, SECONDARY_REDUCTION) for h in heads}
    rec["best_primary"] = primaries[best_name]
    rec["best_head"] = best_name
    rec["metric_rule"] = bar["metric_rule"]
    rec["arm_metric"] = primaries[best_name]["median_spearman"]
    rec["comparator"] = {"atlas_median": inh["comparator"]["atlas_median"],
                         "atlas_per_seed": ctx["committed_atlas"], "recomputed_in_run": False}
    rec["margin"] = inh["margin"]["value"]
    rec["multiplicity"] = bar["multiplicity"]
    rec["honest_limits"] = bar["honest_limits_fixed_in_advance"]
    rec["forecast_recorded_in_bar"] = bar["forecast"]
    rec["secondary_is_ineligible"] = True

    verdict, detail = verdict_from_artifact(rec, bar)
    rec["verdict"], rec["verdict_detail"] = verdict, detail
    rec["verdict_is_rederivable_from_this_artifact"] = True
    rec["wall_seconds"] = round(time.time() - t0, 1)
    return rec


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--block-dir", default=str(BLOCK_DIR))
    ap.add_argument("--limit", type=int, default=None,
                    help="cap the tile count; stamps the block PARTIAL, which the scorer refuses")
    ap.add_argument("--out", default=str(OUT))
    ap.add_argument("--job", default="expression_heads",
                    help="heartbeat job name read by scripts/detached_launch.py status; "
                         "pass empty to disable")
    args = ap.parse_args(argv)

    if not BAR_PATH.is_file():
        print("REFUSED: %s absent; freeze the bar before scoring." % BAR_PATH, flush=True)
        return 2
    bar = json.loads(BAR_PATH.read_text(encoding="utf-8"))

    try:
        rec = run(bar=bar, block_dir=Path(args.block_dir), limit=args.limit,
                  job=args.job or None)
    except ScoreError as e:
        print("REFUSED: %s" % e, flush=True)
        return 2
    except Exception as e:                                     # noqa: BLE001
        print("SCORING FAILED (%s): %s" % (type(e).__name__, e), flush=True)
        return 1

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(rec, indent=2, sort_keys=True), encoding="utf-8")

    for name, g in rec["gates"].items():
        print("GATE %-20s %s  %s" % (name, "PASS" if g["passed"] else "FAIL", g["detail"]))
    for name, row in sorted(rec.get("primaries", {}).items()):
        b = row["gain_bootstrap"]
        print("HEAD %-38s median %.4f  gain %+.4f  CI [%+.4f, %+.4f]  wins %d/%d"
              % (name, row["median_spearman"], b["median"], b["ci_lo"], b["ci_hi"],
                 row["seed_wins_over_atlas"], row["n_seeds"]))
    print()
    print("VERDICT  %s" % rec["verdict"])
    print("         %s" % rec["verdict_detail"])
    print("wall     %.1f s" % rec["wall_seconds"])
    print("-> %s" % out)
    return 0 if rec["verdict"] in OUTCOMES else 1


if __name__ == "__main__":
    raise SystemExit(main())
