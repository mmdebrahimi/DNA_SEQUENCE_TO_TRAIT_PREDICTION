"""The long-horizon fine-tune run for the `finetuned-promoteratlas-encoder-150bp` arm.

ONE RUN, CARRIED LONG, RECORDING EVERY EPOCH. The frozen budget's verdict is a PREFIX of this
run rather than a separate run: same seed, same initialisation, same data order, same
optimiser trajectory. The frozen early-stopping rule is then applied by READING the recorded
curve. That is what makes the longer horizons free instead of a sevenfold cost -- three runs at
E, 2E and 4E epochs would spend 7E epochs per seed for strictly less information, because they
would not share a trajectory.

THE RESUME CONTRACT IS BY CONSTRUCTION, NOT BY BOOKKEEPING. Each epoch's shuffle comes from a
generator seeded by `(seed, epoch)` -- never from the accumulated global RNG. The consequence is
that epoch 20's data order does not depend on whether epochs 0-19 ran in this process or in one
that died last Tuesday, so a resumed run reproduces an uninterrupted one exactly. Threading the
global RNG through a checkpoint would also work on a good day and would silently diverge the
moment anything else drew a random number between epochs.

A DEATH COSTS ONE EPOCH. The per-seed checkpoint carries encoder, head, optimiser state and the
curve so far, and is rewritten after every epoch. Unfinished seeds are retrained; finished seeds
are read from their committed curve and never recomputed. This host has six recorded long-run
deaths, two of them unexplained, so a ten-day run that could not resume would be a ten-day run
that probably produces nothing.

THE WALL CLOCK IS HARD, AND A BOUNDED RUN SAYS SO. Overrunning the window leaves every seed
half-done, which is worse than fewer seeds finished cleanly. When the deadline arrives the run
stops at the next epoch boundary and stamps itself BOUNDED_BY_WALL_CLOCK -- never COMPLETE. A run
that stopped early must not be readable as a finished one.

THE TARGET TRANSFORM IS NOT A FREE PARAMETER. Training is on `log1p(expression)` and validation
Spearman is scored against RAW expression, matching `atlas_embed.ridge_spearman`, whose docstring
states any deviation makes the number incomparable with the committed incumbent tables. An
earlier version of the probe's trainer used MSE on RAW expression, which would have confounded
the OBJECTIVE TRANSFORM with the REPRESENTATION and introduced a THIRD reading on top of the
budget-versus-representation pair the bar names. `_targets` therefore REFUSES any other
transform rather than accepting it as an option.

THE NULL IS BUDGETED INTO THE CEILING, NOT ADDED ON TOP. The shuffled-label null is itself a
fine-tune, so it costs epochs. `plan_run` reserves them before sizing the horizon; a plan that
sized the horizon first and then discovered it had no room for the null would either overrun the
window or ship without the null, and the bar requires the null.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import time
from pathlib import Path

import numpy as np

from scripts.freeze_finetune_arm_bar import ARM
from scripts.positional_arm_probe import load_arm_cohort
from scripts.positional_arm_score import leave_peak_out_indices

REPO = Path(__file__).resolve().parents[1]
BAR_PATH = REPO / "wiki/finetune_arm_acceptance_bar.json"
OUT = REPO / "wiki/finetune_arm_run_2026-11-13.json"
CKPT_DIR = REPO / "data/promoteratlas/finetune_checkpoints"

TARGET_TRANSFORM = "log1p"
"""The ONLY accepted transform. Not a default -- a constraint. See the module note."""

N_NULL_SEEDS = 3
"""Shuffled-label null replicates, at the FROZEN budget each.

Three rather than ten because each one is a full fine-tune and they come out of the same window
the real seeds do. Three is enough to show the fit collapses when the labels are destroyed, which
is the claim; it is not enough to estimate a null DISTRIBUTION, and the artifact says so.
"""

STATUS_COMPLETE = "COMPLETE"
STATUS_BOUNDED = "BOUNDED_BY_WALL_CLOCK"


class RunRefused(ValueError):
    """A run configuration that cannot mean what the caller intends."""


# --------------------------------------------------------------------------- planning


def plan_run(*, ceiling_hours: float, seconds_per_epoch: float, n_seeds: int,
             frozen_epochs: int, n_null_seeds: int = N_NULL_SEEDS) -> dict:
    """Epoch budget for the whole run, with the null's cost RESERVED FIRST.

    The ordering matters: sizing the horizon against the full ceiling and then adding the null
    either overruns the window or drops the null, and the frozen bar requires the null.
    """
    if ceiling_hours <= 0 or seconds_per_epoch <= 0 or n_seeds <= 0 or frozen_epochs <= 0:
        raise RunRefused("plan_run needs positive ceiling, per-epoch cost, seeds and budget")

    total_epochs = int((ceiling_hours * 3600.0) // seconds_per_epoch)
    null_epochs = n_null_seeds * frozen_epochs
    usable = total_epochs - null_epochs
    horizon = int(usable // n_seeds)

    if horizon < frozen_epochs:
        raise RunRefused(
            "after reserving %d epochs for %d null replicates, the %.2f h ceiling leaves %d "
            "epochs per seed across %d seeds -- fewer than the frozen budget of %d. The arm's own "
            "verdict would not fit, so this run must not start."
            % (null_epochs, n_null_seeds, ceiling_hours, horizon, n_seeds, frozen_epochs))

    projected = (n_seeds * horizon + null_epochs) * seconds_per_epoch / 3600.0
    return {
        "horizon_epochs": horizon,
        "frozen_budget_epochs": frozen_epochs,
        "multiple_of_frozen_budget": round(horizon / frozen_epochs, 1),
        "n_seeds": n_seeds,
        "n_null_seeds": n_null_seeds,
        "null_epochs_reserved": null_epochs,
        "total_epochs": n_seeds * horizon + null_epochs,
        "seconds_per_epoch": round(seconds_per_epoch, 1),
        "projected_hours": round(projected, 1),
        "ceiling_hours": ceiling_hours,
        "fits": projected <= ceiling_hours,
        "why_the_null_is_reserved_first": (
            "the shuffled-label null is itself a fine-tune and costs %d epochs. Sizing the "
            "horizon against the whole ceiling and then adding it would overrun the window or "
            "ship without the null, and the frozen bar requires the null." % null_epochs),
    }


def epoch_generator_seed(seed: int, epoch: int) -> int:
    """A deterministic per-(seed, epoch) shuffle seed.

    Explicit arithmetic rather than `hash`, so the value does not depend on interpreter state,
    and so epoch 20's data order is identical whether epochs 0-19 ran here or in a process that
    died days ago. This is the whole of the resume contract.
    """
    return (int(seed) * 1_000_003 + int(epoch) * 7_919 + 90_001) % (2 ** 31 - 1)


# --------------------------------------------------------------------------- targets


def _targets(values, transform: str = TARGET_TRANSFORM):
    """Training targets. REFUSES any transform but log1p.

    Not an option with a default: training on RAW expression would confound the objective
    transform with the representation, and the bar's whole budget-versus-representation
    distinction assumes that confound is absent. The refusal is what stops that regressing
    silently.
    """
    import torch

    if transform != TARGET_TRANSFORM:
        raise RunRefused(
            "target transform %r is refused; this arm trains on %s(expression) and scores "
            "Spearman against RAW expression, matching atlas_embed.ridge_spearman, whose "
            "docstring states any deviation makes the number incomparable with the committed "
            "incumbent tables. Training on raw expression would confound the OBJECTIVE TRANSFORM "
            "with the REPRESENTATION." % (transform, TARGET_TRANSFORM))
    return torch.log1p(torch.tensor([float(v) for v in values]))


# --------------------------------------------------------------------------- the frozen rule


def frozen_budget_value(curve: list[float], *, frozen_epochs: int, patience: int) -> dict:
    """Apply the frozen early-stopping rule to the PREFIX of a long curve.

    This is the function the design rests on. It reads ONLY `curve[:frozen_epochs]`, so a
    spectacular epoch beyond the budget cannot reach the verdict -- which is exactly what makes
    the longer horizons ineligible diagnostics rather than a bigger arm.
    """
    if frozen_epochs <= 0:
        raise RunRefused("frozen_epochs must be positive")
    prefix = list(curve[:frozen_epochs])
    if not prefix:
        return {"value": None, "best_epoch": None, "n_epochs_read": 0,
                "stopped_early": False,
                "detail": "no epoch inside the frozen budget was recorded"}

    best_i, best = 0, prefix[0]
    for i, v in enumerate(prefix):
        if v > best:
            best_i, best = i, v
    stopped = (len(prefix) - 1 - best_i) >= patience
    return {
        "value": best,
        "best_epoch": best_i,
        "n_epochs_read": len(prefix),
        "n_epochs_available": len(curve),
        "stopped_early": bool(stopped),
        "detail": ("best of the first %d epochs is %.4f at epoch %d; %d epochs beyond the budget "
                   "were recorded and are INELIGIBLE"
                   % (len(prefix), best, best_i, max(0, len(curve) - frozen_epochs))),
    }


def diagnostic_tail(curve: list[float], *, frozen_epochs: int) -> dict:
    """The beyond-budget epochs, explicitly marked ineligible.

    Reported so the budget-limited reading can be closed empirically, and labelled so it can
    never become the arm's result.
    """
    tail = list(curve[frozen_epochs:])
    best_tail = max(tail) if tail else None
    inside = frozen_budget_value(curve, frozen_epochs=frozen_epochs, patience=10 ** 9)["value"]
    return {
        "eligible_as_the_arms_result": False,
        "status": "INELIGIBLE diagnostic, pre-declared in the frozen bar",
        "n_epochs": len(tail),
        "best_beyond_budget": best_tail,
        "best_inside_budget": inside,
        "improvement_beyond_budget": (None if (best_tail is None or inside is None)
                                      else round(best_tail - inside, 4)),
        "purpose": ("answering 'would more budget have helped?' and nothing else. A curve flat "
                    "from the frozen budget out to many times it closes the budget-limited "
                    "reading by MEASUREMENT instead of by inference from a subsample curve."),
    }


# --------------------------------------------------------------------------- checkpoints


def _ckpt_path(ckpt_dir: Path, seed: int, *, kind: str = "seed") -> Path:
    return Path(ckpt_dir) / ("%s_%d.pt" % (kind, seed))


def _curve_path(ckpt_dir: Path, seed: int, *, kind: str = "seed") -> Path:
    return Path(ckpt_dir) / ("%s_%d_curve.json" % (kind, seed))


def read_curve(ckpt_dir: Path, seed: int, *, kind: str = "seed") -> dict | None:
    """A finished-or-partial curve record, or None.

    A SEPARATE small JSON beside the heavy checkpoint, so the result artifact can be assembled
    -- and progress inspected -- without loading torch or a 20 MB tensor file.
    """
    p = _curve_path(ckpt_dir, seed, kind=kind)
    if not p.is_file():
        return None
    try:
        rec = json.loads(p.read_text(encoding="utf-8"))
        return rec if isinstance(rec, dict) else None
    except json.JSONDecodeError:
        return None


def _write_curve(ckpt_dir: Path, seed: int, rec: dict, *, kind: str = "seed") -> None:
    p = _curve_path(ckpt_dir, seed, kind=kind)
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(rec, indent=1, sort_keys=True), encoding="utf-8")
    os.replace(tmp, p)


# --------------------------------------------------------------------------- one seed


def train_one_seed(tiles, *, seed: int, horizon_epochs: int, lr: float, batch_size: int,
                   split_frac: float, ckpt_dir: Path, build, deadline: float | None = None,
                   shuffle_labels: bool = False, kind: str = "seed",
                   transform: str = TARGET_TRANSFORM, heartbeat=None,
                   verbose: bool = True) -> dict:
    """Train one seed to the horizon, recording validation Spearman at EVERY epoch.

    Resumes from the per-seed checkpoint when one is present. Returns the curve record.
    """
    import torch

    from dna_decode.glm.genomewide import spearman

    tr, va, _rest = leave_peak_out_indices(tiles, frac=split_frac, seed=seed)
    length = len(tiles[0].seq)

    from scripts.finetune_arm_probe import _onehot

    Xtr = _onehot([tiles[i].seq for i in tr], length)
    raw_tr = [float(tiles[i].expression) for i in tr]
    if shuffle_labels:
        # The null: TRAINING labels permuted, everything else -- split, order, hyperparameters --
        # untouched. Unlike the heads arm's one-feature predictor, this one genuinely fits, so
        # destroying the labels destroys the fit. Seeded off the seed so the null is reproducible.
        rng = np.random.default_rng(770_000 + seed)
        raw_tr = list(np.array(raw_tr)[rng.permutation(len(raw_tr))])
    ytr = _targets(raw_tr, transform)
    Xva = _onehot([tiles[i].seq for i in va], length)
    yva = [float(tiles[i].expression) for i in va]          # RAW, deliberately

    torch.manual_seed(seed)                                  # model init, before any draw
    enc, head = build()
    opt = torch.optim.Adam(list(enc.parameters()) + list(head.parameters()), lr=lr)

    curve: list[float] = []
    untrained = None
    start_epoch = 0
    ck = _ckpt_path(ckpt_dir, seed, kind=kind)
    if ck.is_file():
        state = torch.load(ck, map_location="cpu", weights_only=False)
        enc.load_state_dict(state["enc"])
        head.load_state_dict(state["head"])
        opt.load_state_dict(state["opt"])
        curve = list(state["curve"])
        untrained = state.get("untrained_val_spearman")
        start_epoch = int(state["epoch"]) + 1
        if verbose:
            print("  seed %d: resuming at epoch %d (%d recorded)"
                  % (seed, start_epoch, len(curve)), flush=True)
    else:
        # THE UNTRAINED BASELINE, measured before any gradient step. One forward pass.
        #
        # It exists because a high shuffled-label null has TWO possible causes that look
        # identical in the number: the fine-tune learned something the labels did not justify, or
        # a RANDOM READOUT of an already-informative representation retains signal on its own.
        # Only the second is benign, and only this baseline distinguishes them. Found by a
        # synthetic fixture whose untrained model already scored 0.76 because mean-pooling made
        # its representation essentially base composition and its target WAS base composition --
        # so a null that failed to collapse there said nothing about the fine-tune at all.
        enc.eval()
        with torch.no_grad():
            p0 = head(enc(Xva, return_latent=True)[2].mean(dim=2)).squeeze(-1).numpy()
        if np.isfinite(p0).all():
            untrained = round(float(spearman(list(p0), yva)), 4)

    t0 = time.time()
    bounded = False
    for epoch in range(start_epoch, horizon_epochs):
        if deadline is not None and time.time() >= deadline:
            bounded = True
            break
        enc.train()
        # The resume contract: the permutation is a function of (seed, epoch) alone.
        g = torch.Generator().manual_seed(epoch_generator_seed(seed, epoch))
        perm = torch.randperm(len(tr), generator=g)
        for b0 in range(0, len(tr), batch_size):
            sel = perm[b0:b0 + batch_size]
            opt.zero_grad()
            z = enc(Xtr[sel], return_latent=True)[2]
            loss = ((head(z.mean(dim=2)).squeeze(-1) - ytr[sel]) ** 2).mean()
            loss.backward()
            opt.step()
        enc.eval()
        with torch.no_grad():
            pv = head(enc(Xva, return_latent=True)[2].mean(dim=2)).squeeze(-1).numpy()
        if not np.isfinite(pv).all():
            raise RunRefused(
                "seed %d epoch %d produced non-finite validation predictions. Refusing: the "
                "scorer maps NaN to 0.0, so divergence would be recorded as 'no signal'."
                % (seed, epoch))
        curve.append(round(float(spearman(list(pv), yva)), 4))

        ck.parent.mkdir(parents=True, exist_ok=True)
        tmp = ck.with_suffix(".pt.tmp")
        torch.save({"enc": enc.state_dict(), "head": head.state_dict(),
                    "opt": opt.state_dict(), "curve": curve, "epoch": epoch,
                    "seed": seed, "kind": kind,
                    "untrained_val_spearman": untrained}, tmp)
        os.replace(tmp, ck)

        rec = {"seed": seed, "kind": kind, "n_train": len(tr), "n_val": len(va),
               "n_peaks_train": len({tiles[i].peak for i in tr}),
               "n_peaks_val": len({tiles[i].peak for i in va}),
               "peaks_shared_between_splits": len({tiles[i].peak for i in tr}
                                                  & {tiles[i].peak for i in va}),
               "val_spearman_by_epoch": curve,
               "untrained_val_spearman": untrained,
               "epochs_done": len(curve), "horizon_epochs": horizon_epochs,
               "shuffled_labels": shuffle_labels,
               "target_transform": transform,
               "complete": len(curve) >= horizon_epochs,
               "wall_seconds_this_process": round(time.time() - t0, 1)}
        _write_curve(ckpt_dir, seed, rec, kind=kind)
        if heartbeat is not None:
            heartbeat.beat(heartbeat_units(seed, len(curve), horizon_epochs),
                           force=(epoch == start_epoch), seed=seed, epoch=len(curve),
                           kind=kind, last_val_spearman=curve[-1])
        if verbose:
            print("  seed %d epoch %3d/%d  val rho %+.4f" % (seed, len(curve), horizon_epochs,
                                                             curve[-1]), flush=True)

    out = read_curve(ckpt_dir, seed, kind=kind) or {
        "seed": seed, "kind": kind, "val_spearman_by_epoch": curve,
        "epochs_done": len(curve), "horizon_epochs": horizon_epochs,
        "shuffled_labels": shuffle_labels, "target_transform": transform,
        "complete": len(curve) >= horizon_epochs}
    out["bounded_by_wall_clock"] = bool(bounded)
    return out


def heartbeat_units(seed: int, epochs_done: int, horizon: int) -> int:
    """Total epochs completed across the run so far, for the progress fraction."""
    return seed * horizon + epochs_done


# --------------------------------------------------------------------------- the run


def run(*, bar: dict, tiles, ckpt_dir: Path, build=None, wall_ceiling_hours: float | None = None,
        n_seeds: int | None = None, horizon_epochs: int | None = None,
        n_null_seeds: int = N_NULL_SEEDS, heartbeat=None, verbose: bool = True,
        now: float | None = None) -> dict:
    """The whole arm: every real seed to the horizon, then the nulls at the frozen budget."""
    t0 = now if now is not None else time.time()
    hp = bar["frozen_hyperparameters"]
    inh = bar["inherited_not_chosen"]
    frozen_epochs = int(hp["n_epochs_budget"])
    patience = int(hp["patience"])
    lr = float(hp["learning_rate"])
    batch_size = int(hp["batch_size"])
    split_frac = float(inh["comparator"]["split_frac"])
    n_seeds = int(n_seeds if n_seeds is not None else inh["comparator"]["n_seeds"])

    ceiling = float(wall_ceiling_hours if wall_ceiling_hours is not None
                    else bar["ceiling"]["ceiling_hours"])
    spe = float(bar["horizon"].get("minutes_per_epoch") or 0.0) * 60.0
    if horizon_epochs is None:
        if spe <= 0:
            raise RunRefused("the bar records no per-epoch cost, so no horizon can be planned")
        plan = plan_run(ceiling_hours=ceiling, seconds_per_epoch=spe, n_seeds=n_seeds,
                        frozen_epochs=frozen_epochs, n_null_seeds=n_null_seeds)
        horizon_epochs = plan["horizon_epochs"]
    else:
        plan = {"horizon_epochs": int(horizon_epochs), "frozen_budget_epochs": frozen_epochs,
                "n_seeds": n_seeds, "n_null_seeds": n_null_seeds,
                "note": "horizon supplied by the caller, not planned from the ceiling"}

    if build is None:
        from scripts.finetune_arm_probe import _build_trainable

        build = _build_trainable

    deadline = t0 + ceiling * 3600.0
    seeds_rec, nulls_rec, bounded = {}, {}, False

    for s in range(n_seeds):
        prior = read_curve(ckpt_dir, s)
        if prior and prior.get("complete") and \
                len(prior.get("val_spearman_by_epoch", [])) >= horizon_epochs:
            # A finished seed is READ, never recomputed -- the resume property.
            if verbose:
                print("  seed %d: already complete (%d epochs), reusing"
                      % (s, len(prior["val_spearman_by_epoch"])), flush=True)
            seeds_rec[str(s)] = prior
            continue
        if time.time() >= deadline:
            bounded = True
            break
        rec = train_one_seed(tiles, seed=s, horizon_epochs=horizon_epochs, lr=lr,
                             batch_size=batch_size, split_frac=split_frac, ckpt_dir=ckpt_dir,
                             build=build, deadline=deadline, heartbeat=heartbeat,
                             verbose=verbose)
        seeds_rec[str(s)] = rec
        bounded = bounded or bool(rec.get("bounded_by_wall_clock"))

    if not bounded:
        for s in range(n_null_seeds):
            prior = read_curve(ckpt_dir, s, kind="null")
            if prior and len(prior.get("val_spearman_by_epoch", [])) >= frozen_epochs:
                nulls_rec[str(s)] = prior
                continue
            if time.time() >= deadline:
                bounded = True
                break
            rec = train_one_seed(tiles, seed=s, horizon_epochs=frozen_epochs, lr=lr,
                                 batch_size=batch_size, split_frac=split_frac,
                                 ckpt_dir=ckpt_dir, build=build, deadline=deadline,
                                 shuffle_labels=True, kind="null", verbose=verbose)
            nulls_rec[str(s)] = rec
            bounded = bounded or bool(rec.get("bounded_by_wall_clock"))

    per_seed = {}
    for s, rec in seeds_rec.items():
        curve = rec["val_spearman_by_epoch"]
        per_seed[s] = {
            "val_spearman_by_epoch": curve,
            "frozen_budget": frozen_budget_value(curve, frozen_epochs=frozen_epochs,
                                                 patience=patience),
            "diagnostic_tail": diagnostic_tail(curve, frozen_epochs=frozen_epochs),
            "untrained_val_spearman": rec.get("untrained_val_spearman"),
            "n_train": rec.get("n_train"), "n_val": rec.get("n_val"),
            "peaks_shared_between_splits": rec.get("peaks_shared_between_splits"),
        }

    null_values, null_untrained = [], []
    for s, rec in nulls_rec.items():
        v = frozen_budget_value(rec["val_spearman_by_epoch"], frozen_epochs=frozen_epochs,
                                patience=patience)["value"]
        if v is not None:
            null_values.append(round(float(v), 4))
        if rec.get("untrained_val_spearman") is not None:
            null_untrained.append(rec["untrained_val_spearman"])

    complete = (not bounded
                and len(seeds_rec) == n_seeds
                and all(len(r["val_spearman_by_epoch"]) >= horizon_epochs
                        for r in seeds_rec.values())
                and len(nulls_rec) == n_null_seeds)

    return {
        "schema": "finetune-arm-run-v1",
        "arm": ARM,
        "bar": "wiki/finetune_arm_acceptance_bar.json",
        "status": STATUS_COMPLETE if complete else STATUS_BOUNDED,
        "why_status": ("every seed reached the horizon and every null replicate reached the "
                       "frozen budget" if complete else
                       "the run stopped at an epoch boundary before finishing; a bounded run is "
                       "NOT readable as a complete one"),
        "plan": plan,
        "frozen_budget_epochs": frozen_epochs,
        "horizon_epochs": horizon_epochs,
        "target_transform": TARGET_TRANSFORM,
        "per_seed": per_seed,
        "n_seeds_done": len(seeds_rec),
        "n_seeds_requested": n_seeds,
        "null": {
            "kind": "shuffled_training_labels",
            "n_replicates": len(null_values),
            "n_requested": n_null_seeds,
            "values_at_frozen_budget": null_values,
            "median": (round(float(np.median(null_values)), 4) if null_values else None),
            "max_abs": (round(float(np.max(np.abs(null_values))), 4) if null_values else None),
            "untrained_baseline_per_replicate": null_untrained,
            "how_to_read_a_HIGH_null": (
                "a high null has TWO causes that look identical in the number: the fine-tune "
                "learned something the labels did not justify, or a RANDOM READOUT of an "
                "already-informative representation retains signal on its own. Compare against "
                "untrained_baseline_per_replicate -- a null that merely matches the untrained "
                "baseline is the second, benign case and says nothing about the fine-tune. The "
                "frozen bar's null gate does not make this distinction, so it would return "
                "INDETERMINATE_NULL_FAIL either way, which is the correct pre-committed outcome; "
                "this field is what makes the failure DIAGNOSABLE rather than ambiguous."),
            "scope_limit": ("%d replicates is enough to show the fit collapses when the labels "
                            "are destroyed, which is the claim. It is NOT enough to estimate a "
                            "null distribution." % len(null_values)),
        },
        "ineligible_diagnostics_note": (
            "every epoch beyond %d is an INELIGIBLE diagnostic, pre-declared in the frozen bar. "
            "frozen_budget_value reads only the prefix, so a spectacular late epoch cannot reach "
            "the verdict." % frozen_epochs),
        "wall_seconds": round(time.time() - t0, 1),
        "wall_ceiling_hours": ceiling,
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ckpt-dir", default=str(CKPT_DIR))
    ap.add_argument("--out", default=str(OUT))
    ap.add_argument("--job", default="finetune_arm",
                    help="heartbeat job name read by scripts/detached_launch.py status")
    ap.add_argument("--seeds", type=int, default=None)
    ap.add_argument("--horizon", type=int, default=None,
                    help="override the planned horizon (diagnostic use)")
    ap.add_argument("--wall-ceiling-hours", type=float, default=None)
    args = ap.parse_args(argv)

    if not BAR_PATH.is_file():
        print("REFUSED: %s absent; freeze the bar before running." % BAR_PATH, flush=True)
        return 2
    bar = json.loads(BAR_PATH.read_text(encoding="utf-8"))
    if bar["route"]["chosen"] != "local_cpu":
        print("REFUSED: the frozen bar routes this arm to %r, not local_cpu."
              % bar["route"]["chosen"], flush=True)
        return 2

    tiles, st = load_arm_cohort()
    print("COHORT   %d tiles x %d bp from %s" % (len(tiles), len(tiles[0].seq), st.path),
          flush=True)

    hb = None
    if args.job:
        from scripts.detached_launch import Heartbeat

        n_seeds = int(args.seeds if args.seeds is not None
                      else bar["inherited_not_chosen"]["comparator"]["n_seeds"])
        horizon = args.horizon or int(bar["horizon"]["diagnostic_horizon_epochs"] or 1)
        hb = Heartbeat(args.job, units_total=max(1, n_seeds * horizon),
                       label="epochs completed (%d seeds x %d horizon)" % (n_seeds, horizon))

    try:
        rec = run(bar=bar, tiles=tiles, ckpt_dir=Path(args.ckpt_dir), n_seeds=args.seeds,
                  horizon_epochs=args.horizon, wall_ceiling_hours=args.wall_ceiling_hours,
                  heartbeat=hb)
    except RunRefused as e:
        print("REFUSED: %s" % e, flush=True)
        return 2
    except Exception as e:                                     # noqa: BLE001
        print("RUN FAILED (%s): %s" % (type(e).__name__, e), flush=True)
        return 1

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(rec, indent=2, sort_keys=True), encoding="utf-8")

    print()
    print("STATUS   %s" % rec["status"])
    print("         %s" % rec["why_status"])
    for s in sorted(rec["per_seed"], key=int):
        fb = rec["per_seed"][s]["frozen_budget"]
        dt = rec["per_seed"][s]["diagnostic_tail"]
        print("SEED %-3s frozen %.4f @epoch %s   tail best %s (+%s, INELIGIBLE)"
              % (s, fb["value"] if fb["value"] is not None else float("nan"), fb["best_epoch"],
                 dt["best_beyond_budget"], dt["improvement_beyond_budget"]))
    print("NULL     %s  median %s  max|null| %s"
          % (rec["null"]["values_at_frozen_budget"], rec["null"]["median"],
             rec["null"]["max_abs"]))
    print("wall     %.1f s" % rec["wall_seconds"])
    print("-> %s" % out)
    return 0 if rec["status"] == STATUS_COMPLETE else 1


if __name__ == "__main__":
    raise SystemExit(main())
