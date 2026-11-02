"""Feasibility and compute-route probe for the fine-tuned-encoder arm.

This is the only arm whose cost must be measured BEFORE it is allowed to start, and the reason is
not budgeting -- it is interpretability. A fine-tune that fails to beat the frozen encoder's 0.5076
has two completely different readings:

* the frozen REPRESENTATION was not the limit, so a trainable encoder buys nothing here; or
* the SUBSTRATE or the TRAINING BUDGET was the limit, and the arm simply never converged.

Those look identical in the result artifact and mean opposite things. The only thing that separates
them is an observed convergence curve, which is why this probe measures one and why the arm should
not be run without it.

**Why a fine-tune is tractable here at all.** The encoder is 1,847,428 parameters by
`model.parameters()` -- not a 1.2B language model -- so this is NOT comparable to this project's
blocked large-LM fine-tunes. That is a reason to try it, not a promise that it works.

**Hyperparameters are ANCHORED to the checkpoint's own pretraining recipe, not chosen freely.** The
base checkpoint carries a `config` key recording `learning_rate`, `patience`, `batch_size` and
`n_epochs` as actually used to train it. A fine-tune whose recipe is unrelated to the pretraining
recipe is a weaker test of the same question, so this probe reads those values and records any
deviation WITH its justification rather than silently substituting its own.

**The route is chosen from the measurement, in both directions.** Local is CPU-only here
(`torch.cuda.is_available()` is False, 4 threads). The project's established GPU route is a Kaggle
T4 kernel. The probe projects the local cost and compares it against a NAMED wall-clock ceiling: a
run that would exceed it routes to Kaggle rather than being started and abandoned halfway, which is
the failure mode that produces a `PARTIAL` artifact nobody can read.

**What the convergence curve here does and does NOT establish.** Measuring convergence on the full
cohort is itself the expensive run this probe exists to size, so the curve is observed on a
SUBSAMPLE whose size is derived from the measured throughput to fit a bounded budget. That bounds
the SHAPE of convergence -- whether the validation metric turns, and roughly when -- and it does NOT
establish the full-cohort epoch count. Stated here rather than implied, because treating a subsample
curve as a full-cohort budget is exactly the over-read this probe is meant to prevent.

Read-only with respect to the repo's data. Writes one artifact. Exits non-zero on any refusal.
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import time
from pathlib import Path

import numpy as np

from dna_decode.glm.atlas_embed import EXPECTED_PARAMS
from scripts.expression_heads_probe import EXPECTED_BYTES, EXPECTED_SHA256, ProbeRefused
from scripts.positional_arm_probe import (
    extrapolate_wall_seconds,
    free_ram_bytes,
    load_arm_cohort,
    sha256_of,
    stride_sample,
)
from scripts.positional_arm_score import leave_peak_out_indices

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "wiki/finetune_arm_feasibility_2026-11-02.json"

ARM = "finetuned-promoteratlas-encoder-150bp"

#: The arm's comparator and split, inherited verbatim from the committed reconcile artifact.
N_SEEDS = 10
SPLIT_FRAC = 0.25
FROZEN_INCUMBENT = 0.5076

#: The NAMED local wall-clock ceiling. Above this the route is Kaggle.
#:
#: Six hours is chosen, and the reason is operational rather than statistical: a run that finishes
#: overnight is recoverable if it fails, while a multi-day local run on a 4-thread CPU blocks the
#: machine and tends to be abandoned part-way, which yields a PARTIAL artifact that cannot be read
#: as the arm's result. The ceiling is a constant so the bar can freeze it and a reader can disagree
#: with it explicitly.
LOCAL_WALL_CEILING_SECONDS = 6 * 3600

#: Budget the probe will spend OBSERVING convergence. Bounded on purpose: the convergence run is a
#: scaled-down stand-in, not the arm.
DEFAULT_CONVERGENCE_BUDGET_SECONDS = 900

#: Minimum epochs worth calling a curve. Below this nothing about convergence is established.
MIN_CURVE_EPOCHS = 6

#: Epochs the sizing AIMS for, with the subsample derived to fit rather than the other way round.
#:
#: This direction is deliberate and was corrected after getting it backwards once. Sizing the
#: subsample first and taking whatever epoch count the budget then allows spends a larger budget on
#: MORE DATA at a fixed ~6 epochs, which is the wrong trade: what this probe needs to see is whether
#: the validation metric TURNS, and that is a property of the epoch axis. More data per epoch
#: improves the level of the curve and tells you nothing new about its shape.
TARGET_CURVE_EPOCHS = 15

ROUTE_LOCAL, ROUTE_KAGGLE = "local_cpu", "kaggle_t4"


# ---------------------------------------------------------------------------------------------------
# pure arithmetic -- the parts that decide things, testable without touching a model
# ---------------------------------------------------------------------------------------------------
def seconds_per_epoch(n_train: int, train_seq_per_second: float) -> float:
    """Wall-clock for one pass over `n_train` sequences at a measured throughput."""
    if train_seq_per_second <= 0:
        raise ValueError("throughput must be positive, got %r" % (train_seq_per_second,))
    if int(n_train) < 1:
        raise ValueError("n_train must be >= 1")
    return float(n_train) / float(train_seq_per_second)


def project_arm_cost(*, n_train: int, train_seq_per_second: float, n_epochs: int,
                     n_seeds: int = N_SEEDS) -> dict:
    """Total local cost of the arm: one pass per epoch, per seed."""
    if int(n_epochs) < 1 or int(n_seeds) < 1:
        raise ValueError("n_epochs and n_seeds must both be >= 1")
    per_epoch = seconds_per_epoch(n_train, train_seq_per_second)
    per_seed = per_epoch * int(n_epochs)
    total = per_seed * int(n_seeds)
    return {"seconds_per_epoch": round(per_epoch, 1),
            "minutes_per_epoch": round(per_epoch / 60.0, 1),
            "n_epochs": int(n_epochs), "n_seeds": int(n_seeds),
            "seconds_per_seed": round(per_seed, 1),
            "seconds_total": round(total, 1),
            "hours_total": round(total / 3600.0, 2)}


def choose_route(local_seconds: float, *,
                 ceiling_seconds: float = LOCAL_WALL_CEILING_SECONDS) -> tuple[str, str]:
    """Pick the compute route from a measured projection. Deliberately exercised BOTH ways.

    A rule that can only ever return one answer is not a rule, so the test suite drives this in
    both directions rather than only the one this host happens to produce.
    """
    if local_seconds <= 0:
        raise ValueError("local_seconds must be positive, got %r" % (local_seconds,))
    if local_seconds <= ceiling_seconds:
        return (ROUTE_LOCAL,
                "the projected local cost %.2f h is within the named ceiling of %.2f h, so the arm "
                "runs here and needs no Kaggle round trip."
                % (local_seconds / 3600.0, ceiling_seconds / 3600.0))
    return (ROUTE_KAGGLE,
            "the projected local cost %.2f h EXCEEDS the named ceiling of %.2f h by %.1fx, so the "
            "arm routes to a Kaggle T4 rather than being started locally and abandoned part-way."
            % (local_seconds / 3600.0, ceiling_seconds / 3600.0, local_seconds / ceiling_seconds))


def size_convergence_run(*, train_seq_per_second: float, budget_seconds: float,
                         n_available: int, frac: float = SPLIT_FRAC,
                         target_epochs: int = TARGET_CURVE_EPOCHS,
                         min_epochs: int = MIN_CURVE_EPOCHS) -> dict:
    """Derive a subsample size that buys `target_epochs` within the budget at the MEASURED rate.

    The EPOCH COUNT is the target and the subsample is derived to fit it, not the reverse. Sizing
    the subsample first spends a larger budget on more data at a fixed epoch count, which improves
    the curve's LEVEL and says nothing more about its SHAPE -- and the shape is the only thing this
    probe needs, because the question is whether the validation metric turns.

    Returns `feasible=False` rather than a tiny unusable run when the budget cannot buy
    `min_epochs` over a workable subsample: a two-epoch curve establishes nothing about convergence
    and must not be dressed up as a measurement.
    """
    if train_seq_per_second <= 0 or budget_seconds <= 0:
        raise ValueError("throughput and budget must both be positive")
    if int(target_epochs) < int(min_epochs):
        raise ValueError("target_epochs must be >= min_epochs")
    per_epoch_fraction = 1.0 - float(frac)
    affordable_seq_steps = float(budget_seconds) * float(train_seq_per_second)

    n_sub = int(min(n_available, affordable_seq_steps / (int(target_epochs) * per_epoch_fraction)))
    epochs = int(target_epochs)
    if n_sub < 64:
        # cannot afford the target; fall back toward the floor and see whether even that fits
        n_sub = int(min(n_available, affordable_seq_steps / (int(min_epochs) * per_epoch_fraction)))
        epochs = int(min_epochs)
    if n_sub < 64:
        return {"feasible": False, "n_subsample": max(0, n_sub), "epochs_attempted": epochs,
                "why": ("the budget of %.0f s at %.1f seq/s buys only %d subsample tiles over %d "
                        "epochs, which is too few to say anything about convergence"
                        % (budget_seconds, train_seq_per_second, max(0, n_sub), epochs))}

    n_train = int(round(n_sub * per_epoch_fraction))
    per_epoch = seconds_per_epoch(n_train, train_seq_per_second)
    return {"feasible": True, "n_subsample": n_sub, "n_train_per_epoch": n_train,
            "seconds_per_epoch": round(per_epoch, 1), "max_epochs": epochs,
            "target_epochs": int(target_epochs), "min_epochs": int(min_epochs),
            "hit_target_epochs": epochs == int(target_epochs),
            "projected_seconds": round(per_epoch * epochs, 1),
            "why": ("the EPOCH COUNT is the target and the subsample is derived to fit it, because "
                    "the curve's shape lives on the epoch axis; a bigger subsample at the same "
                    "epoch count would only raise the curve's level")}


def read_pretraining_recipe(config: dict | None) -> dict:
    """The checkpoint's own hyperparameters, read verbatim. Refuses a missing config."""
    if not config:
        raise ProbeRefused(
            "the base checkpoint carries no `config`, so the fine-tune's hyperparameters cannot be "
            "anchored to the pretraining recipe and would have to be invented. Refusing."
        )
    wanted = ("learning_rate", "patience", "batch_size", "n_epochs")
    missing = [k for k in wanted if k not in config]
    if missing:
        raise ProbeRefused("the checkpoint's config is missing %r" % missing)
    return {k: config[k] for k in wanted}


def deviations_from_recipe(recipe: dict, used: dict) -> list[dict]:
    """Every hyperparameter the probe did not inherit, with the reason it did not.

    A deviation is allowed; an UNRECORDED deviation is not. Any key that differs and carries no
    justification raises, so a silent substitution cannot reach the artifact.
    """
    out = []
    for k, v in sorted(used.items()):
        if k not in recipe:
            continue
        if recipe[k] == v["value"]:
            continue
        if not v.get("why"):
            raise ProbeRefused(
                "hyperparameter %r deviates from the pretraining recipe (%r -> %r) with no "
                "justification recorded. A fine-tune whose recipe is unrelated to the pretraining "
                "recipe is a weaker test, so the reason must be stated." % (k, recipe[k], v["value"])
            )
        out.append({"parameter": k, "pretraining_value": recipe[k],
                    "probe_value": v["value"], "why": v["why"]})
    return out


def curve_turned(val_scores: list[float], *, patience: int) -> dict:
    """Did the observed validation curve turn, and where?

    `turned` is True only when the best epoch is followed by `patience` epochs that do not beat it.
    A curve that is still improving at its last epoch has NOT turned, and saying so is the point:
    an untuned epoch budget read off a still-rising curve is the budget-limited null this arm must
    avoid being unable to interpret.
    """
    if not val_scores:
        return {"turned": None, "why": "no epochs were scored"}
    best_i = int(np.argmax(val_scores))
    after = len(val_scores) - 1 - best_i
    return {
        "n_epochs_scored": len(val_scores),
        "best_epoch": best_i,
        "best_value": float(val_scores[best_i]),
        "final_value": float(val_scores[-1]),
        "epochs_after_best": after,
        "patience_required": int(patience),
        "turned": bool(after >= int(patience)),
        "why": ("the best epoch is followed by %d non-improving epochs against a required patience "
                "of %d, so the curve %s"
                % (after, patience,
                   "HAS turned and the epoch budget is bounded by observation"
                   if after >= int(patience) else
                   "has NOT turned; the budget it suggests is a LOWER bound, and a null from a "
                   "budget this short could not be distinguished from a representation-limited one")),
    }


# ---------------------------------------------------------------------------------------------------
# measurement
# ---------------------------------------------------------------------------------------------------
def _onehot(seqs, length):
    import torch

    from dna_decode.glm.atlas_embed import NUC

    X = torch.zeros((len(seqs), 4, length))
    for b, s in enumerate(seqs):
        for j, ch in enumerate(s):
            k = NUC.get(ch)
            if k is not None:
                X[b, k, j] = 1.0
    return X


def _build_trainable(pa_repo=None):
    """Frozen-architecture encoder with gradients ON, plus a mean-pool linear read-out.

    The read-out mirrors the incumbent's SHAPE deliberately: the frozen arm mean-pools the latent to
    128 features and fits a ridge, so a mean-pool linear head makes the contrast a question about
    the REPRESENTATION rather than about the read-out.
    """
    import torch.nn as nn

    from dna_decode.glm.atlas_embed import load_atlas

    enc = load_atlas(pa_repo) if pa_repo is not None else load_atlas()
    enc.train()
    for p in enc.parameters():
        p.requires_grad_(True)
    head = nn.Linear(128, 1)
    return enc, head


def measure_train_throughput(*, batch_sizes=(32, 64, 128), length=150, steps=3,
                             pa_repo=None) -> dict:
    """Forward AND backward, warmup excluded. Reports the best sustained rate."""
    import torch

    enc, head = _build_trainable(pa_repo)
    opt = torch.optim.Adam(list(enc.parameters()) + list(head.parameters()), lr=1e-3)
    rows = {}
    g = torch.Generator().manual_seed(0)
    for bs in batch_sizes:
        X = torch.rand(bs, 4, length, generator=g)
        y = torch.randn(bs, generator=g)
        t0 = None
        for it in range(steps):
            if it == 1:
                t0 = time.time()
            opt.zero_grad()
            z = enc(X, return_latent=True)[2]
            loss = ((head(z.mean(dim=2)).squeeze(-1) - y) ** 2).mean()
            loss.backward()
            opt.step()
        el = (time.time() - t0) / max(1, steps - 1)
        rows[str(bs)] = {"seconds_per_step": round(el, 3),
                         "train_seq_per_second": round(bs / el, 1)}
    best_bs = max(rows, key=lambda k: rows[k]["train_seq_per_second"])
    return {"per_batch_size": rows,
            "best_batch_size": int(best_bs),
            "best_train_seq_per_second": rows[best_bs]["train_seq_per_second"],
            "includes_backward": True,
            "warmup_excluded": True,
            "rate_and_batch_come_from_the_SAME_row": (
                "The projected cost uses the rate measured AT the batch size the arm will use. "
                "Taking the fastest rate across batch sizes while training at a different one "
                "understates the cost -- on this host by about 1.7x, since batch 32 measures ~40 "
                "seq/s and batch 64 ~24 seq/s. So the chosen batch size IS the fastest measured "
                "one, which makes the two consistent by construction rather than by coincidence."),
            "note": ("the first step of each batch size is discarded; a cold step includes lazy "
                     "allocator and kernel-selection cost and would understate the sustained rate")}


def measure_convergence(tiles, sizing: dict, *, lr: float, batch_size: int, patience: int,
                        seed: int = 0, pa_repo=None) -> dict:
    """One seed, peak-blocked inner split, validation Spearman per epoch."""
    import torch

    from dna_decode.glm.genomewide import spearman

    idx = stride_sample(len(tiles), sizing["n_subsample"])
    sub = [tiles[i] for i in idx]
    tr, va, _rest = leave_peak_out_indices(sub, frac=SPLIT_FRAC, seed=seed)

    length = len(sub[0].seq)
    Xtr = _onehot([sub[i].seq for i in tr], length)
    # TRAIN ON log1p(expression), SCORE AGAINST RAW -- matching `ridge_spearman` exactly, whose own
    # docstring says any deviation makes the number incomparable with the committed incumbent
    # tables. The first version of this trainer used MSE on RAW expression, which would have
    # confounded the OBJECTIVE TRANSFORM with the REPRESENTATION: a loss could then have meant the
    # fine-tune optimised the wrong target scale rather than that the frozen representation was not
    # the limit -- a THIRD reading on top of the budget-versus-representation pair the bar names.
    ytr = torch.log1p(torch.tensor([float(sub[i].expression) for i in tr]))
    Xva = _onehot([sub[i].seq for i in va], length)
    yva = [float(sub[i].expression) for i in va]      # RAW, deliberately -- the scoring half

    enc, head = _build_trainable(pa_repo)
    opt = torch.optim.Adam(list(enc.parameters()) + list(head.parameters()), lr=lr)

    curve, t0 = [], time.time()
    for _epoch in range(sizing["max_epochs"]):
        enc.train()
        perm = torch.randperm(len(tr))
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
            raise ProbeRefused(
                "epoch %d produced non-finite validation predictions. Refusing: the scorer maps "
                "NaN to 0.0, so this would be recorded as 'no signal' rather than as divergence."
                % _epoch
            )
        curve.append(round(float(spearman(list(pv), yva)), 4))

    return {"seed": seed, "n_train": len(tr), "n_val": len(va),
            "n_peaks_train": len({sub[i].peak for i in tr}),
            "n_peaks_val": len({sub[i].peak for i in va}),
            "peaks_shared_between_splits": len({sub[i].peak for i in tr}
                                               & {sub[i].peak for i in va}),
            "val_spearman_by_epoch": curve,
            "wall_seconds": round(time.time() - t0, 1),
            "early_stop": curve_turned(curve, patience=patience),
            "scope": ("a SUBSAMPLE curve. It bounds the SHAPE of convergence, NOT the full-cohort "
                      "epoch count, because measuring the latter is the expensive run this probe "
                      "exists to size.")}


def probe(*, pa_repo=None, budget_seconds: float = DEFAULT_CONVERGENCE_BUDGET_SECONDS,
          cohort=None, skip_convergence: bool = False) -> dict:
    import torch

    from dna_decode.glm.atlas_embed import DEFAULT_PA_REPO, audit_checkpoint, load_atlas

    t_start = time.time()
    pa = Path(pa_repo if pa_repo is not None else DEFAULT_PA_REPO)
    ck_path = pa / "trained_weights/base_model/promoteratlas-base.pt"
    if not ck_path.is_file():
        raise ProbeRefused("base checkpoint absent: %s" % ck_path)
    size, digest = ck_path.stat().st_size, sha256_of(ck_path)
    if digest != EXPECTED_SHA256 or size != EXPECTED_BYTES:
        raise ProbeRefused(
            "checkpoint MOVED: sha256 %s (expected %s), %d bytes (expected %d)."
            % (digest[:16], EXPECTED_SHA256[:16], size, EXPECTED_BYTES))
    globals_named = audit_checkpoint(ck_path)

    import pathlib as _pl

    if os.name == "nt":
        _pl.PosixPath = _pl.WindowsPath  # type: ignore[misc, assignment]
    raw = torch.load(ck_path, map_location="cpu", weights_only=False)
    recipe = read_pretraining_recipe(raw.get("config"))

    model = load_atlas(pa)
    n_params = sum(p.numel() for p in model.parameters())
    if n_params != EXPECTED_PARAMS:
        raise ProbeRefused("%d parameters, expected %d" % (n_params, EXPECTED_PARAMS))

    tiles, st = cohort if cohort is not None else load_arm_cohort()
    n_train_full = int(round(len(tiles) * (1.0 - SPLIT_FRAC)))

    thr = measure_train_throughput(pa_repo=pa_repo)
    rate = thr["best_train_seq_per_second"]
    sizing = size_convergence_run(train_seq_per_second=rate, budget_seconds=budget_seconds,
                                  n_available=len(tiles))

    # The batch size is the fastest MEASURED one, so the rate used for the projection and the batch
    # the arm actually trains at are the same row. The pretraining 1024 is a recorded DEVIATION
    # rather than a silent substitution.
    probe_batch = int(thr["best_batch_size"])
    probe_patience = 2
    used = {
        "learning_rate": {"value": recipe["learning_rate"], "why": ""},
        "batch_size": {"value": probe_batch,
                       "why": ("the pretraining batch of %r was used on a CUDA device for a masked "
                               "language-model objective. This host is CPU-only with 4 threads, and "
                               "the throughput was MEASURED across batch sizes rather than assumed: "
                               "%r is the fastest of %r, so it is the batch the arm uses AND the "
                               "rate the projection is computed from. Bigger is not faster here."
                               % (recipe["batch_size"], probe_batch,
                                  {k: v["train_seq_per_second"]
                                   for k, v in thr["per_batch_size"].items()}))},
        "patience": {"value": probe_patience,
                     "why": ("the pretraining patience of %r is calibrated to a %r-epoch budget; a "
                             "subsample curve of %d epochs cannot express it, so the probe uses a "
                             "proportionate value and the artifact reports whether the curve turned "
                             "under it" % (recipe["patience"], recipe["n_epochs"],
                                           sizing.get("max_epochs", 0)))},
        "n_epochs": {"value": sizing.get("max_epochs", 0),
                     "why": ("the pretraining budget of %r epochs is for the masked LM objective on "
                             "the full corpus; this is a bounded convergence OBSERVATION sized from "
                             "the measured throughput" % recipe["n_epochs"])},
    }
    deviations = deviations_from_recipe(recipe, used)

    convergence = None
    if sizing["feasible"] and not skip_convergence:
        convergence = measure_convergence(
            tiles, sizing, lr=float(recipe["learning_rate"]), batch_size=probe_batch,
            patience=probe_patience, pa_repo=pa_repo)

    epochs_for_projection = (convergence["early_stop"]["best_epoch"] + 1
                             if convergence else sizing.get("max_epochs", 1))
    local = project_arm_cost(n_train=n_train_full, train_seq_per_second=rate,
                             n_epochs=max(1, epochs_for_projection), n_seeds=N_SEEDS)
    route, route_why = choose_route(local["seconds_total"])

    verdict, detail = classify(sizing, convergence, route)

    return {
        "schema": "finetune-arm-feasibility-v1",
        "date": "2026-11-02",
        "arm": ARM,
        "host": {"platform": platform.platform(), "python": platform.python_version(),
                 "cpu_count": os.cpu_count(), "cuda_available": bool(torch.cuda.is_available()),
                 "torch_threads": int(torch.get_num_threads()), "torch": torch.__version__,
                 "free_ram_bytes": free_ram_bytes()},
        "checkpoint": {"path": str(ck_path), "bytes": size, "sha256": digest,
                       "pickle_globals": globals_named,
                       "n_parameters_via_model_parameters": n_params,
                       "tractability_note": (
                           "%d parameters, not a 1.2B language model, which is why a fine-tune is "
                           "tractable here at all and why this is NOT comparable to this project's "
                           "blocked large-LM fine-tunes." % n_params)},
        "pretraining_recipe": recipe,
        "hyperparameters_used": {k: v["value"] for k, v in used.items()},
        "deviations_from_pretraining_recipe": deviations,
        "deviation_rule": ("a deviation is allowed; an UNRECORDED deviation raises. Every key that "
                           "differs from the checkpoint's own config carries its reason here."),
        "substrate": {"n_tiles": len(tiles), "n_train_full_cohort": n_train_full,
                      "split_frac": SPLIT_FRAC, "path": st.path},
        "throughput": thr,
        "convergence_sizing": sizing,
        "convergence": convergence,
        "local_projection": local,
        "route": {"chosen": route, "why": route_why,
                  "ceiling_seconds": LOCAL_WALL_CEILING_SECONDS,
                  "ceiling_hours": LOCAL_WALL_CEILING_SECONDS / 3600.0,
                  "ceiling_is_operational_not_statistical": (
                      "six hours is chosen because an overnight run is recoverable while a "
                      "multi-day local run on 4 CPU threads tends to be abandoned part-way, which "
                      "yields a PARTIAL artifact nobody can read. A reader may disagree with the "
                      "number explicitly; it is a constant, not a judgement made per run."),
                  "alternatives_considered": [ROUTE_LOCAL, ROUTE_KAGGLE]},
        "target_transform": {
            "trained_on": "log1p(expression)",
            "scored_against": "raw expression",
            "inherited_not_chosen": True,
            "source": ("dna_decode.glm.atlas_embed.ridge_spearman, which fits on log1p(expression) "
                       "and scores Spearman against RAW expression, and whose docstring states that "
                       "any deviation makes the number incomparable with the committed incumbent "
                       "tables"),
            "why_it_matters": ("training on RAW expression instead would confound the OBJECTIVE "
                               "TRANSFORM with the REPRESENTATION, adding a third reading to the "
                               "budget-versus-representation pair: a loss could mean the fine-tune "
                               "optimised the wrong target scale. Fixed before this measurement, "
                               "not after."),
        },
        "comparator": {"frozen_incumbent_median_spearman": FROZEN_INCUMBENT,
                       "n_seeds": N_SEEDS, "split_frac": SPLIT_FRAC,
                       "comparison": ("the fine-tuned encoder's median Spearman against the FROZEN "
                                      "encoder's %.4f on byte-identical splits" % FROZEN_INCUMBENT)},
        "the_two_readings_that_must_not_be_conflated": (
            "Beating %.4f shows the frozen REPRESENTATION was the limit. Failing to beat it shows "
            "the SUBSTRATE or the TRAINING BUDGET was -- and a budget-limited null is NOT a "
            "representation-limited one. The `convergence.early_stop.turned` field is what "
            "distinguishes them: if the curve never turned, a null is uninterpretable."
            % FROZEN_INCUMBENT),
        "verdict": verdict,
        "verdict_detail": detail,
        "wall_seconds": round(time.time() - t_start, 1),
    }


def classify(sizing: dict, convergence: dict | None, route: str) -> tuple[str, str]:
    """Closed verdict set."""
    if not sizing.get("feasible"):
        return ("INDETERMINATE_BUDGET_TOO_SMALL_TO_OBSERVE_CONVERGENCE", sizing.get("why", ""))
    if convergence is None:
        return ("INDETERMINATE_CONVERGENCE_NOT_MEASURED",
                "the convergence curve was skipped, so a null from this arm could not be "
                "distinguished from a budget-limited one. The arm must not run in this state.")
    es = convergence["early_stop"]
    if convergence["peaks_shared_between_splits"]:
        return ("INDETERMINATE_SPLIT_LEAKED",
                "%d peaks appear in BOTH the inner train and validation splits, so the curve is "
                "not a held-out measurement" % convergence["peaks_shared_between_splits"])
    if not es["turned"]:
        return ("FEASIBLE_BUT_CURVE_DID_NOT_TURN",
                "throughput and route are measured and the arm is runnable on %s, but the observed "
                "curve was still improving at its last epoch, so the epoch budget it suggests is a "
                "LOWER BOUND. A null under this budget could not be read as "
                "representation-limited. %s" % (route, es["why"]))
    return ("FEASIBLE_CURVE_TURNED",
            "throughput, convergence and route are all measured: best epoch %d of %d, followed by "
            "%d non-improving epochs, routed to %s."
            % (es["best_epoch"], es["n_epochs_scored"], es["epochs_after_best"], route))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pa-repo", default=None)
    ap.add_argument("--budget-seconds", type=float, default=DEFAULT_CONVERGENCE_BUDGET_SECONDS)
    ap.add_argument("--skip-convergence", action="store_true",
                    help="throughput and route only; yields an INDETERMINATE verdict by design")
    ap.add_argument("--out", default=str(OUT))
    args = ap.parse_args(argv)

    try:
        rec = probe(pa_repo=args.pa_repo, budget_seconds=args.budget_seconds,
                    skip_convergence=args.skip_convergence)
    except ProbeRefused as e:
        print("REFUSED: %s" % e, flush=True)
        return 2
    except Exception as e:                                     # noqa: BLE001
        print("PROBE FAILED (%s): %s" % (type(e).__name__, e), flush=True)
        return 1

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(rec, indent=2, sort_keys=True), encoding="utf-8")

    print("PARAMS   %d (tractable; not a 1.2B LM)"
          % rec["checkpoint"]["n_parameters_via_model_parameters"])
    print("THRPUT   %.1f train seq/s (fwd+bwd, warmup excluded)"
          % rec["throughput"]["best_train_seq_per_second"])
    if rec["convergence"]:
        c = rec["convergence"]
        print("CURVE    %d epochs on %d train / %d val tiles, %d shared peaks"
              % (len(c["val_spearman_by_epoch"]), c["n_train"], c["n_val"],
                 c["peaks_shared_between_splits"]))
        print("         %s" % c["val_spearman_by_epoch"])
        print("         best epoch %d at %.4f, turned=%s"
              % (c["early_stop"]["best_epoch"], c["early_stop"]["best_value"],
                 c["early_stop"]["turned"]))
    print("LOCAL    %.1f min/epoch -> %.2f h for %d seeds x %d epochs"
          % (rec["local_projection"]["minutes_per_epoch"], rec["local_projection"]["hours_total"],
             rec["local_projection"]["n_seeds"], rec["local_projection"]["n_epochs"]))
    print("ROUTE    %s" % rec["route"]["chosen"])
    print("         %s" % rec["route"]["why"])
    print("VERDICT  %s" % rec["verdict"])
    print("         %s" % rec["verdict_detail"])
    print("wall     %.1f s" % rec["wall_seconds"])
    print("-> %s" % out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
