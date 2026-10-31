"""Does the self-distillation degeneracy PERSIST when the teacher carries a LEARNED component?

    uv run python scripts/glm_selfdistill_learned_teacher.py [--seeds 3] [--rates 0.02,0.05,0.10]

**Why this exists.** `wiki/glm_selfdistill_probe_2026-10-07.json` blocked self-distillation with a
MECHANICAL argument, measured: a ridge teacher's predictions lie in the span of its own features, so a
linear student with the same basis recovered the same function at Spearman agreement **1.0000** -- exactly,
at every mutation rate. The standing rule became *forbidden until a learned representation lands*.

That rule's reason is about the TEACHER'S FUNCTION CLASS, not about which arm wins a ranking contest. So
the decidable question is: swap the linear teacher for a learned one, hold everything else, and does the
1.0000 survive? A persist answer makes the open authority question moot. A dissolve answer is evidence for
one reading of it. **Neither answer ratifies anything** -- the `blocked_by` clause is untouched by this run.

**The bar is frozen at `wiki/glm_selfdistill_learned_teacher_acceptance_bar.json`**, including all five
outcomes, the check order, and a forecast. Thresholds are MIRRORED from the 2026-10-07 probe so the two
results are commensurable; picking new ones here would be the retro-fit that probe was written to avoid.

**The ridge-teacher arm is not redundant and it is not ceremony.** The original probe ran on GSE108535
sigma70 element pairs, whose cache is ABSENT on this host, so its verdict cannot be reproduced on its own
ground. The ridge arm re-measures the degenerate anchor on THIS rung. If it fails to reproduce, a
non-degenerate conv result would be a property of the SUBSTRATE and the conv arm is void -- which is
outcome 2, checked before any headline.

**The teacher must be ONE function, and that is measured rather than assumed.** The conv is refitted once
per mutation rate, relying on CPU determinism to give the same weights each time, and every perturbed call
carries 500 shared real-test rows so the claim is checked rather than trusted. If they diverge the teacher
is not a single function and the arm is void, so the script REFUSES rather than average several models.

**That check's first tolerance was itself the defect, and the fix is a derived number.** It was set at
1e-9, which is BELOW the arithmetic's own floor: the anchor call scores 11,198 rows and a perturbed call
8,500, so they run through different batch shapes and accumulate float32 in a different order. Measured at
full scale, two calls PROVEN to hold the same model (identical `epochs_ran` 11, identical
`best_inner_val_mse` 0.759698, and a repeated anchor bit-identical at 0.0) still differ by **7.15e-07**. A
different model shows ~1.0. `DRIFT_TOL` is therefore 1e-4, with four orders of magnitude of headroom on
each side, and the per-rate drifts are STAMPED into the artifact instead of claimed in prose.

**One unexplained observation, recorded rather than theorised away.** A first full-scale run refused at
rate 0.10 with a drift of **1.230** while rates 0.02 and 0.05 passed at <=1e-9 on identical shapes. Neither
figure reproduces: the controlled repeat gives 7.15e-07 for that same shape pair, and a repeated anchor
gives exactly 0.0. So the 1.230 is a single unreproduced event and the <=1e-9 pair should not have been
achievable either. It is logged here as an open anomaly. It is NOT evidence that `fit_encoder` jitters --
that was checked directly and it does not -- so the published tile-arm numbers are unaffected.

**The non-vacuity control is a SHUFFLED-LABEL student, and the first version of it was wrong.** The bar
originally gated on a GC-only student collapsing 0.20 in rho below its teacher -- copied from the
2026-10-07 probe, where GC fell 0.5963 -> 0.1736. On THIS rung GC scores **0.3000** against the conv's
0.3235, so that gate was arithmetically unreachable and the first run returned
`PROBE_VACUOUS_CONTROL_DID_NOT_COLLAPSE`. The premise "a GC-only student is impoverished" is simply false
on peak tiles -- which is the G-A encoder-gate finding restated. The replacement permutes the teacher's
labels: a permuted label carries no information about the teacher's function however predictive any basis
happens to be here, so the control is substrate-independent. Its threshold was DERIVED BEFORE RUNNING
(expected 0, null SE 1/sqrt(n-1) = 0.0095 at n=11198, gate at 0.20 = ~21 SE) and measured at
agreement 0.0128 / -0.0287. `student_gc` is still fitted and reported, just no longer decisive.
The three headline thresholds were not touched, so no measured verdict could move -- see the bar's
`superseded_v1.outcome_neutrality_of_this_amendment`.

**What this run deliberately cannot do.** The teacher is the conv encoder, NOT the 3-channel stack that
posted 0.4404. Two of that stack's channels come from Promoter Calculator, which is not installed here, and
what exists on disk is a precomputed table over the 44,106 fixed tiles -- structurally unable to score a
perturbed sequence. So this bounds the conv case: the learned half of the stack.
"""
from __future__ import annotations

import argparse
import json
import random
import sys
from datetime import date
from pathlib import Path

from dna_decode.glm.encoder import REGISTERED_PROTOCOL, count_parameters, fit_encoder
from dna_decode.glm.genomewide import (
    DEFAULT_CACHE,
    Tile,
    gc_features,
    leave_peak_out,
    load_peak_tiles,
    positional_onehot,
    spearman,
)

BAR = Path("wiki/glm_selfdistill_learned_teacher_acceptance_bar.json")

#: Mirrored from wiki/glm_selfdistill_probe_2026-10-07.json. See the bar's why_mirrored note.
AGREEMENT_FLOOR = 0.99
DRHO_TOL = 0.01
BEAT_MARGIN = 0.02
CONTROL_MUST_DROP = 0.20          # v1 gate, RETIRED -- still measured and reported, no longer decisive
SHUF_MAX_ABS_AGREEMENT = 0.20     # v2 gate. Derivation in the bar's control_threshold_derivation.
SHUF_MAX_ABS_RHO = 0.20
RECONCILE_CONV_SEED0 = 0.324
RECONCILE_TOL = 0.02

#: Tolerance for the teacher-is-one-function check. DERIVED, and the first value was WRONG.
#:
#: Targets are standardised on the inner-train rows, so predictions are ~unit scale. The forward pass
#: is float32 and the anchor call scores 11,198 rows while a perturbed call scores 8,500, so the two
#: run through different batch shapes and accumulate in a different order. MEASURED on this substrate
#: at full scale: 7.15e-07 between two calls PROVEN to hold the same model (identical epochs_ran 11,
#: identical best_inner_val_mse 0.759698, and a repeated anchor call bit-identical at 0.0).
#:
#: A genuinely DIFFERENT model shows ~1.0 on standardised predictions. So 1e-4 sits ~140x above
#: measured rounding and ~10,000x below a model difference -- four orders of magnitude of headroom on
#: each side. The original 1e-9 was below the arithmetic's own floor and could only ever fail.
DRIFT_TOL = 1e-4

SAMPLE_TILES = 8000
COPIES = 1


class ProbeError(RuntimeError):
    """A condition under which no number may be reported."""


def perturb(seqs: list[str], rate: float, copies: int, seed: int) -> list[str]:
    """Random point substitutions -- AlphaGenome's 'randomly perturb the input', at nucleotide level."""
    rng = random.Random(seed)
    alts = {b: [x for x in "ACGT" if x != b] for b in "ACGT"}
    out: list[str] = []
    for s in seqs:
        for _ in range(copies):
            chars = list(s)
            n_mut = max(1, int(round(len(chars) * rate)))
            for pos in rng.sample(range(len(chars)), n_mut):
                chars[pos] = rng.choice(alts.get(chars[pos], ["A", "C", "G"]))
            out.append("".join(chars))
    return out


def _fit_predict(Xtr, ytr, Xte, kind: str):
    """Student fitting. Mirrors the 2026-10-07 probe's _fit exactly, including the MLP's iteration budget."""
    import numpy as np
    from sklearn.linear_model import Ridge
    from sklearn.preprocessing import StandardScaler

    sc = StandardScaler().fit(Xtr)
    if kind == "mlp":
        from sklearn.neural_network import MLPRegressor

        # Generous max_iter ON PURPOSE: this arm tests the information bound, and an under-trained
        # student failing to beat its teacher would be evidence about the optimiser, not the bound.
        m = MLPRegressor(hidden_layer_sizes=(64,), max_iter=400, early_stopping=True,
                         n_iter_no_change=15, random_state=0).fit(sc.transform(Xtr), np.asarray(ytr))
    else:
        m = Ridge(alpha=1.0).fit(sc.transform(Xtr), np.asarray(ytr))
    return [float(v) for v in m.predict(sc.transform(Xte))]


def _as_tiles(seqs: list[str]) -> list[Tile]:
    """Perturbed sequences as Tile objects. `expression` is a placeholder and is never scored --
    fit_encoder applies group_key to TRAIN only, and this set is passed as test."""
    return [Tile(seq=s, expression=0.0, peak="__perturbed__", strand="+",
                 category="tile", active="NA", rep1=0.0, rep2=0.0) for s in seqs]


def run_seed(tiles, seed: int, rates: list[float]) -> dict:
    import numpy as np

    train, test, held = leave_peak_out(tiles, frac=0.25, seed=seed)
    tr_seqs = [t.seq for t in train]
    te_seqs = [t.seq for t in test]
    measured = [t.expression for t in test]

    rng = random.Random(1000 + seed)
    sample = rng.sample(tr_seqs, min(SAMPLE_TILES, len(tr_seqs)))

    Xtr_oh = np.asarray(positional_onehot(tr_seqs), dtype=float)
    Xte_oh = np.asarray(positional_onehot(te_seqs), dtype=float)
    Xte_gc = np.asarray(gc_features(te_seqs), dtype=float)
    ytr = np.log1p(np.asarray([t.expression for t in train], dtype=float))

    # ---- teacher R: the degenerate anchor, re-measured on THIS rung -------------------------------
    ridge_te = _fit_predict(Xtr_oh, ytr, Xte_oh, "ridge")
    rho_ridge = spearman(ridge_te, measured)

    # ---- teacher C: the learned teacher -----------------------------------------------------------
    # MEMORY LAYOUT IS LOAD-BEARING, not tidiness. fit_encoder scores all of `test` in ONE
    # model(Xte) call, and the conv's intermediate activations are ~32 channels x ~145 positions per
    # sequence across three branches -- about 1 GB transient at 19k sequences. The first version of
    # this passed test + perturbed together (19,198 rows) and was KILLED by allocation pressure on a
    # 16 GB host with 2.9 GB free, losing its own traceback to a buffered stdout. So the anchor call
    # and the perturbed calls are now SEPARATE, which bounds peak at the anchor's 11,198 rows -- a
    # size already proven to run here.
    anchor = fit_encoder(train, test, seed=seed, device="cpu")
    conv_te = anchor["predictions"]
    rho_conv = spearman(conv_te, measured)
    conv_meta = {k: anchor[k] for k in ("epochs_ran", "early_stopped", "n_parameters", "device")}

    # A DETERMINISM TAG rides on every perturbed call so the teacher is proven to be ONE function
    # rather than assumed to be. `fit_encoder` refits per call; CPU determinism should give identical
    # weights, and if it does not, an agreement number averaged over several models is meaningless.
    TAG = 500
    conv_on_pert: dict[float, list[float]] = {}
    drifts: dict[str, float] = {}
    for rate in rates:
        pert = perturb(sample, rate, COPIES, seed)
        res = fit_encoder(train, test[:TAG] + _as_tiles(pert), seed=seed, device="cpu")
        pred = res["predictions"]
        drift = max(abs(a - b) for a, b in zip(conv_te[:TAG], pred[:TAG]))
        if drift > DRIFT_TOL:
            raise ProbeError(
                f"seed {seed} rate {rate}: the conv teacher is NOT deterministic across refits "
                f"(max abs drift {drift:.3e} on {TAG} shared real-test rows). Several different "
                f"models cannot be one teacher; refusing to report an agreement number over them.")
        drifts[str(rate)] = drift
        conv_on_pert[rate] = pred[TAG:]
        del res, pred

    teachers = {
        "ridge_positional_onehot": {"rho": rho_ridge, "te": ridge_te},
        "conv_encoder": {"rho": rho_conv, "te": conv_te},
    }

    rows = []
    for rate in rates:
        pert = perturb(sample, rate, COPIES, seed)
        Xp_oh = np.asarray(positional_onehot(pert), dtype=float)
        Xp_gc = np.asarray(gc_features(pert), dtype=float)

        for tname, t in teachers.items():
            # labels are the TEACHER'S OWN outputs on the perturbed sequences -- no measurement
            if tname == "ridge_positional_onehot":
                yp = np.asarray(_fit_predict(Xtr_oh, ytr, Xp_oh, "ridge"), dtype=float)
            else:
                yp = np.asarray(conv_on_pert[rate], dtype=float)

            # THE v2 NON-VACUITY CONTROL: the same basis and the same perturbed rows, with the
            # teacher's labels PERMUTED. A permuted label carries no information about the teacher's
            # function whatever the feature basis is worth on this rung, which is exactly why it
            # replaced the GC student -- see the bar's superseded_v1 block.
            yp_shuf = np.asarray(yp, dtype=float).copy()
            np.random.default_rng(30_000 + seed).shuffle(yp_shuf)

            arms = {}
            for sname, kind, Xp_f, Xte_f, y_use in (
                ("student_linear", "ridge", Xp_oh, Xte_oh, yp),
                ("student_mlp", "mlp", Xp_oh, Xte_oh, yp),
                ("student_shuffled", "ridge", Xp_oh, Xte_oh, yp_shuf),
                ("student_gc", "ridge", Xp_gc, Xte_gc, yp),
            ):
                pred = _fit_predict(Xp_f, y_use, Xte_f, kind)
                r = spearman(pred, measured)
                arms[sname] = {
                    "rho_vs_measured": round(r, 4),
                    "agreement_with_teacher": round(spearman(pred, t["te"]), 4),
                    "delta_rho_vs_teacher": round(r - t["rho"], 4),
                }
            rows.append({"mutation_rate": rate, "teacher": tname,
                         "n_perturbed": len(pert), "arms": arms})
            bits = "  ".join(
                f"{n.replace('student_', '')}: rho={a['rho_vs_measured']:+.4f} "
                f"agree={a['agreement_with_teacher']:.4f}" for n, a in arms.items())
            print(f"  seed {seed} rate={rate:<5.2f} {tname:<26s} {bits}", flush=True)

    return {
        "seed": seed, "held_out_peaks": held if isinstance(held, int) else len(held),
        "n_train": len(train), "n_test": len(test), "n_sampled_for_perturbation": len(sample),
        "teacher_rho": {"ridge_positional_onehot": round(rho_ridge, 4),
                        "conv_encoder": round(rho_conv, 4)},
        # MEASURED per rate, not asserted in prose. These are the max abs differences between the
        # anchor's and each perturbed call's predictions on TAG shared real-test rows. Expect ~1e-6
        # (float32 rounding across differing batch shapes); ~1.0 would mean a different model.
        "conv_refit_drift_on_shared_rows": {k: f"{v:.3e}" for k, v in drifts.items()},
        "conv_refit_drift_tolerance": DRIFT_TOL,
        "conv_meta": conv_meta,
        "rows": rows,
    }


def verdict(seeds: list[dict]) -> tuple[str, str]:
    """Mechanical, against the frozen check order. First match wins."""
    flat = [(s["seed"], r) for s in seeds for r in s["rows"]]
    # 1 -- non-vacuity FIRST: without a collapsing control nothing is interpretable.
    # v2 gates on the SHUFFLED-LABEL student. The v1 GC gate is retired (unachievable on this rung)
    # but is still computed below and reported, because it is informative that it cannot pass.
    bad_ctrl = []
    for sd, r in flat:
        a = r["arms"]["student_shuffled"]
        if (abs(a["agreement_with_teacher"]) > SHUF_MAX_ABS_AGREEMENT
                or abs(a["rho_vs_measured"]) > SHUF_MAX_ABS_RHO):
            bad_ctrl.append((sd, r["teacher"], r["mutation_rate"],
                             a["agreement_with_teacher"], a["rho_vs_measured"]))
    if bad_ctrl:
        return ("PROBE_VACUOUS_CONTROL_DID_NOT_COLLAPSE",
                f"the SHUFFLED-LABEL student failed to collapse in {len(bad_ctrl)} cell(s) "
                f"(seed, teacher, rate, agreement, rho) = {bad_ctrl[:4]}; the probe cannot distinguish "
                f"a real match from a fabricated one, and per the bar's pre-commitment the control is "
                f"NOT weakened a second time")

    def reproduces(r):
        a = r["arms"]["student_linear"]
        return (a["agreement_with_teacher"] >= AGREEMENT_FLOOR
                and abs(a["delta_rho_vs_teacher"]) <= DRHO_TOL)

    ridge_rows = [(sd, r) for sd, r in flat if r["teacher"] == "ridge_positional_onehot"]
    conv_rows = [(sd, r) for sd, r in flat if r["teacher"] == "conv_encoder"]

    # 2 -- the anchor guard. Without it a conv dissolve could be a property of the substrate.
    anchor_fail = [(sd, r["mutation_rate"], r["arms"]["student_linear"]["agreement_with_teacher"],
                    r["arms"]["student_linear"]["delta_rho_vs_teacher"])
                   for sd, r in ridge_rows if not reproduces(r)]
    if anchor_fail:
        return ("ANCHOR_FAILED_RIDGE_TEACHER_NOT_DEGENERATE_HERE",
                f"the ridge teacher's linear student did NOT reproduce it on this rung in "
                f"{len(anchor_fail)} cell(s): (seed, rate, agreement, d_rho) = {anchor_fail}. This rung "
                f"does not reproduce the 2026-10-07 condition, so THE CONV ARM IS UNINTERPRETABLE and "
                f"must not be read")

    # 3 -- did any student actually beat its own teacher on real measured expression?
    controls = {"student_gc", "student_shuffled"}
    beats = [(sd, r["teacher"], r["mutation_rate"], n, a["delta_rho_vs_teacher"])
             for sd, r in flat for n, a in r["arms"].items()
             if n not in controls and a["delta_rho_vs_teacher"] > BEAT_MARGIN]
    if beats:
        return "DISTILLATION_HELPS_NOW", f"students beating their teacher by > {BEAT_MARGIN}: {beats}"

    conv_repro = [(sd, r["mutation_rate"], r["arms"]["student_linear"]["agreement_with_teacher"],
                   r["arms"]["student_linear"]["delta_rho_vs_teacher"]) for sd, r in conv_rows]
    # 4 / 5 -- the headline, now that the anchor licenses reading it
    if all(reproduces(r) for _, r in conv_rows):
        return ("MECHANISM_PERSISTS_BAN_STANDS",
                f"the conv teacher's linear student ALSO reproduces it at every rate and seed "
                f"(seed, rate, agreement, d_rho) = {conv_repro}; the degeneracy is not a consequence of "
                f"the teacher being linear, so the ban stands on mechanism")
    return ("MECHANISM_DISSOLVES_WITH_LEARNED_TEACHER",
            f"the ridge anchor reproduced but the conv teacher's linear student did NOT "
            f"(seed, rate, agreement, d_rho) = {conv_repro}; the mechanical objection is measurably "
            f"absent with a learned teacher. READ NARROWLY -- the objection is gone, not the bound")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cache-dir", default=str(DEFAULT_CACHE))
    ap.add_argument("--seeds", type=int, default=3)
    ap.add_argument("--rates", default="0.02,0.05,0.10")
    ap.add_argument("--out", default=None)
    # Per-seed process isolation. One seed per process releases the conv's ~1 GB of transient
    # activations between seeds and bounds the cost of a crash to a single seed instead of the run.
    ap.add_argument("--seed-only", type=int, default=None,
                    help="run exactly this seed and checkpoint it; do not aggregate or verdict")
    ap.add_argument("--seed-dir", default="wiki/fragment_null_batches/selfdistill_seeds")
    ap.add_argument("--aggregate", action="store_true",
                    help="read the per-seed checkpoints, reconcile, and write the final artifact")
    args = ap.parse_args()

    if not BAR.exists():
        print(f"REFUSED: the frozen bar {BAR} is absent; this run must not define its own outcomes")
        return 2
    bar = json.loads(BAR.read_text(encoding="utf-8"))

    tiles_all, stats = load_peak_tiles(args.cache_dir, allow_download=False)
    tiles = [t for t in tiles_all if t.category == "tile" and len(t.seq) == 150]
    if len(tiles) != bar["substrate"]["n_tiles"]:
        print(f"REFUSED: substrate drift -- {len(tiles)} tiles, bar pins "
              f"{bar['substrate']['n_tiles']}")
        return 3
    print(f"substrate OK: {len(tiles)} tiles from {stats.path}")

    rates = [float(r) for r in args.rates.split(",")]
    seed_dir = Path(args.seed_dir)
    seed_dir.mkdir(parents=True, exist_ok=True)

    if args.seed_only is not None:
        try:
            rec = run_seed(tiles, args.seed_only, rates)
        except ProbeError as e:
            print(f"REFUSED: {e}", flush=True)
            return 4
        p = seed_dir / f"seed-{args.seed_only:02d}.json"
        p.write_text(json.dumps(rec, indent=2), encoding="utf-8")
        print(f"checkpointed seed {args.seed_only} -> {p}", flush=True)
        return 0

    if args.aggregate:
        seeds = []
        for s in range(args.seeds):
            p = seed_dir / f"seed-{s:02d}.json"
            if not p.exists():
                print(f"REFUSED: seed checkpoint {p} is absent; aggregating a PARTIAL set would "
                      f"report a verdict over fewer seeds than the bar pins ({args.seeds})")
                return 6
            seeds.append(json.loads(p.read_text(encoding="utf-8")))
    else:
        seeds = []
        try:
            for seed in range(args.seeds):
                seeds.append(run_seed(tiles, seed, rates))
        except ProbeError as e:
            print(f"REFUSED: {e}", flush=True)
            return 4

    # reconcile: the conv teacher must reproduce the PUBLISHED seed-0 figure, or this is not the
    # published protocol on the published substrate and no comparison to it is licensed
    conv0 = seeds[0]["teacher_rho"]["conv_encoder"]
    if abs(conv0 - RECONCILE_CONV_SEED0) > RECONCILE_TOL:
        print(f"REFUSED: reconcile failed -- conv seed 0 = {conv0:.4f}, published "
              f"{RECONCILE_CONV_SEED0} +/- {RECONCILE_TOL}")
        return 5
    print(f"reconcile OK: conv seed 0 = {conv0:.4f} vs published {RECONCILE_CONV_SEED0}")

    v, detail = verdict(seeds)
    artifact = {
        "schema": "glm-selfdistill-learned-teacher-v1",
        "date": str(date.today()),
        "bar": str(BAR),
        "question": bar["question"],
        "why_this_is_not_the_authority_call": bar["why_this_is_not_the_authority_call"],
        "substrate": {**bar["substrate"], "n_tiles_loaded": len(tiles), "path": str(stats.path)},
        "thresholds": bar["thresholds_mirrored_from_2026_10_07"],
        "outcomes_defined_before_running": bar["outcomes_defined_before_running"],
        "conv_protocol": REGISTERED_PROTOCOL.as_dict(),
        "conv_n_parameters": count_parameters(),
        "reconcile": {"conv_seed0": conv0, "published": RECONCILE_CONV_SEED0, "ok": True},
        "seeds": seeds,
        "verdict": v,
        "verdict_detail": detail,
        "forecast": bar["forecast"],
        "honest_limits": bar["honest_limits_fixed_in_advance"],
    }
    out = Path(args.out) if args.out else Path(f"wiki/glm_selfdistill_learned_teacher_{date.today()}.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(artifact, indent=2), encoding="utf-8")
    print(f"\nVERDICT: {v}\n  {detail}\nwrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
