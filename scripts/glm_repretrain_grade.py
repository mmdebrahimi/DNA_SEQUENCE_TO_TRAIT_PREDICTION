"""Grade the representative re-pretrain arm against its frozen reading rule.

FROZEN BEFORE THE NUMBERS. This file is committed while both pretrains are still running, so the
gate, the branch logic and every threshold are fixed before any representative or reseed
cross-entropy exists. A grader written after seeing the numbers is a post-hoc rescue; that is the
objection the whole companion discipline exists to pre-empt.

WHAT IT GRADES. `wiki/glm_representative_repretrain_reading_rule_2026-11-24.json`:

  PRIMARY  paired per-window change in MODEL cross-entropy between the representative checkpoint
           and the committed one, on identical eval windows at identical sampled positions, on the
           `orthology_free` set at the frozen 120x30 geometry, selection seed 0, scoring seed 0,
           boot 2000.
  GATE     reported RESOLVED only if |primary| exceeds |reseed floor| AND the primary bootstrap CI
           excludes zero. Otherwise UNRESOLVED_WITHIN_TRAJECTORY_NOISE.
  ESCALATE if |primary| lands within 2x the single floor draw, NO branch is reported and seeds 2
           and 3 must be run first. Pre-registered in the rule, so it is checked FIRST and it
           dominates the gate.

SIGN. The rule contradicts itself: its `quantity` says `representative minus committed` while its
`sign_convention` says `positive = the representative model has LOWER cross-entropy = BETTER`. In
cross-entropy nats lower is better, so those specify opposite subtractions. This grader follows the
SIGN CONVENTION and computes `committed - representative`, which also matches the scorer's own
`d = [mk - ce]` orientation (baseline minus assessed). The reasoning, the evidence and why the
choice is outcome-neutral are recorded in
`wiki/glm_representative_repretrain_sign_correction_2026-11-24.json`. The rule file is NOT edited.

HOW TO PRODUCE THE THREE INPUTS. One `scripts/glm_beyond_kmer.py` run per checkpoint. The corpus is
deliberately NOT parameterised with `--ckpt`, so all three runs score identical windows against an
identical order-3 chain -- which this grader VERIFIES rather than assumes (see REFUSAL 5).

  ALWAYS PASS --out. Its default overwrites the committed wiki/glm_beyond_kmer_2026-11-22.json.

  python scripts/glm_beyond_kmer.py --ckpt data/glm_diffusion_ckpt/latest.pt \
      --eval-set orthology_free --out <dir>/run_committed_orthology_free_120x30.json
  python scripts/glm_beyond_kmer.py --ckpt data/glm_ckpt_representative/latest.pt \
      --eval-set orthology_free --out <dir>/run_representative_orthology_free_120x30.json
  python scripts/glm_beyond_kmer.py --ckpt data/glm_ckpt_reseed_control/latest.pt \
      --eval-set orthology_free --out <dir>/run_reseed_seed1_orthology_free_120x30.json

Read-only apart from `--out`. No network, no GPU, no model: it reads three JSON artifacts.
"""
from __future__ import annotations

import argparse
import json
import math
import random
import statistics
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

RULE = ROOT / "wiki" / "glm_representative_repretrain_reading_rule_2026-11-24.json"
SIGN_CORRECTION = ROOT / "wiki" / "glm_representative_repretrain_sign_correction_2026-11-24.json"
OUT = ROOT / "wiki" / "glm_representative_repretrain_result_2026-11-24.json"

#: The frozen geometry. ORDERED, never a set: {30,120} == {120,30} once accepted a 30-window x
#: 120-position run as "the frozen 120x30", and the two numbers are not interchangeable.
FROZEN_GEOMETRY = (120, 30)
FROZEN_SCORING_SEED = 0
FROZEN_SELECTION_SEED = 0
FROZEN_BOOT = 2000
FROZEN_EVAL_SET = "orthology_free"
FINISHED_STEP = 6000

#: The arm being compared, and the null the published self-test is stated against.
PRIMARY_ARM = "left_only"
PUBLISHED_NULL = "markov_ecoli"

#: Columns that do NOT depend on which checkpoint was scored: two Markov chains and uniform. They
#: are the pairing invariant -- see REFUSAL 5.
CHECKPOINT_INDEPENDENT = ("markov_ecoli", "markov_pooled", "uniform")

#: The committed checkpoint's published primary cell, for the self-test. Source:
#: wiki/glm_beyond_kmer_2026-11-22.json, cells["left_only|markov_ecoli"].
PUBLISHED_DELTA = -0.0191
PUBLISHED_CI = [-0.032, -0.0068]

#: The model-side near-ortholog leak, measured BEFORE this arm was launched. Source:
#: wiki/glm_beyond_kmer_two_sets_2026-11-23.json -- delta_model +0.0054 nats, i.e. the model GAINED
#: about this much from windows that had draft near-orthologs. It is the yardstick separating
#: LEAK_WAS_REAL_AND_SMALL from CROSS_ORGANISM_COST_DOMINATES.
LEAK_SIZE_NATS = 0.0054


#: The leak-identification layer. See
#: wiki/glm_repretrain_leak_identification_prereg_2026-11-24.json -- the PRIMARY is scored on a set
#: whose genes have NO draft counterpart, so the leak cannot act there and two of the rule's four
#: branches are not separable from it. These two DISJOINT halves of the withheld pool differ in
#: exactly that one property, so a difference-in-differences across them isolates the leak.
BEARING_SET = "ortholog_bearing"
FREE_SET = "withheld_orthology_free"
EXPECTED_DID_NATS = -0.0054


class Refusal(Exception):
    """A condition under which no number may be reported. Carries its own exit code."""

    def __init__(self, code: int, msg: str):
        super().__init__(msg)
        self.code = code


# --------------------------------------------------------------------------------------------
# the bootstrap, replicated from the scorer and PROVEN faithful by the self-test
# --------------------------------------------------------------------------------------------
def paired_bootstrap(d: list[float], *, b: int, seed: int) -> dict:
    """Percentile bootstrap over paired per-window differences.

    Replicated call-for-call from `scripts/glm_beyond_kmer.py::summarise` -- same RNG, same
    `r.randrange(len(d))` draw pattern, same `int(0.025 * len(boot))` percentile indexing. It is a
    replication rather than an import because the scorer defines it as a closure over its own
    argparse namespace, so there is nothing importable. `self_test` is what makes the replication
    checkable instead of merely claimed: it must reproduce a PUBLISHED number to 4 decimals.
    """
    if not d:
        raise Refusal(7, "empty difference vector; nothing to bootstrap")
    r = random.Random(seed)
    boot = []
    for _ in range(b):
        s = [d[r.randrange(len(d))] for _ in range(len(d))]
        boot.append(statistics.fmean(s))
    boot.sort()
    lo, hi = boot[int(0.025 * len(boot))], boot[int(0.975 * len(boot))]
    return {
        "n_windows": len(d),
        "mean_nats": statistics.fmean(d),
        "ci95": [lo, hi],
        "excludes_zero": bool(lo > 0 or hi < 0),
        "b": b,
        "seed": seed,
    }


def two_sample_bootstrap(x: list[float], y: list[float], *, b: int, seed: int) -> dict:
    """Percentile bootstrap on mean(x) - mean(y) for DISJOINT samples.

    Mirrors scripts/glm_two_set_contrast.py::two_sample_bootstrap. The two eval halves are
    `a - c` and `a & c`, so they share no window and a PAIRED statistic would be wrong here.
    """
    if not x or not y:
        raise Refusal(7, "a leak-identification half is empty; nothing to bootstrap")
    r = random.Random(seed)
    boot = []
    for _ in range(b):
        mx = statistics.fmean([x[r.randrange(len(x))] for _ in range(len(x))])
        my = statistics.fmean([y[r.randrange(len(y))] for _ in range(len(y))])
        boot.append(mx - my)
    boot.sort()
    lo, hi = boot[int(0.025 * len(boot))], boot[int(0.975 * len(boot))]
    return {
        "observed": statistics.fmean(x) - statistics.fmean(y),
        "ci95": [lo, hi],
        "excludes_zero": bool(lo > 0 or hi < 0),
        "n_x": len(x), "n_y": len(y), "b": b, "seed": seed,
    }


def check_arm_identity(arts: dict[str, dict]) -> dict:
    """REFUSAL 8 -- the arms must be DIFFERENT checkpoints, and the right ones.

    Every artifact passes the step-6000 check, so nothing else distinguishes a representative
    artifact from a reseed one. Passing the committed artifact as `--reseed` would yield a floor of
    exactly zero and hand back a RESOLVED verdict on a comparison against itself. Two guards:
    the recorded checkpoint path when the scorer wrote one, and -- always available -- the fact
    that two genuinely different checkpoints cannot produce bit-identical cross-entropies.
    """
    rec: dict = {"ckpt_paths": {}, "ckpt_is_committed": {}}
    for label, art in arts.items():
        a = art["run_args"]
        rec["ckpt_paths"][label] = a.get("ckpt")
        rec["ckpt_is_committed"][label] = a.get("ckpt_is_committed")

    # (a) distinct columns. This is the guard that works even on artifacts predating `ckpt`.
    cols = {label: [w[PRIMARY_ARM] for w in art["per_window"]] for label, art in arts.items()}
    labels = list(cols)
    for i, li in enumerate(labels):
        for lj in labels[i + 1:]:
            if cols[li] == cols[lj]:
                raise Refusal(8, f"arms {li!r} and {lj!r} have BIT-IDENTICAL {PRIMARY_ARM} "
                                 f"columns, so they are the same checkpoint scored twice. A floor "
                                 f"or primary computed against itself is exactly zero and would "
                                 f"read as a clean result.")
    # (b) recorded paths, when present. Absent on artifacts written before the field existed.
    paths = {k: v for k, v in rec["ckpt_paths"].items() if v}
    if len(set(paths.values())) != len(paths):
        raise Refusal(8, f"two arms name the SAME checkpoint path: {paths}")
    flag = rec["ckpt_is_committed"]
    if flag.get("committed") is False:
        raise Refusal(8, "the artifact passed as --committed records ckpt_is_committed=false, so "
                         "it is not the pinned committed checkpoint")
    for label in [k for k in arts if k != "committed"]:
        if flag.get(label) is True:
            raise Refusal(8, f"the artifact passed as --{label} records ckpt_is_committed=true, so "
                             f"it IS the committed checkpoint and cannot be an arm against it")
    rec["distinct_columns_verified"] = True
    rec["path_fields_present"] = bool(paths)
    rec["note_when_absent"] = (
        "artifacts written before glm_beyond_kmer recorded run_args.ckpt carry no path, so arm "
        "identity rests on the bit-identical-column check alone. That check is sufficient to catch "
        "the same checkpoint passed twice, which is the failure mode that matters."
        if not paths else "")
    return rec


def self_test(committed: dict) -> dict:
    """Reproduce the committed checkpoint's PUBLISHED primary cell, or refuse.

    Guards the whole grader: if this passes, the bootstrap, the percentile convention, the seed
    derivation and the per-window column reads are all confirmed against a number published before
    this file existed. If it fails, every downstream CI would be wrong in the same way and silently.
    """
    pw = committed["per_window"]
    args = committed["run_args"]
    d = [w[PUBLISHED_NULL] - w[PRIMARY_ARM] for w in pw]  # scorer orientation: positive = model better
    got = paired_bootstrap(d, b=args["boot"], seed=args["seed"] + 1)
    delta = round(got["mean_nats"], 4)
    ci = [round(got["ci95"][0], 4), round(got["ci95"][1], 4)]
    ok = (delta == PUBLISHED_DELTA) and (ci == PUBLISHED_CI)
    rec = {
        "what": "reproduce the committed checkpoint's published left_only|markov_ecoli cell",
        "published_delta_nats": PUBLISHED_DELTA,
        "recomputed_delta_nats": delta,
        "published_ci95": PUBLISHED_CI,
        "recomputed_ci95": ci,
        "seed_used": args["seed"] + 1,
        "boot_used": args["boot"],
        "passed": ok,
        "why_it_gates_everything": (
            "the grader's bootstrap is a replication of the scorer's closure, not an import. This "
            "test is the only thing that makes the replication checkable. A mismatch means every "
            "confidence interval below is wrong in the same way, so no number is reported."),
    }
    if not ok:
        raise Refusal(6, f"SELF-TEST FAILED: recomputed {delta} {ci} vs published "
                         f"{PUBLISHED_DELTA} {PUBLISHED_CI}. The bootstrap replication is not "
                         f"faithful; no number is reported.")
    return rec


# --------------------------------------------------------------------------------------------
# refusals
# --------------------------------------------------------------------------------------------
def load(path: Path, label: str) -> dict:
    if not path.exists():
        raise Refusal(2, f"missing {label} artifact: {path}")
    try:
        d = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as e:
        raise Refusal(2, f"{label} artifact is not parseable JSON: {path} ({e})") from e
    for k in ("per_window", "run_args", "cells", "checkpoint_step"):
        if k not in d:
            raise Refusal(7, f"{label} artifact lacks `{k}`: {path}")
    return d


def check_geometry(arts: dict[str, dict], *, expect_set: str = FROZEN_EVAL_SET) -> dict:
    """REFUSAL 3 -- the three runs must share the frozen geometry, exactly and in order.

    `expect_set` is the eval set required of all of them. It defaults to the PRIMARY's frozen set;
    the leak-identification layer passes its own two sets, which must still share every other
    geometry field with the primary.
    """
    geoms = {}
    for label, art in arts.items():
        a = art["run_args"]
        geoms[label] = (a["n_windows"], a["positions"], a["seed"], a["selection_seed"],
                        a["boot"], a["eval_set"])
    if len(set(geoms.values())) != 1:
        raise Refusal(3, f"geometry differs across arms, so any difference is a geometry "
                         f"difference and not a corpus difference: {geoms}")
    n_win, n_pos, seed, sel, boot, eval_set = next(iter(geoms.values()))
    if (n_win, n_pos) != FROZEN_GEOMETRY:
        raise Refusal(3, f"geometry is not the frozen {FROZEN_GEOMETRY} (windows, positions): "
                         f"{(n_win, n_pos)}. Swapping the two is a different experiment with the "
                         f"same total number of predictions.")
    if (seed, sel, boot) != (FROZEN_SCORING_SEED, FROZEN_SELECTION_SEED, FROZEN_BOOT):
        raise Refusal(3, f"seeds/boot are not the frozen "
                         f"{(FROZEN_SCORING_SEED, FROZEN_SELECTION_SEED, FROZEN_BOOT)}: "
                         f"{(seed, sel, boot)}")
    if eval_set != expect_set:
        raise Refusal(3, f"eval set is {eval_set!r}, not the expected {expect_set!r}")
    return {"n_windows": n_win, "positions_per_window": n_pos, "scoring_seed": seed,
            "selection_seed": sel, "boot": boot, "eval_set": eval_set,
            "n_predictions_per_arm": n_win * n_pos}


def check_finished(arts: dict[str, dict]) -> dict:
    """REFUSAL 4 -- every checkpoint must be a FINISHED 6000-step run.

    The first representative attempt was killed at step 3500 and its file is preserved on disk. A
    partially-trained checkpoint scores perfectly cleanly and would be read as the arm, so an
    unfinished step is refused rather than noted.
    """
    steps = {label: art["checkpoint_step"] for label, art in arts.items()}
    bad = {k: v for k, v in steps.items() if v != FINISHED_STEP}
    if bad:
        raise Refusal(4, f"these arms were not scored on a FINISHED {FINISHED_STEP}-step "
                         f"checkpoint: {bad}. A partial checkpoint scores cleanly and would be "
                         f"indistinguishable from the arm.")
    return steps


def check_pairing(arts: dict[str, dict]) -> dict:
    """REFUSAL 5 -- the pairing invariant, and the only check that can detect a swapped window set.

    `per_window` rows carry no window identifier, so pairing is BY INDEX. Index pairing is only
    meaningful if row i is the same window with the same sampled positions in all three artifacts.
    The two Markov columns and uniform do not depend on the checkpoint, so they must be identical
    row by row. Any mismatch means the windows, the sampled positions or the fitted chain differ --
    for instance a run given a different `--corpus`, which would refit the E. coli chain. That
    would make every paired difference a comparison of different things.
    """
    labels = list(arts)
    ref = labels[0]
    n = len(arts[ref]["per_window"])
    for label, art in arts.items():
        if len(art["per_window"]) != n:
            raise Refusal(5, f"{label} has {len(art['per_window'])} windows, {ref} has {n}")
    mismatches = []
    for col in CHECKPOINT_INDEPENDENT:
        base = [w[col] for w in arts[ref]["per_window"]]
        for label in labels[1:]:
            other = [w[col] for w in arts[label]["per_window"]]
            diffs = [i for i, (x, y) in enumerate(zip(base, other)) if x != y]
            if diffs:
                mismatches.append({"column": col, "vs": label, "n_rows_differing": len(diffs),
                                   "first_rows": diffs[:5]})
    if mismatches:
        raise Refusal(5, "PAIRING INVARIANT VIOLATED -- checkpoint-independent columns differ "
                         f"across arms, so row i is not the same window everywhere: {mismatches}")
    return {
        "columns_verified_identical": list(CHECKPOINT_INDEPENDENT),
        "n_windows": n,
        "what_this_proves": (
            "row i is the same eval window at the same sampled positions against the same fitted "
            "order-3 chain in all three arms, so an index-paired difference isolates the MODEL."),
        "what_it_does_not_prove": (
            "it cannot detect a change that moves the model columns alone, which is exactly what "
            "the arm is measuring and therefore cannot be controlled for here."),
    }


# --------------------------------------------------------------------------------------------
# the graded reading
# --------------------------------------------------------------------------------------------
def identify_leak(bearing: dict[str, dict], free: dict[str, dict], *, b: int, seed: int) -> dict:
    """The difference-in-differences that separates the leak from the cross-organism cost.

    `bearing` and `free` each map {committed, representative, reseed_seed1} -> run artifact for that
    eval set. Pre-registered in wiki/glm_repretrain_leak_identification_prereg_2026-11-24.json.
    """
    for label, group in (("bearing", bearing), ("free", free)):
        sets = {g["run_args"]["eval_set"] for g in group.values()}
        want = BEARING_SET if label == "bearing" else FREE_SET
        if sets != {want}:
            raise Refusal(9, f"the {label} group must all be eval_set={want!r}, got {sets}")
        check_geometry(group, expect_set=want)
        check_finished(group)
        check_pairing(group)
        check_arm_identity(group)

    def deltas(group: dict[str, dict], arm: str) -> list[float]:
        c = [w[PRIMARY_ARM] for w in group["committed"]["per_window"]]
        o = [w[PRIMARY_ARM] for w in group[arm]["per_window"]]
        return [ci - oi for ci, oi in zip(c, o)]

    out: dict = {
        "quantity": ("DiD = mean paired delta on ortholog_bearing minus mean paired delta on "
                     "withheld_orthology_free; each delta is committed - assessed, so positive = "
                     "the assessed model is better"),
        "why_it_identifies_the_leak": (
            "the two halves are DISJOINT and differ in exactly one property -- whether the window's "
            "gene has a draft counterpart in training. Removing the draft removes the near-ortholog "
            "from the bearing half only, while cross-organism cost acts on both and cancels."),
        "expected_did_nats": EXPECTED_DID_NATS,
        "bearing_set": BEARING_SET, "free_set": FREE_SET,
    }
    for arm, key in (("representative", "primary"), ("reseed_seed1", "floor")):
        bd, fd = deltas(bearing, arm), deltas(free, arm)
        boot = two_sample_bootstrap(bd, fd, b=b, seed=seed)
        out[key] = {
            "bearing_mean_delta_nats": round(statistics.fmean(bd), 6),
            "free_mean_delta_nats": round(statistics.fmean(fd), 6),
            "did_nats": round(boot["observed"], 6),
            "ci95": [round(boot["ci95"][0], 6), round(boot["ci95"][1], 6)],
            "excludes_zero": boot["excludes_zero"],
            "n_bearing": boot["n_x"], "n_free": boot["n_y"],
        }

    p, f = out["primary"], out["floor"]
    exceeds = abs(p["did_nats"]) > abs(f["did_nats"])
    out["gate"] = {
        "rule": "|DiD_primary| must exceed |DiD_floor| AND the DiD CI must exclude zero",
        "primary_magnitude": abs(p["did_nats"]),
        "floor_magnitude": abs(f["did_nats"]),
        "exceeds_floor": exceeds,
        "ci_excludes_zero": p["excludes_zero"],
    }
    if not (exceeds and p["excludes_zero"]):
        out["verdict"] = "LEAK_NOT_IDENTIFIABLE_AT_THIS_SAMPLE_SIZE"
        out["reading"] = ("the leak-specific term does not clear its own noise floor, so no leak "
                          "size is quoted and the primary's branch keeps "
                          "attribution_identifiable=false.")
        return out
    out["verdict"] = "LEAK_IDENTIFIED"
    if p["did_nats"] < 0:
        out["reading"] = (
            f"the representative model lost {-p['did_nats']:.6f} nats MORE on the ortholog-bearing "
            f"half than on the orthology-free half, CI {p['ci95']}. Cross-organism cost cancels in "
            f"this difference, so this is the leak-specific term. The decomposition predicted "
            f"{EXPECTED_DID_NATS} and that prediction was frozen before any number was read.")
        out["consistent_with_the_decomposition"] = bool(
            p["ci95"][0] <= EXPECTED_DID_NATS <= p["ci95"][1])
    else:
        out["reading"] = (
            f"the representative model lost LESS on the ortholog-bearing half (DiD "
            f"{p['did_nats']:+.6f}, CI {p['ci95']}), which is the OPPOSITE of the pre-registered "
            f"direction. The committed model's bearing-half advantage was therefore not caused by "
            f"the draft near-ortholog, which is a finding about the +0.0054 DECOMPOSITION rather "
            f"than a rescue of the model.")
        out["contradicts_pre_registered_direction"] = True
    return out


def grade(primary: dict, floor: dict, *, n_floor_draws: int = 1,
          attribution: dict | None = None) -> dict:
    """Apply the frozen gate, the pre-registered escalation, and the four named branches.

    ORDER IS PRE-REGISTERED and first match wins: escalation, then the gate, then the branch. The
    escalation criterion is a precondition on reporting ANY branch, so it dominates -- a result can
    clear the gate and still be escalated.

    `n_floor_draws` is how many reseed checkpoints the floor was built from. The rule pre-registered
    "run seeds 2 and 3, then re-grade" WITHOUT saying how several draws should gate, so that is fixed
    here, before any number: the gate uses the LARGEST |floor| across draws (the most conservative
    available reading), and the escalation criterion RETIRES at 3 or more draws because its stated
    purpose -- turning a point estimate into a distribution -- has then been served.
    """
    p_mean, f_mean = primary["mean_nats"], floor["mean_nats"]
    p_mag, f_mag = abs(p_mean), abs(f_mean)
    exceeds_floor = p_mag > f_mag
    excludes_zero = primary["excludes_zero"]

    rec = {
        "primary_mean_nats": round(p_mean, 6),
        "primary_ci95": [round(primary["ci95"][0], 6), round(primary["ci95"][1], 6)],
        "primary_magnitude_nats": round(p_mag, 6),
        "floor_mean_nats": round(f_mean, 6),
        "floor_magnitude_nats": round(f_mag, 6),
        "primary_exceeds_floor": exceeds_floor,
        "primary_ci_excludes_zero": excludes_zero,
        "floor_is_a_point_estimate_not_a_null": (
            "a single reseed is ONE draw of trajectory noise. The rule says so and this grader does "
            "not dress it up as a calibrated null."),
        "direction": ("representative_BETTER" if p_mean > 0 else
                      "representative_WORSE" if p_mean < 0 else "exactly_equal"),
    }

    # ---- 1. escalation, checked FIRST because it gates reporting a branch at all ----
    rec["n_floor_draws"] = n_floor_draws
    rec["floor_selection_rule"] = (
        "the LARGEST |floor| across reseed draws, the most conservative available reading"
        if n_floor_draws > 1 else "the single available floor draw")
    within_2x = (p_mag <= 2.0 * f_mag) and n_floor_draws < 3
    rec["escalation_criterion"] = {
        "rule": ("|primary| <= 2x |floor| requires seeds 2 and 3 before ANY branch; it RETIRES at "
                 ">=3 floor draws, because its purpose is to turn a point estimate into a "
                 "distribution and 3 draws have done that"),
        "threshold_nats": round(2.0 * f_mag, 6),
        "within_2x_of_floor": bool(p_mag <= 2.0 * f_mag),
        "retired_by_draw_count": n_floor_draws >= 3,
        "triggered": bool(within_2x),
    }
    if within_2x:
        rec["verdict"] = "ESCALATION_REQUIRED_ADDITIONAL_RESEEDS"
        rec["branch"] = None
        rec["reading"] = (
            f"|primary| {p_mag:.6f} nats sits within 2x the single floor draw "
            f"({2.0 * f_mag:.6f} nats), so the rule's pre-registered escalation fires and NO branch "
            f"is reported. Seeds 2 and 3 on the COMMITTED corpus must be run to turn the floor from "
            f"a point estimate into a distribution. Deciding to escalate on this criterion was "
            f"pre-registered before any number existed, so it is not a post-hoc rescue.")
        rec["next_action"] = ("two more pretrains on wiki/glm_training_corpus_2026-11-21.txt with "
                             "--seed 2 and --seed 3, everything else identical, then re-grade")
        return rec

    # ---- 2. the gate ----
    if not (exceeds_floor and excludes_zero):
        why = []
        if not exceeds_floor:
            why.append(f"|primary| {p_mag:.6f} does not exceed the floor {f_mag:.6f}")
        if not excludes_zero:
            why.append(f"the primary CI {rec['primary_ci95']} includes zero")
        rec["verdict"] = "UNRESOLVED_WITHIN_TRAJECTORY_NOISE"
        rec["branch"] = "UNRESOLVED_WITHIN_TRAJECTORY_NOISE"
        rec["reading"] = ("the gate is not met (" + "; and ".join(why) +
                          "), so no leak size is quoted. A difference smaller than trajectory "
                          "noise is not evidence about the corpus.")
        return rec

    # ---- 3. the branch ----
    rec["verdict"] = "RESOLVED"

    # The PRIMARY is scored on a set whose genes have NO draft counterpart, so the leak cannot act
    # there and LEAK_WAS_REAL_AND_SMALL vs CROSS_ORGANISM_COST_DOMINATES is NOT decidable from it.
    # Stamped rather than silently chosen. See the leak-identification pre-registration.
    att = {
        "leak_vs_cross_organism_identifiable_from_the_primary_set": False,
        "why": ("the primary's eval set `orthology_free` is built as held-out MG1655 windows whose "
                "gene symbol has NO counterpart in the draft assembly, so there is no "
                "near-ortholog there to lose. A worse model on that set is cross-organism cost "
                "plus trajectory, not the leak."),
        "the_identifying_measurement": ("paired difference-in-differences on ortholog_bearing minus "
                                        "withheld_orthology_free, pre-registered in "
                                        "wiki/glm_repretrain_leak_identification_prereg_2026-11-24.json"),
        "supplied": attribution is not None,
    }
    if attribution is not None:
        att["verdict"] = attribution.get("verdict")
        att["did_nats"] = attribution.get("primary", {}).get("did_nats")
        att["ci95"] = attribution.get("primary", {}).get("ci95")
    rec["attribution"] = att

    if p_mean > 0:
        rec["branch"] = "LEAK_WAS_NOT_HELPING_THE_MODEL"
        rec["reading"] = (
            f"the representative model is BETTER by {p_mean:.6f} nats, outside the noise floor and "
            f"with a CI excluding zero. The near-ortholog windows were not carrying the model's "
            f"performance. This CONTRADICTS the pre-registered direction, which expected a worse "
            f"model, and that is recorded rather than smoothed over.")
        rec["contradicts_pre_registered_direction"] = True
        return rec

    # representative is WORSE. Separate `consistent with the measured leak` from `substantially
    # more` using the primary's own CI against the leak size. Operationalised HERE, before any
    # number, because the rule names the two branches without giving a numeric boundary.
    cost = -p_mean
    cost_ci = [-rec["primary_ci95"][1], -rec["primary_ci95"][0]]  # flip to a positive cost
    rec["cost_nats"] = round(cost, 6)
    rec["cost_ci95"] = [round(cost_ci[0], 6), round(cost_ci[1], 6)]
    rec["leak_size_nats"] = LEAK_SIZE_NATS
    rec["how_consistent_was_operationalised"] = {
        "decision": ("the measured leak is `consistent` iff it lies INSIDE the primary's own "
                     "bootstrap CI expressed as a positive cost; `exceeds` iff the whole CI lies "
                     "above it; `falls_short` iff the whole CI lies below it"),
        "why_the_CI_and_not_a_hand_picked_ratio": (
            "the rule names LEAK_WAS_REAL_AND_SMALL and CROSS_ORGANISM_COST_DOMINATES without a "
            "numeric boundary between them. Using the interval the data already produces avoids "
            "inventing a threshold, and it uses only quantities frozen before the run: the "
            "+0.0054 leak decomposition and the pre-registered bootstrap."),
        "frozen_before_any_number": True,
    }
    if cost_ci[0] > LEAK_SIZE_NATS:
        rec["branch"] = "CROSS_ORGANISM_COST_DOMINATES"
        rec["magnitude_vs_prediction"] = "exceeds"
        rec["branch_fit"] = "CLEAN"
        rec["reading"] = (
            f"the representative model is worse by {cost:.6f} nats, CI {rec['cost_ci95']}, whose "
            f"whole interval lies ABOVE the {LEAK_SIZE_NATS} nat leak the decomposition measured. "
            f"Removing ALL E. coli cost same-genus generalisation, which is confounded with the "
            f"leak, so this run does NOT cleanly price the leak. The disambiguating follow-up is a "
            f"Klebsiella-substituted corpus.")
    elif cost_ci[1] < LEAK_SIZE_NATS:
        # The rule has no branch for `worse, but by significantly LESS than predicted`. Report the
        # nearest frozen branch and mark the fit PARTIAL rather than invent a fifth branch or
        # pretend a mismatch is a match.
        rec["branch"] = "LEAK_WAS_REAL_AND_SMALL"
        rec["magnitude_vs_prediction"] = "falls_short"
        rec["branch_fit"] = "PARTIAL"
        rec["reading"] = (
            f"the representative model is worse by {cost:.6f} nats, CI {rec['cost_ci95']}, whose "
            f"whole interval lies BELOW the {LEAK_SIZE_NATS} nat leak the decomposition measured. "
            f"The leak was real and outside the noise floor, so LEAK_WAS_REAL_AND_SMALL is the "
            f"nearest frozen branch -- but its wording says `consistent with the measured +0.0054`, "
            f"and this is significantly smaller. The fit is marked PARTIAL rather than reported as "
            f"a match, and no fifth branch is invented after the fact.")
        rec["what_this_would_mean"] = (
            "the +0.0054 decomposition over-estimated how much the model was gaining from the "
            "near-ortholog windows, OR the lost 5.4% of corpus was partly compensated. Both are "
            "hypotheses, not findings, and neither is settled by this arm.")
    else:
        rec["branch"] = "LEAK_WAS_REAL_AND_SMALL"
        rec["magnitude_vs_prediction"] = "consistent"
        rec["branch_fit"] = "CLEAN"
        rec["reading"] = (
            f"the representative model is worse by {cost:.6f} nats, CI {rec['cost_ci95']}, which "
            f"CONTAINS the {LEAK_SIZE_NATS} nat leak the decomposition measured independently and "
            f"before this arm ran. The beyond-k-mer negative stands and is now measured on a corpus "
            f"that cannot leak at all.")
    rec["contradicts_pre_registered_direction"] = False

    # The two "worse" branches are the ones the primary set cannot separate. A CLEAN fit is only
    # honest once the identifying measurement has actually resolved the leak.
    if rec["branch"] in ("LEAK_WAS_REAL_AND_SMALL", "CROSS_ORGANISM_COST_DOMINATES"):
        if attribution is None or attribution.get("verdict") != "LEAK_IDENTIFIED":
            rec["branch_fit"] = "PARTIAL_ATTRIBUTION_NOT_IDENTIFIED"
            rec["reading"] += (
                " ATTRIBUTION NOT ESTABLISHED: this branch names the leak as the cause, but the "
                "primary's set excludes draft counterparts by construction, so the leak cannot "
                "act there and this branch is not separable from the other `worse` branch on this "
                "evidence. The branch is reported as the rule defines it and its fit is marked "
                "PARTIAL rather than presented as a settled attribution.")
    return rec


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--committed", type=Path, required=True,
                    help="run artifact for the committed checkpoint (the baseline arm)")
    ap.add_argument("--representative", type=Path, required=True,
                    help="run artifact for the representative checkpoint (the arm)")
    ap.add_argument("--reseed", type=Path, required=True, nargs="+",
                    help="run artifact(s) for the committed corpus at seed 1 (and optionally seeds "
                         "2, 3 ...) -- the noise floor. Several draws gate via the LARGEST |floor|.")
    ap.add_argument("--out", type=Path, default=OUT)
    for grp in ("bearing", "free"):
        for arm in ("committed", "representative", "reseed"):
            ap.add_argument(f"--{grp}-{arm}", type=Path, default=None,
                            help=f"OPTIONAL leak-identification layer: {arm} artifact on the "
                                 f"{BEARING_SET if grp == 'bearing' else FREE_SET} set. All SIX "
                                 f"must be given together or none at all.")
    a = ap.parse_args(argv)

    leak_args = {f"{g}_{m}": getattr(a, f"{g}_{m}")
                 for g in ("bearing", "free") for m in ("committed", "representative", "reseed")}
    n_leak = sum(v is not None for v in leak_args.values())
    if n_leak not in (0, 6):
        missing = [k for k, v in leak_args.items() if v is None]
        print(f"REFUSING: the leak-identification layer needs all SIX artifacts or none; "
              f"{n_leak} given, missing {missing}", file=sys.stderr)
        return 9

    if not RULE.exists():
        print(f"REFUSING: the frozen rule is missing: {RULE}", file=sys.stderr)
        return 2

    try:
        arts = {
            "committed": load(a.committed, "committed"),
            "representative": load(a.representative, "representative"),
        }
        for i, p in enumerate(a.reseed):
            arts[f"reseed_{i + 1}"] = load(p, f"reseed_{i + 1}")
        geometry = check_geometry(arts)
        steps = check_finished(arts)
        pairing = check_pairing(arts)
        identity = check_arm_identity(arts)
        st = self_test(arts["committed"])

        c = [w[PRIMARY_ARM] for w in arts["committed"]["per_window"]]
        r = [w[PRIMARY_ARM] for w in arts["representative"]["per_window"]]

        # SIGN: committed - assessed, so positive = the assessed model has LOWER cross-entropy
        # = BETTER, exactly as the frozen sign convention states. See the sign-correction artifact.
        boot_seed = arts["committed"]["run_args"]["seed"] + 1
        boot_n = arts["committed"]["run_args"]["boot"]
        primary = paired_bootstrap([ci - ri for ci, ri in zip(c, r)], b=boot_n, seed=boot_seed)

        floors = {}
        for label in [k for k in arts if k.startswith("reseed_")]:
            s = [w[PRIMARY_ARM] for w in arts[label]["per_window"]]
            floors[label] = paired_bootstrap([ci - si for ci, si in zip(c, s)],
                                             b=boot_n, seed=boot_seed)
        # MOST CONSERVATIVE available floor, fixed before any number was read.
        worst = max(floors, key=lambda k: abs(floors[k]["mean_nats"]))
        floor = floors[worst]

        attribution = None
        if n_leak == 6:
            bearing = {m if m != "reseed" else "reseed_seed1":
                       load(leak_args[f"bearing_{m}"], f"bearing_{m}")
                       for m in ("committed", "representative", "reseed")}
            free = {m if m != "reseed" else "reseed_seed1":
                    load(leak_args[f"free_{m}"], f"free_{m}")
                    for m in ("committed", "representative", "reseed")}
            attribution = identify_leak(bearing, free, b=boot_n, seed=boot_seed)

        graded = grade(primary, floor, n_floor_draws=len(floors), attribution=attribution)
    except Refusal as e:
        print(f"REFUSING: {e}", file=sys.stderr)
        return e.code

    rec = {
        "record": "glm-representative-repretrain-result-v1",
        "date": __import__("time").strftime("%Y-%m-%d"),
        "graded_against": str(RULE.relative_to(ROOT)).replace("\\", "/"),
        "sign_convention_corrected_by": str(SIGN_CORRECTION.relative_to(ROOT)).replace("\\", "/"),
        "grader_frozen_before_any_number": (
            "scripts/glm_repretrain_grade.py was committed while both pretrains were still "
            "running, so the gate, the escalation threshold and the branch boundary were all fixed "
            "before any representative or reseed cross-entropy existed. Verifiable by commit order."),
        "inputs": {k: str(p).replace("\\", "/") for k, p in
                   (("committed", a.committed), ("representative", a.representative),
                    ("reseed_seed1", a.reseed))},
        "checkpoint_steps": steps,
        "geometry": geometry,
        "REFUSALS_THAT_PASSED": {
            "geometry_identical_and_frozen": True,
            "all_checkpoints_finished_at_6000": True,
            "pairing_invariant": pairing,
            "arm_identity": identity,
            "bootstrap_self_test": st,
        },
        "SIGN": {
            "implemented": "committed.left_only - assessed.left_only",
            "positive_means": "the assessed model has LOWER cross-entropy = BETTER",
            "why_not_the_rule_s_literal_subtraction_order": (
                "the rule's `quantity` and its `sign_convention` specify opposite subtractions. The "
                "convention states the MEANING and is followed; the correction artifact records the "
                "evidence and shows the choice cannot move the gate or the branch."),
        },
        "PRIMARY": {
            "quantity": ("paired per-window change in model cross-entropy, representative vs "
                         "committed, on identical windows at identical positions"),
            "bootstrap": primary,
        },
        "NOISE_FLOOR": {
            "quantity": "the same quantity with a committed-corpus reseed as the second arm",
            "n_draws": len(floors),
            "selected": worst,
            "selection_rule": "largest |mean| across draws, the most conservative available",
            "bootstrap": floor,
            "all_draws": floors,
        },
        "LEAK_IDENTIFICATION": attribution if attribution is not None else {
            "supplied": False,
            "why_it_matters": (
                "without it, LEAK_WAS_REAL_AND_SMALL and CROSS_ORGANISM_COST_DOMINATES are not "
                "separable, because the primary's eval set excludes draft counterparts by "
                "construction and the leak cannot act there."),
            "how_to_supply_it": ("score ortholog_bearing and withheld_orthology_free for the "
                                 "representative and reseed checkpoints and pass the six "
                                 "--bearing-* / --free-* artifacts; the committed pair is already "
                                 "measured at wiki/glm_beyond_kmer_runs_committed_2026-11-24/"),
            "prereg": "wiki/glm_repretrain_leak_identification_prereg_2026-11-24.json",
        },
        "GRADED": graded,
        "honest_limits": [
            "The noise floor is ONE draw of trajectory noise, not a distribution. The rule says so.",
            "The representative corpus has NO E. coli, so this model is a cross-ORGANISM generator "
            "and this arm produces no same-strain beyond-k-mer verdict. That needs an additional "
            "fine-tune on MG1655's own 991 trainable windows, which is not part of this arm.",
            "Removing the draft removes BOTH the leak and the only same-genus genome, so the leak "
            "and the cross-organism cost are confounded by construction. Only the "
            "CROSS_ORGANISM_COST_DOMINATES branch is diagnostic of that confound; the other "
            "branches do not separate them.",
            "A closed verbatim channel alongside a model much worse than +0.0054 is COHERENT, not "
            "contradictory: the companion froze homology-mediated transfer as a second channel that "
            "a longest-substring measure is structurally blind to.",
        ],
    }
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(rec, indent=2), encoding="utf-8")

    g = graded
    print(f"self-test: reproduced published {PUBLISHED_DELTA} {PUBLISHED_CI} -> PASS")
    print(f"pairing:   {pairing['n_windows']} windows, "
          f"{len(CHECKPOINT_INDEPENDENT)} checkpoint-independent columns identical")
    print(f"primary:   {g['primary_mean_nats']:+.6f} nats  CI {g['primary_ci95']}  "
          f"({g['direction']})")
    print(f"floor:     {g['floor_mean_nats']:+.6f} nats  |floor| {g['floor_magnitude_nats']:.6f}"
          f"  ({len(floors)} draw(s), using {worst})")
    if attribution is not None:
        ap_ = attribution["primary"]
        print(f"leak DiD:  {ap_['did_nats']:+.6f} nats  CI {ap_['ci95']}  -> "
              f"{attribution['verdict']}")
    else:
        print("leak DiD:  not supplied -> leak vs cross-organism cost NOT identifiable")
    print(f"verdict:   {g['verdict']}")
    print(f"branch:    {g['branch']}"
          + (f"  [fit {g['branch_fit']}]" if g.get("branch_fit") else ""))
    print(f"-> {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
