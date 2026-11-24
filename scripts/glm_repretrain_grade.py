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


def check_geometry(arts: dict[str, dict]) -> dict:
    """REFUSAL 3 -- the three runs must share the frozen geometry, exactly and in order."""
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
    if eval_set != FROZEN_EVAL_SET:
        raise Refusal(3, f"eval set is {eval_set!r}, not the frozen {FROZEN_EVAL_SET!r}")
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
def grade(primary: dict, floor: dict) -> dict:
    """Apply the frozen gate, the pre-registered escalation, and the four named branches.

    ORDER IS PRE-REGISTERED and first match wins: escalation, then the gate, then the branch. The
    escalation criterion is a precondition on reporting ANY branch, so it dominates -- a result can
    clear the gate and still be escalated.
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
    within_2x = p_mag <= 2.0 * f_mag
    rec["escalation_criterion"] = {
        "rule": "|primary| <= 2x |single floor draw| requires seeds 2 and 3 before ANY branch",
        "threshold_nats": round(2.0 * f_mag, 6),
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
    return rec


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--committed", type=Path, required=True,
                    help="run artifact for the committed checkpoint (the baseline arm)")
    ap.add_argument("--representative", type=Path, required=True,
                    help="run artifact for the representative checkpoint (the arm)")
    ap.add_argument("--reseed", type=Path, required=True,
                    help="run artifact for the committed corpus at seed 1 (the noise floor)")
    ap.add_argument("--out", type=Path, default=OUT)
    a = ap.parse_args(argv)

    if not RULE.exists():
        print(f"REFUSING: the frozen rule is missing: {RULE}", file=sys.stderr)
        return 2

    try:
        arts = {
            "committed": load(a.committed, "committed"),
            "representative": load(a.representative, "representative"),
            "reseed_seed1": load(a.reseed, "reseed_seed1"),
        }
        geometry = check_geometry(arts)
        steps = check_finished(arts)
        pairing = check_pairing(arts)
        st = self_test(arts["committed"])

        c = [w[PRIMARY_ARM] for w in arts["committed"]["per_window"]]
        r = [w[PRIMARY_ARM] for w in arts["representative"]["per_window"]]
        s1 = [w[PRIMARY_ARM] for w in arts["reseed_seed1"]["per_window"]]

        # SIGN: committed - assessed, so positive = the assessed model has LOWER cross-entropy
        # = BETTER, exactly as the frozen sign convention states. See the sign-correction artifact.
        boot_seed = arts["committed"]["run_args"]["seed"] + 1
        boot_n = arts["committed"]["run_args"]["boot"]
        primary = paired_bootstrap([ci - ri for ci, ri in zip(c, r)], b=boot_n, seed=boot_seed)
        floor = paired_bootstrap([ci - si for ci, si in zip(c, s1)], b=boot_n, seed=boot_seed)
        graded = grade(primary, floor)
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
            "quantity": "the same quantity with the committed corpus at seed 1 as the second arm",
            "bootstrap": floor,
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
    print(f"floor:     {g['floor_mean_nats']:+.6f} nats  |floor| {g['floor_magnitude_nats']:.6f}")
    print(f"verdict:   {g['verdict']}")
    print(f"branch:    {g['branch']}"
          + (f"  [fit {g['branch_fit']}]" if g.get("branch_fit") else ""))
    print(f"-> {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
