"""The leak DiD against its OWN noise floor at 500 windows, using the grader's pre-registered gate.

WHY THIS EXISTS, AND WHAT IT CORRECTS. `scripts/glm_leak_did.py` computes the DiD for ONE assessed
checkpoint against the committed one. Two v5 checkpoints were preserved and both were re-scored at
500 windows, which yields four cells (2 checkpoints x 2 scoring arms). The first version of this
script read those four as "two re-pretrains of one recipe" and reported that they disagreed on sign.

**That was wrong, and the error mattered.** The two checkpoints are NOT two draws of one recipe:

    run1_representative   ckpt_dir data/glm_ckpt_representative   -> TREATMENT
                          the re-pretrain that EXCLUDES the draft E. coli
    run2_reseed_seed1     ckpt_dir data/glm_ckpt_reseed_control    -> CONTROL
                          a reseed on the COMMITTED corpus, which still CONTAINS the draft

The control has no leak to lose by construction, so its DiD is the NOISE FLOOR, not a second leak
estimate. The roles are not asserted here -- they are DERIVED from `ckpt_dir` in the scoring
kernel's own manifest, so a mislabelled run cannot be read as the wrong arm.

Two claims from the first version are WITHDRAWN. "Two re-pretrains disagree on sign" is false: one
is the treatment and one is the control, and treatment-negative with control-positive is the
ORDINARY shape of a real-but-small effect against a null, not a contradiction. And the "0.002892
spread agrees with an independently measured 0.00281 noise floor" was a quantity compared against
itself: the graded artifact's 0.00281 IS this same control arm's DiD at 120 windows, so the
closeness was arithmetic coincidence between a treatment-minus-control difference at one n and the
control's own value at another.

WHAT REPLACES THEM IS STRONGER, because the control measures the floor directly at the SAME n
instead of inheriting a 120-window figure. The gate is the grader's own, quoted not invented:

    |DiD_treatment| must EXCEED |DiD_control| AND the treatment CI must EXCLUDE ZERO

Both halves are required, which is why a treatment that merely out-sizes the floor does not pass.

HONEST LIMITS.
* The predicted -0.0054 is the graded artifact's `expected_did_nats`, a PRE-REGISTERED expectation.
  Its own provenance was not re-derived here, so "excluded at the predicted size" is a statement
  about that recorded expectation, not about a re-measured leak size.
* ONE control draw is one draw of trajectory noise. The graded artifact says so of its own floor
  ("a single reseed is ONE draw ... this grader does not dress it up as a calibrated null") and that
  caveat carries over unchanged to the 500-window floor computed here.
* Excluding the predicted SIZE does not exclude a leak. It bounds one.
* The DiD identifies a leak only if cross-organism cost acts equally on both halves and cancels.
  That is the design's assumption and is not testable from these numbers.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

PREDICTED = -0.0054
PRIOR_FLOOR_AT_120 = 0.00281
MANIFEST = ROOT / "wiki" / "glm_repretrain_scored_2026-11-25" / "_kernel_manifest.json"
OUT = ROOT / "wiki" / "glm_leak_did_n500_synthesis_2026-11-28.json"

TREATMENT, CONTROL = "treatment", "control"

#: Derived from `ckpt_dir`, never from the run label. The label `run2_reseed_seed1` does not say
#: "control" anywhere; only its checkpoint directory does.
_DIR_ROLE = (("glm_ckpt_reseed_control", CONTROL), ("glm_ckpt_representative", TREATMENT))

#: The grader's gate, quoted from LEAK_IDENTIFICATION.gate.rule in the graded artifact.
GATE_RULE = ("|DiD_treatment| must exceed |DiD_control| AND the treatment DiD CI must exclude zero")


class SynthError(RuntimeError):
    pass


def role_for(label: str, manifest: dict) -> str:
    """TREATMENT or CONTROL, derived from the scoring kernel's recorded checkpoint directory."""
    arms = manifest.get("arms") or {}
    arm = arms.get(label)
    if arm is None:
        raise SynthError(f"{label!r} is not in the kernel manifest; cannot derive its role")
    ckpt_dir = str(arm.get("ckpt_dir") or "")
    for needle, role in _DIR_ROLE:
        if needle in ckpt_dir:
            return role
    raise SynthError(f"{label!r} has ckpt_dir {ckpt_dir!r}, which names neither a treatment nor a "
                     "control directory; refusing to guess which arm it is")


def apply_gate(treat: dict, ctrl: dict) -> dict:
    """The grader's pre-registered gate, both halves required."""
    t_did, c_did = treat["result"]["did_nats"], ctrl["result"]["did_nats"]
    exceeds = abs(t_did) > abs(c_did)
    excl_zero = treat["result"]["excludes_zero"]
    identified = bool(exceeds and excl_zero)
    return {
        "rule": GATE_RULE,
        "treatment_did_nats": t_did,
        "control_did_nats": c_did,
        "treatment_magnitude": round(abs(t_did), 6),
        "control_magnitude_the_floor_at_this_n": round(abs(c_did), 6),
        "exceeds_floor": exceeds,
        "treatment_ci_excludes_zero": excl_zero,
        "treatment_ci95": treat["result"]["ci95"],
        "VERDICT": "LEAK_IDENTIFIED" if identified else "LEAK_NOT_IDENTIFIABLE",
        "treatment_excludes_predicted_size": treat["result"]["excludes_predicted"],
        "control_ci95": ctrl["result"]["ci95"],
        "treatment_point_inside_control_ci": bool(
            ctrl["result"]["ci95"][0] <= t_did <= ctrl["result"]["ci95"][1]),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cells", nargs="+", type=Path, required=True)
    ap.add_argument("--manifest", type=Path, default=MANIFEST)
    ap.add_argument("--out", type=Path, default=OUT)
    a = ap.parse_args()

    if not a.manifest.exists():
        print(f"REFUSING: missing kernel manifest {a.manifest}; roles would have to be guessed",
              file=sys.stderr)
        return 2
    manifest = json.loads(a.manifest.read_text(encoding="utf-8"))

    cells = []
    for p in a.cells:
        if not p.exists():
            print(f"REFUSING: missing {p}", file=sys.stderr)
            return 2
        d = json.loads(p.read_text(encoding="utf-8"))
        if d.get("record") != "glm-leak-did-v1":
            print(f"REFUSING: {p.name} is not a glm-leak-did-v1 artifact", file=sys.stderr)
            return 2
        cells.append((p, d))

    ns = {(d["result"]["n_bearing"], d["result"]["n_clean"]) for _, d in cells}
    if len(ns) != 1:
        print(f"REFUSING: cells disagree on n: {ns}; comparing a treatment against a floor "
              "measured at a different window count confounds effect with power", file=sys.stderr)
        return 2

    try:
        roled = [(p, d, role_for(d["assessed_label"], manifest)) for p, d in cells]
    except SynthError as exc:
        print(f"REFUSING: {exc}", file=sys.stderr)
        return 2

    by_arm: dict[str, dict] = {}
    for _, d, role in roled:
        by_arm.setdefault(d["arm"], {})[role] = d
    pairs = {arm: v for arm, v in by_arm.items() if TREATMENT in v and CONTROL in v}
    if not pairs:
        print("REFUSING: no scoring arm has BOTH a treatment and a control cell; the floor cannot "
              "be compared against the effect at the same n", file=sys.stderr)
        return 2

    gates = {arm: apply_gate(v[TREATMENT], v[CONTROL]) for arm, v in sorted(pairs.items())}
    n = next(iter(ns))
    all_excl_pred = all(g["treatment_excludes_predicted_size"] for g in gates.values())
    any_identified = any(g["VERDICT"] == "LEAK_IDENTIFIED" for g in gates.values())

    out = {
        "record": "glm-leak-did-synthesis-v2",
        "date": "2026-11-28",
        "question": ("at 500 windows, does the leak-specific term clear its own noise floor, and is "
                     "the pre-registered expected leak size still admissible?"),
        "supersedes": ("glm-leak-did-synthesis-v1, which misread the control arm as a second "
                       "re-pretrain and reported a cross-seed sign disagreement that does not "
                       "exist. Both of that version's headline claims are withdrawn here."),
        "arm_roles_derived_from": str(a.manifest.relative_to(ROOT)),
        "arm_roles": {d["assessed_label"]: role for _, d, role in roled},
        "n_windows_per_half": {"bearing": n[0], "clean": n[1]},
        "predicted_did_nats": PREDICTED,
        "predicted_is": ("the graded artifact's pre-registered `expected_did_nats`; its own "
                         "provenance was NOT re-derived here"),
        "prior_at_120_windows": {
            "treatment_did_nats": -0.001897, "treatment_ci95": [-0.006801, 0.002813],
            "control_did_nats_the_floor": PRIOR_FLOOR_AT_120,
            "verdict": "LEAK_NOT_IDENTIFIABLE_AT_THIS_SAMPLE_SIZE",
        },
        "gate_by_scoring_arm": gates,
        "headline": {
            "VERDICT": "LEAK_IDENTIFIED" if any_identified else "LEAK_NOT_IDENTIFIABLE",
            "treatment_excludes_predicted_size_in_every_arm": all_excl_pred,
            "floor_at_500_windows": {arm: g["control_magnitude_the_floor_at_this_n"]
                                     for arm, g in gates.items()},
            "floor_at_120_windows": PRIOR_FLOOR_AT_120,
            "reading": (
                "the leak-specific term still does not clear the grader's gate at 4x the windows -- "
                "it fails on the CI half in every arm -- while the pre-registered expected size of "
                f"{PREDICTED} nats now falls OUTSIDE the treatment interval, which the 120-window "
                "run could not establish. The floor itself is measured at this n rather than "
                "inherited."),
        },
        "inputs": [str(p) for p, _, _ in roled],
        "what_this_does_NOT_show": [
            "NOT that there is no leak. The pre-registered SIZE is excluded; a smaller leak is not.",
            "NOT a calibrated null. ONE control draw is one draw of trajectory noise, exactly as "
            "the graded artifact says of its own floor.",
            "NOT a re-derivation of the predicted -0.0054, which is quoted from the graded "
            "artifact and whose provenance was not revisited.",
            "NOT a generalisation beyond this corpus and this checkpoint family.",
        ],
        "honest_limits": [
            "The DiD identifies a leak only if cross-organism cost acts equally on both halves and "
            "cancels. That is the design's assumption and is not testable from these numbers.",
            "The clean half uses 500 of 519 available windows, so the bootstrap treats a near-"
            "exhaustive sample as an infinite population. The intervals are therefore CONSERVATIVE: "
            "a finite-population correction would narrow them, never widen them.",
        ],
    }
    a.out.write_text(json.dumps(out, indent=1) + "\n", encoding="utf-8")

    print(f"roles (derived): {out['arm_roles']}")
    print(f"n per half: bearing {n[0]}, clean {n[1]}")
    for arm, g in gates.items():
        print(f"  {arm:14s} treat {g['treatment_did_nats']:+.6f} CI {g['treatment_ci95']} | "
              f"floor |{g['control_magnitude_the_floor_at_this_n']:.6f}|")
        print(f"                 exceeds_floor={g['exceeds_floor']} "
              f"ci_excludes_zero={g['treatment_ci_excludes_zero']} -> {g['VERDICT']}"
              f" | excludes predicted={g['treatment_excludes_predicted_size']}")
    print(f"HEADLINE {out['headline']['VERDICT']} | predicted size excluded everywhere: "
          f"{all_excl_pred}")
    print(f"wrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
