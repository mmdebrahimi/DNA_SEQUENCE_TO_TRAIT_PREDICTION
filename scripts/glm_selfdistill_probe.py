"""Is AlphaGenome's SELF-DISTILLATION implementable on our oracle today? A kill-test, not an argument.

    uv run python scripts/glm_selfdistill_probe.py [--axis Minus10] [--rates 0.02,0.05,0.10]

**The method, and why it is the most transferable idea in the AlphaGenome programme.** Their phase 2 takes
the pre-trained model, RANDOMLY PERTURBS the input sequence, and trains on the pre-trained model's OWN
outputs as labels. Because the labels are predicted rather than measured, the held-out-region restriction
lifts and they can train across the whole genome. The problem it solves is exactly ours: real measured data
exists only for sequences someone assayed, so there is no direct training signal for what an EDIT does --
and our generative loop will query the oracle on generated sequences off the measured distribution.

**THE CLAIM UNDER TEST (mine, and reasoned rather than measured, which is why it needs killing):**
self-distillation is DEGENERATE on the current oracle. A ridge teacher's predictions lie exactly in the span
of its own features, so a linear student with the SAME feature basis should recover the same function no
matter where in sequence space the perturbed samples are drawn. If that holds, there is nothing to gain
until the student has capacity the perturbation can regularise -- and the sequencing is FORCED: learned
representation first, distillation second.

**The arms.**

- `teacher`          -- ridge on the real measured training data. The incumbent.
- `student_linear`   -- ridge, SAME basis, fitted on perturbed sequences labelled by the teacher.
                        Prediction: indistinguishable from the teacher.
- `student_mlp`      -- a small MLP, RICHER basis, same distillation targets. Prediction: cannot EXCEED the
                        teacher, because distillation cannot create information the teacher lacks. This is
                        the caveat that must travel with the method and it is tested, not asserted.
- `student_gc`       -- **the non-vacuity control.** Distil the same teacher into a GC-only student. This
                        MUST collapse. A probe that cannot detect an impoverished student detects nothing,
                        and "the student matched the teacher" would otherwise be unfalsifiable.

Everything is scored on the SAME leave-element-out test set the oracle is scored on, against REAL measured
expression -- so a student that merely agrees with the teacher cannot hide behind the teacher's own metric.
"""
from __future__ import annotations

import argparse
import json
import random
from datetime import date
from pathlib import Path

from dna_decode.glm.expression import (
    DEFAULT_CACHE,
    FEATURES,
    Pair,
    leave_element_out,
    load_pairs,
    spearman,
)

#: What the probe is allowed to conclude. Fixed before running.
OUTCOMES = {
    "DEGENERATE_DISTILLATION_BLOCKED_UNTIL_LEARNED_REPRESENTATION":
        "the linear student reproduces the teacher (agreement >= 0.99 and |d rho| <= 0.01) AND the MLP "
        "student does not exceed it -- so distillation buys nothing on a linear oracle and must wait",
    "DISTILLATION_HELPS_BUILD_IT_NOW":
        "some student BEATS the teacher on held-out measured expression by > 0.02 rho",
    "PROBE_VACUOUS_CONTROL_DID_NOT_COLLAPSE":
        "the GC-only student did NOT collapse, so the probe cannot distinguish a real match from a "
        "degenerate one and no conclusion is licensed",
}

AGREEMENT_FLOOR = 0.99
DRHO_TOL = 0.01
BEAT_MARGIN = 0.02
CONTROL_MUST_DROP = 0.20


def _perturb(seqs: list[str], rate: float, copies: int, seed: int) -> list[str]:
    """Random point substitutions -- AlphaGenome's 'randomly perturb the input', at nucleotide level."""
    rng = random.Random(seed)
    out: list[str] = []
    alts = {b: [x for x in "ACGT" if x != b] for b in "ACGT"}
    for s in seqs:
        for _ in range(copies):
            chars = list(s)
            n_mut = max(1, int(round(len(chars) * rate)))
            for pos in rng.sample(range(len(chars)), n_mut):
                chars[pos] = rng.choice(alts[chars[pos]])
            out.append("".join(chars))
    return out


def _fit(Xtr, ytr, Xte, kind: str):
    import numpy as np
    from sklearn.linear_model import Ridge
    from sklearn.preprocessing import StandardScaler
    sc = StandardScaler().fit(Xtr)
    if kind == "mlp":
        from sklearn.neural_network import MLPRegressor
        # max_iter is generous ON PURPOSE. The MLP arm exists to test "distillation cannot create
        # information the teacher lacks", and an UNDER-TRAINED student that fails to beat the teacher is a
        # weak negative -- it would be evidence about the optimiser, not about the information bound. At 80
        # iterations it stopped on max_iter without converging, so it gets 400 with early stopping to end it
        # as soon as held-out loss plateaus.
        m = MLPRegressor(hidden_layer_sizes=(64,), max_iter=400, early_stopping=True,
                         n_iter_no_change=15, random_state=0).fit(sc.transform(Xtr), ytr)
    else:
        m = Ridge(alpha=1.0).fit(sc.transform(Xtr), ytr)
    return [float(v) for v in m.predict(sc.transform(Xte))]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cache-dir", default=str(DEFAULT_CACHE))
    ap.add_argument("--axis", default="Minus10")
    ap.add_argument("--rates", default="0.02,0.05,0.10")
    ap.add_argument("--copies", type=int, default=2)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--no-download", action="store_true")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    import numpy as np

    pairs, rep = load_pairs(args.cache_dir, allow_download=not args.no_download)
    train, test, held = leave_element_out(pairs, args.axis, seed=args.seed)
    print(f"{args.axis}: held {held} -> train {len(train)} / test {len(test)}")

    onehot, gc = FEATURES["onehot"], FEATURES["gc"]
    tr_seqs = [p.seq for p in train]
    te_seqs = [p.seq for p in test]
    measured = [p.expression for p in test]

    # --- the teacher: ridge on the REAL measured data -------------------------------------------------
    Xtr = np.asarray(onehot(tr_seqs), dtype=float)
    Xte = np.asarray(onehot(te_seqs), dtype=float)
    ytr = np.log1p(np.asarray([p.expression for p in train], dtype=float))
    teacher_te = _fit(Xtr, ytr, Xte, "ridge")
    rho_teacher = spearman(teacher_te, measured)
    print(f"  teacher (real labels)           rho={rho_teacher:+.4f}")

    rows = []
    for rate in [float(r) for r in args.rates.split(",")]:
        pert = _perturb(tr_seqs, rate, args.copies, args.seed)
        Xp = np.asarray(onehot(pert), dtype=float)
        # labels are the TEACHER's own outputs on the perturbed sequences -- no measurement involved
        teacher_on_pert = _fit(Xtr, ytr, Xp, "ridge")
        yp = np.asarray(teacher_on_pert, dtype=float)

        arms = {}
        for name, kind, feats in (("student_linear", "ridge", onehot),
                                  ("student_mlp", "mlp", onehot),
                                  ("student_gc", "ridge", gc)):
            Xp_f = Xp if feats is onehot else np.asarray(feats(pert), dtype=float)
            Xte_f = Xte if feats is onehot else np.asarray(feats(te_seqs), dtype=float)
            pred = _fit(Xp_f, yp, Xte_f, kind)
            arms[name] = {
                "rho_vs_measured": round(spearman(pred, measured), 4),
                "agreement_with_teacher": round(spearman(pred, teacher_te), 4),
                "delta_rho_vs_teacher": round(spearman(pred, measured) - rho_teacher, 4),
            }
        rows.append({"mutation_rate": rate, "n_perturbed": len(pert), "arms": arms})
        bits = "  ".join(f"{n}: rho={a['rho_vs_measured']:+.4f} agree={a['agreement_with_teacher']:.4f}"
                         for n, a in arms.items())
        print(f"  rate={rate:<5.2f} n={len(pert):<6d} {bits}")

    verdict, detail = _verdict(rows, rho_teacher)
    artifact = {
        "schema": "glm-selfdistill-probe-v1",
        "date": str(date.today()),
        "claim_under_test": ("self-distillation is DEGENERATE on a linear oracle with a fixed feature "
                             "basis, so it must wait for a learned representation"),
        "method_source": ("AlphaGenome phase 2: perturb the input, train on the pre-trained model's own "
                          "outputs as labels (technical talk, 2026)"),
        "axis": args.axis, "held_out_values": held,
        "n_train": len(train), "n_test": len(test),
        "teacher_rho_vs_measured": round(rho_teacher, 4),
        "thresholds": {"agreement_floor": AGREEMENT_FLOOR, "delta_rho_tol": DRHO_TOL,
                       "beat_margin": BEAT_MARGIN, "control_must_drop_by": CONTROL_MUST_DROP},
        "outcomes_defined_before_running": OUTCOMES,
        "rows": rows,
        "verdict": verdict,
        "verdict_detail": detail,
        "honest_limits": [
            "ONE element axis and one teacher family; this bounds the LINEAR case, which is the case that "
            "matters because it is what currently ships.",
            "A negative here does NOT say distillation is useless -- it says it is useless on THIS oracle. "
            "The whole point is that it becomes testable once a learned representation exists.",
            "Perturbation is uniform random point substitution; a generator's candidates are not uniform, "
            "so this does not simulate the real off-distribution gap, only the method's mechanics.",
            "The MLP arm is a capacity probe, not a tuned model -- it is not evidence about what a "
            "well-trained CNN could do on REAL labels.",
        ],
    }
    out = Path(args.out) if args.out else Path(f"wiki/glm_selfdistill_probe_{date.today()}.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(artifact, indent=2), encoding="utf-8")
    print(f"\nVERDICT: {verdict}\n  {detail}\nwrote {out}")
    return 0


def _verdict(rows: list[dict], rho_teacher: float) -> tuple[str, str]:
    """Mechanical, against thresholds fixed before the run."""
    if not rows:
        return "INDETERMINATE_NO_ROWS", "no mutation rate produced a result"
    # non-vacuity FIRST: if the control did not collapse, nothing else is interpretable
    ctrl = [r["arms"]["student_gc"]["rho_vs_measured"] for r in rows]
    if not all(rho_teacher - c >= CONTROL_MUST_DROP for c in ctrl):
        return ("PROBE_VACUOUS_CONTROL_DID_NOT_COLLAPSE",
                f"GC-only student reached {ctrl} against teacher {rho_teacher:.4f}; the probe cannot "
                f"distinguish a real match from a degenerate one")
    beats = [(r["mutation_rate"], n, a["delta_rho_vs_teacher"])
             for r in rows for n, a in r["arms"].items()
             if n != "student_gc" and a["delta_rho_vs_teacher"] > BEAT_MARGIN]
    if beats:
        return "DISTILLATION_HELPS_BUILD_IT_NOW", f"students beating the teacher: {beats}"
    lin = [(r["mutation_rate"], r["arms"]["student_linear"]["agreement_with_teacher"],
            r["arms"]["student_linear"]["delta_rho_vs_teacher"]) for r in rows]
    reproduced = all(ag >= AGREEMENT_FLOOR and abs(d) <= DRHO_TOL for _, ag, d in lin)
    if reproduced:
        return ("DEGENERATE_DISTILLATION_BLOCKED_UNTIL_LEARNED_REPRESENTATION",
                f"linear student reproduces the teacher at every rate (rate, agreement, d_rho) = {lin}; "
                f"control collapsed as required")
    return ("INDETERMINATE_LINEAR_STUDENT_DIVERGED_WITHOUT_HELPING",
            f"the linear student neither reproduced the teacher nor beat it: {lin} -- the degeneracy "
            f"argument does not hold as stated, and that is a finding, not a pass")


if __name__ == "__main__":
    raise SystemExit(main())
