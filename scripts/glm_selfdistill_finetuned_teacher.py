"""Does the self-distillation degeneracy SURVIVE swapping in the fine-tuned encoder as teacher?

    uv run python scripts/glm_selfdistill_finetuned_teacher.py \
        --cache-dir data/mpra_gse108535 --ckpt data/promoteratlas/finetune_checkpoints/seed_0.pt

**WHY THIS EXISTS.** `glm_selfdistill_probe.py` measured self-distillation DEGENERATE on 2026-10-07:
a ridge teacher's predictions lie exactly in the span of its own features, so a linear student with
the SAME basis recovered it at agreement **1.0000**. The umbrella therefore blocks
`glm-self-distillation` until "a learned representation LANDS", and records the blocker as
"a measured degeneracy, not a preference" -- so lifting it needs a MEASUREMENT, not an argument.

On 2026-11-17 the `finetuned-promoteratlas-encoder-150bp` arm returned FINETUNE_ADDS (median 0.6405
over 10 seeds vs the frozen 0.5076). The question this script answers is narrow and mechanical:
**with that network as the teacher instead of a ridge, does the linear student still reproduce it?**

**THE MECHANISM BEING TESTED, STATED PLAINLY.** The 1.0000 is an IDENTITY, not a biological finding:
a linear function on basis B, resampled anywhere in sequence space and refit on basis B, recovers
itself. A 1,847,428-parameter conv network is not in the span of the one-hot basis, so the identity
cannot hold by construction. That is a PREDICTION, and predicting an outcome is not measuring it.

**THE POSITIVE CONTROL IS THE LOAD-BEARING PART.** This script runs BOTH teachers through the SAME
harness -- same substrate, same splits, same perturbations, same students, same thresholds -- and the
ONLY thing that differs is the teacher. The ridge arm must reproduce ~1.0000. If it does not, this
harness is broken and NOTHING about the neural teacher is interpretable; that is a named outcome, not
a footnote. Without it, a low agreement for the neural teacher is indistinguishable from a bug in my
own scoring code.

**THE COMPETENCE GATE, and why its threshold is borrowed rather than invented.** A student failing to
copy a teacher that carries no signal proves nothing -- it is copying noise. So the new teacher must
itself be competent on THIS substrate before any agreement number is read. The floor is not a fresh
number I chose: it reuses the original probe's own `CONTROL_MUST_DROP` (0.20, fixed 2026-10-07).
The probe already demands that a GC-only student fall at least that far BELOW the teacher, so a
teacher sitting at GC level makes the probe's own control incoherent. Hence:

    teacher rho  >=  standalone GC rho  +  CONTROL_MUST_DROP

**SUBSTRATE, AND THE CIRCULARITY CHECK THAT PERMITS IT.** The probe's substrate is GEO GSE108535
(Urtecho/Kosuri sigma70 variants, 10,898 pairs, every sequence exactly 150 bp). The encoder was
FINE-TUNED on GSE144621 peak tiles -- a different dataset -- so this is a cross-dataset application
of the teacher and its rho here is not its rho there. `wiki/promoteratlas_circularity_2026-10-31.md`
resolves the one hazard that would void it: the `promoteratlas-trspred-urtecho2018` head IS fit on
GSE108535 and "must never be a comparator on that rung", but this arm descends from
`trained_weights/base_model/promoteratlas-base.pt`, whose corpus is genomic sequence with no
expression file -- `cleared_vs_GSE108535_for_the_BASE_model`. Verified in code, not assumed:
`_build_trainable` -> `load_atlas` loads exactly that base checkpoint.

Everything is scored on the SAME leave-element-out test set against REAL measured expression, so a
student that merely agrees with its teacher cannot hide behind the teacher's own metric.
"""
from __future__ import annotations

import argparse
import json
from datetime import date
from pathlib import Path

import numpy as np

from dna_decode.glm.expression import (
    FEATURES,
    leave_element_out,
    load_pairs,
    spearman,
)

# Imported, never re-typed: the thresholds and the verdict ladder were fixed on 2026-10-07 and a
# copy here could drift into a bar chosen after seeing a number.
from scripts.glm_selfdistill_probe import (
    AGREEMENT_FLOOR,
    BEAT_MARGIN,
    CONTROL_MUST_DROP,
    DRHO_TOL,
    _fit,
    _perturb,
)

SEQ_LENGTH = 150

#: What this script is allowed to conclude. Fixed before running, and the first two are failures of
#: the instrument rather than findings about the world -- which is the distinction that keeps a
#: broken run from reading as a result.
OUTCOMES = {
    "HARNESS_BROKEN_RIDGE_CONTROL_DID_NOT_REPRODUCE":
        "the ridge teacher did NOT recover agreement >= the floor in this harness, so the harness "
        "fails to reproduce the 2026-10-07 result and no statement about any other teacher is "
        "licensed",
    "TEACHER_NOT_COMPETENT_ON_THIS_SUBSTRATE":
        "the fine-tuned teacher does not clear standalone-GC + CONTROL_MUST_DROP on GSE108535, so a "
        "student failing to copy it would be copying noise and the agreement is uninterpretable",
    "DEGENERACY_PERSISTS_WITH_A_LEARNED_TEACHER":
        "the linear student reproduces the LEARNED teacher too (agreement >= floor and |d rho| <= "
        "tol), so the measured degeneracy is NOT specific to a linear oracle and the block stands",
    "DEGENERACY_DOES_NOT_TRANSFER_TO_THE_LEARNED_TEACHER":
        "the ridge control reproduced as required while the learned teacher's linear student did "
        "NOT, so the identity that produced 1.0000 does not hold for this teacher",
    "DISTILLATION_HELPS_ON_THE_LEARNED_TEACHER":
        "some student BEATS the learned teacher on held-out measured expression by > beat margin",
}


def load_finetuned_teacher(ckpt: Path, pa_repo: Path | None):
    """The fine-tuned encoder as a scoring function. Mirrors the arm's own forward pass exactly."""
    import torch

    from scripts.finetune_arm_probe import _build_trainable, _onehot

    enc, head = _build_trainable(pa_repo)
    state = torch.load(ckpt, map_location="cpu", weights_only=False)
    enc.load_state_dict(state["enc"])
    head.load_state_dict(state["head"])
    enc.eval()
    head.eval()

    def score(seqs: list[str], chunk: int = 256) -> np.ndarray:
        outs = []
        with torch.no_grad():
            for i in range(0, len(seqs), chunk):
                X = _onehot(seqs[i:i + chunk], SEQ_LENGTH)
                z = enc(X, return_latent=True)[2]
                outs.append(head(z.mean(dim=2)).squeeze(-1).numpy())
                del z, X
        return np.concatenate(outs) if outs else np.zeros(0)

    meta = {"checkpoint": str(ckpt), "epoch_in_checkpoint": int(state.get("epoch", -1))}
    return score, meta


def students_against(teacher_te, teacher_on_pert, *, pert, te_seqs, measured, rho_teacher):
    """The three distillation arms, labelled by whatever teacher was handed in."""
    onehot, gc = FEATURES["onehot"], FEATURES["gc"]
    Xp = np.asarray(onehot(pert), dtype=float)
    Xte = np.asarray(onehot(te_seqs), dtype=float)
    yp = np.asarray(teacher_on_pert, dtype=float)

    arms = {}
    for name, kind, feats in (("student_linear", "ridge", onehot),
                              ("student_mlp", "mlp", onehot),
                              ("student_gc", "ridge", gc)):
        Xp_f = Xp if feats is onehot else np.asarray(feats(pert), dtype=float)
        Xte_f = Xte if feats is onehot else np.asarray(feats(te_seqs), dtype=float)
        pred = _fit(Xp_f, yp, Xte_f, kind)
        r = spearman(pred, measured)
        arms[name] = {
            "rho_vs_measured": round(r, 4),
            "agreement_with_teacher": round(spearman(pred, teacher_te), 4),
            "delta_rho_vs_teacher": round(r - rho_teacher, 4),
        }
    return arms


def verdict(ridge_rows, atlas_rows, *, rho_ridge, rho_atlas, rho_gc_standalone):
    """Mechanical, against thresholds fixed before the run. Instrument failures are checked FIRST."""
    # 1. positive control. If the harness cannot reproduce the known result, stop.
    ridge_lin = [r["arms"]["student_linear"]["agreement_with_teacher"] for r in ridge_rows]
    if not ridge_lin or not all(a >= AGREEMENT_FLOOR for a in ridge_lin):
        return ("HARNESS_BROKEN_RIDGE_CONTROL_DID_NOT_REPRODUCE",
                "ridge teacher's linear-student agreement %s did not clear the %.2f floor that the "
                "2026-10-07 probe measured at 1.0000" % (ridge_lin, AGREEMENT_FLOOR))

    # 2. competence gate, borrowed threshold.
    floor = rho_gc_standalone + CONTROL_MUST_DROP
    if rho_atlas < floor:
        return ("TEACHER_NOT_COMPETENT_ON_THIS_SUBSTRATE",
                "fine-tuned teacher rho %+.4f on GSE108535 is below standalone-GC %+.4f + %.2f = "
                "%+.4f, so its students would be copying noise" % (
                    rho_atlas, rho_gc_standalone, CONTROL_MUST_DROP, floor))

    # 3. did anything beat the learned teacher outright?
    beats = [(r["mutation_rate"], n, a["delta_rho_vs_teacher"])
             for r in atlas_rows for n, a in r["arms"].items()
             if n != "student_gc" and a["delta_rho_vs_teacher"] > BEAT_MARGIN]
    if beats:
        return "DISTILLATION_HELPS_ON_THE_LEARNED_TEACHER", "students beating the teacher: %s" % beats

    # 4. the question itself.
    lin = [(r["mutation_rate"], r["arms"]["student_linear"]["agreement_with_teacher"],
            r["arms"]["student_linear"]["delta_rho_vs_teacher"]) for r in atlas_rows]
    if all(ag >= AGREEMENT_FLOOR and abs(d) <= DRHO_TOL for _, ag, d in lin):
        return ("DEGENERACY_PERSISTS_WITH_A_LEARNED_TEACHER",
                "linear student reproduces the LEARNED teacher at every rate "
                "(rate, agreement, d_rho) = %s" % lin)
    return ("DEGENERACY_DOES_NOT_TRANSFER_TO_THE_LEARNED_TEACHER",
            "ridge control reproduced at %s while the learned teacher's linear student reached "
            "(rate, agreement, d_rho) = %s -- below the %.2f floor" % (
                ridge_lin, lin, AGREEMENT_FLOOR))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cache-dir", default="data/mpra_gse108535")
    ap.add_argument("--ckpt", default="data/promoteratlas/finetune_checkpoints/seed_0.pt")
    ap.add_argument("--pa-repo", default=None)
    ap.add_argument("--axis", default="Minus10")
    ap.add_argument("--rates", default="0.02,0.05,0.10")
    ap.add_argument("--copies", type=int, default=2)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    pairs, _rep = load_pairs(args.cache_dir, allow_download=False)
    train, test, held = leave_element_out(pairs, args.axis, seed=args.seed)
    print("%s: held %s -> train %d / test %d" % (args.axis, held, len(train), len(test)))

    onehot, gc = FEATURES["onehot"], FEATURES["gc"]
    tr_seqs = [p.seq for p in train]
    te_seqs = [p.seq for p in test]
    measured = [p.expression for p in test]

    Xtr = np.asarray(onehot(tr_seqs), dtype=float)
    Xte = np.asarray(onehot(te_seqs), dtype=float)
    ytr = np.log1p(np.asarray([p.expression for p in train], dtype=float))

    # --- teacher A: the incumbent ridge. THE POSITIVE CONTROL. ------------------------------------
    ridge_te = _fit(Xtr, ytr, Xte, "ridge")
    rho_ridge = spearman(ridge_te, measured)
    print("  teacher A  ridge (real labels)        rho=%+.4f" % rho_ridge)

    # --- teacher B: the fine-tuned encoder --------------------------------------------------------
    pa_repo = Path(args.pa_repo) if args.pa_repo else None
    score, ck_meta = load_finetuned_teacher(Path(args.ckpt), pa_repo)
    atlas_te = score(te_seqs)
    rho_atlas = spearman([float(v) for v in atlas_te], measured)
    print("  teacher B  fine-tuned encoder         rho=%+.4f  (%s)" % (rho_atlas, args.ckpt))

    # --- the competence reference: GC alone on real labels ----------------------------------------
    gc_tr = np.asarray(gc(tr_seqs), dtype=float)
    gc_te = np.asarray(gc(te_seqs), dtype=float)
    rho_gc_standalone = spearman(_fit(gc_tr, ytr, gc_te, "ridge"), measured)
    print("  reference  GC-only (real labels)      rho=%+.4f   -> teacher floor %+.4f" % (
        rho_gc_standalone, rho_gc_standalone + CONTROL_MUST_DROP))

    ridge_rows, atlas_rows = [], []
    for rate in [float(r) for r in args.rates.split(",")]:
        pert = _perturb(tr_seqs, rate, args.copies, args.seed)

        ridge_on_pert = _fit(Xtr, ytr, np.asarray(onehot(pert), dtype=float), "ridge")
        a_ridge = students_against(ridge_te, ridge_on_pert, pert=pert, te_seqs=te_seqs,
                                  measured=measured, rho_teacher=rho_ridge)
        ridge_rows.append({"mutation_rate": rate, "n_perturbed": len(pert), "arms": a_ridge})

        atlas_on_pert = [float(v) for v in score(pert)]
        a_atlas = students_against([float(v) for v in atlas_te], atlas_on_pert, pert=pert,
                                  te_seqs=te_seqs, measured=measured, rho_teacher=rho_atlas)
        atlas_rows.append({"mutation_rate": rate, "n_perturbed": len(pert), "arms": a_atlas})

        print("  rate=%-5.2f n=%-6d RIDGE lin_agree=%.4f | ATLAS lin_agree=%.4f  mlp_agree=%.4f  "
              "gc_rho=%+.4f" % (
                  rate, len(pert), a_ridge["student_linear"]["agreement_with_teacher"],
                  a_atlas["student_linear"]["agreement_with_teacher"],
                  a_atlas["student_mlp"]["agreement_with_teacher"],
                  a_atlas["student_gc"]["rho_vs_measured"]))

    v, detail = verdict(ridge_rows, atlas_rows, rho_ridge=rho_ridge, rho_atlas=rho_atlas,
                        rho_gc_standalone=rho_gc_standalone)

    artifact = {
        "schema": "glm-selfdistill-finetuned-teacher-v1",
        "date": str(date.today()),
        "question": ("does the 2026-10-07 self-distillation degeneracy survive replacing the ridge "
                     "teacher with the fine-tuned PromoterAtlas encoder?"),
        "why_the_ridge_arm_is_here": ("positive control -- same substrate, splits, perturbations, "
                                      "students and thresholds; only the teacher differs. It must "
                                      "reproduce ~1.0000 or this harness proves nothing."),
        "substrate": {
            "dataset": "GEO GSE108535 (Urtecho/Kosuri sigma70 variants)",
            "n_pairs": len(pairs), "axis": args.axis, "held_out_values": held,
            "n_train": len(train), "n_test": len(test), "seq_length": SEQ_LENGTH,
            "teacher_was_finetuned_on": "GSE144621 peak tiles -- a DIFFERENT dataset",
            "circularity": ("base model cleared_vs_GSE108535 per "
                            "wiki/promoteratlas_circularity_2026-10-31.md; the contaminated "
                            "trspred-urtecho2018 head is NOT used -- verified in code"),
        },
        "teacher_finetuned": ck_meta,
        "teacher_rho_ridge_control": round(rho_ridge, 4),
        "teacher_rho_finetuned": round(rho_atlas, 4),
        "gc_standalone_rho": round(rho_gc_standalone, 4),
        "competence_floor": round(rho_gc_standalone + CONTROL_MUST_DROP, 4),
        "thresholds_imported_not_retyped": {
            "agreement_floor": AGREEMENT_FLOOR, "delta_rho_tol": DRHO_TOL,
            "beat_margin": BEAT_MARGIN, "control_must_drop_by": CONTROL_MUST_DROP,
            "source": "scripts/glm_selfdistill_probe.py, fixed 2026-10-07",
        },
        "outcomes_defined_before_running": OUTCOMES,
        "rows_ridge_control": ridge_rows,
        "rows_finetuned_teacher": atlas_rows,
        "verdict": v,
        "verdict_detail": detail,
        "honest_limits": [
            "ONE element axis, ONE checkpoint (one of ten seeds), ONE substrate. This bounds the "
            "MECHANISM, not the eventual value of distillation.",
            "CROSS-DATASET: the teacher was fine-tuned on GSE144621 tiles and is scored here on "
            "GSE108535, so its rho here is NOT the 0.6405 the arm reported and the two must never "
            "be quoted interchangeably.",
            "A non-transfer verdict says the 1.0000 IDENTITY does not hold for this teacher. It does "
            "NOT show self-distillation WORKS -- that needs the real-data validation half the ban's "
            "own wording requires.",
            "Perturbation is uniform random point substitution; a generator's candidates are not "
            "uniform, so this exercises the method's mechanics and not the real off-distribution gap.",
            "The students use the one-hot basis. A student with a RICHER basis than tested here might "
            "copy the network more closely, so this is a lower bound on copyability.",
        ],
    }
    out = Path(args.out) if args.out else Path(
        "wiki/glm_selfdistill_finetuned_teacher_%s.json" % date.today())
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(artifact, indent=2), encoding="utf-8")
    print("\nVERDICT: %s\n  %s\nwrote %s" % (v, detail, out))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
