"""Does the one-hot blind-spot model have a blind spot OF ITS OWN? (2026-10-01)

THE SHARP HYPOTHESIS. One-hot encodes each (position, amino-acid) pair as an INDEPENDENT column, so a
substitution whose column carries no training signal gets ~zero weight: the model is structurally blind to
it. Positions are islands — one-hot knows nothing about 181 being near 188, or valine resembling isoleucine.
A position-level PRIOR (structure, conservation, a geometric channel) CAN score such a substitution, which
is the orthogonality argument for adding a modality rather than more parameters.

SO MEASURE THE ROOM BEFORE BUILDING A CHANNEL. If blind-spot isolates carry essentially no unseen
substitutions under leave-one-STUDY-out, there is nothing for a prior to add and the honest output is that
finding, not a model. This is the repo's own recorded pattern twice over: an over-caller empties its own
blind spot, and "exhausted lever headroom tracks coverage" — so check whether the lever has room first.

WHAT IS REUSED, NOT RE-DERIVED. The fold structure (GroupKFold by RefID), the one-hot featurizer, the
1%-consensus drift detection, the blind-spot definition (the DEPLOYED `call_hiv_observed`) and the AUROC all
come from `hiv_deployability_per_gene.py` / `hiv_supervised_deployability.py`. The floor model is the same
LogisticRegression, so the OOF scores here reproduce the committed 0.8102 — that is the ANCHOR.

THE VERDICT RULE IS FROZEN BELOW, BEFORE ANY NUMBER WAS SEEN, and a stratum too small to score is REFUSED
rather than reported.

Run: uv run python scripts/hiv_onehot_unseen_position_headroom.py [--self-check]
Exits 2 without the gitignored Stanford data, 3 if the floor does not reproduce.
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import sys
from datetime import date
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
for p in (str(REPO), str(REPO / "scripts")):
    if p not in sys.path:
        sys.path.insert(0, p)

FLOOR_ARTIFACT = REPO / "wiki" / "hiv_deployability_per_gene_2026-09-30.json"
ANCHOR_TOL = 0.02
MIN_STRATUM = 10          # a stratum with fewer than this of either class is not scoreable

# ---------------------------------------------------------------------------
# FROZEN VERDICT RULE (fixed before the first run; deltas are reported either way)
FROZEN_RULE = {
    "frozen": "2026-10-01",
    "frozen_before_any_number_was_seen": True,
    "HEADROOM_EXISTS": ("at least MIN_EXPOSED_FRAC of blind-spot isolates carry >=1 substitution with NO "
                        "training signal in their own fold, AND the model's AUROC on that exposed stratum "
                        "is at least AUROC_GAP below its AUROC on the unexposed stratum"),
    "NO_HEADROOM_ONEHOT_ALREADY_COVERS_IT": ("exposure is below MIN_EXPOSED_FRAC, or the exposed stratum "
                                             "scores as well as the unexposed one -- a position-level "
                                             "prior has nothing to add and the channel should NOT be built"),
    "UNDERPOWERED": "either stratum holds fewer than MIN_STRATUM of either class",
    "MIN_EXPOSED_FRAC": 0.20,
    "AUROC_GAP": 0.05,
}
MIN_EXPOSED_FRAC = FROZEN_RULE["MIN_EXPOSED_FRAC"]
AUROC_GAP = FROZEN_RULE["AUROC_GAP"]
# ---------------------------------------------------------------------------


class ProbeRefused(RuntimeError):
    """Raised instead of returning a number that cannot be trusted."""


def _mod(name: str):
    spec = importlib.util.spec_from_file_location(name, REPO / "scripts" / f"{name}.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def probe_gene(cls: str, D, PG) -> dict:
    """One gene: per-fold unseen-column exposure, then stratified AUROC of the floor model."""
    import numpy as np
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import GroupKFold

    import dna_decode.data.hiv_amr as H
    ev, S = D.ev, D.S
    cfg = PG.GENES[cls]

    path = PG.RAW / cfg["dataset"]
    if not path.exists():
        raise ProbeRefused(f"dataset absent: {path}")
    prot = PG._protein(cfg["protein"], ev)
    rows = list(csv.DictReader(open(path, encoding="utf-8"), delimiter="\t"))
    pcols = [c for c in rows[0] if c.startswith("P") and c[1:].isdigit()]
    col, gene = cfg["col"], cfg["gene"]
    have = [r for r in rows if r.get(col) not in ("NA", "", "-", None)]

    drifted = set()
    for c in pcols:
        p = int(c[1:])
        if p <= len(prot) and sum(1 for r in have if (r[c] or "").strip() == prot[p - 1]) > 0.01 * len(have):
            drifted.add(p)

    def subs_of(r):
        return {f"{prot[p-1]}{p}{aa}" for p, aa in ev.isolate_muts(r, pcols)
                if p <= len(prot) and p not in drifted and prot[p - 1] != aa}

    cat = np.array([1 if H.call_hiv_observed(cfg["drug"], {gene: subs_of(r)}).prediction == "R" else 0
                    for r in have])
    neg = cat == 0
    y = np.array([1 if float(r[col]) >= PG.CUTOFF else 0 for r in have])
    X, feats = S.build_onehot(have, pcols, drifted, prot)
    X = np.asarray(X, dtype=float)
    refid = np.array([r["RefID"] for r in have])
    folds = GroupKFold(n_splits=min(5, len(set(refid.tolist()))))

    n = len(have)
    oof = np.zeros(n)
    unseen = np.zeros(n, dtype=int)          # substitutions with NO training signal, per isolate
    zero_wt_coef_absmax = []                 # VERIFY the blindness instead of assuming it

    for tr, te in folds.split(X, y, groups=refid):
        clf = LogisticRegression(max_iter=2000, C=1.0, solver="liblinear").fit(X[tr], y[tr])
        oof[te] = clf.predict_proba(X[te])[:, 1]
        train_carriers = X[tr].sum(axis=0)          # per-column training carrier count
        zero_cols = np.flatnonzero(train_carriers == 0)
        # the isolate's own active columns that had zero training carriers
        for i in te:
            active = np.flatnonzero(X[i] > 0)
            unseen[i] = int(np.intersect1d(active, zero_cols, assume_unique=False).size)
        if zero_cols.size:
            zero_wt_coef_absmax.append(float(np.abs(clf.coef_[0][zero_cols]).max()))

    out = {"class": cls, "gene": gene, "drug": cfg["drug"], "n": n,
           "n_features": X.shape[1], "n_studies": len(set(refid.tolist())),
           "floor_blindspot_auroc": round(ev.auroc(y[neg].tolist(), oof[neg].tolist()), 4),
           "blind_spot_n": int(neg.sum()), "blind_spot_R": int(y[neg].sum()),
           "zero_training_weight_absmax": (round(max(zero_wt_coef_absmax), 6)
                                           if zero_wt_coef_absmax else None)}

    # --- exposure on the blind spot ---
    bs = np.flatnonzero(neg)
    exp_mask = unseen[bs] > 0
    ybs, sbs = y[bs], oof[bs]
    out["exposure"] = {
        "frac_blindspot_with_unseen": round(float(exp_mask.mean()), 4),
        "mean_unseen_blindspot": round(float(unseen[bs].mean()), 3),
        "mean_unseen_R": round(float(unseen[bs][ybs == 1].mean()), 3) if (ybs == 1).any() else None,
        "mean_unseen_S": round(float(unseen[bs][ybs == 0].mean()), 3) if (ybs == 0).any() else None,
        "frac_R_with_unseen": (round(float((unseen[bs][ybs == 1] > 0).mean()), 4)
                               if (ybs == 1).any() else None),
        "frac_S_with_unseen": (round(float((unseen[bs][ybs == 0] > 0).mean()), 4)
                               if (ybs == 0).any() else None),
    }

    # --- WHY: is the blindness absorbed by per-isolate REDUNDANCY? ---
    # These clinical isolates carry many substitutions each, so one zero-weight column may simply be
    # outvoted by the rest. The decisive descriptive statistic is the FRACTION of an isolate's active
    # columns that carry no training signal: if no isolate is mostly-blind, the blindness cannot bite.
    # Measured rather than asserted -- this repo has had three asserted causes refuted by measurement.
    active_counts = (X[bs] > 0).sum(axis=1).astype(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        frac_blind = np.where(active_counts > 0, unseen[bs] / active_counts, 0.0)
    expf = frac_blind[exp_mask]
    out["redundancy"] = {
        "mean_active_columns_per_isolate": round(float(active_counts.mean()), 2),
        "median_active_columns_per_isolate": float(np.median(active_counts)),
        "among_exposed_mean_blind_fraction": round(float(expf.mean()), 4) if expf.size else None,
        "among_exposed_median_blind_fraction": round(float(np.median(expf)), 4) if expf.size else None,
        "among_exposed_MAX_blind_fraction": round(float(expf.max()), 4) if expf.size else None,
        "n_isolates_majority_blind": int((frac_blind > 0.5).sum()),
        "reading": ("if MAX_blind_fraction is small and n_isolates_majority_blind is 0, then no isolate's "
                    "signal is mostly unseen, so the surviving seen substitutions outvote the blind ones "
                    "-- redundancy ABSORBS the blindness. That also predicts the blindness WOULD bite on "
                    "SINGLE-mutant prediction, which is a different task."),
    }

    # --- stratified AUROC: does the model do WORSE where it is blind? ---
    strata = {}
    for name, m in (("exposed", exp_mask), ("unexposed", ~exp_mask)):
        yy, ss = ybs[m], sbs[m]
        n_r, n_s = int(yy.sum()), int(len(yy) - yy.sum())
        rec = {"n": int(len(yy)), "n_R": n_r, "n_S": n_s}
        rec["auroc"] = (round(ev.auroc(yy.tolist(), ss.tolist()), 4)
                        if n_r >= MIN_STRATUM and n_s >= MIN_STRATUM else None)
        if rec["auroc"] is None:
            rec["refused"] = (f"{n_r} R / {n_s} S -- fewer than {MIN_STRATUM} of a class cannot support "
                              "an AUROC, so NO number is reported")
        strata[name] = rec
    out["strata"] = strata
    a_e, a_u = strata["exposed"]["auroc"], strata["unexposed"]["auroc"]
    out["auroc_gap_unexposed_minus_exposed"] = (round(a_u - a_e, 4)
                                                if (a_e is not None and a_u is not None) else None)
    return out


def verdict(cell: dict) -> dict:
    """Mechanical application of the FROZEN rule."""
    ex = cell["exposure"]["frac_blindspot_with_unseen"]
    a_e = cell["strata"]["exposed"]["auroc"]
    a_u = cell["strata"]["unexposed"]["auroc"]
    if a_e is None or a_u is None:
        return {"verdict": "UNDERPOWERED",
                "why": "a stratum could not be scored; see strata[*].refused"}
    gap = a_u - a_e
    if ex >= MIN_EXPOSED_FRAC and gap >= AUROC_GAP:
        v = "HEADROOM_EXISTS"
    else:
        v = "NO_HEADROOM_ONEHOT_ALREADY_COVERS_IT"
    return {"verdict": v, "exposed_frac": ex, "auroc_gap": round(gap, 4),
            "needed_exposed_frac": MIN_EXPOSED_FRAC, "needed_gap": AUROC_GAP,
            "rule": FROZEN_RULE[v]}


def self_check() -> None:
    assert FROZEN_RULE["frozen_before_any_number_was_seen"] is True
    assert MIN_EXPOSED_FRAC == 0.20 and AUROC_GAP == 0.05 and MIN_STRATUM == 10
    # the rule must be able to return BOTH substantive verdicts, or it is not a test
    hi = {"exposure": {"frac_blindspot_with_unseen": 0.5},
          "strata": {"exposed": {"auroc": 0.70}, "unexposed": {"auroc": 0.85}}}
    lo = {"exposure": {"frac_blindspot_with_unseen": 0.5},
          "strata": {"exposed": {"auroc": 0.84}, "unexposed": {"auroc": 0.85}}}
    assert verdict(hi)["verdict"] == "HEADROOM_EXISTS"
    assert verdict(lo)["verdict"] == "NO_HEADROOM_ONEHOT_ALREADY_COVERS_IT"
    none = {"exposure": {"frac_blindspot_with_unseen": 0.5},
            "strata": {"exposed": {"auroc": None}, "unexposed": {"auroc": 0.8}}}
    assert verdict(none)["verdict"] == "UNDERPOWERED"
    print("self-check OK -- frozen rule reachable in all three directions")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default=str(REPO / "wiki" / f"hiv_onehot_unseen_headroom_{date.today()}.json"))
    ap.add_argument("--self-check", action="store_true")
    a = ap.parse_args(argv)
    if a.self_check:
        self_check()
        return 0

    D, PG = _mod("hiv_supervised_deployability"), _mod("hiv_deployability_per_gene")
    floors = json.loads(FLOOR_ARTIFACT.read_text(encoding="utf-8"))["cells"] \
        if FLOOR_ARTIFACT.exists() else {}

    cells, refused = {}, {}
    for cls in ("NNRTI", "INSTI"):
        try:
            cells[cls] = probe_gene(cls, D, PG)
        except ProbeRefused as exc:
            refused[cls] = str(exc)
            print(f"{cls}: REFUSED -- {exc}")
    if not cells:
        print("no gene could be probed (gitignored Stanford data absent?)")
        return 2

    # ANCHOR: the floor must reproduce, or the strata are measured against the wrong model.
    anchors = {}
    for cls, c in cells.items():
        want = (floors.get(cls) or {}).get("leave_study_out_blindspot_auroc_MEASURED")
        got = c["floor_blindspot_auroc"]
        anchors[cls] = {"committed": want, "here": got,
                        "reproduces": bool(want is not None and abs(got - want) <= ANCHOR_TOL)}
    print("\nANCHOR (the floor model must reproduce the committed deployability AUROC):")
    for cls, an in anchors.items():
        print(f"  {cls:6s} committed {an['committed']} vs here {an['here']} -> "
              f"{'REPRODUCES' if an['reproduces'] else 'MISMATCH'}")
    if not all(an["reproduces"] for an in anchors.values()):
        print("REFUSING: strata measured against a floor that does not reproduce are not interpretable.")
        return 3

    for cls, c in cells.items():
        e = c["exposure"]
        print(f"\n=== {cls} ({c['gene']}) blind spot n={c['blind_spot_n']} R={c['blind_spot_R']} ===")
        print(f"  max |coef| over zero-training-signal columns: {c['zero_training_weight_absmax']} "
              f"(the blindness, MEASURED not assumed)")
        print(f"  carry >=1 unseen substitution: {e['frac_blindspot_with_unseen']:.1%}  "
              f"(R {e['frac_R_with_unseen']} vs S {e['frac_S_with_unseen']})")
        print(f"  mean unseen count: R {e['mean_unseen_R']} vs S {e['mean_unseen_S']}")
        for nm, st in c["strata"].items():
            extra = f"AUROC {st['auroc']}" if st["auroc"] is not None else f"WITHHELD -- {st['refused']}"
            print(f"  {nm:10s} n={st['n']:5d} R={st['n_R']:4d} S={st['n_S']:4d}  {extra}")
        v = verdict(c)
        c["verdict_block"] = v
        print(f"  VERDICT (frozen rule): {v['verdict']}  (exposed {v.get('exposed_frac')}, "
              f"gap {v.get('auroc_gap')}; need >={MIN_EXPOSED_FRAC} and >={AUROC_GAP})")

    art = {"schema": "hiv-onehot-unseen-headroom-v1", "date": str(date.today()),
           "question": ("one-hot treats every (position, amino-acid) pair as an INDEPENDENT column, so a "
                        "substitution with no training signal gets ~zero weight. Does that blindness "
                        "actually bite on the blind spot -- i.e. is there room for an orthogonal "
                        "position-level prior (structure/conservation) to add anything?"),
           "frozen_rule": FROZEN_RULE, "anchor": anchors,
           "floor_artifact": str(FLOOR_ARTIFACT.relative_to(REPO)).replace("\\", "/"),
           "honest_scope": [
               "IN-DISTRIBUTION to the Stanford knowledge base; ~96% subtype B.",
               ("The one-hot FEATURE SET is built globally, so a zero-training-signal column still EXISTS; "
                "what is absent is training signal, and the reported max |coef| over those columns "
                "measures the resulting blindness directly rather than assuming it."),
               ("Exposure is measured per FOLD against that fold's own training rows, which is the only "
                "honest denominator -- a globally-rare substitution may be well-covered in some folds."),
               ("This bounds the room for a position-level prior. It does NOT say such a prior would "
                "capture the room, only whether any room exists to capture."),
           ],
           "cells": cells, "refused": refused}
    Path(a.out).write_text(json.dumps(art, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
