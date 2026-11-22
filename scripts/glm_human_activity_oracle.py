"""The HUMAN activity oracle — and its OWN held-out accuracy, reported FIRST as a gate.

WHAT THIS IS. A supervised sequence -> measured-activity predictor for human regulatory elements,
trained on Agarwal lentiMPRA (`LentiRecord.activity` = observed log2(RNA/DNA)). It is the functional
scorer this project has never had on a substrate where it can be made INDEPENDENT of the generator.

WHY IT IS THE FIRST THING RUN, AND A GATE. A weak oracle cannot validate a generator's designs: its
verdict on generated sequence is uninterpretable in exactly the way a sub-chance eval set was on the
E. coli arm, where the model and both Markov nulls all scored worse than uniform and the model-vs-null
delta carried no information. So oracle validity is CONTROL ZERO. If the oracle does not clearly beat
a crude baseline on held-out folds, no generator number should be computed at all.

THE FOLD CONTRACT — three disjoint roles, taken from the PUBLISHERS' own 10-fold assignment, never
invented here:

    folds 1-6   GENERATOR training          (reserved; this script does not touch them)
    folds 7-8   ORACLE training
    folds 9-10  HELD OUT -- the oracle's test set AND, later, the natural comparison set

The oracle is therefore trained on folds the generator will never see, which makes it an INDEPENDENT
judge. That was structurally impossible on the E. coli arm, where the only available scorer was
trained on the same genome family, and it is the single reason the human substrate can answer a
design question at all.

GEOMETRY. The publishers' sequence is `15nt 5' adaptor - 200nt element - 15nt 3' adaptor` (the column
name says so verbatim). The PRIMARY geometry is the 200 nt ELEMENT, adaptors stripped, fixed in
advance on a priori grounds: the flanks are measurably CONSTANT -- one distinct value each, share
1.0000, identical across both cell types, and a Markov chain reproduced them 0 times in 120 -- so
they carry zero natural variability and cannot help characterise real sequence. Leaving them in
inflated a distinguishability margin by up to 0.43 on this very substrate
(`wiki/human_adaptor_probe_2026-11-21.json`). The full 230 nt is reported as a secondary row so the
choice is auditable, not hidden.

CELLS ARE NEVER POOLED. HepG2 and K562 are different assays in different cells, and this project has
already measured them CONTRADICTING each other on a related verdict. Pooling would average a
disagreement into a number that describes neither.

WHAT A GOOD NUMBER HERE DOES NOT ESTABLISH. The oracle predicting activity well does not mean a
generator can design to it. It means the judge is competent. The design question is separate and
gated on this one.
"""
from __future__ import annotations

import argparse
import json
import math
import statistics
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from dna_decode.glm.lentimpra import load_cell  # noqa: E402

#: The data is on C: on this host. `lentimpra.DEFAULT_DIR` points at D:/dna_decode_cache/mpra_human,
#: which does not exist here, so a bare load_cell() reports "absent, fetch it" for data already on
#: disk. Passing the directory explicitly is mandatory, not stylistic.
LOCAL_DIR = ROOT / "data" / "mpra_human"
OUT = ROOT / "wiki" / "glm_human_activity_oracle_2026-11-22.json"

ORACLE_TRAIN_FOLDS = {7, 8}
HELD_OUT_FOLDS = {9, 10}
GENERATOR_FOLDS = {1, 2, 3, 4, 5, 6}        # reserved, deliberately untouched here

ADAPTOR = 15                                 # 5' and 3', per the publishers' column name
BASES = "ACGT"

#: Pre-registered gate. The oracle must beat GC by this margin on held-out folds, or its verdict on a
#: generator is uninterpretable. Mirrors the bar the E. coli genome-wide oracle was graded against
#: (`beats_crude_baseline_by: 0.05`), imported rather than re-invented.
MIN_MARGIN_OVER_GC = 0.05
#: And it must reach at least this absolute Spearman, or "beats GC" is a race between two weak things.
MIN_ABS_SPEARMAN = 0.20


class DegenerateSpearman(ValueError):
    """Raised when one side has zero rank variance, so no correlation is defined.

    WHY THIS RAISES INSTEAD OF RETURNING A NUMBER. This function used to return 0.0 here, and 0.0 is a
    perfectly ordinary Spearman value, so THREE different conditions arrived as one indistinguishable
    measurement: a dead conditioning knob (constant predictions), a broken run where no target varies
    (constant targets), and total plumbing failure (both). The middle one is the dangerous one -- a
    generator run that never varied its target would have reported rho 0.0 and read as the bar's
    NO_TARGET_CONTROL outcome, i.e. a broken run wearing the costume of a null result.

    The repo's own convention already agreed: five other copies of this same guard
    (`hiv_esm_likelihood_probe`, `humsavar_am_vs_esm`, `kaggle_proteingym_sweep`, `tem1_forward_cell`,
    `protein_effect_facevalidity`) return NaN or None -- a value nobody can mistake for a result. This
    arm was the lone copy returning 0.0. Raising is strictly stronger than NaN because `side` names
    WHICH input is constant, and the bar needs a dead knob told apart from a broken harness.
    """

    def __init__(self, side: str, n: int):
        self.side = side                 # "x" (targets) | "y" (predictions) | "both"
        self.n = n
        super().__init__(
            "Spearman undefined: %s has zero rank variance over n=%d. "
            "side='y' means the scorer returned a CONSTANT prediction (check the conditioning signal "
            "actually moves the output before reading this as a null); side='x' means the requested "
            "targets never varied, which is a harness defect, not a result." % (side, n))


def spearman(x, y) -> float:
    def ranks(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(order):               # average ties, or sort order becomes a finding
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1.0
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r
    rx, ry = ranks(x), ranks(y)
    mx, my = statistics.fmean(rx), statistics.fmean(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    dx = math.sqrt(sum((a - mx) ** 2 for a in rx))
    dy = math.sqrt(sum((b - my) ** 2 for b in ry))
    if not dx or not dy:
        raise DegenerateSpearman("both" if (not dx and not dy) else ("x" if not dx else "y"), len(rx))
    return float(num / (dx * dy))


def gc(seq: str) -> list[float]:
    return [(seq.count("G") + seq.count("C")) / max(len(seq), 1)]


_KMER_IDX: dict[int, dict[str, int]] = {}


def _kmer_index(k: int) -> dict[str, int]:
    """Cached. Rebuilding a 4**6 = 4096-entry index per sequence would cost ~260M dict builds
    across the run -- minutes of pure overhead for a table that never changes."""
    if k not in _KMER_IDX:
        import itertools
        _KMER_IDX[k] = {"".join(p): i
                        for i, p in enumerate(itertools.product(BASES, repeat=k))}
    return _KMER_IDX[k]


def kmer_feats(seq: str, k: int) -> list[float]:
    import numpy as np
    idx = _kmer_index(k)
    v = np.zeros(len(idx), dtype=np.float32)
    n = 0
    for i in range(len(seq) - k + 1):
        j = idx.get(seq[i:i + k])
        if j is not None:
            v[j] += 1.0
            n += 1
    if n:
        v /= n
    return v.tolist()


def positional_onehot(seq: str, bins: int) -> list[float]:
    """Mean base composition per positional bin. REFUSES a non-dividing bin count rather than
    truncating or padding, because a ragged final bin changes what a bin index means."""
    import numpy as np
    L = len(seq)
    if L % bins:
        raise ValueError(f"{L} positions do not divide into {bins} bins")
    w = L // bins
    out = np.zeros((bins, 4), dtype=np.float32)
    for b in range(bins):
        chunk = seq[b * w:(b + 1) * w]
        for ci, base in enumerate(BASES):
            out[b, ci] = chunk.count(base) / w
    return out.reshape(-1).tolist()


def build(seqs, mode: str, bins: int):
    if mode == "gc":
        return [gc(s) for s in seqs]
    if mode == "kmer4":
        return [kmer_feats(s, 4) for s in seqs]
    if mode == "kmer6":
        return [kmer_feats(s, 6) for s in seqs]
    if mode == "positional":
        return [positional_onehot(s, bins) for s in seqs]
    raise ValueError(mode)


def ridge_spearman(Xtr, ytr, Xte, yte, alpha: float = 1.0) -> float:
    import numpy as np
    A = np.asarray(Xtr, dtype=np.float64)
    B = np.asarray(Xte, dtype=np.float64)
    mu, sd = A.mean(0), A.std(0)
    sd[sd == 0] = 1.0
    A = (A - mu) / sd
    B = (B - mu) / sd
    A = np.hstack([A, np.ones((A.shape[0], 1))])
    B = np.hstack([B, np.ones((B.shape[0], 1))])
    y = np.asarray(ytr, dtype=np.float64)
    reg = alpha * np.eye(A.shape[1])
    reg[-1, -1] = 0.0                      # never penalise the intercept
    w = np.linalg.solve(A.T @ A + reg, A.T @ y)
    return spearman((B @ w).tolist(), list(yte))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cells", default="HepG2,K562")
    ap.add_argument("--bins", type=int, default=10)
    ap.add_argument("--alpha", type=float, default=1.0)
    a = ap.parse_args()

    if not LOCAL_DIR.exists():
        print(f"REFUSE: {LOCAL_DIR} absent", file=sys.stderr)
        return 2

    res = {
        "record": "glm-human-activity-oracle-v1",
        "date": "2026-11-22",
        "substrate": "Agarwal lentiMPRA (Zenodo 10558183, CC-BY-4.0), observed log2(RNA/DNA)",
        "data_dir": str(LOCAL_DIR),
        "default_dir_note": ("lentimpra.DEFAULT_DIR is D:/dna_decode_cache/mpra_human and does NOT "
                             "exist on this host; the data is on C: at data/mpra_human, so the "
                             "directory must be passed explicitly or the loader reports absent data"),
        "fold_contract": {"generator_reserved": sorted(GENERATOR_FOLDS),
                          "oracle_train": sorted(ORACLE_TRAIN_FOLDS),
                          "held_out": sorted(HELD_OUT_FOLDS),
                          "why": "publishers' own 10-fold assignment; the oracle never sees the "
                                 "generator's folds, which makes it an INDEPENDENT judge"},
        "primary_geometry": "element_200nt_adaptors_stripped",
        "gate": {"min_margin_over_gc": MIN_MARGIN_OVER_GC,
                 "min_abs_spearman": MIN_ABS_SPEARMAN,
                 "why": "a weak oracle cannot validate a generator's designs; its verdict would be "
                        "uninterpretable the way a sub-chance eval set was on the E. coli arm"},
        "cells": {},
    }
    t0 = time.time()

    for cell in [c.strip() for c in a.cells.split(",") if c.strip()]:
        recs, stats = load_cell(cell, LOCAL_DIR)
        tr = [r for r in recs if r.fold in ORACLE_TRAIN_FOLDS]
        te = [r for r in recs if r.fold in HELD_OUT_FOLDS]
        print(f"{cell}: {len(recs)} records | oracle-train {len(tr)} | held-out {len(te)}")
        assert not ({r.fold for r in tr} & {r.fold for r in te}), "fold roles overlap"

        cell_out = {"n_records": len(recs), "n_oracle_train": len(tr), "n_held_out": len(te),
                    "sequence_column": stats.sequence_column, "geometries": {}}

        for geo, cut in (("element_200nt_adaptors_stripped", ADAPTOR), ("full_230nt_as_shipped", 0)):
            strip = (lambda s: s[cut:len(s) - cut]) if cut else (lambda s: s)
            Str, Ste = [strip(r.seq) for r in tr], [strip(r.seq) for r in te]
            ytr, yte = [r.activity for r in tr], [r.activity for r in te]
            L = len(Str[0])
            bins = a.bins if L % a.bins == 0 else (8 if L % 8 == 0 else 5)
            arms = {}
            for mode in ("gc", "kmer4", "kmer6", "positional"):
                s = ridge_spearman(build(Str, mode, bins), ytr, build(Ste, mode, bins), yte, a.alpha)
                arms[mode] = round(s, 4)
                print(f"   {geo:36s} {mode:11s} rho {s:+.4f}")
            best_learned = max(arms[m] for m in ("kmer4", "kmer6", "positional"))
            margin = round(best_learned - arms["gc"], 4)
            cell_out["geometries"][geo] = {
                "length_bp": L, "positional_bins": bins, "arms": arms,
                "best_learned": best_learned, "margin_over_gc": margin,
                "passes_gate": bool(margin >= MIN_MARGIN_OVER_GC
                                    and best_learned >= MIN_ABS_SPEARMAN),
            }
        res["cells"][cell] = cell_out

    pri = [res["cells"][c]["geometries"]["element_200nt_adaptors_stripped"] for c in res["cells"]]
    all_pass = all(g["passes_gate"] for g in pri)
    any_pass = any(g["passes_gate"] for g in pri)
    if all_pass:
        res["verdict"] = "ORACLE_VALID_PROCEED"
        res["verdict_detail"] = (
            "Every cell clears the pre-registered gate on the primary geometry, so the oracle is a "
            "competent judge and a generator design question is interpretable against it. This does "
            "NOT show a generator can design to it -- only that the judge works.")
    elif any_pass:
        res["verdict"] = "ORACLE_VALID_IN_SOME_CELLS_ONLY"
        res["verdict_detail"] = (
            "The gate passes in some cells and fails in others. Cells are NOT pooled, so proceed "
            "only in the passing cells and name them; a generator verdict in a failing cell would "
            "be uninterpretable.")
    else:
        res["verdict"] = "ORACLE_TOO_WEAK_DO_NOT_SCORE_A_GENERATOR"
        res["verdict_detail"] = (
            "No cell clears the gate on the primary geometry. A generator number scored against "
            "this oracle would be uninterpretable. REFUSING to proceed to the design question; the "
            "oracle is the thing to fix first.")

    res["seconds"] = round(time.time() - t0, 1)
    print("\n" + "=" * 74)
    print(f"VERDICT: {res['verdict']}\n  {res['verdict_detail']}")
    OUT.write_text(json.dumps(res, indent=2), encoding="utf-8")
    print(f"\n-> {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
