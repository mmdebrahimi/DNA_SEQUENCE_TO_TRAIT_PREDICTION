"""The mandated baseline nobody ran: a one-hot sequence ridge on the 44,106 peak tiles.

WHY THIS EXISTS. The project's mandatory baseline panel names `one-hot ridge` for ANY learned
number. Measured 2026-11-28 against `soraya_runs/arms-ledger.json`: the string `one-hot` appears
**ZERO** times across all sixteen arms of the activity ceiling. The four `ridge` hits are
`promcalc-full-decomposition-ridge-150bp`, a ridge over 55 BIOPHYSICAL features -- a different
object. So every arm on that criterion is benchmarked against `gc-single-feature` (0.3112) as its
floor, and a 600-feature linear sequence model -- cheaper than every arm above it -- was never in
the comparison.

That matters in one direction only. A missing WEAK baseline cannot make an arm look worse; it can
only have let arms claim gains over a floor that was too low. So this is a back-run on the existing
claims, not a new arm competing for the ceiling.

COMPARABILITY IS THE WHOLE POINT, so nothing here is re-chosen:
  * tiles come from the OWNER, `dna_decode.glm.genomewide.load_peak_tiles`, with the arm's own
    `category == 'tile' and len(seq) == 150` filter
  * splits are the arm's own `leave_peak_out(frac=0.25, seed=0..9)` -- peak-level, disjoint BY
    CONSTRUCTION because membership is a pure function of the peak id
  * the estimator is the arm's own `ridge_spearman`, called bare at its default `alpha=1.0`, which
    is exactly how both `gc-single-feature` and the 128-dim atlas embedding were run

THE GC RECONCILE IS A GATE, NOT A FORMALITY. Before any new number is read, seed 0 must reproduce
0.3000 and the 10-seed median 0.3112, both within the arm's own 0.01 tolerance. Reproducing a known
number to four decimals validates the loader, the filter, the split, the log1p target, the ridge and
the scorer TOGETHER. A pipeline that cannot reproduce a known number cannot establish a new one, so
this REFUSES (exit 2) rather than printing a caveat. The pattern is lifted from
`scripts/promoteratlas_local_reconcile.py`, which owns it.

THE BAR IS FROZEN HERE, BEFORE ANY NUMBER EXISTS (see PREREGISTERED), with all three outcomes named
and the consequence of each spelled out. Alpha is NOT tuned for the primary: a swept alpha is
reported as a clearly-labelled SECONDARY diagnostic, tuned on a split of TRAIN only so it never sees
test, and it is explicitly not the comparable number.

WHAT THIS CANNOT SHOW. It bounds what a LINEAR ONE-HOT model reaches on this substrate under this
protocol. It is not a claim about any arm's architecture, and beating an arm here does not mean that
arm is worthless -- it means that arm's reported GAIN was measured against a floor that was lower
than the panel requires.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

ARM = REPO / "wiki/promoteratlas_tile_arm_2026-11-01.json"
LEDGER = REPO / "soraya_runs/arms-ledger.json"
CACHE = REPO / "data/mpra"

EXPECTED_TILES = 44106
N_SEEDS = 10
FRAC = 0.25
GC_TOL = 0.01
ALPHA_PRIMARY = 1.0            # the arm's own default -- NOT re-chosen
#: Extended upward after the first run: EVERY seed chose the grid's MAXIMUM (1000), which means the
#: tuning wanted more regularization than was offered and the grid was saturated at its ceiling. A
#: saturated grid does not bound the model, so the conclusion "one-hot does not reach GC" would have
#: rested on an un-bounded sweep. Extended rather than asserted.
ALPHA_SWEEP = (1.0, 10.0, 100.0, 1_000.0, 10_000.0, 100_000.0, 1_000_000.0)

#: Arms on a DIFFERENT substrate cannot be compared to a 150 bp tile Spearman. The ledger mixes
#: substrates under one criterion key and says so itself -- `fragment-rung-321k-rows` carries the
#: verbatim warning "this criterion quotes the 0.941 TILE ceiling; this arm measures the FRAGMENT
#: rung, whose derived ceiling is 0.7928 - do not compare the two fractions without saying which
#: ceiling." The first version of this script walked straight past that and reported three
#: "exceeded" arms, two of which were fragment-substrate. Excluded here WITH the exclusions recorded,
#: so the filter is auditable rather than silent.
OFF_SUBSTRATE_MARKERS = ("fragment", "499bp")
WINDOW = 150
BASES = "ACGT"

#: Frozen before any number existed. All three outcomes named, each with its consequence.
PREREGISTERED = {
    "margin_nats_na": "metric is Spearman, not nats",
    "margin": 0.02,
    "margin_derivation": ("INHERITED verbatim from the arm acceptance bars so this stays "
                          "commensurable with every arm it is compared against. Not re-chosen."),
    "outcomes": {
        "ADDS_NOTHING_OVER_GC": {
            "when": "one-hot median <= gc median + 0.02",
            "consequence": ("the existing baseline ladder stands; no arm's claim changes, and the "
                            "panel's one-hot requirement is satisfied by a number that happens to "
                            "add nothing. This is the outcome that CLEARS the back-run."),
        },
        "RE_READS_CLAIMS_BELOW_IT": {
            "when": ("one-hot median > gc median + 0.02 AND it exceeds at least one arm that "
                     "reported a GAIN"),
            "consequence": ("every arm it exceeds had its gain measured against a floor lower than "
                            "the panel requires; those claims must be restated against this "
                            "baseline before they mean what they say. The arms are NOT invalidated "
                            "-- their reported MARGIN is."),
        },
        "DOMINATES_FROZEN_EMBEDDING": {
            "when": "one-hot median >= 0.5076 - 0.02",
            "consequence": ("a 600-feature linear model matches the frozen 128-dim atlas embedding, "
                            "which would make the embedding arm's contribution the open question "
                            "rather than its headline. The largest possible finding here."),
        },
    },
    "gate": "GC seed0 0.3000 and median 0.3112, both within 0.01, or REFUSE with exit 2",
    "alpha_primary_is_not_tuned": True,
}


def _fail(code: int, msg: str) -> int:
    print("REFUSED: " + msg, flush=True)
    return code


def load_tiles():
    """The OWNER's loader plus the arm's own filter. Not a second loader."""
    from dna_decode.glm import genomewide as gw
    tiles_all, st = gw.load_peak_tiles(CACHE, allow_download=False)
    tiles = [t for t in tiles_all if t.category == "tile" and len(t.seq) == WINDOW]
    print("substrate: %d rows parsed (delim=%s), %d after filter(category=='tile', len==%d)"
          % (st.n_parsed, st.delimiter, len(tiles), WINDOW), flush=True)
    return tiles


def onehot(seqs: list[str]) -> np.ndarray:
    """150 positions x 4 bases = 600 binary features. Any non-ACGT base contributes nothing."""
    ix = {b: i for i, b in enumerate(BASES)}
    X = np.zeros((len(seqs), WINDOW * 4), dtype=np.float32)
    for r, s in enumerate(seqs):
        for p, ch in enumerate(s[:WINDOW]):
            j = ix.get(ch)
            if j is not None:
                X[r, p * 4 + j] = 1.0
    return X


def committed_arms() -> tuple[dict[str, float], dict[str, str]]:
    """Measured arms, split into same-substrate (comparable) and off-substrate (excluded + reason)."""
    d = json.loads(LEDGER.read_text(encoding="utf-8"))
    crit = [k for k in d["criteria"] if k.startswith("ceiling: best achievable")][0]
    same: dict[str, float] = {}
    excluded: dict[str, str] = {}
    for name, a in d["criteria"][crit].items():
        if a.get("status") != "measured" or a.get("metric") is None:
            continue
        # NAME only, deliberately. Matching the EVIDENCE blob over-fired: it excluded
        # `finetuned-promoteratlas-encoder-150bp` -- a genuine 150 bp tile arm -- because its
        # evidence text mentions a fragment arm in passing. An arm's substrate is declared in its
        # NAME (`...-150bp`, `...-499bp-fragments`); the evidence is prose that may discuss others.
        hay = name.lower()
        hit = next((m for m in OFF_SUBSTRATE_MARKERS if m in hay), None)
        if hit:
            excluded[name] = (
                "off-substrate (%r in name/evidence): a 150 bp TILE Spearman is not comparable to a "
                "FRAGMENT-ceiling fraction; the ledger's own SCOPE NOTE forbids it" % hit)
        else:
            same[name] = float(a["metric"])
    return same, excluded


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sweep-alpha", action="store_true",
                    help="also report the SECONDARY swept-alpha diagnostic (train-internal tuning)")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    if not ARM.is_file():
        return _fail(3, "committed arm artifact absent at %s" % ARM)
    arm = json.loads(ARM.read_text(encoding="utf-8"))
    want_gc0 = float(arm["gc_anchor"]["seed0"])
    want_gcmed = float(arm["gc_anchor"]["median"])
    print("committed GC targets: seed0=%.4f median=%.4f (tol %.2f)" % (want_gc0, want_gcmed, GC_TOL),
          flush=True)

    from dna_decode.glm import genomewide as gw
    from dna_decode.glm.atlas_embed import ridge_spearman

    try:
        tiles = load_tiles()
    except Exception as exc:
        return _fail(3, "%s: %s" % (type(exc).__name__, exc))
    if len(tiles) != EXPECTED_TILES:
        return _fail(2, "%d tiles, the arm scored %d -- substrate differs, nothing is comparable"
                     % (len(tiles), EXPECTED_TILES))

    seqs = [t.seq for t in tiles]
    gc = np.asarray(gw.gc_features(seqs), dtype=float)
    t0 = time.time()
    OH = onehot(seqs)
    print("one-hot: %s built in %.1fs" % (OH.shape, time.time() - t0), flush=True)

    idx = {id(t): i for i, t in enumerate(tiles)}
    gc_per, oh_per, sweep_per, n_test = {}, {}, {}, {}

    for seed in range(N_SEEDS):
        train, test, held = gw.leave_peak_out(tiles, frac=FRAC, seed=seed)
        tr = [idx[id(t)] for t in train]
        te = [idx[id(t)] for t in test]
        ytr = [tiles[i].expression for i in tr]
        yte = [tiles[i].expression for i in te]
        n_test[seed] = len(te)

        g, _ = ridge_spearman(gc[tr], ytr, gc[te], yte, alpha=ALPHA_PRIMARY)
        o, _ = ridge_spearman(OH[tr], ytr, OH[te], yte, alpha=ALPHA_PRIMARY)
        gc_per[seed], oh_per[seed] = round(float(g), 4), round(float(o), 4)

        line = "   seed %d: n_test=%-6d held_peaks=%-4d gc=%.4f onehot=%.4f" % (
            seed, len(te), held, g, o)

        if a.sweep_alpha:
            # TRAIN-INTERNAL tuning only. The test split is never consulted, so this cannot leak.
            cut = int(len(tr) * 0.75)
            itr, iva = tr[:cut], tr[cut:]
            ytr2 = [tiles[i].expression for i in itr]
            yva = [tiles[i].expression for i in iva]
            best_al, best_s = None, None
            for al in ALPHA_SWEEP:
                s, _ = ridge_spearman(OH[itr], ytr2, OH[iva], yva, alpha=al)
                if best_s is None or s > best_s:
                    best_al, best_s = al, s
            sv, _ = ridge_spearman(OH[tr], ytr, OH[te], yte, alpha=best_al)
            sweep_per[seed] = {"alpha": best_al, "spearman": round(float(sv), 4)}
            line += " onehot_swept=%.4f(a=%g)" % (sv, best_al)

        print(line, flush=True)

    gc_med = float(np.median([gc_per[s] for s in range(N_SEEDS)]))
    oh_med = float(np.median([oh_per[s] for s in range(N_SEEDS)]))

    # ---- THE GATE. A pipeline that cannot reproduce a known number cannot establish a new one.
    gate = []
    for label, got, want in (("gc.seed0", gc_per[0], want_gc0), ("gc.median", gc_med, want_gcmed)):
        ok = abs(got - want) <= GC_TOL
        gate.append({"check": label, "got": round(got, 4), "want": want, "tol": GC_TOL, "pass": ok})
        print("gate %-10s got %.4f want %.4f  %s" % (label, got, want, "PASS" if ok else "FAIL"),
              flush=True)
    if not all(c["pass"] for c in gate):
        return _fail(2, "GC reconcile FAILED -- this pipeline cannot reproduce a committed number, "
                        "so its one-hot figure of %.4f means nothing and is withheld" % oh_med)

    # ---- the back-run comparison, only now that the gate passed
    arms, off_substrate = committed_arms()
    exceeded = sorted([(n, m) for n, m in arms.items() if oh_med > m + PREREGISTERED["margin"]],
                      key=lambda kv: -kv[1])
    beaten_by = sorted([(n, m) for n, m in arms.items() if m > oh_med + PREREGISTERED["margin"]],
                       key=lambda kv: -kv[1])
    if off_substrate:
        print("excluded %d off-substrate arm(s) from the comparison:" % len(off_substrate),
              flush=True)
        for n, why in off_substrate.items():
            print("    %s -- %s" % (n, why), flush=True)

    frozen_emb = arms.get("promoteratlas-base-embeddings-frozen-150bp-padded-to-200", 0.5076)
    if oh_med >= frozen_emb - PREREGISTERED["margin"]:
        verdict = "DOMINATES_FROZEN_EMBEDDING"
    elif oh_med > gc_med + PREREGISTERED["margin"] and exceeded:
        verdict = "RE_READS_CLAIMS_BELOW_IT"
    else:
        verdict = "ADDS_NOTHING_OVER_GC"

    res = {
        "record": "glm-onehot-ridge-backrun-v1",
        "date": "2026-11-28",
        "question": ("the mandatory panel names one-hot ridge and the arms ledger contains it ZERO "
                     "times -- what happens to the existing claims when it is actually run?"),
        "PREREGISTERED": PREREGISTERED,
        "protocol": {
            "loader": "dna_decode.glm.genomewide.load_peak_tiles (the OWNER)",
            "filter": "category == 'tile' and len(seq) == 150",
            "split": "dna_decode.glm.genomewide.leave_peak_out(frac=0.25, seed=0..9)",
            "estimator": "dna_decode.glm.atlas_embed.ridge_spearman, alpha=1.0 (the arm's default)",
            "features": "150 positions x 4 bases = 600 binary",
            "n_tiles": len(tiles),
        },
        "gc_reconcile_gate": gate,
        "gc": {"per_seed": gc_per, "median": round(gc_med, 4)},
        "onehot": {"per_seed": oh_per, "median": round(oh_med, 4), "alpha": ALPHA_PRIMARY},
        "onehot_swept_SECONDARY": (
            {"per_seed": sweep_per,
             "median": round(float(np.median([v["spearman"] for v in sweep_per.values()])), 4),
             "grid": list(ALPHA_SWEEP),
             "alphas_chosen": sorted({v["alpha"] for v in sweep_per.values()}),
             "grid_saturated_at_max": all(v["alpha"] == max(ALPHA_SWEEP)
                                          for v in sweep_per.values()),
             "why_saturation_matters": ("if every seed picks the grid MAXIMUM the tuning wanted more "
                                        "regularization than offered, so the sweep does NOT bound "
                                        "the model and cannot support 'this is not a regularization "
                                        "artifact'. The grid was extended once for exactly this."),
             "note": ("tuned on a TRAIN-INTERNAL 75/25 split, test never consulted. NOT the "
                      "comparable number -- every arm was run at alpha=1.0.")}
            if a.sweep_alpha else {"note": "not run; pass --sweep-alpha"}),
        "n_test_per_seed": n_test,
        "committed_arms_same_substrate": arms,
        "excluded_off_substrate": off_substrate,
        "onehot_exceeds_these_arms": [{"arm": n, "metric": m} for n, m in exceeded],
        "these_arms_exceed_onehot": [{"arm": n, "metric": m} for n, m in beaten_by],
        "VERDICT": verdict,
        "consequence": PREREGISTERED["outcomes"][verdict]["consequence"],
        "honest_limits": [
            "Bounds a LINEAR ONE-HOT model on this substrate under this protocol. Not a claim about "
            "any arm's architecture.",
            "Exceeding an arm does NOT invalidate that arm -- it invalidates the MARGIN that arm "
            "reported, because the margin was measured against a floor lower than the panel requires.",
            "alpha=1.0 is the arm's own default and is deliberately NOT tuned for the primary. A "
            "600-feature model and a 1-feature model do not share an optimal penalty, so the swept "
            "number is reported separately and is not comparable to the arms.",
            "One substrate (44,106 peak tiles), one geometry (native 150 bp), one metric (Spearman).",
        ],
    }

    out = Path(a.out) if a.out else REPO / "wiki" / "glm_onehot_ridge_backrun_2026-11-28.json"
    out.write_text(json.dumps(res, indent=1) + "\n", encoding="utf-8")
    print("\nVERDICT %s | gc %.4f -> onehot %.4f (margin %.2f)"
          % (verdict, gc_med, oh_med, PREREGISTERED["margin"]))
    if exceeded:
        print("one-hot EXCEEDS %d committed arm(s):" % len(exceeded))
        for n, m in exceeded:
            print("    %.4f  %s" % (m, n))
    print("wrote %s" % out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
