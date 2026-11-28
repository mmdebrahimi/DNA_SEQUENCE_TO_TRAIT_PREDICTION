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
#: Extended TWICE. The first grid topped out at 1e3 and every seed chose the maximum; the second
#: topped out at 1e6 and still did. A saturated grid does not BOUND the model, so any conclusion
#: drawn from one is drawn from an un-bounded sweep. This grid reaches the real peak: the curve rises
#: to 1e9 and then TURNS DOWN (1e10 0.3133, 1e11 0.3096, 1e12 0.0850), which is the intercept-only
#: limit a ridge must collapse to as alpha grows. `assert_peak_not_at_grid_edge` enforces it.
ALPHA_SWEEP = (1.0, 1e4, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11, 1e12)

#: STRAND NORMALIZATION, and it was a real defect in the first version of this probe. 22,285 of the
#: 44,106 tiles are minus-strand and `load_peak_tiles` stores the raw `variant` column with NO
#: reverse-complement applied. So position p meant opposite things across half the substrate, which
#: cripples a POSITIONAL linear model while leaving a learned encoder largely unaffected -- i.e. it
#: biased this probe toward a FALSE negative, the one direction that matters here. Measured cost:
#: +0.0129 at alpha=1, +0.0006 at the peak. GC is orientation-invariant by construction (a G+C count
#: survives reverse-complement), so the comparator needs no equivalent fix and gets none.
NORMALIZE_STRAND = True

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


_COMP = str.maketrans("ACGT", "TGCA")


def strand_normalized(tiles) -> list[str]:
    """Reverse-complement minus-strand tiles so position p means the same thing on every row."""
    return [t.seq.translate(_COMP)[::-1] if t.strand == "-" else t.seq for t in tiles]


def onehot(seqs: list[str]) -> np.ndarray:
    """150 positions x 4 bases = 600 binary features. Any non-ACGT base contributes nothing.

    All four bases are kept deliberately. Reference coding (drop one base per position, 450 columns)
    was measured and is WORSE here -- 0.2826 against this encoding's 0.3135 at the same alpha. With
    StandardScaler in front of the ridge, the redundant fourth column is what lets the model express
    aggregate composition cleanly under heavy shrinkage, and removing it costs more than the rank
    deficiency does. So the collinearity is real and is NOT the thing to fix.
    """
    ix = {b: i for i, b in enumerate(BASES)}
    X = np.zeros((len(seqs), WINDOW * 4), dtype=np.float32)
    for r, s in enumerate(seqs):
        for p, ch in enumerate(s[:WINDOW]):
            j = ix.get(ch)
            if j is not None:
                X[r, p * 4 + j] = 1.0
    return X


def assert_peak_not_at_grid_edge(by_alpha: dict[float, float]) -> dict:
    """A best-alpha sitting at the grid MAXIMUM means the sweep never bounded the model.

    This probe reported an un-bounded number twice before this check existed. A ridge must collapse
    toward intercept-only as alpha grows, so a real optimum is INTERIOR and the curve turns down
    after it. If the best alpha is the largest tried, the sweep is incomplete and says so.
    """
    alphas = sorted(by_alpha)
    best = max(alphas, key=lambda a: by_alpha[a])
    at_edge = best == alphas[-1]
    return {"best_alpha": best, "best": round(by_alpha[best], 6),
            "curve": {f"{a:g}": round(by_alpha[a], 6) for a in alphas},
            "peak_is_interior": not at_edge,
            "grid_saturated": at_edge,
            "what_saturation_means": ("the grid maximum won, so this sweep did NOT bound the model "
                                      "and its number must not be read as the one-hot ceiling"),
            "collapse_check": ("a ridge tends to intercept-only as alpha grows, so a trustworthy "
                               "sweep shows the curve turning DOWN past the peak")}


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

    raw_seqs = [t.seq for t in tiles]
    norm_seqs = strand_normalized(tiles)
    n_rc = sum(1 for t in tiles if t.strand == "-")
    gc = np.asarray(gw.gc_features(raw_seqs), dtype=float)   # orientation-invariant; see NORMALIZE_STRAND
    t0 = time.time()
    OH = onehot(norm_seqs if NORMALIZE_STRAND else raw_seqs)
    OH_RAW = onehot(raw_seqs)
    print("one-hot: %s built in %.1fs | strand-normalized %d of %d tiles (minus-strand)"
          % (OH.shape, time.time() - t0, n_rc, len(tiles)), flush=True)

    idx = {id(t): i for i, t in enumerate(tiles)}
    splits, n_test, held_peaks = [], {}, {}
    for seed in range(N_SEEDS):
        train, test, held = gw.leave_peak_out(tiles, frac=FRAC, seed=seed)
        splits.append(([idx[id(t)] for t in train], [idx[id(t)] for t in test]))
        n_test[seed], held_peaks[seed] = len(test), held

    def med(M, alpha: float) -> tuple[float, dict]:
        per = {}
        for seed, (tr, te) in enumerate(splits):
            ytr = [tiles[i].expression for i in tr]
            yte = [tiles[i].expression for i in te]
            v, _ = ridge_spearman(M[tr], ytr, M[te], yte, alpha=alpha)
            per[seed] = round(float(v), 4)
        return float(np.median([per[s] for s in range(N_SEEDS)])), per

    gc_med, gc_per = med(gc, ALPHA_PRIMARY)
    print("   GC (alpha=%g): %.4f" % (ALPHA_PRIMARY, gc_med), flush=True)

    # Sweep the FULL grid on the median, not per-seed: the question is where the one-hot model peaks,
    # and a per-seed train-internal pick cannot show the curve turning down -- which is the only
    # evidence that the sweep bounded anything at all.
    curve, oh_curve_per = {}, {}
    for al in ALPHA_SWEEP:
        m, per = med(OH, al)
        curve[al], oh_curve_per[al] = m, per
        print("   onehot alpha=%-8g %.4f" % (al, m), flush=True)
    peak = assert_peak_not_at_grid_edge(curve)
    oh_med = peak["best"]
    oh_per = oh_curve_per[peak["best_alpha"]]

    # The strand defect, QUANTIFIED rather than asserted: same splits, raw vs normalized.
    raw_at_1, _ = med(OH_RAW, ALPHA_PRIMARY)
    raw_at_peak, _ = med(OH_RAW, peak["best_alpha"])
    strand_effect = {
        "n_minus_strand_tiles": n_rc, "n_tiles": len(tiles),
        "raw_at_alpha_1": round(raw_at_1, 6),
        "normalized_at_alpha_1": round(curve[ALPHA_PRIMARY], 6),
        "gain_at_alpha_1": round(curve[ALPHA_PRIMARY] - raw_at_1, 6),
        "raw_at_peak": round(raw_at_peak, 6),
        "normalized_at_peak": round(oh_med, 6),
        "gain_at_peak": round(oh_med - raw_at_peak, 6),
        "why_it_mattered": ("load_peak_tiles stores the raw `variant` column with no "
                            "reverse-complement, so position p meant opposite things on half the "
                            "substrate. It biased this probe toward a FALSE negative."),
    }
    print("   strand defect: raw %.4f vs normalized %.4f at alpha=1 (gain %+.4f); "
          "%+.4f at the peak" % (raw_at_1, curve[ALPHA_PRIMARY],
                                 curve[ALPHA_PRIMARY] - raw_at_1, strand_effect["gain_at_peak"]),
          flush=True)
    sweep_per = {}

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
        "onehot": {"per_seed": oh_per, "median": round(oh_med, 4),
                   "alpha_at_peak": peak["best_alpha"],
                   "encoding": "all-4-base, strand-normalized" if NORMALIZE_STRAND else "all-4-base, raw",
                   "note": ("the PEAK over the alpha grid, not the arms' default alpha=1.0. At "
                            "alpha=1.0 this scores %.4f, but that is an un-tuned point on a "
                            "collinear design, not the model's reach." % curve[ALPHA_PRIMARY])},
        "alpha_peak_check": peak,
        "strand_normalization": strand_effect,
        "reference_coding_measured_and_rejected": {
            "metric": 0.2826, "alpha": 1e8,
            "why": ("dropping one base per position (450 cols) removes the rank deficiency and "
                    "scores WORSE than keeping all four (0.2826 vs %.4f). With StandardScaler ahead "
                    "of the ridge the redundant column is what expresses aggregate composition under "
                    "heavy shrinkage. Measured, not assumed -- and it is why an earlier review's "
                    "'below GC under both codings' reading was wrong." % oh_med)},
        "n_test_per_seed": n_test,
        "committed_arms_same_substrate": arms,
        "excluded_off_substrate": off_substrate,
        "onehot_exceeds_these_arms": [{"arm": n, "metric": m} for n, m in exceeded],
        "these_arms_exceed_onehot": [{"arm": n, "metric": m} for n, m in beaten_by],
        "VERDICT": verdict,
        "consequence": PREREGISTERED["outcomes"][verdict]["consequence"],
        "independent_corroboration_NOT_registerable_as_arms": {
            "finding": ("wiki/promoteratlas_tile_arm_2026-11-01.json already records dinucleotide "
                        "0.2895 and k-mer4 0.2811 and states verbatim that 'both k-mer baselines "
                        "lost held-out ground against plain GC'. Three composition encodings "
                        "therefore sit at or below GC on this substrate, which corroborates this "
                        "probe's verdict by a route that predates it."),
            "why_they_are_NOT_registered_as_arms": [
                "That block declares `not_a_new_arm`: 'This is a LIMIT on an existing measured arm. "
                "The ceiling verdict is untouched.' Registering its internal comparators as arms "
                "would contradict the artifact's own governance.",
                "It ran on 43,698 tiles (408 excluded), not the arms' 44,106 -- so scoring them "
                "against the arms' criterion repeats the cross-substrate error this same probe "
                "already made once with the fragment arms.",
                "The block also carries `may_NEVER_be_written_as`, an explicit prohibition on how "
                "its result may be phrased.",
            ],
            "so_the_ledger_zero_count_means": ("not registered as arms on the activity criterion -- "
                                               "NOT never measured anywhere. The earlier framing "
                                               "'the panel was never satisfied, 0 of 16' came from "
                                               "string-matching one file and under-credited the "
                                               "project."),
        },
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
