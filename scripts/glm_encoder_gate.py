"""Does a learned conv encoder beat a SINGLE GC NUMBER on real genomic promoter tiles?

    uv run python scripts/glm_encoder_gate.py [--seeds 10] [--no-download] [--with-grid]

**The gate that Step 3 predicted would fail, run on the user's explicit decision (2026-10-08).**
`wiki/glm_tile_headroom_2026-10-08.json` measured a position gain of **-0.0068** against a +0.02 bar and
recommended ratification before this build. The user ratified building it. The encoder is therefore the
test of its own gate: Step 3's proxy was positional GC in 10 bins under a LINEAR ridge, and both of those
are limits an encoder can escape (a motif at a precise offset; a composition x position interaction). If
the encoder also fails, the negative stops being a prediction from a weak proxy and becomes a measurement
from the model class the proxy stood in for.

**FOUR ARMS on tiles (primary):**

| arm | what it is | role |
|---|---|---|
| `gc` | ridge on ONE GC number | the incumbent, and the reconcile anchor |
| `encoder-w20` | conv, single width 20 | is one resolution enough? |
| `encoder-multi` | conv, widths 6/20/75 | the full candidate |
| `encoder-shuffled` | `encoder-multi`, TRAIN labels permuted | **validity** -- must show no signal |

**RECONCILE GATE FIRST.** The `gc` arm on seed 0 must reproduce the committed **0.3000 +/- 0.01**
(`wiki/glm_genomewide_oracle_2026-10-07.json`). If it does not, this script is not reading the same data
or split as the published baseline and NOTHING below is interpretable. The committed single-seed artifact
is left untouched and POINTED AT; the GC baseline is re-derived over the registered seed list so the
comparison is seed-matched rather than comparing a 10-seed median against a 1-seed number.

**TWO FRAMINGS REPORTED, because regression may be the wrong one.** The tile set is 14,913 active /
29,193 inactive, so a Spearman over everything can be substantially ranking noise among inactive tiles.
So the gate also reports **AUROC on the `active` column** per arm. If the encoder loses on Spearman and
wins on AUROC, the finding is about the framing, not the architecture.

**The grid arm (`--with-grid`) is a POSITIVE CONTROL, not evidence of transferable grammar.** On the
designed grid, positional one-hot scores 0.5914 against GC 0.1414 -- position massively matters there --
so an encoder that works on the grid and fails on tiles locates the failure in the SUBSTRATE rather than
the architecture. It must NOT be read as the encoder having learned promoter grammar: the grid is
separable from genomic sequence at AUROC 1.0000, and with ~3x8x8x8x8 element variants a model scoring
well under leave-element-out can be exploiting the four axes it did see.
"""
from __future__ import annotations

import argparse
import json
import statistics
from datetime import date
from pathlib import Path

from dna_decode.glm.encoder import (
    DEFAULT_WIDTHS,
    REGISTERED_PROTOCOL,
    fit_encoder,
)
from dna_decode.glm.genomewide import (
    DEFAULT_CACHE,
    gc_features,
    leave_peak_out,
    load_peak_tiles,
    noise_ceiling,
    spearman,
)

#: Frozen BEFORE the run. Every threshold here is ASSERTED; the ceiling that contextualises them is DERIVED.
PREREGISTERED = {
    "question": "does a learned convolutional encoder beat a single GC number on real genomic promoter "
                "tiles, under leave-peak-out, over a registered seed list?",
    "primary_metric": "spearman",
    "secondary_metric": "auroc_active",
    "comparator": "gc",
    "candidate_arms": ["encoder-w20", "encoder-multi"],
    "null_arm": "encoder-shuffled",
    "generalises_margin": 0.05,
    "null_max_abs": 0.05,
    "registered_seeds": list(range(10)),
    "reconcile_anchor": {"arm": "gc", "seed": 0, "expected": 0.3000, "tolerance": 0.01,
                         "source": "wiki/glm_genomewide_oracle_2026-10-07.json"},
    "verdicts": {
        "GENERALISES": "best encoder arm beats gc by >= 0.05 at the median -- the learned representation "
                       "earns its keep and G-A is a positive",
        "WEAK": "best encoder arm beats gc by > 0 but < 0.05 -- a real but small edge; not enough to "
                "carry a learned representation on its own",
        "DOES_NOT_GENERALISE": "best encoder arm does NOT beat gc at the median -- GA-H3 confirmed and "
                               "Step 3's predicted-negative becomes a measured one",
        "INDETERMINATE_NULL_NOT_CLEAN": "the shuffled-label arm itself scores, so the pipeline fabricates "
                                        "signal and NOTHING is graded",
        "INDETERMINATE_RECONCILE_FAILED": "the gc arm did not reproduce the committed 0.3000 on seed 0, "
                                          "so this script is not reading the same data as the published "
                                          "number and nothing is interpretable",
        "INDETERMINATE_PROTOCOL_DRIFT": "an arm ran under a protocol other than the registered one",
    },
    "thresholds_are": "ASSERTED; the 0.941 tile ceiling that contextualises them is DERIVED from the "
                      "tiles' own replicates",
}


def ridge_gc(train, test) -> tuple[float, list[float]]:
    """The incumbent arm, the same shape as the oracle's and Step 3's `_fit`.

    Re-implemented rather than imported because the published number came from a script, not a module --
    and the reconcile gate below is what proves the re-implementation is faithful instead of assuming it.

    Returns `(spearman, predictions)`; the predictions are what let the SECOND framing (AUROC on the
    active column) cover the incumbent too. Scoring only the encoder arms on AUROC would make the second
    framing unable to answer the question it exists for -- whether the encoder wins under a different
    frame -- because there would be nothing to compare it against.
    """
    import numpy as np
    from sklearn.linear_model import Ridge
    from sklearn.preprocessing import StandardScaler
    Xtr = np.asarray(gc_features([t.seq for t in train]), dtype=np.float64)
    Xte = np.asarray(gc_features([t.seq for t in test]), dtype=np.float64)
    ytr = np.log1p(np.asarray([t.expression for t in train], dtype=np.float64))
    sc = StandardScaler().fit(Xtr)
    m = Ridge(alpha=1.0).fit(sc.transform(Xtr), ytr)
    pred = [float(v) for v in m.predict(sc.transform(Xte))]
    return spearman(pred, [t.expression for t in test]), pred


def auroc(scores: list[float], labels: list[int]) -> float | None:
    """Rank-based AUROC with tie handling. `None` when one class is absent (never 0.5, which would read
    as a measured coin-flip rather than as an unmeasurable cell)."""
    pos = sum(labels)
    neg = len(labels) - pos
    if pos == 0 or neg == 0:
        return None
    order = sorted(range(len(scores)), key=lambda i: scores[i])
    ranks = [0.0] * len(scores)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and scores[order[j + 1]] == scores[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[order[k]] = avg
        i = j + 1
    rsum = sum(r for r, l in zip(ranks, labels) if l)
    return (rsum - pos * (pos + 1) / 2.0) / (pos * neg)


def _active_labels(test) -> list[int]:
    return [1 if getattr(t, "active", "") == "active" else 0 for t in test]


def verdict(arms: dict, *, seeds_used: list[int]) -> tuple[str, str, dict]:
    """Frozen six-branch rule. VALIDITY AND PROVENANCE ARE CHECKED BEFORE THE MARGIN (decision D1).

    Pure: takes the aggregated per-arm medians and returns the verdict, so every branch is unit-testable
    without touching torch or the data.
    """
    anc = PREREGISTERED["reconcile_anchor"]
    flags = {"registered_protocol_used": True, "null_clean": True, "reconcile_ok": True,
             "seeds_match_registered": seeds_used == PREREGISTERED["registered_seeds"]}

    if not flags["seeds_match_registered"]:
        flags["registered_protocol_used"] = False
        return ("INDETERMINATE_PROTOCOL_DRIFT",
                f"seeds used {seeds_used} differ from the registered list "
                f"{PREREGISTERED['registered_seeds']}; a median over a different seed set is not the "
                f"registered gate.", flags)

    for name, a in arms.items():
        if a.get("protocol") is not None and a["protocol"] != REGISTERED_PROTOCOL.as_dict():
            flags["registered_protocol_used"] = False
            return ("INDETERMINATE_PROTOCOL_DRIFT",
                    f"arm {name} ran under a protocol other than the registered one; the artifact must "
                    f"not claim the registered protocol was used.", flags)

    gc_seed0 = arms[anc["arm"]]["per_seed"][anc["seed"]]
    if abs(gc_seed0 - anc["expected"]) > anc["tolerance"]:
        flags["reconcile_ok"] = False
        return ("INDETERMINATE_RECONCILE_FAILED",
                f"gc arm on seed {anc['seed']} gave {gc_seed0:.4f}, expected {anc['expected']} "
                f"+/-{anc['tolerance']} from {anc['source']}. Not reading the same split or data as the "
                f"published number; nothing below is interpretable.", flags)

    base = arms["gc"]["median"]
    null = arms[PREREGISTERED["null_arm"]]["median"]
    best_name = max(PREREGISTERED["candidate_arms"], key=lambda a: arms[a]["median"])
    best = arms[best_name]["median"]
    margin = best - base
    d = (f"gc {base:+.4f} (re-derived over the registered seeds; committed seed0 = "
         f"{anc['expected']}) | best encoder = {best_name} {best:+.4f} | MARGIN {margin:+.4f} against "
         f"{PREREGISTERED['generalises_margin']} | null {null:+.4f}; reconcile OK (gc seed0 "
         f"{gc_seed0:.4f})")

    if abs(null) >= PREREGISTERED["null_max_abs"]:
        flags["null_clean"] = False
        return ("INDETERMINATE_NULL_NOT_CLEAN",
                f"the shuffled-label arm scored {null:+.4f}, at or above the "
                f"{PREREGISTERED['null_max_abs']} bar. The pipeline fabricates signal; nothing is "
                f"graded. {d}", flags)
    if margin >= PREREGISTERED["generalises_margin"]:
        return "GENERALISES", d, flags
    if margin > 0:
        return "WEAK", d, flags
    return "DOES_NOT_GENERALISE", d, flags


def _agg(runs: list[dict], key: str = "spearman") -> dict:
    vals = [r[key] for r in runs]
    out = {
        "median": round(statistics.median(vals), 4),
        "stdev": round(statistics.stdev(vals), 4) if len(vals) > 1 else 0.0,
        "min": round(min(vals), 4), "max": round(max(vals), 4),
        "per_seed": [round(v, 4) for v in vals],
    }
    au = [r["auroc_active"] for r in runs if r.get("auroc_active") is not None]
    out["auroc_active_median"] = round(statistics.median(au), 4) if au else None
    out["protocol"] = runs[0].get("protocol")
    if runs[0].get("n_parameters") is not None:
        out["n_parameters"] = runs[0]["n_parameters"]
        out["epochs_ran_per_seed"] = [r.get("epochs_ran") for r in runs]
        out["device"] = runs[0].get("device")
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cache-dir", default=str(DEFAULT_CACHE))
    ap.add_argument("--seeds", type=int, default=10)
    ap.add_argument("--device", default="cpu", choices=["cpu", "cuda", "auto"],
                    help="cpu is the DETERMINISM REFERENCE; cuda convolutions may be nondeterministic")
    ap.add_argument("--with-grid", action="store_true",
                    help="also run the designed-grid POSITIVE CONTROL (not evidence of grammar)")
    ap.add_argument("--no-download", action="store_true")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    tiles_all, st = load_peak_tiles(args.cache_dir, allow_download=not args.no_download)
    # the SAME two filters the published baseline used, or the reconcile gate cannot pass
    tiles = [t for t in tiles_all if t.category == "tile" and len(t.seq) == 150]
    ceil = noise_ceiling(tiles)["ceiling_for_a_sequence_model"]
    n_active = sum(_active_labels(tiles))
    print(f"tiles {len(tiles)} of {st.n_parsed} | ceiling {ceil} | "
          f"{n_active} active / {len(tiles) - n_active} inactive | peaks {len({t.peak for t in tiles})}")

    seeds = list(range(args.seeds))
    runs: dict[str, list[dict]] = {}
    for s in seeds:
        tr, te, held = leave_peak_out(tiles, seed=s)
        lab = _active_labels(te)

        g, g_pred = ridge_gc(tr, te)
        runs.setdefault("gc", []).append(
            {"spearman": g, "auroc_active": auroc(g_pred, lab), "protocol": None})
        print(f"  seed {s}: gc {g:+.4f} auroc {auroc(g_pred, lab):.4f} "
              f"(held {held} peaks, {len(te)} test)")

        for name, kw in (("encoder-w20", {"widths": (20,)}),
                         ("encoder-multi", {"widths": DEFAULT_WIDTHS}),
                         ("encoder-shuffled", {"widths": DEFAULT_WIDTHS, "shuffle_labels": True})):
            r = fit_encoder(tr, te, seed=s, device=args.device, **kw)
            # compute the second framing from the SAME fit, then DROP the raw predictions -- one float
            # per test row x 30 fits has no business in a committed artifact
            r["auroc_active"] = auroc(r.pop("predictions"), lab)
            runs.setdefault(name, []).append(r)
            print(f"           {name:18s} {r['spearman']:+.4f} auroc {r['auroc_active']:.4f} "
                  f"({r['n_parameters']} params, {r['epochs_ran']} epochs, {r['device']})")

    arms = {k: _agg(v) for k, v in runs.items()}
    v, detail, flags = verdict(arms, seeds_used=seeds)

    artifact = {
        "schema": "glm-encoder-gate-v1",
        "date": str(date.today()),
        "substrate": "GEO GSE144621 peak tiles, category=tile, len=150, leave-peak-out",
        "n_tiles": len(tiles),
        "n_active": n_active,
        "noise_ceiling": ceil,
        "seeds": seeds,
        "device": args.device,
        "preregistered": PREREGISTERED,
        "registered_protocol": REGISTERED_PROTOCOL.as_dict(),
        "arms": arms,
        "verdict": v,
        "verdict_detail": detail,
        "registered_protocol_used": flags["registered_protocol_used"],
        "null_clean": flags["null_clean"],
        "reconcile_ok": flags["reconcile_ok"],
        "built_against_a_predicted_negative": {
            "gate": "wiki/glm_tile_headroom_2026-10-08.json",
            "gate_verdict": "NO_POSITIONAL_HEADROOM",
            "gate_position_gain": -0.0068,
            "user_decision": "build the conv encoder anyway (2026-10-08)",
            "why_that_is_informative": "the gate's proxy was positional GC in 10 bins under a LINEAR "
                                       "ridge; an encoder escapes both limits, so its result upgrades a "
                                       "predicted-negative to a measured one either way",
        },
        "honest_limits": [
            "A negative here is about THIS encoder at THIS size (32 channels/branch, ~13k parameters) on "
            "THIS substrate. It is strong evidence against the conv-encoder lever on genomic promoter "
            "tiles and is NOT a statement about learned representations in general.",
            "The encoder is deliberately small because capacity has already HURT on this substrate "
            "(6-mers, 4,096 features, scored worst of every arm at 0.2022). A larger encoder is a "
            "different experiment, and the shuffled-label null is what would detect it over-fitting.",
            "Spearman over all tiles may be substantially ranking noise among the 29,193 inactive tiles; "
            "AUROC on the active column is reported as the second framing for exactly that reason.",
            "CPU is the determinism reference. A --device cuda run is faster but CUDA convolutions are "
            "documented as potentially nondeterministic, and the device used is recorded per arm.",
            "The margin thresholds (0.05 generalises, 0.05 null bar) are ASSERTED, not derived.",
        ],
    }
    out = Path(args.out) if args.out else Path(f"wiki/glm_encoder_gate_{date.today()}.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(artifact, indent=2), encoding="utf-8")

    print(f"\nVERDICT: {v}\n  {detail}")
    print(f"registered_protocol_used={flags['registered_protocol_used']} "
          f"null_clean={flags['null_clean']} reconcile_ok={flags['reconcile_ok']}")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
