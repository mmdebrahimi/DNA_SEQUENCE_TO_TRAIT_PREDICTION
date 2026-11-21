"""Measure the human within-substrate margin: how far a 3-mer Markov null sits above its own floor.

Graded by `wiki/human_falsifier_bar_2026-11-21.json`, frozen before any number here existed.

THE QUESTION. The generator's one robust failure is losing to a 3-mer Markov chain (worse than its
own seed-matched null at 10 of 10 seeds, mean gap +0.1186), and every measurement of that sits on
ONE substrate: 4,340 E. coli upstream-of-CDS PROXY windows at 150 bp, which the project has already
measured as the weaker one (SUBSTRATE_LIMITED, 2026-11-01). Is the defeat a property of that
SUBSTRATE or of GENERATORS?

WHAT IS COMPARED, AND WHAT IS NOT. Per the bar, the decision quantity is each substrate's
WITHIN-substrate margin, `null_distinguishability - self_control_ceiling`. The raw human level is
NOT compared to the bacterial 0.6679: the panels differ in length (230 vs 150 bp), assay, kingdom
and library construction, and the predecessor memo states that only a within-substrate contrast
crosses a substrate boundary.

NO GENERATOR RUNS HERE, by design. The current generator is GENERator-v2-prokaryote-1.2b-base, so a
human generation arm would need a new base model and would confound substrate with model family.
This characterises the OPPONENT only: CPU, no GPU, no training, no network.

DISCRIMINATOR CAPACITY IS HELD FIXED at the bacterial run's settings (k=4, bins=6, n=120 generated,
1200 natural, seed 7). The margin involves the self-control, so changing the discriminator's
capacity would move both terms and make the two substrates' margins incomparable.

Run: uv run python scripts/human_falsifier_margin.py
"""

from __future__ import annotations

import json
import time
from pathlib import Path

from dna_decode.glm.baselines import MarkovGenerator, UniformGenerator
from dna_decode.glm.falsifier import MIN_PER_CLASS, heldout_auroc, self_control_ceiling
from dna_decode.glm.human_natural import ELEMENT_BP, FULL_BP, fit_eval_by_fold

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "mpra_human"
OUT = ROOT / "wiki" / "human_falsifier_margin_2026-11-21.json"

#: Held fixed at the bacterial run's geometry so the two margins stay comparable.
K = 4
BINS = 6
N_GEN = 120
N_NATURAL = 1200
SEED = 7
MARKOV_ORDER = 3
CEILING_SEEDS = (0, 1, 2, 3, 4)
#: The publisher's fold held out for scoring; the rest fits the null.
EVAL_FOLDS = {1}

CELLS = ("HepG2", "K562")
#: (label, strip_adaptors, bp). The bar fixes full-230 as PRIMARY.
GEOMETRIES = (("full_230nt_as_shipped", False, FULL_BP),
              ("element_200nt_adaptors_stripped", True, ELEMENT_BP))
MODES = ("kmer", "both")

#: Committed bacterial anchors, from wiki/glm_generator_paired_sweep_2026-11-20.json. NOT recomputed.
BACTERIAL = {"both": {"null": 0.6679, "self_control": 0.5865, "margin": 0.0814},
             "kmer": {"null": 0.5990, "self_control": 0.5447, "margin": 0.0543}}


def one_cell(cell: str, geom_label: str, strip: bool, bp: int, mode: str) -> dict:
    fit, ev, split = fit_eval_by_fold(cell, DATA, eval_folds=EVAL_FOLDS,
                                      limit_eval=N_NATURAL, strip_adaptors=strip)
    rec = {"cell": cell, "geometry": geom_label, "bp": bp, "mode": mode,
           "n_fit": len(fit), "n_eval": len(ev),
           "fit_folds": split["fit_folds"], "eval_folds": split["eval_folds"],
           "content_overlap": split["content_overlap"]}
    if len(ev) < MIN_PER_CLASS:
        rec.update(status="INSUFFICIENT_NATURAL", min_per_class=MIN_PER_CLASS)
        return rec

    markov = MarkovGenerator(order=MARKOV_ORDER, seed=SEED).fit(fit)
    null_seqs = markov.generate(N_GEN, bp)
    uni_seqs = UniformGenerator(seed=SEED).generate(N_GEN, bp)

    # natural FIRST, candidate SECOND -- the convention the falsifier uses at lines 310/325/328.
    null_auc = heldout_auroc(ev, null_seqs, k=K, seed=SEED, mode=mode, bins=BINS)
    uni_auc = heldout_auroc(ev, uni_seqs, k=K, seed=SEED, mode=mode, bins=BINS)
    ceiling, ceil_vals = self_control_ceiling(ev, k=K, seeds=CEILING_SEEDS, mode=mode, bins=BINS)

    rec.update(null_distinguishability=None if null_auc is None else round(null_auc, 4),
               uniform_control=None if uni_auc is None else round(uni_auc, 4),
               self_control_ceiling=None if ceiling is None else round(ceiling, 4),
               self_control_values=[round(v, 4) for v in ceil_vals])

    if null_auc is None or ceiling is None:
        rec["status"] = "INSTRUMENT_INVALID"
        rec["reason"] = ("the null or the self-control could not be computed; no margin exists"
                         if ceiling is None else "the null produced no score")
        return rec

    rec["margin"] = round(null_auc - ceiling, 4)
    rec["self_control_seed_spread"] = round(max(ceil_vals) - min(ceil_vals), 4) if ceil_vals else None
    rec["bacterial_margin"] = BACTERIAL[mode]["margin"]
    rec["margin_minus_bacterial"] = round(rec["margin"] - BACTERIAL[mode]["margin"], 4)
    rec["status"] = "SCORED"
    return rec


def main() -> int:
    t0 = time.time()
    rows: list[dict] = []
    for cell in CELLS:
        for geom_label, strip, bp in GEOMETRIES:
            for mode in MODES:
                r = one_cell(cell, geom_label, strip, bp, mode)
                rows.append(r)
                # Written after EVERY cell: a truncated run then leaves a readable partial table
                # rather than an average over an unknown denominator.
                OUT.write_text(json.dumps(_wrap(rows, t0, complete=False), indent=2),
                               encoding="utf-8")
                print("  %-7s %-32s %-5s null=%-7s ceil=%-7s margin=%-8s uni=%s"
                      % (cell, geom_label, mode,
                         r.get("null_distinguishability"), r.get("self_control_ceiling"),
                         r.get("margin"), r.get("uniform_control")), flush=True)

    OUT.write_text(json.dumps(_wrap(rows, t0, complete=True), indent=2), encoding="utf-8")
    print("\nwrote %s in %.0f s" % (OUT.relative_to(ROOT).as_posix(), time.time() - t0))
    return 0


def _wrap(rows: list[dict], t0: float, *, complete: bool) -> dict:
    return {
        "record": "human-falsifier-margin-v1",
        "date": "2026-11-21",
        "bar": "wiki/human_falsifier_bar_2026-11-21.json",
        "substrate_artifact": "wiki/human_falsifier_substrate_2026-11-21.json",
        "complete": complete,
        "no_generator_was_run": ("By design. The opponent (3-mer Markov null) and the uniform "
                                 "control are characterised; no neural generator is involved, so "
                                 "substrate is not confounded with model family."),
        "discriminator": {"k": K, "bins": BINS, "n_generated": N_GEN,
                          "n_natural_eval": N_NATURAL, "seed": SEED,
                          "markov_order": MARKOV_ORDER, "ceiling_seeds": list(CEILING_SEEDS),
                          "held_fixed_because": ("the margin involves the self-control, so changing "
                                                 "discriminator capacity would move both terms and "
                                                 "make the two substrates' margins incomparable")},
        "split": {"eval_folds": sorted(EVAL_FOLDS),
                  "source": "the publisher's own 10-fold assignment, carried through load_cell"},
        "bacterial_anchors": BACTERIAL,
        "PROHIBITED_COMPARISON": ("The raw human null level may NOT be compared to the bacterial "
                                  "0.6679. Only the within-substrate MARGIN crosses the boundary."),
        "rows": rows,
        "wall_seconds": round(time.time() - t0, 1),
    }


# ---------------------------------------------------------------------------------------------
# THE FROZEN READING, applied to the committed grid. Rules are fixed in
# wiki/human_falsifier_bar_2026-11-21.json and are NOT recomputed or relaxed here.
#
# The tolerance is the HUMAN self-control's own seed spread, per the bar -- not a constant, because
# the magnitude is what the run measures and freezing one beforehand would be guessing at the scale.
# ---------------------------------------------------------------------------------------------

VERDICT_OUT = ROOT / "wiki" / "human_falsifier_verdict_2026-11-21.json"


def row_verdict(row: dict) -> tuple[str, dict]:
    """One cell's verdict: margin-vs-margin against the measured seed spread."""
    if row.get("status") != "SCORED":
        return "INSTRUMENT_INVALID", {"reason": row.get("status")}
    hm = row["margin"]
    bm = BACTERIAL[row["mode"]]["margin"]
    spread = row["self_control_seed_spread"]
    diff = round(hm - bm, 4)
    if diff > spread:
        v = "CHAIN_IS_WEAKER_ON_HUMAN"
    elif diff < -spread:
        v = "CHAIN_IS_STRONGER_ON_HUMAN"
    else:
        v = "CHAIN_IS_COMPARABLE"
    return v, {"human_margin": hm, "bacterial_margin": bm, "diff": diff,
               "tolerance_seed_spread": spread, "exceeds_tolerance": abs(diff) > spread}


def combine_cells(verdicts: list[str]) -> tuple[str, bool]:
    """The bar's CONSERVATIVE disagreement rule: if the cells do not agree, report COMPARABLE.

    WEAKER + COMPARABLE -> COMPARABLE. COMPARABLE + STRONGER -> COMPARABLE.
    WEAKER + STRONGER -> COMPARABLE, and the contradiction is recorded.
    Fixed in advance so a split result cannot be resolved by picking the agreeable cell.
    """
    uniq = set(verdicts)
    if "INSTRUMENT_INVALID" in uniq:
        return "INSTRUMENT_INVALID", False
    if len(uniq) == 1:
        return uniq.pop(), False
    contradiction = {"CHAIN_IS_WEAKER_ON_HUMAN", "CHAIN_IS_STRONGER_ON_HUMAN"} <= uniq
    return "CHAIN_IS_COMPARABLE", contradiction


def compute_verdict(rows: list[dict]) -> dict:
    per_row = []
    for r in rows:
        v, detail = row_verdict(r)
        per_row.append({"cell": r["cell"], "geometry": r["geometry"], "mode": r["mode"],
                        "verdict": v, **detail})

    by_geom_mode = {}
    for geom_label, _strip, _bp in GEOMETRIES:
        for mode in MODES:
            sel = [p for p in per_row if p["geometry"] == geom_label and p["mode"] == mode]
            v, contra = combine_cells([p["verdict"] for p in sel])
            by_geom_mode["%s|%s" % (geom_label, mode)] = {
                "verdict": v, "cells_contradict": contra,
                "per_cell": {p["cell"]: p["verdict"] for p in sel}}

    # The adaptor diagnostic: how much of the full-230 distinguishability survives stripping.
    adaptor = {}
    for cell in CELLS:
        for mode in MODES:
            full = next((r for r in rows if r["cell"] == cell and r["mode"] == mode
                         and r["geometry"] == GEOMETRIES[0][0]), None)
            elem = next((r for r in rows if r["cell"] == cell and r["mode"] == mode
                         and r["geometry"] == GEOMETRIES[1][0]), None)
            if full and elem and full.get("status") == "SCORED" == elem.get("status"):
                adaptor["%s|%s" % (cell, mode)] = {
                    "null_full_230": full["null_distinguishability"],
                    "null_element_200": elem["null_distinguishability"],
                    "null_drop": round(full["null_distinguishability"]
                                       - elem["null_distinguishability"], 4),
                    "ceiling_full_230": full["self_control_ceiling"],
                    "ceiling_element_200": elem["self_control_ceiling"],
                    "ceiling_move": round(elem["self_control_ceiling"]
                                          - full["self_control_ceiling"], 4)}
    return {"per_row": per_row, "by_geometry_and_mode": by_geom_mode,
            "adaptor_diagnostic": adaptor}


def verdict_main() -> int:
    grid = json.loads(OUT.read_text(encoding="utf-8"))
    v = compute_verdict(grid["rows"])
    rec = {"record": "human-falsifier-verdict-v1", "date": "2026-11-21",
           "bar": "wiki/human_falsifier_bar_2026-11-21.json",
           "grid": "wiki/human_falsifier_margin_2026-11-21.json",
           "primary_geometry": GEOMETRIES[0][0],
           "secondary_geometry": GEOMETRIES[1][0],
           "primary_was_fixed_in_advance": ("The bar named full-230 PRIMARY before any number "
                                            "existed, for comparability with the activity labels "
                                            "and the 2026-11-01 encoder arm. It is NOT swapped "
                                            "after the fact; see the memo for why that choice is "
                                            "nonetheless uninformative about this question."),
           **v}
    VERDICT_OUT.write_text(json.dumps(rec, indent=2), encoding="utf-8")

    print("PER (geometry, mode):")
    for key, g in v["by_geometry_and_mode"].items():
        print("  %-40s %-28s contradict=%-5s %s"
              % (key, g["verdict"], g["cells_contradict"], g["per_cell"]))
    print()
    print("ADAPTOR DIAGNOSTIC (null distinguishability lost by stripping 15 nt each end):")
    for key, a in v["adaptor_diagnostic"].items():
        print("  %-14s null %.4f -> %.4f  (drop %+.4f)  ceiling moves %+.4f"
              % (key, a["null_full_230"], a["null_element_200"], -a["null_drop"],
                 a["ceiling_move"]))
    print()
    print("wrote", VERDICT_OUT.relative_to(ROOT).as_posix())
    return 0


if __name__ == "__main__":
    import sys
    raise SystemExit(verdict_main() if "--verdict-only" in sys.argv else main())
