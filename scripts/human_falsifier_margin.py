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


if __name__ == "__main__":
    raise SystemExit(main())
