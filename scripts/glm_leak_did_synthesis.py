"""Two re-pretrains of ONE recipe put the leak DiD on OPPOSITE sides of zero. Read them together.

WHY THIS EXISTS. `scripts/glm_leak_did.py` computes the DiD for ONE assessed checkpoint against the
committed one. Run for both preserved v5 arms it produces four cells (2 checkpoints x 2 scoring
arms), and the single most informative fact about the quantity is visible only ACROSS them: the two
arms are re-pretrains of the SAME recipe differing only in seed, so the spread between their DiDs is
a direct, assumption-free read on how much of a single-seed DiD is noise.

At 500 windows per half the four cells disagree on SIGN -- one arm negative, the other positive --
while every interval contains zero and every interval EXCLUDES the predicted -0.0054. A
between-seed spread comparable to the point estimates is what "noise-dominated" looks like, and it
corroborates, by a completely different route, the 0.00281 noise floor the graded artifact measured
from a committed-corpus reseed.

WHAT THIS IS NOT. Sign disagreement between two seeds does NOT prove the true DiD is zero; it shows
a single-seed DiD of this magnitude is uninterpretable. And excluding the PREDICTED size does not
exclude a leak -- it bounds one. Both readings are stated in the output rather than left to a
reader, because "the leak is excluded" is the tempting over-read and it is wrong.

Pure JSON synthesis over already-written artifacts: no scoring, no torch, nothing re-derived.
"""
from __future__ import annotations

import argparse
import json
import statistics
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

PREDICTED = -0.0054
NOISE_FLOOR = 0.00281
OUT = ROOT / "wiki" / "glm_leak_did_n500_synthesis_2026-11-28.json"


class SynthError(RuntimeError):
    pass


def synthesise(cells: list[dict]) -> dict:
    """Read the four cells together. Every claim below is computed, not asserted."""
    if len(cells) < 2:
        raise SynthError("a cross-seed reading needs at least two cells")
    dids = [c["result"]["did_nats"] for c in cells]
    verdicts = sorted({c["VERDICT"]["VERDICT"] for c in cells})
    signs = {("negative" if d < 0 else "positive" if d > 0 else "zero") for d in dids}

    # The cross-SEED spread, held at a fixed scoring arm, is the cleanest noise read: same recipe,
    # same windows, same arm, different seed only.
    by_arm: dict[str, list[float]] = {}
    for c in cells:
        by_arm.setdefault(c["arm"], []).append(c["result"]["did_nats"])
    seed_spread = {arm: round(max(v) - min(v), 6) for arm, v in by_arm.items() if len(v) > 1}

    all_exclude_predicted = all(c["result"]["excludes_predicted"] for c in cells)
    any_excludes_zero = any(c["result"]["excludes_zero"] for c in cells)
    widest_lower = min(c["result"]["ci95"][0] for c in cells)
    widest_upper = max(c["result"]["ci95"][1] for c in cells)

    return {
        "n_cells": len(cells),
        "did_nats": [round(d, 6) for d in dids],
        "mean_did_nats": round(statistics.fmean(dids), 6),
        "signs_observed": sorted(signs),
        "signs_disagree": len(signs) > 1,
        "cross_seed_spread_nats_by_arm": seed_spread,
        "noise_floor_nats": NOISE_FLOOR,
        "spread_is_comparable_to_the_estimates": bool(
            seed_spread and max(seed_spread.values()) >= min(abs(d) for d in dids)),
        "verdicts": verdicts,
        "all_cells_exclude_predicted_size": all_exclude_predicted,
        "any_cell_excludes_zero": any_excludes_zero,
        "envelope_over_all_cells": [round(widest_lower, 6), round(widest_upper, 6)],
        "reading": (
            "the predicted leak of {p} nats is excluded in every cell, and no cell separates the "
            "DiD from zero; the two seeds of one recipe disagree on sign, so a single-seed DiD of "
            "this magnitude is not interpretable".format(p=PREDICTED)
            if all_exclude_predicted and not any_excludes_zero and len(signs) > 1 else
            "mixed; read the per-cell artifacts rather than this summary"),
        "what_this_does_NOT_show": [
            "NOT that there is no leak. The predicted SIZE is excluded; the envelope still admits "
            f"a leak anywhere inside {[round(widest_lower, 6), round(widest_upper, 6)]} nats.",
            "NOT that the true DiD is zero. Sign disagreement between seeds shows a single-seed "
            "estimate of this magnitude is uninterpretable, which is a statement about the "
            "measurement, not about the parameter.",
            "NOT a generalisation beyond this corpus and this checkpoint family.",
        ],
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cells", nargs="+", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=OUT)
    a = ap.parse_args()

    cells = []
    for p in a.cells:
        if not p.exists():
            print(f"REFUSING: missing {p}", file=sys.stderr)
            return 2
        d = json.loads(p.read_text(encoding="utf-8"))
        if d.get("record") != "glm-leak-did-v1":
            print(f"REFUSING: {p.name} is not a glm-leak-did-v1 artifact", file=sys.stderr)
            return 2
        cells.append(d)

    # Every cell must describe the same n, or the spread mixes power levels with seeds.
    ns = {(c["result"]["n_bearing"], c["result"]["n_clean"]) for c in cells}
    if len(ns) != 1:
        print(f"REFUSING: cells disagree on n: {ns}; a cross-seed spread computed across "
              "different window counts confounds noise with power", file=sys.stderr)
        return 2

    try:
        syn = synthesise(cells)
    except SynthError as exc:
        print(f"REFUSING: {exc}", file=sys.stderr)
        return 2

    out = {
        "record": "glm-leak-did-synthesis-v1",
        "date": "2026-11-28",
        "question": "how much of a single-seed leak DiD is noise, and is the predicted leak real?",
        "predicted_did_nats": PREDICTED,
        "prior_at_120_windows": {"did_nats": -0.001897, "ci95": [-0.006801, 0.002813],
                                 "verdict": "contained BOTH zero and the predicted size"},
        "inputs": [str(p) for p in a.cells],
        "synthesis": syn,
        "honest_limits": [
            "Two seeds is the smallest set that can show disagreement; it bounds nothing tightly.",
            "The DiD identifies a leak only under the design's assumption that cross-organism cost "
            "acts equally on both halves and cancels. Not testable from these numbers.",
            "The clean half uses 500 of 519 available windows, so the bootstrap treats a near-"
            "exhaustive sample as an infinite population. That makes the intervals CONSERVATIVE: a "
            "finite-population correction would narrow them, never widen them.",
        ],
    }
    a.out.write_text(json.dumps(out, indent=1) + "\n", encoding="utf-8")

    print(f"cells {syn['n_cells']} | DiDs {syn['did_nats']}")
    print(f"  signs {syn['signs_observed']} disagree={syn['signs_disagree']}")
    print(f"  cross-seed spread by arm {syn['cross_seed_spread_nats_by_arm']} "
          f"(independently measured noise floor {NOISE_FLOOR})")
    print(f"  all exclude predicted {PREDICTED}: {syn['all_cells_exclude_predicted_size']} | "
          f"any excludes zero: {syn['any_cell_excludes_zero']}")
    print(f"  envelope {syn['envelope_over_all_cells']}")
    print(f"wrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
