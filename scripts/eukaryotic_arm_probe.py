"""Feasibility probe for the eukaryotic arm: does the substrate assemble, and is the target scorable?

Runs BEFORE the arm so the arm's own run cannot be the thing that discovers its substrate is
unusable. Three questions, in the order the frozen bar checks them (validity first):

  1. join rate -- does sequence join to activity above the loader's MIN_JOIN_RATE?
  2. folds     -- does the published fold assignment resolve onto the loaded rows?
  3. degeneracy -- is the activity readout degenerate under the project's OWN shipped bar?

(3) imports `assay_degeneracy` from the forward/inverse work rather than restating its thresholds.
That screen exists because CcdB posted this project's best-looking number purely because 79.3% of
its values sat at an assay ceiling, and a correlation read off a degenerate target is a figure, not
a measurement. It runs before any correlation is computed.

Read-only. Writes one JSON artifact and no model is fit here.
"""

from __future__ import annotations

import json
import time
from collections import Counter
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "wiki/eukaryotic_arm_feasibility_2026-11-01.json"


def gc_fraction(seq: str) -> float:
    s = seq.upper()
    return (s.count("G") + s.count("C")) / len(s) if s else 0.0


def probe_cell(cell: str) -> dict:
    from dna_decode.glm.lentimpra import load_cell

    t0 = time.time()
    records, stats = load_cell(cell)
    secs = time.time() - t0

    acts = [r.activity for r in records]
    folds = Counter(r.fold for r in records)
    lengths = Counter(len(r.seq) for r in records)

    # The shipped degeneracy bar, imported rather than restated.
    from scripts.forward_inverse_roundtrip import assay_degeneracy
    deg = assay_degeneracy(acts)

    return {
        "cell": cell,
        "n_records": len(records),
        "load_seconds": round(secs, 1),
        "stats": {k: v for k, v in vars(stats).items()} if hasattr(stats, "__dict__") else str(stats),
        "activity": {
            "min": round(min(acts), 4),
            "max": round(max(acts), 4),
            "distinct": len(set(acts)),
            "mode_share": round(max(Counter(acts).values()) / len(acts), 4),
        },
        "degeneracy_screen": deg if isinstance(deg, dict) else {"result": str(deg)},
        "folds": {
            "n_distinct": len(folds),
            "sizes": {str(k): v for k, v in sorted(folds.items(), key=lambda kv: str(kv[0]))},
            "unassigned": folds.get(None, 0) + folds.get("", 0),
        },
        "length": {
            "n_distinct": len(lengths),
            "uniform": len(lengths) == 1,
            "values": {str(k): v for k, v in sorted(lengths.items())[:6]},
        },
        "gc": {
            "min": round(min(gc_fraction(r.seq) for r in records), 4),
            "max": round(max(gc_fraction(r.seq) for r in records), 4),
        },
    }


def main() -> int:
    from dna_decode.glm.lentimpra import published_incumbents, usable_cells

    rec = {
        "schema": "eukaryotic-arm-feasibility-v1",
        "date": "2026-11-01",
        "bar": "wiki/eukaryotic_arm_acceptance_bar.json",
        "purpose": ("validity-first feasibility, so the arm's own run is not what discovers its "
                    "substrate is unusable"),
    }

    t0 = time.time()
    rec["usable_cells"] = usable_cells()
    inc = published_incumbents()
    rec["published_incumbents"] = {
        "metric_unstated": True,
        "refused_as_verdict_input": True,
        "why": ("Table S6 names no metric. These are context for what the panel supports, not a "
                "comparator; the frozen bar refuses the head-to-head."),
        "value": inc,
    }
    rec["cells"] = [probe_cell(c) for c in ("HepG2", "K562")]
    rec["wall_seconds"] = round(time.time() - t0, 1)

    OUT.write_text(json.dumps(rec, indent=2, sort_keys=True, default=str), encoding="utf-8")
    for c in rec["cells"]:
        a, f = c["activity"], c["folds"]
        print("%-6s n=%-7d folds=%-3d mode_share=%.4f distinct=%-6d uniform_len=%s"
              % (c["cell"], c["n_records"], f["n_distinct"], a["mode_share"], a["distinct"],
                 c["length"]["uniform"]))
        print("       degeneracy: %s" % json.dumps(c["degeneracy_screen"])[:160])
    print("-> %s  (%.1fs)" % (OUT, rec["wall_seconds"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
