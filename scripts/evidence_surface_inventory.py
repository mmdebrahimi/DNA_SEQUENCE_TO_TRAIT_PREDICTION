"""Which disclosure layers exist, and can a user of the TOOL actually see each one?

Read-only; exit 0 always (a report, not a gate).

WHY. The evaluation machinery is this project's differentiator, and most of it was legible only by
reading `wiki/`. Four per-cell disclosure layers rendered on the standing report card while only two
reached a decoder call, and one of the missing pair — source concentration — is the caveat that
explains the project's own worst-known metric gap (`escherichia_coli_shigella x gentamicin`: sens
0.893 from a cohort 95% one BioProject with no `rmt` carriers; source-diverse measurements of the same
cell report 0.523).

Reachability is checked at THREE levels, because they fail independently:
  card      the layer renders on wiki/decoder_validation_report_card.json
  record    it reaches `trust_block` -> the JSON a decoder call emits
  human     it reaches the printed CLI output (a JSON-only disclosure is not one for most users)

Everything is DERIVED from the artifact — a hand-listed layer set beside the data that defines it is
the drift bug this repo has hit four times.

Run: uv run python scripts/evidence_surface_inventory.py
"""
from __future__ import annotations

import json
import sys
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
WIKI = ROOT / "wiki"

# Layers with a dedicated human-readable renderer in the CLI output path. Kept explicit BECAUSE it is
# the level most easily forgotten: `doubt` shipped in the record a day before it printed anywhere.
_HUMAN_RENDERED = {
    "doubt_layer": "doubt_one_line (target-site paths) + render_note (bacterial uncounted-determinant path)",
    "source_concentration": "concentration_one_line, printed only for a SINGLE-SOURCE cell",
    "prospective": "folded into the validation caveat via trust_block's prospective-regression note",
    "lineage": "lineage_one_line — effective lineage N + cluster-weighted metrics, CI mandatory",
}


def main() -> int:
    from dna_decode.data.trust_surface import DISCLOSURE_LAYERS, trust_block

    card_p = WIKI / "decoder_validation_report_card.json"
    if not card_p.exists():
        print("report card artifact absent — nothing to inventory")
        return 0
    cells = json.loads(card_p.read_text(encoding="utf-8"))["cells"]

    observed = sorted({k for c in cells for k, v in c.items() if isinstance(v, dict)})
    rows = []
    for layer in observed:
        on_card = [c for c in cells if isinstance(c.get(layer), dict)]
        reaching = [c for c in cells if layer in trust_block(c["drug"], c["organism"])]
        rows.append({
            "layer": layer,
            "cells_on_card": len(on_card),
            "cells_reaching_record": len(reaching),
            "record_reachable": bool(reaching),
            "human_readable": _HUMAN_RENDERED.get(layer),
            "human_reachable": layer in _HUMAN_RENDERED,
        })

    undeclared = sorted(set(observed) - set(DISCLOSURE_LAYERS))
    out = {
        "schema": "evidence-surface-layer-inventory-v1", "generated": date.today().isoformat(),
        "contract": ("Reachability, not correctness. A layer that renders on the card but reaches no "
                     "decoder call is not a disclosure to anyone using the tool."),
        "n_cells": len(cells), "layers": rows,
        "layers_undeclared_in_DISCLOSURE_LAYERS": undeclared,
        "all_record_reachable": all(r["record_reachable"] for r in rows),
        "layers_without_a_human_renderer": [r["layer"] for r in rows if not r["human_reachable"]],
    }

    print(f"\n{len(cells)} cells | {len(observed)} disclosure layers\n")
    print(f"  {'layer':24} {'on card':>8} {'in record':>10}  human-readable")
    for r in rows:
        human = "yes" if r["human_reachable"] else "NO — JSON only"
        print(f"  {r['layer']:24} {r['cells_on_card']:>8} {r['cells_reaching_record']:>10}  {human}")
    if undeclared:
        print(f"\n  UNDECLARED in DISCLOSURE_LAYERS: {undeclared}")
    missing = out["layers_without_a_human_renderer"]
    if missing:
        print(f"\n  no human-readable renderer: {missing}")
        print("  (visible in --json-only output and in the record file, but not in printed output)")

    dest = WIKI / "evidence_surface_layer_inventory.json"
    dest.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"\nwrote wiki/{dest.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
