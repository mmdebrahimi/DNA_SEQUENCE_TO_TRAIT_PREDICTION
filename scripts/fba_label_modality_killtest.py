"""KILL-TEST: is the flatness a competitive-vs-monoculture LABEL-MODALITY mismatch (candidate C-G5-1)?

The claim relaxes an assumption: RB-TnSeq measures POOLED COMPETITIVE fitness while FBA computes
MONOCULTURE growth, so a gene could be competitively unfit yet monoculture-viable -- making the flatness
a label artifact rather than a model defect.

The disproof is already in hand and needs no new computation. The Orth 2011 4-media screen scores
**individual deletion strains** for growth/no-growth -- MONOCULTURE, the same modality FBA computes. If
the model also fails to reproduce the switch on THAT substrate, the mismatch cannot be what causes the
flatness.

POLARITY: exit 0 = DISPROVED = claim KILLED.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ART = Path(__file__).resolve().parent.parent / "wiki/fba_conditional_essentiality_ecoli_2026-08-12.json"
# The model reproduces <= this many of the 4-media switches -> monoculture fails too -> claim disproved.
MAX_EXACT_SET_FOR_DISPROOF = 10


def main() -> int:
    if not ART.exists():
        print(f"INDETERMINATE: {ART.name} absent", file=sys.stderr)
        return 2
    d = json.loads(ART.read_text(encoding="utf-8"))
    blob = json.dumps(d)
    if "monoculture" not in blob and "growth/no-growth" not in blob and "Orth" not in blob:
        print("INDETERMINATE: cannot confirm the substrate is monoculture", file=sys.stderr)
        return 2

    exact = None
    for key in ("model_scored", "our_model", "ours", "switch", "our_switch"):
        node = d.get(key)
        if isinstance(node, dict):
            inner = node.get("switch", node)
            if isinstance(inner, dict) and "exact_set_match" in inner:
                exact = inner["exact_set_match"]
                break
    if exact is None:
        for v in d.values():
            if isinstance(v, dict) and "exact_set_match" in v:
                exact = v["exact_set_match"]
                break
    if exact is None:
        print("INDETERMINATE: no exact_set_match found in the 4-media artifact", file=sys.stderr)
        return 2

    print(f"4-media (MONOCULTURE) substrate exact-set match: {exact}")
    if exact <= MAX_EXACT_SET_FOR_DISPROOF:
        print("CLAIM KILLED: the model fails to reproduce the switch on a MONOCULTURE substrate too, "
              "so competitive-vs-monoculture modality cannot be what causes the flatness.")
        return 0
    print("CLAIM SURVIVES: the monoculture substrate does NOT show the same failure.")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
