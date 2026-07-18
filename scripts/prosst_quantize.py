"""CLI over `dna_decode.forward.prosst_scorer.quantize_structure` — the novel-protein ProSST quantize (local).

The plan tagged the ProSST `PdbQuantizer` step as Kaggle-only (torch_geometric wall). It runs LOCALLY on this
Windows/CPU host — the quantize logic + its three shims (pure-python `torch_scatter`, biotite 1.x
`filter_backbone` rename, serial pathos + `num_workers=0`) live in `prosst_scorer.quantize_structure`.
Validated 2026-07-18: self-quantized GRB2 == ProteinGym's pre-quantized tokens 217/217. Needs the cloned repo
at `$PROSST_REPO` (default D:/prosst_repo; bundles AE.pt + {vocab}.joblib).

Run:  uv run python scripts/prosst_quantize.py <pdb_path> [--vocab 2048] [--expect t0,t1,...]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from dna_decode.forward.prosst_scorer import quantize_structure  # noqa: E402


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("pdb")
    ap.add_argument("--vocab", type=int, default=2048)
    ap.add_argument("--expect", default=None, help="comma-separated expected first-N tokens to verify")
    args = ap.parse_args(argv)
    toks = quantize_structure(args.pdb, args.vocab)
    print(f"quantized {Path(args.pdb).name}: {len(toks)} tokens; first 8 = {toks[:8]}")
    if args.expect:
        exp = [int(x) for x in args.expect.split(",")]
        ok = toks[: len(exp)] == exp
        print(f"MATCH expected {exp}: {ok}")
        return 0 if ok else 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
