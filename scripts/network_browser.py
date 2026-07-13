"""CLI: render an AMR-determinant co-occurrence network to a standalone HTML (2026-07-13).

Reads the committed co-occurrence + (optional) cross-axis-lineage artifacts and writes ONE self-contained
HTML file (force-directed, offline, no CDN). Read-only over committed artifacts; frozen surface untouched.

  uv run python scripts/network_browser.py \
    --cooc wiki/determinant_cooccurrence_result_2026-07-11.json \
    --organism escherichia_coli_shigella \
    --crossaxis wiki/crossaxis_lineage_deconfound_determinant_2026-07-12.json \
    --out wiki/network_cooccurrence_ecoli.html

The de-confound (leave-one-clade-out) is rendered as node-border + edge style: solid = generalizes,
dashed = lineage-mediated (possible clonal structure), dotted = untested.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from dna_decode.viz.network_adapter import build_graph  # noqa: E402
from dna_decode.viz.network_browser import build_network_html  # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cooc", required=True, help="determinant_cooccurrence_result_*.json")
    ap.add_argument("--organism", required=True, help="e.g. escherichia_coli_shigella / klebsiella")
    ap.add_argument("--crossaxis", default=None,
                    help="crossaxis_lineage_deconfound_determinant*.json for the SAME organism (node border)")
    ap.add_argument("--min-cooc", type=int, default=8, help="prune edges below this co-occurrence count")
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)

    graph = build_graph(a.cooc, a.organism, a.crossaxis, min_cooc=a.min_cooc)
    a.out.write_text(build_network_html(graph), encoding="utf-8")
    m = graph["meta"]
    print(f"[network] {a.organism}: {m['n_nodes']} nodes / {m['n_edges']} edges "
          f"({m['n_lineage_mediated_nodes']} lineage-mediated) -> {a.out}")
    print(f"[provenance] {m['cooc_artifact']}"
          + (f" + {m['crossaxis_artifact']}" if m["crossaxis_artifact"] else "")
          + f" · verdict={m['verdict']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
