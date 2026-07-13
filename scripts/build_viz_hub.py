"""CLI: build the unified visual-decoder hub by discovering generated view HTMLs (2026-07-13).

Scans a directory (default wiki/) for the standalone view HTMLs the other viz CLIs emit
(network_cooccurrence_*, circular_genome_*, genome_map_*, heatmap_*) and writes ONE hub landing page
(offline, no CDN) that ties them into the genome -> network -> protein semantic-zoom story. Read-only;
frozen surface untouched. Hrefs are basenames (hub + views live in the same dir).

  uv run python scripts/build_viz_hub.py --dir wiki --out wiki/index_visual_decoder.html
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from dna_decode.viz.hub import build_hub_html, classify_view_file  # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dir", type=Path, default=REPO / "wiki", help="directory to scan for view HTMLs")
    ap.add_argument("--out", type=Path, default=None, help="default <dir>/index_visual_decoder.html")
    a = ap.parse_args(argv)
    out = a.out or (a.dir / "index_visual_decoder.html")

    views = []
    for p in sorted(a.dir.glob("*.html")):
        if p.resolve() == out.resolve():
            continue                                  # never link the hub to itself
        v = classify_view_file(p.name)
        if v:
            views.append(v)
    out.write_text(build_hub_html(views), encoding="utf-8")
    by = {}
    for v in views:
        by[v["scale"]] = by.get(v["scale"], 0) + 1
    print(f"[hub] {len(views)} views ({', '.join(f'{k}={n}' for k, n in sorted(by.items()))}) -> {out}")
    for v in views:
        print(f"   [{v['scale']:>7}] {v['href']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
