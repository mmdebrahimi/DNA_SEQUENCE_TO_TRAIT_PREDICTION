"""Render an EXISTING genome-map JSON to a standalone HTML browser.

The live `scripts/genome_map.py` now emits the HTML automatically alongside the .json/.md. This standalone
renderer lets you browse a map that was ALREADY produced (e.g. the committed spike maps under
`wiki/genome_map_spike_*/`) without re-running Bakta/AMRFinder or touching the D: cache — the map JSON is
the only input.

Run:
  uv run python scripts/genome_map_browser.py wiki/genome_map_spike_2026-06-19/genome_map_GCA_002180195.1.json
  # -> writes <same-stem>.html next to the JSON (or pass --out PATH)
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from dna_decode.genome_map.browser import build_genome_map_html  # noqa: E402


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("map_json", type=Path, help="a genome_map_<sample>.json produced by scripts/genome_map.py")
    ap.add_argument("--out", type=Path, default=None, help="output HTML path (default: <map_json>.html)")
    a = ap.parse_args(argv)

    if not a.map_json.exists():
        print(f"ERROR: {a.map_json} not found", file=sys.stderr)
        return 2
    gm = json.loads(a.map_json.read_text(encoding="utf-8"))
    if "features" not in gm:
        print(f"ERROR: {a.map_json} is not a genome-map JSON (no 'features' key)", file=sys.stderr)
        return 2

    out = a.out or a.map_json.with_suffix(".html")
    out.write_text(build_genome_map_html(gm), encoding="utf-8")
    m = gm.get("metrics", {})
    print(f"Wrote {out} | features={m.get('total_features')} "
          f"determinant-phenotype={m.get('determinant_phenotype_feature_count')} "
          f"contigs={len(set(f.get('seqid') for f in gm['features']))}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
