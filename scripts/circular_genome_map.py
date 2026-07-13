"""CLI: render a genome-map JSON as a CIRCULAR (Circos/Proksee-style) ring to standalone HTML (2026-07-13).

Reads a committed genome-map JSON (the same artifact scripts/genome_map.py emits + the linear browser
renders) and writes ONE self-contained circular-ring HTML (offline, no CDN). Read-only; frozen surface
untouched.

  uv run python scripts/circular_genome_map.py \
    wiki/genome_map_spike_2026-06-19/genome_map_GCA_002180195.1.json \
    --out wiki/circular_genome_GCA_002180195.1.html
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from dna_decode.genome_map.circular_browser import build_circular_html  # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("map_json", help="genome-map JSON (genome_map_<acc>.json)")
    ap.add_argument("--size", type=int, default=820)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)

    gm = json.loads(Path(a.map_json).read_text(encoding="utf-8"))
    a.out.write_text(build_circular_html(gm, size=a.size), encoding="utf-8")
    m = gm.get("metrics", {}) or {}
    n_contigs = len({f.get("seqid") for f in gm.get("features", []) if f.get("seqid")})
    print(f"[circular] {gm.get('genome_accession')}: {len(gm.get('features', []))} features / "
          f"{n_contigs} contigs / {m.get('determinant_phenotype_feature_count', '?')} determinant-phenotype "
          f"-> {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
