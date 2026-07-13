"""CLI: render a DMS-style protein damage heatmap to standalone HTML (2026-07-13).

Reads a protein sequence + its cached ESM masked-marginals (the SAME cache the rung-2 predictor writes)
and writes ONE self-contained HTML heatmap (offline, no CDN). Read-only; frozen surface untouched.

  # HIV RT, with catalogued DRMs pinned on a separate axis
  uv run python scripts/protein_heatmap.py \
    --cache data/processed/hiv_rt_esm650m_masked_marginals.json \
    --title "HIV-1 RT" --markers K103N,Y181C,M184V \
    --out wiki/heatmap_hiv_rt.html

The colour is the zero-shot ESM molecular rank (damage_llr), NOT a resistance call; catalog markers
are shown separately (the molecular-vs-phenotype distinction is explicit).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from dna_decode.protein_effect.integration import load_logp  # noqa: E402  (handles raw + wrapped + int keys)
from dna_decode.viz.protein_heatmap import build_heatmap_html, build_matrix  # noqa: E402


def _read_seq(cache_path, seq, seq_file) -> str:
    if seq:
        return seq.strip().upper()
    if seq_file:
        lines = [l.strip() for l in open(seq_file, encoding="utf-8") if l.strip() and not l.startswith(">")]
        return "".join(lines).upper()
    # else pull the sequence embedded in the predictor cache ({sequence, model, logp})
    import json
    doc = json.loads(Path(cache_path).read_text(encoding="utf-8"))
    if isinstance(doc, dict) and doc.get("sequence"):
        return doc["sequence"].upper()
    raise SystemExit("no sequence: pass --seq/--seq-file, or use a cache that embeds 'sequence'")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cache", required=True, help="masked-marginal cache JSON (predictor format or raw)")
    ap.add_argument("--seq", default=None, help="protein sequence (else taken from the cache's 'sequence')")
    ap.add_argument("--seq-file", default=None)
    ap.add_argument("--title", default="protein")
    ap.add_argument("--markers", default=None, help="comma-sep catalogued mutations to pin, e.g. K103N,Y181C")
    ap.add_argument("--max-positions", type=int, default=400)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)

    seq = _read_seq(a.cache, a.seq, a.seq_file)
    logp = load_logp(a.cache)
    markers = [m for m in (a.markers.split(",") if a.markers else []) if m.strip()]
    a.out.write_text(build_heatmap_html(seq, logp, title=a.title, markers=markers,
                                        max_positions=a.max_positions), encoding="utf-8")
    mtx = build_matrix(seq, logp)
    print(f"[heatmap] {a.title}: {mtx['length']} aa, {mtx['n_scored_positions']} scored positions, "
          f"scale ±{mtx['scale']} llr -> {a.out}")
    if markers:
        print(f"[markers] {', '.join(markers)} (catalog phenotype axis — separate from the molecular colour)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
