"""Stage A of the mechanism-decomposition arm: export the peak tiles the sweep will decompose.

**Why this is a separate script from the decomposition.** `promoter_calculator` pulls sklearn 1.9.1
while this project is pinned at 1.8.0, so it lives in an ISOLATED interpreter at
`D:/tools/promcalc_env` that cannot import `dna_decode` at all. The substrate therefore has to cross
the process boundary as a file. Stage A (this script, project env) writes it; Stage B
(`promcalc_decompose.py`, isolated env) reads it and never imports anything from this package.

**The tiles come from `genomewide.load_peak_tiles`, which OWNS GSE144621 — this script never opens the
substrate files itself.** That is not incidental tidiness: on 2026-10-13 a second reader of a different
GEO accession silently reproduced a solved join trap, and `tests/test_substrate_loader_ownership.py`
now fails any module that names another owner's accession files. Calling the owner is what that guard
asks for, and the loader's delimiter sniff (the peak-tile file is TAB-delimited while its siblings in
the same series are space-delimited) is exactly the kind of thing a second reader gets wrong.

**The two filters are load-bearing, not cosmetic.** `category == "tile"` and `len(seq) == 150` are the
same two the published baseline used; without them the reconcile anchor (`gc` seed 0 == 0.3000) cannot
pass, and a sweep run over a different row set is not comparable to the incumbent arms.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

from dna_decode.glm.genomewide import DEFAULT_CACHE, load_peak_tiles, noise_ceiling

#: Where stage B looks by default. On D: because the per-site table is ~3M rows (the standing D: rule).
DEFAULT_OUT = Path("D:/dna_decode_cache/promcalc")

#: Asserted by the published baseline; re-asserted here so a drifting loader fails loudly.
EXPECTED_TILES = 44106
EXPECTED_CEILING = 0.941


def export(cache_dir: Path | str = DEFAULT_CACHE, out_dir: Path | str = DEFAULT_OUT,
           *, allow_download: bool = True, limit: int | None = None) -> dict:
    """Write `tiles.tsv` (tile_id, seq) plus a `substrate.json` provenance sidecar."""
    tiles_all, st = load_peak_tiles(cache_dir, allow_download=allow_download)
    # The SAME two filters the published baseline used, or the reconcile anchor cannot pass.
    tiles = [t for t in tiles_all if t.category == "tile" and len(t.seq) == 150]

    if limit is None and len(tiles) != EXPECTED_TILES:
        raise SystemExit(
            f"expected {EXPECTED_TILES} filtered tiles, got {len(tiles)}. The substrate or the filters "
            f"have drifted; a sweep over a different row set is NOT comparable to the incumbent arms."
        )
    ceil = noise_ceiling(tiles)["ceiling_for_a_sequence_model"]

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = tiles if limit is None else tiles[:limit]

    tsv = out_dir / ("tiles_smoke.tsv" if limit else "tiles.tsv")
    with tsv.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["tile_id", "seq"])
        for i, t in enumerate(rows):
            # The index IS the identity: `leave_peak_out` partitions this same filtered list, so
            # stage C rejoins on position. A name would be a second source of truth for the same thing.
            w.writerow([i, t.seq])

    meta = {
        "substrate": "GSE144621 peak tiles",
        "owner_module": "dna_decode/glm/genomewide.py",
        "n_filtered_tiles": len(tiles),
        "n_exported": len(rows),
        "n_parsed": st.n_parsed,
        "filters": ['category == "tile"', "len(seq) == 150"],
        "replicate_ceiling": ceil,
        "tile_length_bp": 150,
        "tile_id_is": "index into the filtered tile list, in load order",
        "smoke": bool(limit),
    }
    if limit is None and abs(ceil - EXPECTED_CEILING) > 1e-9:
        raise SystemExit(f"ceiling {ceil} != published {EXPECTED_CEILING}; substrate drift")
    (out_dir / ("substrate_smoke.json" if limit else "substrate.json")).write_text(
        json.dumps(meta, indent=1), encoding="utf-8")
    return {"tsv": str(tsv), **meta}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cache-dir", default=str(DEFAULT_CACHE))
    ap.add_argument("--out-dir", default=str(DEFAULT_OUT))
    ap.add_argument("--no-download", action="store_true")
    ap.add_argument("--limit", type=int, default=None,
                    help="export only the first N tiles, to a *_smoke.* name (skips the count assert)")
    a = ap.parse_args()
    r = export(a.cache_dir, a.out_dir, allow_download=not a.no_download, limit=a.limit)
    print(json.dumps(r, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
