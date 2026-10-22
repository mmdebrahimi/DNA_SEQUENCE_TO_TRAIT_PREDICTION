"""Stage A of the mechanism-decomposition arm: export the peak tiles the sweep will decompose.

**Why this is a separate script from the decomposition.** Stage B runs under an ISOLATED
interpreter that cannot import `dna_decode` at all, so the substrate has to cross the process
boundary as a file. Stage A (this script, project env) writes it; Stage B
(`promcalc_decompose.py`, isolated env) reads it and never imports anything from this package.

**The reason for that isolation was MEASURED 2026-10-22, and the two halves of the original claim
have DIFFERENT truth values -- do not collapse them.** The old wording said `promoter_calculator`
pulls sklearn 1.9.1 against this project's pinned 1.8.0, and concluded that installing it here would
move sklearn under the whole suite.

  * TRUE: the distribution declares no hard sklearn pin, so an UNCONSTRAINED resolve takes the newest.
    A fresh `uv venv` + `uv pip install promotercalculator` really did install scikit-learn 1.9.1
    (with numpy 2.4.6, scipy 1.17.1, biopython 1.88).
  * FALSE: that does NOT make it unsatisfiable here. `uv pip install --dry-run promotercalculator` in
    the PROJECT env resolves 7 packages, would install exactly ONE, and leaves sklearn at 1.8.0 --
    because 1.8.0 already satisfies it. Nothing would move under the suite.
  * So the stated CONSEQUENCE is refuted while its antecedent holds. A first pass at this correction
    called the whole premise "REFUTED", which over-corrected in the other direction; the precise
    statement is the three bullets above.

The isolation is still right, for a different and now-primary reason: **`uv.lock` is a frozen
reproducibility unit**, so `uv add` is unavailable and an imperative `uv pip install` into the project
env is silently reverted by the next `uv sync`. A separate venv is near-free and keeps the boundary
`tests/test_promcalc_decompose.py` now enforces.

**The interpreter path moved.** `D:/tools/promcalc_env` was on the volume that detached permanently.
The reachable one is `C:/Users/Laptop Portable/tools/promcalc_env/Scripts/python.exe`. Commands,
measured cost and honest limits: `wiki/promcalc_regeneration_runbook_2026-10-22.md`.

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


def _require_creatable_out_dir(out_dir: Path | str, *, flag: str = "--out-dir") -> Path:
    """`out_dir`, or a refusal naming the flag instead of an OSError from `mkdir`.

    **Why not `resolve_external_root` here.** That resolver reports `configured` only when the path
    ALREADY is a directory, which is right for an INPUT root and wrong for this one: on a first export
    the leaf legitimately does not exist yet, so asking the resolver about it would refuse a perfectly
    good run. What actually has to be true for `mkdir(parents=True)` to work is that the nearest
    EXISTING ancestor is reachable -- on a detached volume no ancestor resolves at all, and `mkdir`
    then raises a bare WinError that names neither the flag nor where the data really is.

    The D: default is unchanged, so a host that has the volume takes the identical path.
    """
    d = Path(out_dir)
    for anc in [d, *d.parents]:
        if anc.exists():
            if anc.is_dir():
                return d
            raise SystemExit(
                f"{flag} {d}\n\n"
                f"  {anc} exists but is a FILE, so it cannot hold this directory.\n\n"
                f"Pass a directory path."
            )
    from dna_decode.data.external_root import describe

    raise SystemExit(
        f"{flag} {d} is not reachable -- no ancestor of it exists, which is what a detached volume "
        f"looks like.\n\n{describe(d)}\n\n"
        f"Pass {flag} explicitly. The root is NOT substituted automatically: a directory of the right "
        f"NAME is not evidence of the right CONTENTS."
    )


def export(cache_dir: Path | str = DEFAULT_CACHE, out_dir: Path | str = DEFAULT_OUT,
           *, allow_download: bool = True, limit: int | None = None) -> dict:
    """Write `tiles.tsv` (tile_id, seq) plus a `substrate.json` provenance sidecar."""
    # Checked BEFORE the ~40 s substrate load, not at the later mkdir: an unreachable output volume
    # is knowable immediately, and spending the load first only to die on mkdir wastes the operator's
    # time and buries the real cause under a traceback.
    out_dir = _require_creatable_out_dir(out_dir)
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
