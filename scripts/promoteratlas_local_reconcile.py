"""Reproduce the PromoterAtlas arm's committed numbers LOCALLY, on CPU, or refuse.

This is a GATE, not a measurement. Before any new number is read off this encoder on this host, the
local harness must reproduce two committed figures from `wiki/promoteratlas_tile_arm_2026-11-01.json`:

  * the GC pipeline anchor -- seed 0 at 0.3000 +/- 0.01 and the 10-seed median 0.3112
  * the arm itself -- `arms.native150` median 0.5076

**The GC anchor is the load-bearing one.** Reproducing a known number to four decimals validates the
loader, the `category == 'tile' and len == 150` filter, `leave_peak_out`, the `log1p` target, the
ridge protocol and the scorer TOGETHER. A pipeline that cannot reproduce a known number cannot
establish a new one, which is why this script exits non-zero rather than printing a caveat.

Offline: substrate from `data/mpra/`, encoder from the local clone, no network, no GPU, no Kaggle.
Measured on this host: ~77 seq/s on 8 CPU threads, so the single embedding pass over 44,106 tiles
takes about 10 minutes. The pass is cached to `data/promoteratlas/` (gitignored) so the downstream
experiments never re-embed -- and so a float dump never lands in `wiki/`, which is the decoy pool
`scripts/contract_number_audit.py` uses and which a 102,755-decimal artifact broke on 2026-11-01.

Exit codes: 0 reconciled; 2 a committed figure did not reproduce; 3 the substrate or clone is absent.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
ARM = REPO / "wiki/promoteratlas_tile_arm_2026-11-01.json"
CACHE = REPO / "data/mpra"
EMB_DIR = REPO / "data/promoteratlas"

#: Tolerances. The GC anchor's +/- 0.01 is the arm's own `gc_anchor.tol`. The atlas tolerance is
#: tighter because both runs are deterministic on fixed splits and the arm reproduced itself
#: EXACTLY across two independent Kaggle runs -- so a drift here is a harness difference, not noise.
GC_TOL = 0.01
ATLAS_TOL = 0.02

EXPECTED_TILES = 44106
N_SEEDS = 10
FRAC = 0.25


def _fail(code: int, msg: str) -> int:
    print("REFUSED: " + msg, flush=True)
    return code


def load_tiles():
    from dna_decode.glm import genomewide as gw

    tiles_all, st = gw.load_peak_tiles(CACHE, allow_download=False)
    tiles = [t for t in tiles_all if t.category == "tile" and len(t.seq) == 150]
    print(
        "substrate: %d rows parsed (delim=%s, cols=%s), %d after filter(category=='tile', len==150)"
        % (st.n_parsed, st.delimiter, st.n_header_columns, len(tiles)),
        flush=True,
    )
    return tiles


def embed_tiles(seqs, *, force=False):
    """One embedding pass over the tiles, cached to a gitignored .npy."""
    from dna_decode.glm.atlas_embed import load_atlas, masked_mean_pool

    EMB_DIR.mkdir(parents=True, exist_ok=True)
    out = EMB_DIR / ("tile_embeddings_native150_n%d.npy" % len(seqs))
    if out.is_file() and not force:
        X = np.load(out)
        if X.shape == (len(seqs), 128):
            print("embeddings: reusing cache %s %s" % (out.name, X.shape), flush=True)
            return X, out
        print("embeddings: cache shape %s != expected, re-embedding" % (X.shape,), flush=True)

    model = load_atlas()
    t0 = time.time()

    def prog(done, total):
        if done % 6400 == 0 or done == total:
            el = time.time() - t0
            rate = done / el if el else 0.0
            print(
                "   embedded %d/%d (%.0fs, %.1f seq/s, eta %.1f min)"
                % (done, total, el, rate, (total - done) / rate / 60 if rate else 0.0),
                flush=True,
            )

    X = masked_mean_pool(model, seqs, device="cpu", batch=64, progress=prog)
    np.save(out, X)
    print("embeddings: %s saved in %.1f min" % (X.shape, (time.time() - t0) / 60), flush=True)
    return X, out


def run(force_embed=False):
    if not ARM.is_file():
        return _fail(3, "committed arm artifact absent at %s" % ARM)
    arm = json.loads(ARM.read_text(encoding="utf-8"))
    want_gc_seed0 = float(arm["gc_anchor"]["seed0"])
    want_gc_med = float(arm["gc_anchor"]["median"])
    want_atlas = float(arm["arms"]["native150"]["median"])
    print(
        "committed targets: gc.seed0=%.4f gc.median=%.4f atlas.native150=%.4f"
        % (want_gc_seed0, want_gc_med, want_atlas),
        flush=True,
    )

    from dna_decode.glm import genomewide as gw
    from dna_decode.glm.atlas_embed import AtlasEmbedError, ridge_spearman

    try:
        tiles = load_tiles()
    except Exception as exc:
        return _fail(3, "%s: %s" % (type(exc).__name__, exc))
    if len(tiles) != EXPECTED_TILES:
        return _fail(2, "%d tiles, the arm scored %d" % (len(tiles), EXPECTED_TILES))

    seqs = [t.seq for t in tiles]
    try:
        X, emb_path = embed_tiles(seqs, force=force_embed)
    except AtlasEmbedError as exc:
        return _fail(3, str(exc))

    gc = np.asarray(gw.gc_features(seqs), dtype=float)
    idx = {id(t): i for i, t in enumerate(tiles)}

    gc_per, atlas_per, n_test = {}, {}, {}
    for seed in range(N_SEEDS):
        train, test, _ = gw.leave_peak_out(tiles, frac=FRAC, seed=seed)
        tr = [idx[id(t)] for t in train]
        te = [idx[id(t)] for t in test]
        ytr = [tiles[i].expression for i in tr]
        yte = [tiles[i].expression for i in te]
        n_test[seed] = len(te)
        g, _ = ridge_spearman(gc[tr], ytr, gc[te], yte)
        a, _ = ridge_spearman(X[tr], ytr, X[te], yte)
        gc_per[seed], atlas_per[seed] = round(g, 4), round(a, 4)
        print("   seed %d: n_test=%d gc=%.4f atlas=%.4f" % (seed, len(te), g, a), flush=True)

    gc_med = float(np.median([gc_per[s] for s in range(N_SEEDS)]))
    atlas_med = float(np.median([atlas_per[s] for s in range(N_SEEDS)]))

    checks = [
        ("gc.seed0", gc_per[0], want_gc_seed0, GC_TOL),
        ("gc.median", gc_med, want_gc_med, GC_TOL),
        ("atlas.native150.median", atlas_med, want_atlas, ATLAS_TOL),
    ]
    print("", flush=True)
    bad = []
    for name, got, want, tol in checks:
        ok = abs(got - want) <= tol
        print(
            "%-24s got %.4f  want %.4f  tol %.2f  %s"
            % (name, got, want, tol, "OK" if ok else "MISMATCH"),
            flush=True,
        )
        if not ok:
            bad.append(name)

    report = {
        "schema": "promoteratlas-local-reconcile-v1",
        "date": "2026-11-01",
        "device": "cpu",
        "n_tiles": len(tiles),
        "n_seeds": N_SEEDS,
        "frac": FRAC,
        "embedding_cache": str(emb_path.relative_to(REPO)).replace("\\", "/"),
        "n_test_per_seed": n_test,
        "gc_per_seed": gc_per,
        "atlas_per_seed": atlas_per,
        "gc_median": round(gc_med, 4),
        "atlas_median": round(atlas_med, 4),
        "targets": {
            "gc_seed0": want_gc_seed0,
            "gc_median": want_gc_med,
            "atlas_native150_median": want_atlas,
        },
        "mismatched": bad,
        "verdict": "RECONCILED" if not bad else "HARNESS_DRIFT",
    }
    out = EMB_DIR / "local_reconcile_2026-11-01.json"
    out.write_text(json.dumps(report, indent=2, sort_keys=True), encoding="utf-8")
    print("\nreport -> %s" % out, flush=True)

    if bad:
        return _fail(
            2,
            "the local harness does NOT reproduce %s. No exposure claim may rest on a harness that "
            "cannot reproduce the published number -- diagnose the harness, do not widen the "
            "tolerance." % bad,
        )
    print("RECONCILED: the local CPU harness reproduces the committed figures.", flush=True)
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--force-embed", action="store_true", help="ignore the cached embedding pass")
    ap.add_argument(
        "--perturb-target",
        type=float,
        default=None,
        help="NON-VACUITY CHECK: pretend the committed atlas median is this value. The gate must "
        "then FAIL. Proves the gate can reject, rather than merely being present.",
    )
    a = ap.parse_args(argv)
    if a.perturb_target is not None:
        global ATLAS_TOL
        arm = json.loads(ARM.read_text(encoding="utf-8"))
        real = float(arm["arms"]["native150"]["median"])
        print(
            "NON-VACUITY: committed median is %.4f; a target of %.4f must make this gate FAIL"
            % (real, a.perturb_target),
            flush=True,
        )
        ATLAS_TOL = abs(real - a.perturb_target) / 2.0
    return run(force_embed=a.force_embed)


if __name__ == "__main__":
    sys.exit(main())
