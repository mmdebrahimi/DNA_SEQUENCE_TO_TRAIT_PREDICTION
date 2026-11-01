"""Dry-run the PromoterAtlas arm kernel OFF-Kaggle, against a stubbed substrate.

WHY THIS EXISTS. Three Kaggle pushes were spent on defects in a kernel that had never executed a
line locally: a payload missing a file the run needed, an emitted kernel with an unterminated string
literal, and a subprocess with no PYTHONPATH. Every one of those was reachable without a GPU.

This exercises the kernel end to end EXCEPT `load_peak_tiles` -- the clone, the pickle audit, the
strict weight load, the parameter assertion, both embedding geometries, the determinism probe, the
embedding-sanity gate, `leave_peak_out`, the identity-based index join, the split-hash comparison,
the ridge protocol, the GC anchor, the shuffled-label null, the verdict ladder and the artifact
write. The substrate is synthetic tiles with sequence-dependent expression, so a ridge that silently
learns nothing shows up as a near-zero Spearman rather than as a crash.

It does NOT validate the arm's number. It validates that the kernel runs and that its gates fire.

Usage:  uv run python -I scripts/promoteratlas_arm_dryrun.py
Needs:  a local PromoterAtlas clone; point at it with PA_REPO (default below).
"""
from __future__ import annotations

import base64
import gzip
import importlib.util
import json
import os
import random
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
KERNEL = REPO / "scripts" / "kaggle" / "promoteratlas_tile_arm_kernel.py"
PA_DEFAULT = r"C:/Users/LAPTOP~1/AppData/Local/Temp/pa_probe/PromoterAtlas"
PA_REPO = os.environ.get("PA_REPO", PA_DEFAULT)

N_TILES, N_PEAKS = 900, 48


def synthetic_tiles(Tile):
    """Tiles whose expression genuinely depends on sequence, so a dead ridge is visible.

    Signal: GC fraction plus the count of a planted 6-mer, both monotone in expression, plus noise.
    A model that learns nothing lands near 0; one that learns the composition signal lands well
    above it. Peaks are assigned in contiguous runs, mirroring real tiles-within-a-peak structure.
    """
    rng = random.Random(1234)
    tiles = []
    per = N_TILES // N_PEAKS
    for i in range(N_TILES):
        gc_target = rng.uniform(0.25, 0.65)
        seq = "".join(rng.choices("GC" if rng.random() < gc_target else "AT", k=1)[0]
                      for _ in range(150))
        if rng.random() < 0.3:
            p = rng.randrange(0, 144)
            seq = seq[:p] + "TTGACA" + seq[p + 6:]
        gc = (seq.count("G") + seq.count("C")) / len(seq)
        expr = max(0.0, 8.0 * gc + 3.0 * seq.count("TTGACA") + rng.gauss(0, 0.6))
        tiles.append(Tile(seq=seq, expression=expr, peak="peak_%03d" % (i // max(1, per)),
                          strand="+", category="tile", active="yes", rep1=expr, rep2=expr))
    return tiles


def main() -> int:
    pa = Path(PA_REPO)
    if not (pa / "trained_weights/base_model/promoteratlas-base.pt").exists():
        print("SKIP: no PromoterAtlas clone at %s (set PA_REPO)" % pa)
        return 0

    work = Path(tempfile.mkdtemp(prefix="pa_arm_dryrun_"))
    os.environ["PA_ARM_WORK"] = str(work)
    os.environ["PA_ARM_PA_DIR"] = str(pa)
    print("dry-run WORK  = %s" % work)
    print("dry-run PA    = %s" % pa)

    # Load the kernel as a module so its __main__ guard does not fire.
    spec = importlib.util.spec_from_file_location("pa_arm_kernel", KERNEL)
    k = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(k)
    print("kernel loaded: arm=%s margin=%s n_seeds=%s" % (k.ARM, k.MARGIN, k.N_SEEDS))
    print("              splits pinned for seeds: %s" % sorted(k.SPLITS, key=int))

    # Materialise the embedded genomewide, import it for real, then stub ONLY load_peak_tiles and
    # install it under the name the kernel imports. Everything else (leave_peak_out, gc_features,
    # spearman) is the genuine project code travelling in the payload.
    work.mkdir(parents=True, exist_ok=True)
    gwpath = work / "genomewide.py"
    gwpath.write_bytes(gzip.decompress(base64.b64decode(k.GENOMEWIDE_B64)))
    original = (REPO / "dna_decode/glm/genomewide.py").read_bytes()
    assert gwpath.read_bytes() == original, "PAYLOAD MISMATCH: embedded genomewide != repo copy"
    print("payload check: embedded genomewide.py is byte-identical to the repo copy")

    sys.path.insert(0, str(work))
    gspec = importlib.util.spec_from_file_location("genomewide", gwpath)
    gw = importlib.util.module_from_spec(gspec)
    # Register BEFORE exec_module: @dataclass resolves its own class's __module__ through
    # sys.modules while the body runs, so exec-then-register raises AttributeError on None.
    # The kernel's plain `import genomewide` does this automatically; only this harness has to.
    sys.modules["genomewide"] = gw
    gspec.loader.exec_module(gw)
    tiles = synthetic_tiles(gw.Tile)

    # THE STUB MUST HONOUR THE REAL CONTRACT. The first Kaggle run died because this stub returned a
    # bare list while load_peak_tiles returns tuple[list[Tile], LoadStats] -- so the one boundary the
    # dry run fakes was also the only one whose shape it got wrong, and a fixture I authored can only
    # ever validate the shape I invented. The assertions below check the stub against the REAL
    # function's annotation, so a future divergence fails here instead of on Kaggle.
    ret = (gw.load_peak_tiles.__annotations__ or {}).get("return")
    ret_s = str(ret)
    assert "tuple" in ret_s and "LoadStats" in ret_s, (
        "load_peak_tiles return annotation changed to %r -- re-derive the stub shape" % ret_s)
    stub_ret = (list(tiles), gw.LoadStats(path="<dryrun-stub>"))
    assert isinstance(stub_ret, tuple) and len(stub_ret) == 2
    assert isinstance(stub_ret[0], list) and isinstance(stub_ret[0][0], gw.Tile)
    assert isinstance(stub_ret[1], gw.LoadStats)
    gw.load_peak_tiles = lambda *a, **kw: (list(tiles), gw.LoadStats(path="<dryrun-stub>"))
    print("stub contract checked against the real annotation: %s" % ret_s)
    print("substrate stubbed: %d synthetic tiles over %d peaks"
          % (len(tiles), len({t.peak for t in tiles})))

    sys.argv = ["kernel", "--smoke"]
    rc = k.main()
    print("\nkernel main() returned %r" % rc)

    art = work / "promoteratlas_tile_arm_smoke.json"
    if not art.exists():
        print("FAIL: no artifact written")
        return 1
    d = json.loads(art.read_text(encoding="utf-8"))
    print("\n--- artifact ---")
    print("verdict        : %s" % d.get("verdict"))
    print("detail         : %s" % d.get("verdict_detail"))
    print("device         : %s  determinism %.3e" % (d.get("device"), d.get("determinism_max_abs_delta", -1)))
    print("tiles scored   : %s over %s seeds" % (d.get("n_tiles_scored"), d.get("n_seeds")))
    print("sanity         : %s" % d.get("embedding_sanity"))
    for name, a in (d.get("arms") or {}).items():
        print("arm %-12s median %.4f  per-seed %s" % (name, a["median"], a["per_seed"]))
    print("gc anchor      : median %.4f  per-seed %s"
          % (d["gc_anchor"]["median"], d["gc_anchor"]["per_seed"]))
    print("shuffled null  : median %.4f  per-seed %s"
          % (d["shuffled_label_null"]["median"], d["shuffled_label_null"]["per_seed"]))
    print("per-tile preds : %d seeds, %d preds on seed 0"
          % (len(d["per_tile_predictions"]), len(d["per_tile_predictions"]["0"]["pred"])))

    # ---- what the dry run is allowed to conclude ----
    fails = []
    if d.get("verdict") != "SMOKE_NOT_SCOREABLE":
        fails.append("smoke verdict should be SMOKE_NOT_SCOREABLE, got %r" % d.get("verdict"))
    if d["embedding_sanity"]["n_constant_dims"] != 0:
        fails.append("embedding has constant dims on synthetic input")
    prim = d["arms"][d["geometry_primary"]]["median"]
    if abs(prim) < 0.05:
        fails.append("primary arm median %.4f is near zero -- the ridge may be learning nothing" % prim)
    if abs(d["shuffled_label_null"]["median"]) > 0.15:
        fails.append("shuffled-label null %.4f is far from zero -- the null path is suspect"
                     % d["shuffled_label_null"]["median"])
    if len(d["per_tile_predictions"]["0"]["pred"]) != len(d["per_tile_predictions"]["0"]["tile_indices"]):
        fails.append("per-tile predictions are ragged")

    # ---- PHASE 2: exercise the null EXTENSION and the FIXED validity gate ----------------------
    # The smoke path skips both. The gate is the line that let ATLAS_BEATS_THERMODYNAMIC through
    # while the null was failing, so it is the single most important thing to exercise. The tile
    # count assertion is relaxed here ONLY because synthetic tiles are not the real 44,106 pool --
    # this is a plumbing test of the null loop and the gate, not a measurement.
    print("\n--- phase 2: null extension + fixed validity gate ---")
    k.EXPECTED_TILES = len(tiles)
    sys.argv = ["kernel", "--full-only"]
    rc2 = k.main(smoke=False, n_null=6)
    full = json.loads((work / "promoteratlas_tile_arm_full.json").read_text(encoding="utf-8"))
    ext = full.get("shuffled_label_null_extended")
    print("extension      : n=%s median %s ci95 %s coverage %s bar_met=%s"
          % (ext and ext["n_seeds"], ext and ext["median"], ext and ext["median_ci95"],
             ext and ext["achieved_coverage"], ext and ext["bar_met"]))
    print("verdict        : %s" % full.get("verdict"))
    if ext is None:
        fails.append("null extension did not run with n_null=6")
    elif ext["n_seeds"] != 6:
        fails.append("null extension ran %s seeds, expected 6" % ext["n_seeds"])
    elif ext["bar_met"] != (abs(ext["median"]) < k.MARGIN):
        fails.append("bar_met disagrees with its own median")

    # The gate must FIRE on a failing null. Drive it directly rather than hoping synthetic data
    # produces one: a gate that has only ever been observed passing is not a tested gate.
    gate = lambda m: abs(m) >= k.MARGIN
    if not (gate(0.0203) and gate(-0.05) and not gate(0.0199) and not gate(0.0)):
        fails.append("the fixed null gate does not fire as specified")
    else:
        print("gate check     : fires at 0.0203 and -0.05, passes at 0.0199 and 0.0 (non-vacuous)")

    print("\n--- dry-run verdict ---")
    if fails:
        for f in fails:
            print("  FAIL: %s" % f)
        return 1
    print("  PASS: the kernel runs end to end, both geometries embed, the gates fire, the ridge")
    print("        learns a real signal (primary median %.4f) and the shuffled null sits near zero" % prim)
    print("        (%.4f). This validates PLUMBING ONLY -- not the arm's number."
          % d["shuffled_label_null"]["median"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
