"""ONE aligned per-tile prediction table that every tile-rung measurement can consume.

    uv run python scripts/glm_tile_predictions.py [--seeds 10] [--scorers gc,max_tx,decomposition,conv]

**Why this exists, and it is not a convenience.** Two published tile-rung arms compute per-tile
predictions and then THROW THEM AWAY. `scripts/glm_encoder_gate.py` is explicit about it --
`r["auroc_active"] = auroc(r.pop("predictions"), lab)`, commented *"then DROP the raw predictions --
one float"* -- and `scripts/promcalc_decomposition_arm.py` aggregates to a per-seed Spearman. That was
the right call for those arms, which only ever needed a correlation. But two new questions need the
VECTORS:

* what the scorer is worth to a CHOICE (propose k, score k, keep the best) -- `dna_decode/glm/selection.py`
* whether a learned and a biophysical scorer are ORTHOGONAL enough to combine -- `rank_average_hybrid`

Neither is derivable from a correlation. Re-fitting twice would also risk the two measurements running
on different splits, which would make any comparison between them meaningless. So one table, one split
per seed, every scorer aligned by construction.

**The join key is POSITION in the filtered tile list, and that is asserted rather than assumed.**
`leave_peak_out` returns the very objects it was handed, not copies, so `id(tile) -> index` is a valid
key -- and stage A (`scripts/promcalc_export_tiles.py`) exported that same index. A splitter that
started copying would mis-join every feature row to the wrong tile AND STILL PRODUCE A PLAUSIBLE
NUMBER, which is why the overlap assert is here and not in a comment.

**Targets are asymmetric and that is deliberate, copied from the published arms:** fit on
`log1p(expression)`, score against RAW `expression`. Mixing those up changes every number quietly.

**The reconcile gate is the whole proof of faithfulness.** Anyone can write a script that produces
Spearman values on this substrate. What makes this table the PUBLISHED substrate rather than something
that merely resembles it is that `gc` seed 0 recomputes to 0.3000 and `max_tx` to 0.3968 -- from the
PERSISTED vectors, not from the fit's own return value. On drift it refuses and writes nothing: a
wrong table that looks right is worse than no table.

A `--smoke` run reads the smoke decomposition shard, so most tiles have no sites and the anchors cannot
hold. It therefore stamps `smoke: true` and SKIPS the gate, and every consumer refuses a smoke table.
A smoke artifact that downstream could mistake for a real one would defeat the gate it bypassed.
"""
from __future__ import annotations

import argparse
import gzip
import json
import math
import sys
from datetime import date
from pathlib import Path

DEFAULT_OUT = Path("D:/dna_decode_cache/promcalc/tile_predictions")

#: Registered scorers. `conv` is the long pole (a real training run per seed); the other three are
#: ridge fits and finish in seconds. Order is the reporting order.
SCORERS = ("gc", "max_tx", "decomposition", "conv")

#: Published anchors, imported in spirit from `promcalc_decomposition_arm.PREREGISTERED["reconcile"]`
#: and restated here with their tolerances because this script must refuse WITHOUT importing the arm.
RECONCILE = {
    "gc": {"seed": 0, "expected": 0.3000, "tolerance": 0.01},
    "max_tx": {"seed": 0, "expected": 0.3968, "tolerance": 0.02},
}

#: The published filter. 44,106 at a fixed 150 bp.
EXPECTED_TILES = 44106

SCHEMA = "glm-tile-predictions-v1"


class TilePredictionError(RuntimeError):
    """The table could not be built in a state any consumer should trust."""


# --------------------------------------------------------------------------------------------------
# pure helpers (no I/O, no model, no D: -- every one of these is unit-tested offline)
# --------------------------------------------------------------------------------------------------
def split_indices(tiles: list, train: list, test: list) -> tuple[list[int], list[int]]:
    """Positions of `train` and `test` within `tiles`, by object identity.

    Raises if a member is absent (the splitter copied, so identity no longer joins), if the halves
    overlap, or if they do not exhaust `tiles` -- each of which would silently mis-align features.
    """
    pos = {id(t): i for i, t in enumerate(tiles)}
    try:
        itr = [pos[id(t)] for t in train]
        ite = [pos[id(t)] for t in test]
    except KeyError as exc:                                      # pragma: no cover - defensive
        raise TilePredictionError(
            "a split member is not identical to any tile in the filtered list; `leave_peak_out` "
            "appears to be returning COPIES, which would mis-join every feature row while still "
            "producing a plausible number"
        ) from exc
    overlap = set(itr) & set(ite)
    if overlap:
        raise TilePredictionError(f"train/test index overlap on {len(overlap)} tiles, e.g. "
                                  f"{sorted(overlap)[:5]} -- the split leaks")
    if len(itr) + len(ite) != len(tiles):
        raise TilePredictionError(
            f"split is not exhaustive: {len(itr)} + {len(ite)} != {len(tiles)}; some tile is in "
            f"neither half and would be silently unscored")
    return itr, ite


def targets(train: list, test: list) -> tuple[list[float], list[float]]:
    """`(log1p(train expression), raw test expression)` -- the published arms' asymmetry, verbatim."""
    return ([math.log1p(t.expression) for t in train], [float(t.expression) for t in test])


def max_tx_scalar(sites: list[dict]) -> list[float]:
    """The single biophysical feature: `log1p(max Tx_rate)`, 0.0 when a tile has no candidate site.

    A site-free tile is a real observation, not an error -- the Promoter Calculator finding no promoter
    is information. It is NOT the same as a tile that was never swept, which is what the coverage gate
    is for.
    """
    best = max((s["Tx_rate"] for s in sites), default=0.0)
    return [math.log1p(max(best, 0.0))]


def reconcile_problems(seed: int, spearmans: dict[str, float],
                       anchors: dict | None = None) -> list[str]:
    """Which published anchors this seed's recomputed Spearmans violate. Empty list == clean.

    Only checks anchors whose scorer was actually run, so a `--scorers gc` invocation is gated on the
    gc anchor alone rather than refusing for a scorer nobody asked for.
    """
    anchors = RECONCILE if anchors is None else anchors
    out: list[str] = []
    for scorer, a in anchors.items():
        if a["seed"] != seed or scorer not in spearmans:
            continue
        got = spearmans[scorer]
        if abs(got - a["expected"]) > a["tolerance"]:
            out.append(f"{scorer} seed {seed} = {got:.4f}, not {a['expected']} "
                       f"+-{a['tolerance']} -- this is not the published substrate")
    return out


def table_payload(*, seed: int, held: object, tile_indices: list[int], measured: list[float],
                  predictions: dict[str, list[float]], spearmans: dict[str, float],
                  smoke: bool) -> dict:
    """The per-seed record. Raises on any ragged vector rather than writing a mis-aligned table."""
    n = len(tile_indices)
    if len(measured) != n:
        raise TilePredictionError(f"measured has {len(measured)} values for {n} tiles")
    for name, vec in predictions.items():
        if len(vec) != n:
            raise TilePredictionError(f"scorer {name!r} has {len(vec)} predictions for {n} tiles")
    if not predictions:
        raise TilePredictionError("no scorer produced predictions; refusing to write an empty table")
    return {
        "schema": SCHEMA, "seed": seed, "smoke": smoke,
        "held_out_peaks": held, "n_test": n,
        "tile_indices": tile_indices, "measured": measured,
        "predictions": {k: predictions[k] for k in sorted(predictions)},
        "spearman_recomputed": {k: round(v, 6) for k, v in sorted(spearmans.items())},
    }


def merge_payload(existing: dict, new: dict) -> dict:
    """Union the scorer columns of two tables for the SAME seed and split.

    **Found by running the script, not by reading it.** The plan's own staged workflow is "cheap
    scorers first, verify the anchors, then commit the expensive conv run" -- and a second invocation
    with different `--scorers` rewrote the seed file from scratch, silently destroying the columns the
    first run had computed. The destroyed table is still perfectly well-formed, so nothing downstream
    would have complained; it would simply have had fewer scorers than the manifest claimed.

    The split must be IDENTICAL or the merge is refused: two tables whose `tile_indices` differ came
    from different splits, and uniting their columns would align one scorer's predictions against
    another scorer's tiles -- a mis-join that still produces a plausible number.
    """
    for key in ("seed", "smoke"):
        if existing.get(key) != new.get(key):
            raise TilePredictionError(
                f"refusing to merge: {key} differs ({existing.get(key)!r} vs {new.get(key)!r})")
    if existing["tile_indices"] != new["tile_indices"]:
        raise TilePredictionError(
            f"refusing to merge seed {new['seed']}: the two tables have DIFFERENT splits "
            f"({len(existing['tile_indices'])} vs {len(new['tile_indices'])} test tiles, or a "
            f"different ordering). Uniting their columns would align one scorer's predictions against "
            f"another scorer's tiles. Delete the old table and rebuild all scorers together.")
    if existing["measured"] != new["measured"]:
        raise TilePredictionError(
            f"refusing to merge seed {new['seed']}: same tile_indices but different measured values "
            f"-- the substrate changed underneath")
    merged = dict(new)
    merged["predictions"] = {**existing["predictions"], **new["predictions"]}
    merged["spearman_recomputed"] = {**existing["spearman_recomputed"],
                                     **new["spearman_recomputed"]}
    merged["predictions"] = {k: merged["predictions"][k] for k in sorted(merged["predictions"])}
    merged["spearman_recomputed"] = {k: merged["spearman_recomputed"][k]
                                     for k in sorted(merged["spearman_recomputed"])}
    return merged


def load_table(path: Path | str) -> dict:
    """Read one per-seed table and refuse a shape no consumer should trust."""
    p = Path(path)
    with gzip.open(p, "rt", encoding="utf-8") as fh:
        d = json.load(fh)
    if d.get("schema") != SCHEMA:
        raise TilePredictionError(f"{p}: schema {d.get('schema')!r}, expected {SCHEMA!r}")
    n = d["n_test"]
    if len(d["tile_indices"]) != n or len(d["measured"]) != n:
        raise TilePredictionError(f"{p}: ragged table")
    for name, vec in d["predictions"].items():
        if len(vec) != n:
            raise TilePredictionError(f"{p}: scorer {name!r} ragged ({len(vec)} vs {n})")
    return d


def require_real(tables: list[dict]) -> None:
    """Every consumer calls this. A smoke table skipped the reconcile gate and cannot be scored."""
    smoky = [t["seed"] for t in tables if t.get("smoke")]
    if smoky:
        raise TilePredictionError(
            f"seeds {smoky} are SMOKE tables: they were built from the smoke decomposition shard and "
            f"skipped the reconcile gate, so their numbers are not the published substrate. Rebuild "
            f"without --smoke before scoring.")


def manifest_payload(*, seeds: list[int], scorers: list[str], n_tiles: int, coverage: float,
                     shards: list, conv_params: int | None, protocol: dict | None,
                     smoke: bool) -> dict:
    """Sidecar identifying WHICH table this is, so two cannot be confused for one another."""
    return {
        "schema": "glm-tile-predictions-manifest-v1",
        "date": str(date.today()),
        "substrate": "GEO GSE144621 peak tiles (E. coli), category=='tile' and len(seq)==150",
        "n_tiles": n_tiles, "expected_tiles": EXPECTED_TILES,
        "split": "leave_peak_out(frac=0.25)", "seeds": seeds, "scorers": scorers,
        "targets": "fit on log1p(expression), score against raw expression",
        "join_key": "position in the filtered tile list (asserted via id() identity + disjointness)",
        "decomposition_coverage": round(coverage, 6), "decomposition_shards": [str(s) for s in shards],
        "conv_parameters": conv_params, "conv_protocol": protocol,
        "smoke": smoke,
        "reconcile_anchors": RECONCILE,
        "honest_limits": [
            "Predictions are out-of-sample on each seed's held-out peaks, but the pool is MEASURED "
            "peak tiles -- any selection number computed from it is an upper bound on a generator's "
            "candidates, which sit off this distribution.",
            "The conv arm's numbers depend on REGISTERED_PROTOCOL; a protocol change makes this "
            "table incomparable with the published encoder-gate figures.",
        ],
    }


# --------------------------------------------------------------------------------------------------
# impure: load, fit, persist
# --------------------------------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--in-dir", default="D:/dna_decode_cache/promcalc",
                    help="stage-B decomposition directory (main file + resume shards)")
    ap.add_argument("--out-dir", default=str(DEFAULT_OUT))
    ap.add_argument("--seeds", type=int, default=10)
    ap.add_argument("--scorers", default=",".join(SCORERS))
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--smoke", action="store_true",
                    help="read the smoke decomposition shard; SKIPS the reconcile gate and stamps "
                         "smoke:true so no consumer can score it")
    ap.add_argument("--no-merge", action="store_true",
                    help="overwrite an existing seed table instead of merging scorer columns into "
                         "it. The DEFAULT is merge, because a staged run (cheap scorers first, then "
                         "conv) would otherwise silently destroy the earlier columns.")
    a = ap.parse_args()

    scorers = [s.strip() for s in a.scorers.split(",") if s.strip()]
    unknown = [s for s in scorers if s not in SCORERS]
    if unknown:
        raise SystemExit(f"unknown scorer(s) {unknown}; known: {list(SCORERS)}")

    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from glm_encoder_gate import ridge_gc
    from promcalc_decomposition_arm import (
        _fit_ridge,
        feature_names,
        load_decomposition,
        tile_features,
    )

    from dna_decode.glm.genomewide import leave_peak_out, load_peak_tiles, noise_ceiling, spearman

    tiles_all, _ = load_peak_tiles()
    tiles = [t for t in tiles_all if t.category == "tile" and len(t.seq) == 150]
    if not a.smoke and len(tiles) != EXPECTED_TILES:
        raise SystemExit(f"{len(tiles)} tiles after the published filter, expected {EXPECTED_TILES} "
                         f"-- the substrate moved; refusing rather than scoring a different one")
    ceil = noise_ceiling(tiles)["ceiling_for_a_sequence_model"]
    print(f"tiles {len(tiles)} | ceiling {ceil:.4f}", flush=True)

    dec, shards = load_decomposition(a.in_dir, smoke=a.smoke)
    cov = len(dec) / len(tiles)
    print(f"decomposition: {len(dec)} tiles | coverage {cov:.4f} | shards {shards}", flush=True)
    if not a.smoke and cov < 0.999:
        raise SystemExit(
            f"coverage {cov:.4f} over {len(tiles)} tiles - the sweep is INCOMPLETE "
            f"({len(tiles) - len(dec)} missing). Every missing tile would enter as a no-promoter "
            f"observation and the table would look clean over a partial substrate. Resume stage B: "
            f"D:/tools/promcalc_env/Scripts/python.exe scripts/promcalc_decompose.py --workers 7 "
            f"--resume")

    need_feats = "decomposition" in scorers
    need_scalar = "max_tx" in scorers
    feats = [tile_features(dec.get(i, [])) for i in range(len(tiles))] if need_feats else []
    scalar = [max_tx_scalar(dec.get(i, [])) for i in range(len(tiles))] if need_scalar else []
    if need_feats:
        nF = len(feature_names())
        if not all(len(f) == nF for f in feats):
            raise SystemExit("ragged feature matrix")

    conv_params = protocol = None
    if "conv" in scorers:
        from dna_decode.glm.encoder import REGISTERED_PROTOCOL, count_parameters
        conv_params = count_parameters()
        protocol = REGISTERED_PROTOCOL.as_dict()
        print(f"conv arm: {conv_params} parameters, protocol {protocol}", flush=True)

    out_dir = Path(a.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []

    for s in range(a.seeds):
        tr, te, held = leave_peak_out(tiles, seed=s)
        itr, ite = split_indices(tiles, tr, te)
        ytr, yte = targets(tr, te)

        preds: dict[str, list[float]] = {}
        if "gc" in scorers:
            _, preds["gc"] = ridge_gc(tr, te)
        if "max_tx" in scorers:
            preds["max_tx"] = _fit_ridge([scalar[i] for i in itr], ytr, [scalar[i] for i in ite])
        if "decomposition" in scorers:
            preds["decomposition"] = _fit_ridge([feats[i] for i in itr], ytr,
                                                [feats[i] for i in ite])
        if "conv" in scorers:
            from dna_decode.glm.encoder import fit_encoder
            preds["conv"] = fit_encoder(tr, te, seed=s, device=a.device)["predictions"]

        # Recomputed FROM the vectors we are about to persist -- not from any fit's own return value.
        sp = {k: spearman(v, yte) for k, v in preds.items()}

        problems = [] if a.smoke else reconcile_problems(s, sp)
        if problems:
            for p in problems:
                print(f"  RECONCILE FAILED: {p}", flush=True)
            raise SystemExit(
                "refusing to write: a table whose anchors drift is not the published substrate, and "
                "a wrong table that looks right is worse than no table")

        payload = table_payload(seed=s, held=held, tile_indices=ite, measured=yte,
                                predictions=preds, spearmans=sp, smoke=a.smoke)
        path = out_dir / f"seed-{s:02d}.json.gz"
        # MERGE rather than clobber: a staged run (cheap scorers, then conv) must not destroy the
        # columns an earlier invocation already computed. The merge refuses on a split mismatch.
        carried: list[str] = []
        if path.is_file() and not a.no_merge:
            prior = load_table(path)
            carried = [k for k in prior["predictions"] if k not in payload["predictions"]]
            payload = merge_payload(prior, payload)
        with gzip.open(path, "wt", encoding="utf-8") as fh:
            json.dump(payload, fh)
        written.append(path)
        shown = payload["spearman_recomputed"]
        bits = "  ".join(f"{k} {shown[k]:+.4f}" for k in sorted(shown))
        extra = f"  [carried over: {','.join(carried)}]" if carried else ""
        print(f"  seed {s}: {bits}  (held {held} peaks, {len(ite)} test){extra}", flush=True)

    # The manifest must describe the table AS WRITTEN, not as requested: after a merge the file holds
    # the union of this run's scorers and any earlier run's. Reporting only `scorers` would under-claim
    # the table, which is the same class of falsehood as over-claiming it.
    present: set[str] = set()
    for p in written:
        present |= set(load_table(p)["predictions"])
    man = manifest_payload(seeds=list(range(a.seeds)), scorers=sorted(present), n_tiles=len(tiles),
                           coverage=cov, shards=shards, conv_params=conv_params,
                           protocol=protocol, smoke=a.smoke)
    man["scorers_this_run"] = scorers
    (out_dir / "manifest.json").write_text(json.dumps(man, indent=2), encoding="utf-8")
    print(f"\nwrote {len(written)} seed tables + manifest to {out_dir}")
    if a.smoke:
        print("SMOKE: reconcile gate skipped; consumers will refuse these tables")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
