"""Score the as-given arm. One rung, one interval, one verdict -- and the same three gates first.

The arm reuses the COMMITTED stored block on `D:`, so it costs minutes rather than a forward pass:
`bins10_mean_asgiven` is `bins10_mean` sliced from that block and flattened WITHOUT `orient`, which
is the entire difference from the positional arm's primary.

**The three validity gates are IMPORTED, not re-implemented.** `gate_split_equality`, `gate_parity`
and `gate_shuffled_label_null` come from `positional_arm_score`, driven by this arm's own
`CHECK_ORDER`, with `ScoreCtx.primary_rung` pointing the null at THIS arm's primary. A second copy
of those gates would be a second chance to drift from the parity check that makes any of these
numbers comparable to 0.5076.

**What this arm cannot establish, restated at the top of the runner because it is easy to lose.**
The number is already visible -- 0.5086 was the positional run's ineligible secondary -- so this is
not a blind pre-registration. The bar carries that hazard with its figures and predicts a TIE. The
verdict here is still worth computing, because clearing the margin as a declared PRIMARY through the
full gate chain is a different claim from appearing well-placed as an ineligible secondary.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

from scripts.freeze_asgiven_arm_bar import (
    ARM,
    CHECK_ORDER,
    OUTCOMES,
    PRIMARY,
    classify,
)
from scripts.freeze_positional_arm_bar import MARGIN
from scripts.positional_arm_extract import DEFAULT_BLOCK_DIR, load_block
from scripts.positional_arm_probe import (
    RAM_HEADROOM,
    design_matrix_bytes,
    fits_in_ram,
    free_ram_bytes,
    load_arm_cohort,
)
from scripts.positional_arm_score import (
    GATE_FUNCS,
    ScoreCtx,
    ScoreError,
    build_rung,
    leave_peak_out_indices,
    median_of,
    paired_bootstrap,
)

REPO = Path(__file__).resolve().parents[1]
BAR_PATH = REPO / "wiki/asgiven_arm_acceptance_bar.json"
OUT = REPO / "wiki/asgiven_arm_result_2026-11-02.json"
#: per-tile predictions go under data/, never wiki/ -- a large decimal dump there has previously
#: polluted the contract-number audit's decoy pool and broken two of its guards
PERTILE_OUT = REPO / "data/promoteratlas/asgiven_arm_pertile_2026-11-02.json"

N_BOOT = 10_000

#: Gate name -> the outcome a failure returns, taken from THIS arm's frozen bar.
GATE_SEQUENCE = tuple((g, o) for g, o in CHECK_ORDER if isinstance(o, str))


def verdict_from_artifact(rec: dict, bar: dict) -> tuple[str, str]:
    """Re-derive the verdict from the emitted numbers plus the frozen anchors, and nothing else."""
    gates = rec.get("gates", {})
    margin = bar["inherited_not_chosen"]["margin"]["value"]
    for gate, outcome in GATE_SEQUENCE:
        g = gates.get(gate)
        if g is None:
            return (outcome, "gate %r was never evaluated, so the run stopped at or before it"
                    % gate)
        if not g.get("passed"):
            return (outcome, g.get("detail", "gate %r failed" % gate))
    r = rec.get("rung")
    if r is None or r.get("skipped"):
        return ("ASGIVEN_TIES",
                "the primary rung was not scored (skipped=%r), so neither a gain of at least the "
                "margin nor a loss is demonstrated" % (r or {}).get("skipped"))
    b = r["gain_bootstrap"]
    return classify(b["ci_lo"], b["ci_hi"], margin=margin)


def run(*, block_dir: Path, bar: dict, null_extension: int, n_boot: int,
        scorer=None, cohort=None, block=None, emb_cache: Path | None = None,
        verbose: bool = True) -> dict:
    from dna_decode.glm.atlas_embed import ridge_spearman

    t0 = time.time()
    inh = bar["inherited_not_chosen"]
    tiles, st = cohort if cohort is not None else load_arm_cohort()
    arr, side, layout = load_block(block_dir) if block is None else block
    if side["geometry"]["n_rows"] != len(tiles):
        raise ScoreError(
            "the block holds %d rows and the cohort is %d tiles; the block is indexed positionally "
            "against the cohort, so this pairing is impossible."
            % (side["geometry"]["n_rows"], len(tiles))
        )

    seeds = list(range(inh["split"]["n_seeds"]))
    frac = inh["split"]["frac"]
    pct = bar["multiplicity"]["ci_percentiles"]
    ctx = ScoreCtx(
        tiles=tiles, st=st, arr=arr, side=side, layout=layout,
        strands=[t.strand for t in tiles],
        splits={s: leave_peak_out_indices(tiles, frac=frac, seed=s) for s in seeds},
        seeds=seeds, score=scorer or ridge_spearman, bar=bar, inh=inh,
        committed_atlas=inh["comparator"]["atlas_per_seed"],
        committed_n_test=inh["split"]["n_test_per_seed"],
        parity_tol=bar["validity_controls"]["parity"]["tolerance"],
        co_pct=pct, sec_pct=pct,
        headroom=RAM_HEADROOM,
        frac=frac, null_extension=null_extension, n_boot=n_boot, emb_cache=emb_cache,
        verbose=verbose, t0=t0,
        primary_rung=PRIMARY,          # <- the one thing that differs from the positional run
    )

    rec: dict = {
        "schema": "asgiven-arm-result-v1",
        "date": "2026-11-02",
        "bar": "wiki/asgiven_arm_acceptance_bar.json",
        "arm": ARM,
        "criterion": bar["criterion"],
        "primary_rung": PRIMARY,
        "block": {"path": side["block"], "sidecar_status": side["status"],
                  "checkpoint_sha256": side.get("checkpoint_sha256"),
                  "geometry": side["geometry"]},
        "substrate": {"n_tiles": len(tiles), "n_peaks": len({t.peak for t in tiles}),
                      "path": st.path},
        "free_ram_bytes_at_start": free_ram_bytes(),
        "gates": {},
    }

    for gate_name, outcome in GATE_SEQUENCE:
        g = GATE_FUNCS[gate_name](ctx)
        rec["gates"][gate_name] = g
        if not g["passed"]:
            rec["verdict"], rec["verdict_detail"] = outcome, g["detail"]
            rec["wall_seconds"] = round(time.time() - t0, 1)
            return rec

    # ---- the single primary --------------------------------------------------------------------
    X = ctx.primary_matrix if ctx.primary_matrix is not None else build_rung(
        arr, layout, PRIMARY, ctx.strands)
    need = design_matrix_bytes(len(tiles), X.shape[1])
    free_now = free_ram_bytes()
    affordable = fits_in_ram(need, free_now, headroom=RAM_HEADROOM)
    if affordable is not True:
        rec["rung"] = {
            "skipped": True, "n_features": int(X.shape[1]), "peak_design_bytes": need,
            "free_ram_bytes": free_now, "fits_in_ram": affordable,
            "why": ("the RAM pre-check did not return True (None means UNMEASURED, which is not "
                    "the same as affordable). A rung nobody could afford to score is not evidence "
                    "of no gain, which is why this yields TIES rather than NO_GAIN."),
        }
    else:
        per_seed, gains, wins, pertile = {}, [], 0, {}
        for s in seeds:
            tr, te, _ = ctx.splits[s]
            sp, pred = ctx.score(X[tr], ctx.ys(tr), X[te], ctx.ys(te))
            sp = round(float(sp), 4)
            want = float(ctx.committed_atlas[str(s)])
            per_seed[str(s)] = {"spearman": sp, "atlas_committed": want,
                                "gain": round(sp - want, 4)}
            gains.append(sp - want)
            wins += int(sp > want)
            pertile[str(s)] = {"tile_indices": te, "pred": [round(float(v), 6) for v in pred]}
        boot = paired_bootstrap(gains, percentiles=pct, n_boot=n_boot)
        rec["rung"] = {
            "skipped": False,
            "rung": PRIMARY,
            "role": "the arm's ONLY primary; no secondaries are declared",
            "n_features": int(X.shape[1]),
            "peak_design_bytes": need,
            "median_spearman": median_of(v["spearman"] for v in per_seed.values()),
            "per_seed": per_seed,
            "seed_wins_over_atlas": wins,
            "n_seeds": len(seeds),
            "gain_bootstrap": boot,
        }
        rec["_pertile"] = pertile
        if verbose:
            print("RUNG %-22s median %.4f  gain %+.4f  CI [%+.4f, %+.4f]  wins %d/%d"
                  % (PRIMARY, rec["rung"]["median_spearman"], boot["median"],
                     boot["ci_lo"], boot["ci_hi"], wins, len(seeds)), flush=True)

    rec["comparator"] = {
        "atlas_median": inh["comparator"]["atlas_median"],
        "atlas_per_seed": ctx.committed_atlas,
        "source": inh["comparator"]["source"],
        "recomputed_in_run": False,
    }
    rec["margin"] = MARGIN
    rec["multiplicity"] = bar["multiplicity"]
    rec["not_a_blind_pre_registration"] = bar["NOT_A_BLIND_PRE_REGISTRATION"]
    rec["honest_limits"] = bar["honest_limits_fixed_in_advance"]
    rec["forecast_recorded_in_bar"] = bar["forecast"]
    rec["ram"] = {"free_at_start_bytes": rec["free_ram_bytes_at_start"],
                  "free_before_rung_bytes": free_now, "headroom_factor": RAM_HEADROOM}

    verdict, detail = verdict_from_artifact(rec, bar)
    rec["verdict"], rec["verdict_detail"] = verdict, detail
    rec["verdict_is_rederivable_from_this_artifact"] = True
    rec["wall_seconds"] = round(time.time() - t0, 1)
    return rec


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--block-dir", default=str(DEFAULT_BLOCK_DIR))
    ap.add_argument("--null-seeds", type=int, default=300)
    ap.add_argument("--n-boot", type=int, default=N_BOOT)
    args = ap.parse_args(argv)

    if not BAR_PATH.is_file():
        print("REFUSED: %s absent; freeze the bar before scoring." % BAR_PATH, flush=True)
        return 2
    bar = json.loads(BAR_PATH.read_text(encoding="utf-8"))

    try:
        rec = run(block_dir=Path(args.block_dir), bar=bar, null_extension=args.null_seeds,
                  n_boot=args.n_boot)
    except Exception as e:
        print("SCORING FAILED (%s): %s" % (type(e).__name__, e), flush=True)
        return 1

    pertile = rec.pop("_pertile", {})
    if pertile:
        PERTILE_OUT.parent.mkdir(parents=True, exist_ok=True)
        PERTILE_OUT.write_text(json.dumps(
            {"schema": "asgiven-arm-pertile-v1", "rung": PRIMARY,
             "why_not_in_wiki": ("a large decimal dump in wiki/ has previously polluted the "
                                 "contract-number audit's decoy pool"),
             "per_seed": pertile}, sort_keys=True), encoding="utf-8")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(rec, indent=2, sort_keys=True), encoding="utf-8")

    if rec["verdict"] not in OUTCOMES:
        print("INTERNAL: verdict %r is not in the bar's closed set" % rec["verdict"])
        return 1

    print()
    print("VERDICT  %s" % rec["verdict"])
    print("         %s" % rec["verdict_detail"])
    h = rec.get("not_a_blind_pre_registration", {})
    if h:
        print("HAZARD   the number was already visible as the positional run's INELIGIBLE "
              "secondary (%.4f); this is not a blind pre-registration"
              % h.get("already_observed_median", float("nan")))
    print("wall     %.1f s" % rec["wall_seconds"])
    print("-> %s" % OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
