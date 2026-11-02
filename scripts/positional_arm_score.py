"""Score the positional ladder against the frozen bar. GATES FIRST, so a broken run cannot read as
a finding.

Executes `wiki/positional_arm_acceptance_bar.json`'s check order exactly, in order, stopping at the
first validity failure:

    1. split_equality      -> INDETERMINATE_SPLIT_MISMATCH   (checked first of all)
    2. parity              -> INDETERMINATE_PARITY_FAIL
    3. shuffled_label_null -> INDETERMINATE_NULL_FAIL        (pre-committed: spend MORE SEEDS)
    4. co_primary_decision -> POSITIONAL_ADDS / _NO_GAIN / _AMBIGUOUS

**Why split equality is first.** Every other number is meaningless if the test sets differ: the arm
would not be scored on the incumbent's held-out peaks, so no comparison could be reported at all.

**Why parity is second, and what it buys.** The stored 10-bin means are collapsed back to a global
mean and scored through the SAME `ridge_spearman` on the SAME splits. If that does not reproduce the
committed atlas per-seed values, the forward pass has drifted and no later number is comparable to
0.5076. This is the check that makes the whole ladder commensurable rather than merely adjacent.

**The null statistic is INHERITED from the incumbent's implementation, not from its prose.** The
incumbent bar's text reads `|median gain over GC| < 0.02`, but its kernel computes
`ridge_spearman(primary[train], shuffled_train_labels, primary[test], real_test_labels)` -- the
shuffled-label model's OWN Spearman, with no GC term in it. That is the statistic that produced the
fragment arm's +0.0517 and the atlas arm's 0.0203; a literal `arm_shuffled - gc` would sit near
-0.31 and could never clear a 0.02 bar. So the implementation is inherited and the prose discrepancy
is recorded in the artifact rather than silently resolved either way.

**The null auto-extends on failure, because the bar pre-committed to that.** The incumbent's own
10-seed null failed at 0.0203 against the 0.02 bar on per-seed values spanning -0.0909..+0.1293, and
the median's standard error was ~0.024 -- LARGER than the bar, so ten seeds cannot decide it either
way. The response is more seeds of the SAME quantity, never a weaker control. Only the null extends:
the ladder stays on the ten committed splits, because those are the only ones byte-identical to the
incumbents'.

**Multiplicity is read out of the bar, never recomputed here.** Co-primaries are judged at the
two-comparison percentiles, secondaries at the all-rungs percentiles, and `classify` is imported
from the freeze script so the verdict cannot drift from the bar that defines it.

**A rung that fails its RAM pre-check is SKIPPED and recorded with its byte estimate**, never
silently dropped and never quietly retried at a smaller width. Free RAM is re-measured before every
rung, because it moves during a run.

Per-tile predictions for the primary rung land under `data/`, NEVER `wiki/`: a large decimal dump in
`wiki/` has previously polluted the contract-number audit's decoy pool and broken two of its guards.
"""

from __future__ import annotations

import argparse
import json
import time
from math import comb
from pathlib import Path

import numpy as np

from dna_decode.glm.positional_readout import flatten_bins, global_mean_from_bin_means, orient
from scripts.freeze_positional_arm_bar import (
    CO_PRIMARIES,
    MARGIN,
    OUTCOMES,
    SECONDARIES,
    classify,
)
from scripts.positional_arm_extract import (
    DEFAULT_BLOCK_DIR,
    as_bins,
    load_block,
    rebin,
    slice_block,
)
from scripts.positional_arm_probe import (
    design_matrix_bytes,
    fits_in_ram,
    free_ram_bytes,
    load_arm_cohort,
    parity_report,
)

REPO = Path(__file__).resolve().parents[1]
BAR_PATH = REPO / "wiki/positional_arm_acceptance_bar.json"
EMB_CACHE = REPO / "data/promoteratlas/tile_embeddings_native150_n44106.npy"
OUT = REPO / "wiki/positional_arm_result_2026-11-01.json"
#: per-tile predictions go under data/, never wiki/ -- see the module docstring
PERTILE_OUT = REPO / "data/promoteratlas/positional_arm_pertile_2026-11-01.json"

N_BOOT = 10_000
BOOT_SEED = 707

#: The committed atlas per-seed values are stored rounded to 4 dp, so the floor on any reproduction
#: check is 5e-5. 1e-3 leaves a 20x margin over that rounding while still catching any drift large
#: enough to move a reported median.
SPEARMAN_PARITY_TOL = 1e-3

#: Default size of the pooled null extension if the 10-seed null fails. The incumbent needed 270.
DEFAULT_NULL_EXTENSION = 300


class ScoreError(RuntimeError):
    """The run cannot proceed, as distinct from a verdict the bar legitimately returns."""


# ---------------------------------------------------------------------------------------------------
# splits
# ---------------------------------------------------------------------------------------------------
def leave_peak_out_indices(tiles, *, frac: float, seed: int):
    """`leave_peak_out` membership as ROW INDICES into the cohort order the block is stored in."""
    from dna_decode.glm import genomewide as gw

    _tr, te, n_held = gw.leave_peak_out(tiles, frac=frac, seed=seed)
    held = {id(t) for t in te}
    test_idx = [i for i, t in enumerate(tiles) if id(t) in held]
    train_idx = [i for i, t in enumerate(tiles) if id(t) not in held]
    return train_idx, test_idx, n_held


# ---------------------------------------------------------------------------------------------------
# rung construction -- every rung out of the SAME stored bytes
# ---------------------------------------------------------------------------------------------------
def build_rung(block, layout: dict, name: str, strands: list[str]) -> np.ndarray:
    """Materialise one rung's design matrix from the stored block.

    `rebin` and `orient` COMMUTE for bin means -- the mean within a group is order-invariant, so
    reversing ten bins then grouping into two gives the same two values as grouping then reversing.
    The cheaper order (rebin, then orient) is used, and a test pins the commutation.
    """
    if name == "global_mean":                      # the parity rung, not a ladder rung
        bm = as_bins(slice_block(block, layout, "bins10_mean"), layout)
        return np.ascontiguousarray(global_mean_from_bin_means(bm))
    if name == "autocorr":
        return np.ascontiguousarray(np.asarray(slice_block(block, layout, "autocorr"),
                                               dtype=np.float64))
    if name == "bins10_mean_asgiven":
        return np.ascontiguousarray(np.asarray(slice_block(block, layout, "bins10_mean"),
                                               dtype=np.float64))
    if name == "bins10_mean_oriented":
        bm = as_bins(slice_block(block, layout, "bins10_mean"), layout)
        return np.ascontiguousarray(flatten_bins(orient(bm, strands)))
    if name == "bins10_max_oriented":
        bx = as_bins(slice_block(block, layout, "bins10_max"), layout)
        return np.ascontiguousarray(flatten_bins(orient(bx, strands)))
    if name in ("bins2_mean_oriented", "bins5_mean_oriented"):
        n = 2 if name.startswith("bins2") else 5
        bm = as_bins(slice_block(block, layout, "bins10_mean"), layout)
        return np.ascontiguousarray(flatten_bins(orient(rebin(bm, n), strands)))
    raise ScoreError(
        "unknown rung %r; the bar declares co-primaries %r and secondaries %r"
        % (name, list(CO_PRIMARIES), list(SECONDARIES))
    )


# ---------------------------------------------------------------------------------------------------
# statistics
# ---------------------------------------------------------------------------------------------------
def median_of(values) -> float:
    v = sorted(float(x) for x in values)
    n = len(v)
    if not n:
        raise ScoreError("median of an empty sequence")
    return round(v[n // 2] if n % 2 else 0.5 * (v[n // 2 - 1] + v[n // 2]), 4)


def paired_bootstrap(gains, *, percentiles, n_boot: int = N_BOOT, seed: int = BOOT_SEED) -> dict:
    """Bootstrap the PER-SEED paired gains -- the pairing the frozen bar names."""
    g = np.asarray(list(gains), dtype=float)
    if g.size == 0:
        raise ScoreError("paired_bootstrap got no gains")
    if not np.all(np.isfinite(g)):
        raise ScoreError("paired gains hold non-finite values; refusing to bootstrap them")
    lo_p, hi_p = percentiles
    rng = np.random.default_rng(seed)
    draws = rng.choice(g, size=(n_boot, g.size), replace=True).mean(axis=1)
    return {
        "mean": round(float(g.mean()), 4),
        "median": median_of(g),
        "ci_lo": round(float(np.percentile(draws, lo_p)), 4),
        "ci_hi": round(float(np.percentile(draws, hi_p)), 4),
        "ci_percentiles": [lo_p, hi_p],
        "n_seeds": int(g.size),
        "n_boot": n_boot,
        "seed": seed,
    }


def order_statistic_median_ci(values):
    """Exact distribution-free interval for the MEDIAN, with its ACHIEVED coverage.

    Deliberately not a normal approximation around the MEAN -- a different estimand, and one this
    project has misread before. Copied in form from the incumbent's null extension.
    """
    v = sorted(float(x) for x in values)
    n = len(v)
    if n < 2:
        return None, None, None
    lo = hi = cov = None
    for k in range(n // 2 + 1):
        c = 1.0 - 2.0 * sum(comb(n, i) for i in range(k)) / (2.0 ** n)
        if c >= 0.95:
            lo, hi, cov = v[k], v[n - 1 - k], c
    return lo, hi, cov


# ---------------------------------------------------------------------------------------------------
# the verdict, re-derivable from the artifact alone
# ---------------------------------------------------------------------------------------------------
def verdict_from_artifact(rec: dict, bar: dict) -> tuple[str, str]:
    """Recompute the verdict from the emitted numbers plus the frozen anchors, and nothing else.

    This exists so a reader can check the verdict without re-running anything, and so a test can
    assert the emitted verdict IS a pure function of the measurements rather than of control flow.
    """
    gates = rec.get("gates", {})
    margin = bar["inherited_not_chosen"]["margin"]["value"]
    for gate, outcome in (("split_equality", "INDETERMINATE_SPLIT_MISMATCH"),
                          ("parity", "INDETERMINATE_PARITY_FAIL"),
                          ("shuffled_label_null", "INDETERMINATE_NULL_FAIL")):
        g = gates.get(gate)
        if g is None:
            return (outcome, "gate %r was never evaluated, so the run stopped at or before it"
                    % gate)
        if not g.get("passed"):
            return (outcome, g.get("detail", "gate %r failed" % gate))
    cis = {}
    for name in CO_PRIMARIES:
        r = rec.get("rungs", {}).get(name)
        if r is None or r.get("skipped"):
            return ("POSITIONAL_AMBIGUOUS",
                    "co-primary %r was not scored (skipped=%r), so no difference is detectable at "
                    "the corrected level" % (name, (r or {}).get("skipped")))
        cis[name] = (r["gain_bootstrap"]["ci_lo"], r["gain_bootstrap"]["ci_hi"])
    return classify(cis, margin=margin)


# ---------------------------------------------------------------------------------------------------
# the run
# ---------------------------------------------------------------------------------------------------
def run(*, block_dir: Path, bar: dict, null_extension: int, n_boot: int,
        scorer=None, cohort=None, block=None, emb_cache: Path | None = None,
        verbose: bool = True) -> dict:
    from dna_decode.glm import genomewide as gw
    from dna_decode.glm.atlas_embed import ridge_spearman

    score = scorer or ridge_spearman
    t0 = time.time()

    inh = bar["inherited_not_chosen"]
    frac = inh["split"]["frac"]
    seeds = list(range(inh["split"]["n_seeds"]))
    committed_atlas = inh["comparator"]["atlas_per_seed"]
    committed_n_test = inh["split"]["n_test_per_seed"]
    parity_tol = bar["measured_and_inherited_from_the_probe"]["parity_tolerance"]
    co_pct = bar["multiplicity"]["co_primary_ci_percentiles"]
    sec_pct = bar["multiplicity"]["secondary_ci_percentiles"]
    headroom = bar["measured_and_inherited_from_the_probe"]["ram_headroom_factor"]

    tiles, st = cohort if cohort is not None else load_arm_cohort()
    if block is None:
        arr, side, layout = load_block(block_dir)
    else:
        arr, side, layout = block
    if side["geometry"]["n_rows"] != len(tiles):
        raise ScoreError(
            "the block holds %d rows and the cohort is %d tiles. The block is indexed positionally "
            "against the cohort, so this pairing is impossible."
            % (side["geometry"]["n_rows"], len(tiles))
        )

    strands = [t.strand for t in tiles]
    rec: dict = {
        "schema": "positional-arm-result-v1",
        "date": "2026-11-01",
        "bar": "wiki/positional_arm_acceptance_bar.json",
        "arm": bar["arm"],
        "criterion": bar["criterion"],
        "block": {"path": side["block"], "sidecar_status": side["status"],
                  "checkpoint_sha256": side.get("checkpoint_sha256"),
                  "geometry": side["geometry"]},
        "substrate": {"n_tiles": len(tiles), "n_peaks": len({t.peak for t in tiles}),
                      "path": st.path},
        "free_ram_bytes_at_start": free_ram_bytes(),
        "gates": {},
        "rungs": {},
    }

    # ---- splits, once, reused byte-identically by every rung -----------------------------------
    splits = {}
    for s in seeds:
        tr, te, n_held = leave_peak_out_indices(tiles, frac=frac, seed=s)
        splits[s] = (tr, te, n_held)

    # ---- GATE 1: split equality (first of all) --------------------------------------------------
    per_seed_split = {}
    mismatched = []
    for s in seeds:
        tr, te, n_held = splits[s]
        want = int(committed_n_test[str(s)])
        ok = len(te) == want
        per_seed_split[str(s)] = {"n_test": len(te), "n_test_committed": want,
                                  "n_train": len(tr), "held_out_peaks": n_held, "match": ok}
        if not ok:
            mismatched.append(s)
    rec["gates"]["split_equality"] = {
        "passed": not mismatched,
        "per_seed": per_seed_split,
        "mismatched_seeds": mismatched,
        "detail": ("all %d recomputed splits match the committed n_test_per_seed" % len(seeds))
        if not mismatched else
        ("seeds %r differ from the committed splits, so the arm is not being scored on the "
         "incumbent's held-out peaks and no comparison may be reported" % mismatched),
    }
    if verbose:
        print("GATE split_equality  %s (%d/%d seeds match)"
              % ("PASS" if not mismatched else "FAIL",
                 len(seeds) - len(mismatched), len(seeds)), flush=True)
    if mismatched:
        rec["verdict"], rec["verdict_detail"] = (
            "INDETERMINATE_SPLIT_MISMATCH", rec["gates"]["split_equality"]["detail"])
        rec["wall_seconds"] = round(time.time() - t0, 1)
        return rec

    # ---- GATE 2: parity ------------------------------------------------------------------------
    gm = build_rung(arr, layout, "global_mean", strands)
    emb_par = None
    cache_path = EMB_CACHE if emb_cache is None else Path(emb_cache)
    if cache_path.is_file():
        committed_emb = np.asarray(np.load(cache_path, mmap_mode="r"), dtype=np.float64)
        if committed_emb.shape == gm.shape:
            emb_par = parity_report(gm, committed_emb, atol=parity_tol)
        elif verbose:
            # a shape mismatch here is not a parity failure -- it means this run is not scoring the
            # cohort the cache was built for, which the block/cohort check above already enforces
            print("   (embedding cache is %r, derived is %r -- embedding-level parity skipped, "
                  "the Spearman-level check still runs)" % (committed_emb.shape, gm.shape),
                  flush=True)
    per_seed_parity, parity_bad = {}, []
    for s in seeds:
        tr, te, _ = splits[s]
        ytr = [tiles[i].expression for i in tr]
        yte = [tiles[i].expression for i in te]
        sp, _pred = score(gm[tr], ytr, gm[te], yte)
        want = float(committed_atlas[str(s)])
        dev = abs(round(float(sp), 4) - want)
        per_seed_parity[str(s)] = {"derived": round(float(sp), 4), "committed": want,
                                   "abs_deviation": round(dev, 6)}
        if dev > SPEARMAN_PARITY_TOL:
            parity_bad.append(s)
    parity_ok = (not parity_bad) and (emb_par is None or emb_par["passed"])
    rec["gates"]["parity"] = {
        "passed": bool(parity_ok),
        "embedding_level": emb_par,
        "embedding_tolerance": parity_tol,
        "spearman_level": {"per_seed": per_seed_parity, "tolerance": SPEARMAN_PARITY_TOL,
                           "deviant_seeds": parity_bad},
        "derived_median": median_of(v["derived"] for v in per_seed_parity.values()),
        "committed_median": inh["comparator"]["atlas_median"],
        "detail": ("the global mean derived from the stored bin means reproduces the committed "
                   "embedding and its per-seed Spearmans, so every rung below is comparable to "
                   "%.4f" % inh["comparator"]["atlas_median"]) if parity_ok else
                  ("derived values deviate from the committed atlas per-seed table on seeds %r "
                   "(embedding parity passed=%r); the forward pass has drifted and no number from "
                   "this arm would be comparable to %.4f"
                   % (parity_bad, None if emb_par is None else emb_par["passed"],
                      inh["comparator"]["atlas_median"])),
    }
    if verbose:
        print("GATE parity          %s (emb max|dev| %s, spearman median %.4f vs committed %.4f)"
              % ("PASS" if parity_ok else "FAIL",
                 "n/a" if emb_par is None else "%.3e" % emb_par["max_abs_deviation"],
                 rec["gates"]["parity"]["derived_median"], inh["comparator"]["atlas_median"]),
              flush=True)
    if not parity_ok:
        rec["verdict"], rec["verdict_detail"] = (
            "INDETERMINATE_PARITY_FAIL", rec["gates"]["parity"]["detail"])
        rec["wall_seconds"] = round(time.time() - t0, 1)
        return rec

    # ---- GATE 3: shuffled-label null on the primary rung ----------------------------------------
    primary = CO_PRIMARIES[0]
    prim = build_rung(arr, layout, primary, strands)
    null_by_seed = {}
    for s in seeds:
        tr, te, _ = splits[s]
        ytr = [tiles[i].expression for i in tr]
        yte = [tiles[i].expression for i in te]
        sp, _ = score(prim[tr], ytr, prim[te], yte, shuffle_seed=s)
        null_by_seed[str(s)] = round(float(sp), 4)
    null_median = median_of(null_by_seed.values())
    null_ext, pooled_median, n_pooled = {}, null_median, len(null_by_seed)
    extended = False
    if abs(null_median) >= MARGIN and null_extension > 0:
        extended = True
        if verbose:
            print("   null at %d seeds is %+.4f against the |median|<%.2f bar -- EXTENDING to %d "
                  "seeds, which is the bar's pre-committed response (more seeds of the SAME "
                  "quantity, never a weaker control)"
                  % (len(null_by_seed), null_median, MARGIN, null_extension), flush=True)
        for s in range(null_extension):
            tr, te, _ = leave_peak_out_indices(tiles, frac=frac, seed=s)
            ytr = [tiles[i].expression for i in tr]
            yte = [tiles[i].expression for i in te]
            sp, _ = score(prim[tr], ytr, prim[te], yte, shuffle_seed=s)
            null_ext[str(s)] = round(float(sp), 4)
            if verbose and s % 25 == 0:
                print("      null seed %d/%d -> %+.4f (%.0fs)"
                      % (s, null_extension, sp, time.time() - t0), flush=True)
        pooled_median = median_of(null_ext.values())
        n_pooled = len(null_ext)
    lo, hi, cov = order_statistic_median_ci((null_ext or null_by_seed).values())
    null_ok = abs(pooled_median) < MARGIN
    rec["gates"]["shuffled_label_null"] = {
        "passed": bool(null_ok),
        "statistic": ("the shuffled-label model's OWN Spearman on the primary rung: "
                      "ridge_spearman(primary[train], SHUFFLED train labels, primary[test], REAL "
                      "test labels)"),
        "statistic_provenance": (
            "INHERITED from the incumbent arm's IMPLEMENTATION, not from its prose. That bar's text "
            "reads '|median gain over GC| < 0.02' but its kernel computes no GC term at all. The "
            "implementation is the statistic that produced the fragment arm's +0.0517 and the atlas "
            "arm's 0.0203; a literal arm_shuffled-minus-GC would sit near -0.31 and could never "
            "clear a 0.02 bar. Discrepancy recorded rather than silently resolved."
        ),
        "rung": primary,
        "bar": "|median| < %.2f" % MARGIN,
        "per_seed_10": null_by_seed,
        "median_10": null_median,
        "extended": extended,
        "n_extension_seeds": len(null_ext),
        "per_seed_extension": null_ext,
        "pooled_median": pooled_median,
        "n_pooled": n_pooled,
        "median_ci95": [lo, hi],
        "median_ci_achieved_coverage": cov,
        "ci_method": ("exact distribution-free order-statistic interval for the MEDIAN, with its "
                      "achieved coverage -- not a normal approximation around the MEAN"),
        "detail": ("the shuffled-label control is %+.4f over %d seeds, inside the |median|<%.2f bar"
                   % (pooled_median, n_pooled, MARGIN)) if null_ok else
                  ("the shuffled-label control is %+.4f over %d seeds, outside the |median|<%.2f "
                   "bar; nothing may be graded and the pre-committed response is MORE SEEDS, never "
                   "a weaker control" % (pooled_median, n_pooled, MARGIN)),
    }
    if verbose:
        print("GATE null            %s (median %+.4f over %d seeds, bar |median|<%.2f)"
              % ("PASS" if null_ok else "FAIL", pooled_median, n_pooled, MARGIN), flush=True)
    if not null_ok:
        rec["verdict"], rec["verdict_detail"] = (
            "INDETERMINATE_NULL_FAIL", rec["gates"]["shuffled_label_null"]["detail"])
        rec["wall_seconds"] = round(time.time() - t0, 1)
        return rec

    # ---- the ladder -----------------------------------------------------------------------------
    pertile = {}
    min_free = rec["free_ram_bytes_at_start"]
    for name in list(CO_PRIMARIES) + list(SECONDARIES):
        is_primary = name in CO_PRIMARIES
        pct = co_pct if is_primary else sec_pct
        X = prim if name == primary else build_rung(arr, layout, name, strands)
        need = design_matrix_bytes(len(tiles), X.shape[1])
        free_now = free_ram_bytes()
        if free_now is not None:
            min_free = min(min_free, free_now) if min_free is not None else free_now
        ok = fits_in_ram(need, free_now, headroom=headroom)
        if ok is not True:
            rec["rungs"][name] = {
                "skipped": True,
                "role": "co-primary" if is_primary else "secondary",
                "n_features": int(X.shape[1]),
                "peak_design_bytes": need,
                "free_ram_bytes": free_now,
                "fits_in_ram": ok,
                "why": ("the RAM pre-check did not return True (None means UNMEASURED, which is "
                        "not the same as affordable). Recorded as skipped with its byte estimate "
                        "rather than silently dropped or quietly retried at a smaller width."),
            }
            if verbose:
                print("RUNG %-22s SKIPPED (needs %.3f GB, free %s)"
                      % (name, need / 1e9, "unknown" if free_now is None else
                         "%.2f GB" % (free_now / 1e9)), flush=True)
            continue

        per_seed, gains, wins = {}, [], 0
        for s in seeds:
            tr, te, _ = splits[s]
            ytr = [tiles[i].expression for i in tr]
            yte = [tiles[i].expression for i in te]
            sp, pred = score(X[tr], ytr, X[te], yte)
            sp = round(float(sp), 4)
            want = float(committed_atlas[str(s)])
            per_seed[str(s)] = {"spearman": sp, "atlas_committed": want,
                                "gain": round(sp - want, 4)}
            gains.append(sp - want)
            wins += int(sp > want)
            if name == primary:
                pertile[str(s)] = {"tile_indices": te, "pred": [round(float(v), 6) for v in pred]}
        boot = paired_bootstrap(gains, percentiles=pct, n_boot=n_boot)
        rec["rungs"][name] = {
            "skipped": False,
            "role": "co-primary" if is_primary else "secondary",
            "eligible_as_the_arms_result": is_primary,
            "n_features": int(X.shape[1]),
            "peak_design_bytes": need,
            "median_spearman": median_of(v["spearman"] for v in per_seed.values()),
            "per_seed": per_seed,
            "seed_wins_over_atlas": wins,
            "n_seeds": len(seeds),
            "gain_bootstrap": boot,
        }
        if verbose:
            print("RUNG %-22s median %.4f  gain %+.4f  CI [%+.4f, %+.4f]  wins %d/%d  (%s)"
                  % (name, rec["rungs"][name]["median_spearman"], boot["median"],
                     boot["ci_lo"], boot["ci_hi"], wins, len(seeds),
                     "co-primary" if is_primary else "secondary"), flush=True)
        del X

    rec["ram"] = {
        "free_at_start_bytes": rec["free_ram_bytes_at_start"],
        "min_free_observed_bytes": min_free,
        "headroom_factor": headroom,
        "note": ("free RAM is re-measured before every rung because it moves during a run; the "
                 "minimum observed is what a reader should size against"),
    }
    rec["comparator"] = {
        "atlas_median": inh["comparator"]["atlas_median"],
        "atlas_per_seed": committed_atlas,
        "source": inh["comparator"]["source"],
        "recomputed_in_run": False,
        "why": "the incumbent is never re-run by this arm, so this arm cannot move the number it "
               "is graded against",
    }
    rec["margin"] = MARGIN
    rec["multiplicity"] = bar["multiplicity"]
    rec["honest_limits"] = bar["honest_limits_fixed_in_advance"]
    rec["scoping_corrections"] = bar["scoping_corrections_fixed_in_advance"]
    rec["forecast_recorded_in_bar"] = bar["forecast"]

    verdict, detail = verdict_from_artifact(rec, bar)
    rec["verdict"], rec["verdict_detail"] = verdict, detail
    rec["verdict_is_rederivable_from_this_artifact"] = True
    rec["wall_seconds"] = round(time.time() - t0, 1)
    rec["_pertile"] = pertile
    return rec


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--block-dir", default=str(DEFAULT_BLOCK_DIR))
    ap.add_argument("--null-seeds", type=int, default=DEFAULT_NULL_EXTENSION,
                    help="pooled null extension size used ONLY if the 10-seed null fails; 0 "
                         "disables the extension and the run then reports INDETERMINATE_NULL_FAIL")
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
            {"schema": "positional-arm-pertile-v1", "rung": CO_PRIMARIES[0],
             "why_not_in_wiki": ("a large decimal dump in wiki/ has previously polluted the "
                                 "contract-number audit's decoy pool and broken two of its guards"),
             "per_seed": pertile}, sort_keys=True), encoding="utf-8")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(rec, indent=2, sort_keys=True), encoding="utf-8")

    if rec["verdict"] not in OUTCOMES:
        print("INTERNAL: verdict %r is not in the bar's closed set %r" % (rec["verdict"], OUTCOMES))
        return 1

    print()
    print("VERDICT  %s" % rec["verdict"])
    print("         %s" % rec["verdict_detail"])
    print("wall     %.1f s (%.1f min)" % (rec["wall_seconds"], rec["wall_seconds"] / 60.0))
    print("-> %s" % OUT)
    if pertile:
        print("-> %s (per-tile predictions, NOT in wiki/)" % PERTILE_OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
