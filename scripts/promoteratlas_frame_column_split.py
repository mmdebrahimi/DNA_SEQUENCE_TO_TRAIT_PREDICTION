"""Which HALF of the frame table carries the signal -- the window columns, or the CDS columns?

This is a PREREQUISITE measurement, and it exists because the proposed successor control may be
invalid by construction. The frozen bar's matched-informativeness control failed its own premise
(the non-positional comparators are WEAKER standalone than the frame), and the proposed replacement
was a SHUFFLED-WINDOW control: identical feature construction with the window set relocated along
the genome.

But `frame_features(located, mask, cds_rows, *, tile_length)` takes the window information ONLY as
`mask`. Of its five columns, three are derived from that mask and TWO are derived from `cds_rows`:

    in_window_any                      <- mask      (window-derived)
    in_window_full                     <- mask      (window-derived)
    overlap_fraction                   <- mask      (window-derived)
    signed_distance_to_next_cds_start  <- cds_rows  (CDS-derived, INVARIANT under a window shuffle)
    inside_cds                         <- cds_rows  (CDS-derived, INVARIANT under a window shuffle)

So a window shuffle leaves two of five columns BYTE-IDENTICAL between treatment and control. Worse,
corpus windows are DEFINED as the 200 bp upstream of annotated gene starts, which makes
`signed_distance_to_next_cds_start` a continuous restatement of `in_window_any` -- meaning the
proposed control would retain the strongest positional column in the treatment table.

The question this run answers is therefore NOT whether the frame is informative (already measured,
0.3170 standalone). It is: how much of that is reachable by the CDS columns ALONE? If the CDS half
carries most of it, a window-shuffle control is a near-copy of the treatment, its absorption lands
close to the real absorption, `attributable_absorption` collapses toward zero, and the grader fires
INDETERMINATE_ABSORPTION_NOT_ATTRIBUTABLE for a reason that is an artifact of shared columns rather
than a fact about position. That is a SPURIOUS indeterminate, and it is worth one cheap run to find
out BEFORE building the control.

It is deliberately NOT a verdict on the frame arm. Nothing here re-grades anything, no frozen file
is read for a threshold, and no acceptance bar is consulted -- the output is four numbers and their
seed spreads.

**The anchor runs FIRST**, the same discipline as the conditioning run: the cached embeddings must
reproduce the committed `arms.native150` before any new number is read.

**No k-mer baselines are built.** The conditioning run's slowest step is the 4-mer feature matrix,
and nothing here needs it, so it is skipped rather than computed and discarded.

Run:  uv run python -m scripts.promoteratlas_frame_column_split
"""

from __future__ import annotations

import dataclasses
import json
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
ARM = REPO / "wiki/promoteratlas_tile_arm_2026-11-01.json"
CACHE = REPO / "data/mpra"
EMB_DIR = REPO / "data/promoteratlas"
FNA = REPO / "data/cache/refseq/GCF_000005845.2/genome.fna"
GFF = REPO / "data/cache/refseq/GCF_000005845.2/annotations.gff3"
LOCATE_CACHE = EMB_DIR / "frame_locate_cache.json"
OUT = REPO / "wiki/promoteratlas_frame_column_split_2026-11-05.json"

SCHEMA = "promoteratlas-frame-column-split-v1"
N_SEEDS = 10
FRAC = 0.25
TILE_LEN = 150
ANCHOR_TOL = 0.02

# Names, not indices. Asserted against the live FEATURE_NAMES tuple at run time rather than
# trusted -- a reordering there would silently invert this entire measurement.
WINDOW_COLS = ("in_window_any", "in_window_full", "overlap_fraction")
CDS_COLS = ("signed_distance_to_next_cds_start", "inside_cds")


def _spread(xs) -> dict:
    a = np.asarray([float(x) for x in xs], dtype=float)
    a = a[~np.isnan(a)]
    if a.size == 0:
        return {"n": 0, "note": "no usable seeds"}
    return {
        "n": int(a.size),
        "median": round(float(np.median(a)), 4),
        "mean": round(float(a.mean()), 4),
        "min": round(float(a.min()), 4),
        "max": round(float(a.max()), 4),
        "sd": round(float(a.std(ddof=1)) if a.size > 1 else 0.0, 4),
    }


def run(verbose: bool = True) -> dict:
    from dna_decode.data.annotations import parse_gff3
    from dna_decode.glm import genomewide as gw
    from dna_decode.glm.atlas_embed import ridge_spearman
    from dna_decode.glm.corpus_window import (
        cds_rows_from_annotation_table,
        corpus_windows,
        window_coverage,
    )
    from dna_decode.glm.frame_features import FEATURE_NAMES, frame_features
    from dna_decode.glm.tile_locate import (
        load_locate_cache,
        locate_report,
        locate_sequences,
        read_fasta_single,
    )

    t0 = time.time()
    arm = json.loads(ARM.read_text(encoding="utf-8"))
    committed_atlas = float(arm["arms"]["native150"]["median"])

    names = list(FEATURE_NAMES)
    if set(WINDOW_COLS) | set(CDS_COLS) != set(names):
        raise SystemExit(
            "FEATURE_NAMES is %r; this script's window/CDS partition names %r and is stale"
            % (names, sorted(set(WINDOW_COLS) | set(CDS_COLS)))
        )
    wj = [names.index(c) for c in WINDOW_COLS]
    cj = [names.index(c) for c in CDS_COLS]

    rows, _ = gw.load_peak_tiles(CACHE, allow_download=False)
    tiles = [t for t in rows if t.category == "tile" and len(t.seq) == TILE_LEN]
    n = len(tiles)
    Xa = np.load(EMB_DIR / ("tile_embeddings_native150_n%d.npy" % n))
    if Xa.shape != (n, 128):
        raise SystemExit("embedding cache is %r, expected (%d, 128)" % (Xa.shape, n))

    _contig, genome = read_fasta_single(FNA)
    cds = cds_rows_from_annotation_table(parse_gff3(GFF))
    mask = window_coverage(corpus_windows(cds, len(genome)), len(genome))

    if LOCATE_CACHE.is_file():
        located = load_locate_cache(LOCATE_CACHE)
        if {t.seq for t in tiles} - set(located):
            located = locate_sequences([t.seq for t in tiles], genome)
    else:
        located = locate_sequences([t.seq for t in tiles], genome)
    lrep = locate_report(located, tile_length=TILE_LEN, genome_length=len(genome))
    table = frame_features(located, mask, cds, tile_length=TILE_LEN)
    by_seq = {s: i for i, s in enumerate(table.seqs)}

    y = np.asarray([t.expression for t in tiles], dtype=float)
    gc = np.asarray(gw.gc_features([t.seq for t in tiles]), dtype=float)

    F = np.full((n, len(names)), np.nan, dtype=float)
    has_frame = np.zeros(n, dtype=bool)
    for i, t in enumerate(tiles):
        j = by_seq.get(t.seq)
        if j is not None:
            F[i] = table.matrix[j]
            has_frame[i] = True

    W = F[:, wj]
    C = F[:, cj]
    gc_frame = np.hstack([gc, F])
    gc_window = np.hstack([gc, W])
    gc_cds = np.hstack([gc, C])

    if verbose:
        print("tiles %d, located %d, window cols %s, cds cols %s (%.0fs)"
              % (n, int(has_frame.sum()), wj, cj, time.time() - t0), flush=True)

    keys = ("anchor_atlas_full", "atlas", "gc",
            "frame_alone", "window_alone", "cds_alone",
            "gc_frame", "gc_window", "gc_cds")
    per_seed: dict[str, list] = {k: [] for k in keys}

    for seed in range(N_SEEDS):
        train, test, _ = gw.leave_peak_out(tiles, frac=FRAC, seed=seed)
        pos = {id(t): i for i, t in enumerate(tiles)}
        tr = np.asarray([pos[id(t)] for t in train])
        te = np.asarray([pos[id(t)] for t in test])

        # the anchor, on the FULL cohort, before any new number is read
        a_full, _ = ridge_spearman(Xa[tr], y[tr], Xa[te], y[te])
        per_seed["anchor_atlas_full"].append(round(float(a_full), 4))

        trL = tr[has_frame[tr]]
        teL = te[has_frame[te]]
        for key, X in (("atlas", Xa), ("gc", gc),
                       ("frame_alone", F), ("window_alone", W), ("cds_alone", C),
                       ("gc_frame", gc_frame), ("gc_window", gc_window), ("gc_cds", gc_cds)):
            v, _ = ridge_spearman(X[trL], y[trL], X[teL], y[teL])
            per_seed[key].append(round(float(v), 4))

        if verbose:
            print("  seed %d  frame %.4f | window %.4f | cds %.4f  (%.0fs)"
                  % (seed, per_seed["frame_alone"][-1], per_seed["window_alone"][-1],
                     per_seed["cds_alone"][-1], time.time() - t0), flush=True)

    med = {k: float(np.median(np.asarray(v, dtype=float))) for k, v in per_seed.items()}
    anchor_delta = abs(med["anchor_atlas_full"] - committed_atlas)
    if anchor_delta > ANCHOR_TOL:
        raise SystemExit(
            "ANCHOR FAILED: full-cohort atlas median %.4f vs committed %.4f (delta %.4f > %.4f)"
            % (med["anchor_atlas_full"], committed_atlas, anchor_delta, ANCHOR_TOL)
        )

    # The decisive ratio. `cds_alone / frame_alone` is the fraction of the frame table's standalone
    # signal that a WINDOW SHUFFLE cannot touch, because those two columns are byte-identical under
    # it. A high value means the proposed control shares the treatment's signal and is not a null.
    cds_share = med["cds_alone"] / med["frame_alone"] if med["frame_alone"] else float("nan")
    window_share = med["window_alone"] / med["frame_alone"] if med["frame_alone"] else float("nan")
    absorption_frame = med["gc_frame"] - med["gc"]
    absorption_window = med["gc_window"] - med["gc"]
    absorption_cds = med["gc_cds"] - med["gc"]
    cds_absorption_share = (
        absorption_cds / absorption_frame if absorption_frame else float("nan")
    )

    # A ratio of medians cannot say whether the window columns add anything UNIQUE -- for that the
    # comparison must be PAIRED within seed, because the seed-to-seed spread (sd ~0.02) is larger
    # than the gap being tested. n_seeds_cds_at_least_frame near half the seeds means the two are
    # indistinguishable, i.e. the three window columns are redundant given the two CDS columns.
    fa = np.asarray(per_seed["frame_alone"], dtype=float)
    ca = np.asarray(per_seed["cds_alone"], dtype=float)
    wa = np.asarray(per_seed["window_alone"], dtype=float)
    paired_cds_minus_frame = ca - fa
    paired_window_minus_frame = wa - fa
    paired = {
        "n_seeds": int(fa.size),
        "cds_minus_frame_median": round(float(np.median(paired_cds_minus_frame)), 4),
        "cds_minus_frame_mean": round(float(paired_cds_minus_frame.mean()), 4),
        "n_seeds_cds_at_least_frame": int((ca >= fa).sum()),
        "window_minus_frame_median": round(float(np.median(paired_window_minus_frame)), 4),
        "n_seeds_window_at_least_frame": int((wa >= fa).sum()),
        "note": (
            "Paired within seed. If cds_minus_frame straddles zero with roughly half the seeds on "
            "each side, the two CDS columns ALONE match the full five-column table, so the three "
            "window columns carry no unique standalone signal and a window shuffle cannot remove "
            "the frame table's predictive content."
        ),
    }

    if cds_share >= 0.75:
        verdict = "WINDOW_SHUFFLE_INVALID_CDS_COLUMNS_CARRY_THE_SIGNAL"
    elif cds_share >= 0.40:
        verdict = "WINDOW_SHUFFLE_CONFOUNDED_CDS_COLUMNS_CARRY_A_SUBSTANTIAL_SHARE"
    elif cds_share >= 0.10:
        verdict = "WINDOW_SHUFFLE_WEAKENED_CDS_COLUMNS_CARRY_A_MINORITY_SHARE"
    else:
        verdict = "WINDOW_SHUFFLE_VIABLE_ON_THIS_AXIS"

    rec = {
        "schema": SCHEMA,
        "date": "2026-11-05",
        "question": (
            "How much of the frame table's standalone signal is carried by the two CDS-derived "
            "columns, which a WINDOW SHUFFLE leaves byte-identical between treatment and control?"
        ),
        "why_this_runs_before_the_control_is_built": (
            "The proposed shuffled-window control relocates windows, which changes only the three "
            "mask-derived columns. If the two cds_rows-derived columns carry the signal, the "
            "control is a near-copy of the treatment, attributable_absorption collapses toward "
            "zero, and the grader fires INDETERMINATE_ABSORPTION_NOT_ATTRIBUTABLE for a reason "
            "that is an artifact of shared columns rather than a fact about position."
        ),
        "anchor": {
            "committed_native150_median": round(committed_atlas, 4),
            "reproduced_full_cohort_median": round(med["anchor_atlas_full"], 4),
            "abs_delta": round(anchor_delta, 4),
            "tolerance": ANCHOR_TOL,
            "ok": True,
        },
        "locate": dataclasses.asdict(lrep) if dataclasses.is_dataclass(lrep) else lrep,
        "column_partition": {
            "window_derived": list(WINDOW_COLS),
            "cds_derived": list(CDS_COLS),
            "window_indices": wj,
            "cds_indices": cj,
            "asserted_against_live_feature_names": names,
            "note": (
                "frame_features() receives window information ONLY as the union mask. The "
                "cds_derived columns are computed from cds_rows and are INVARIANT under any "
                "relocation of the window set."
            ),
        },
        "standalone_medians": {
            "frame_alone_all_five": round(med["frame_alone"], 4),
            "window_alone_three_cols": round(med["window_alone"], 4),
            "cds_alone_two_cols": round(med["cds_alone"], 4),
        },
        "handed_to_gc_medians": {
            "gc": round(med["gc"], 4),
            "gc_plus_frame": round(med["gc_frame"], 4),
            "gc_plus_window_cols": round(med["gc_window"], 4),
            "gc_plus_cds_cols": round(med["gc_cds"], 4),
        },
        "absorption": {
            "frame": round(absorption_frame, 4),
            "window_cols_only": round(absorption_window, 4),
            "cds_cols_only": round(absorption_cds, 4),
            "cds_share_of_frame_absorption": round(cds_absorption_share, 4),
            "note": (
                "Absorption here is the GAIN a baseline receives from the columns (gc+X minus gc), "
                "NOT the margin-shrinkage figure the frozen bar grades. It is reported because a "
                "shuffled-window control's absorption is what the grader would compare against."
            ),
        },
        "decisive_ratio": {
            "cds_share_of_standalone": round(cds_share, 4),
            "window_share_of_standalone": round(window_share, 4),
            "note": (
                "Shares do not sum to 1: the two halves overlap, because corpus windows are "
                "defined as the 200 bp upstream of annotated gene starts, so "
                "signed_distance_to_next_cds_start is a continuous restatement of in_window_any."
            ),
        },
        "paired_within_seed": paired,
        "verdict": verdict,
        "per_seed": per_seed,
        "spreads": {k: _spread(v) for k, v in per_seed.items()},
        "honest_limits": [
            "This measures the COLUMN partition only. It does not measure a shuffled control; it "
            "measures what such a control could not change.",
            "It does NOT re-grade the frame arm and consults no acceptance bar. The frame verdict "
            "FRAME_PARTIAL and its resolved=false status are untouched.",
            "The share is a ratio of medians over 10 seeds on 25 percent peak-disjoint held-out "
            "splits; the two halves overlap, so the shares are not a variance decomposition.",
            "A viable verdict on this axis would NOT make a window shuffle matched on gene "
            "density. Independent per-window relocation still breaks density; a circular shift of "
            "the whole annotation preserves it.",
        ],
        "seconds": round(time.time() - t0, 1),
    }
    OUT.write_text(json.dumps(rec, indent=2) + "\n", encoding="utf-8")
    return rec


def main() -> int:
    rec = run()
    print("")
    print("ANCHOR   %.4f vs committed %.4f (delta %.4f) OK"
          % (rec["anchor"]["reproduced_full_cohort_median"],
             rec["anchor"]["committed_native150_median"], rec["anchor"]["abs_delta"]))
    s = rec["standalone_medians"]
    print("STANDALONE  frame(5) %.4f | window(3) %.4f | cds(2) %.4f"
          % (s["frame_alone_all_five"], s["window_alone_three_cols"], s["cds_alone_two_cols"]))
    g = rec["handed_to_gc_medians"]
    print("TO GC       gc %.4f | +frame %.4f | +window %.4f | +cds %.4f"
          % (g["gc"], g["gc_plus_frame"], g["gc_plus_window_cols"], g["gc_plus_cds_cols"]))
    a = rec["absorption"]
    print("ABSORB      frame %+.4f | window %+.4f | cds %+.4f (cds share %.4f)"
          % (a["frame"], a["window_cols_only"], a["cds_cols_only"],
             a["cds_share_of_frame_absorption"]))
    d = rec["decisive_ratio"]
    print("SHARE       cds %.4f of standalone | window %.4f"
          % (d["cds_share_of_standalone"], d["window_share_of_standalone"]))
    pr = rec["paired_within_seed"]
    print("PAIRED      cds-frame median %+.4f (cds >= frame on %d of %d seeds) | "
          "window-frame %+.4f (%d of %d)"
          % (pr["cds_minus_frame_median"], pr["n_seeds_cds_at_least_frame"], pr["n_seeds"],
             pr["window_minus_frame_median"], pr["n_seeds_window_at_least_frame"], pr["n_seeds"]))
    print("VERDICT     %s" % rec["verdict"])
    print("wrote %s (%.0fs)" % (OUT, rec["seconds"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
