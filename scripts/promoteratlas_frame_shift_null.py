"""Is the frame table's absorption specific to the TRUE annotation, or to annotation-SHAPED features?

The frame arm established that handing five positional features to the GC baseline absorbs 0.0687
of its 0.1950 margin on the located subset. The bar's required control for that number -- a
matched-informativeness series of non-positional comparators -- FAILED its own premise, because
dinucleotide (0.2895) and 4-mer (0.2811) are weaker standalone than the frame (0.3170). So
`attributable_absorption = +0.0949` is only weakly controlled: part of it may measure the frame
being a BETTER feature rather than a POSITIONAL one.

This run supplies the control that is matched by CONSTRUCTION instead of by hoping two feature
families happen to be equally strong. The whole CDS annotation is rotated around the circular
chromosome and the entire feature table is rebuilt from the rotated rows -- windows and CDS columns
alike. The rotated annotation has the same gene count, the same gene-length multiset, the same
strand composition, the same intergenic spacing distribution and the same coding density as the
real one, exactly rather than approximately. The only thing destroyed is the PHASE between the
annotation and the measured tiles.

So it answers precisely one question: would an annotation of this STRUCTURE, placed anywhere,
absorb as much of the atlas margin as the real annotation does?

  * If the rotated annotation absorbs about as much, then the absorption is a property of adding
    five annotation-shaped columns to a ridge, NOT evidence that real corpus-window position
    explains the atlas margin.
  * If the rotated annotation absorbs little or nothing, the absorption IS specific to true
    positions, and the frame explanation survives a control it has not yet faced.

The same rotation supplies a STRUCTURALLY matched null for the stratified cell too. The frozen run
already carries a SIZE-matched random stratification null (minimum margin 0.1860, which cleared the
0.02 bar widely). A rotated-annotation stratification is strictly stronger: it matches not only
stratum SIZE but the spatial clustering of the strata, which a random split destroys. If the
in-window/out-of-window collapse (0.2689 against 0.0230) reproduces under rotation, the collapse is
an artifact of stratifying by anything with this spatial structure rather than a fact about corpus
windows.

**NOTHING HERE IS GRADED, and that is deliberate.** The frozen acceptance bar
(`wiki/promoteratlas_frame_acceptance_bar.json`) does not name this control. Adding a control after
seeing the result it bears on, and then letting it select a verdict, is the move the bar exists to
prevent -- so no threshold is read, no verdict is computed, and the published FRAME_PARTIAL stands
untouched with its `resolved: false`. These numbers are REPORTED for the user to decide on, in the
same posture as the run's own `secondary_granularity_in_window_full` block.

**True and rotated values are measured in the SAME run, on the same splits**, so the comparison is
paired rather than taken against a remembered figure. The true figures reproducing the committed
0.3170 and 0.0687 is itself the check that this harness is the one that produced them.

**The anchor runs FIRST**, before any new number is read.

Run:  uv run python -m scripts.promoteratlas_frame_shift_null
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
OUT = REPO / "wiki/promoteratlas_frame_shift_null_2026-11-14.json"

SCHEMA = "promoteratlas-frame-shift-null-v1"
N_SEEDS = 10
N_OFFSETS = 5
FRAC = 0.25
TILE_LEN = 150
ANCHOR_TOL = 0.02
MIN_STRATUM = 50


def _spread(xs) -> dict:
    a = np.asarray([float(x) for x in xs], dtype=float)
    a = a[~np.isnan(a)]
    if a.size == 0:
        return {"n": 0, "note": "no usable values"}
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
    from dna_decode.glm.annotation_shift import (
        circular_shift_cds,
        sample_shift_offset,
        shift_match_report,
    )
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

    rows, _ = gw.load_peak_tiles(CACHE, allow_download=False)
    tiles = [t for t in rows if t.category == "tile" and len(t.seq) == TILE_LEN]
    n = len(tiles)
    Xa = np.load(EMB_DIR / ("tile_embeddings_native150_n%d.npy" % n))
    if Xa.shape != (n, 128):
        raise SystemExit("embedding cache is %r, expected (%d, 128)" % (Xa.shape, n))

    _contig, genome = read_fasta_single(FNA)
    glen = len(genome)
    cds_true = cds_rows_from_annotation_table(parse_gff3(GFF))

    if LOCATE_CACHE.is_file():
        located = load_locate_cache(LOCATE_CACHE)
        if {t.seq for t in tiles} - set(located):
            located = locate_sequences([t.seq for t in tiles], genome)
    else:
        located = locate_sequences([t.seq for t in tiles], genome)
    lrep = locate_report(located, tile_length=TILE_LEN, genome_length=glen)

    y = np.asarray([t.expression for t in tiles], dtype=float)
    gc = np.asarray(gw.gc_features([t.seq for t in tiles]), dtype=float)
    pos = {id(t): i for i, t in enumerate(tiles)}
    primary_j = list(FEATURE_NAMES).index("in_window_any")

    def build(cds_rows):
        """Frame matrix and in-window indicator for one annotation. Returns (F, has_frame, inwin)."""
        mask = window_coverage(corpus_windows(cds_rows, glen), glen)
        table = frame_features(located, mask, cds_rows, tile_length=TILE_LEN)
        by_seq = {s: i for i, s in enumerate(table.seqs)}
        F = np.full((n, len(FEATURE_NAMES)), np.nan, dtype=float)
        hf = np.zeros(n, dtype=bool)
        for i, t in enumerate(tiles):
            j = by_seq.get(t.seq)
            if j is not None:
                F[i] = table.matrix[j]
                hf[i] = True
        inwin = np.zeros(n, dtype=float)
        inwin[hf] = F[hf][:, primary_j]
        return F, hf, inwin, int(mask.count(1) if isinstance(mask, bytes) else sum(mask)), \
            len(corpus_windows(cds_rows, glen))

    # ---- splits, fixed once so true and rotated arms see identical folds ---------------------
    splits = []
    for seed in range(N_SEEDS):
        train, test, _ = gw.leave_peak_out(tiles, frac=FRAC, seed=seed)
        splits.append((np.asarray([pos[id(t)] for t in train]),
                       np.asarray([pos[id(t)] for t in test])))

    # ---- the anchor, before any new number is read -------------------------------------------
    anchor = []
    for tr, te in splits:
        a, _ = ridge_spearman(Xa[tr], y[tr], Xa[te], y[te])
        anchor.append(round(float(a), 4))
    anchor_med = float(np.median(np.asarray(anchor, dtype=float)))
    anchor_delta = abs(anchor_med - committed_atlas)
    if anchor_delta > ANCHOR_TOL:
        raise SystemExit(
            "ANCHOR FAILED: full-cohort atlas median %.4f vs committed %.4f (delta %.4f > %.4f)"
            % (anchor_med, committed_atlas, anchor_delta, ANCHOR_TOL)
        )
    if verbose:
        print("anchor %.4f vs committed %.4f (delta %.4f) OK  (%.0fs)"
              % (anchor_med, committed_atlas, anchor_delta, time.time() - t0), flush=True)

    def measure(F, hf, inwin, label):
        """Absorption and stratum margins for one annotation, over the fixed splits."""
        gcF = np.hstack([gc, F])
        out = {k: [] for k in ("atlas", "gc", "frame_alone", "gc_frame",
                               "margin_in_window", "margin_out_window",
                               "n_in_window", "n_out_window")}
        for tr, te in splits:
            trL, teL = tr[hf[tr]], te[hf[te]]
            for key, X in (("atlas", Xa), ("gc", gc), ("frame_alone", F), ("gc_frame", gcF)):
                v, _ = ridge_spearman(X[trL], y[trL], X[teL], y[teL])
                out[key].append(round(float(v), 4))
            for tag in ("in_window", "out_window"):
                s_tr = inwin[trL] > 0.5 if tag == "in_window" else inwin[trL] <= 0.5
                s_te = inwin[teL] > 0.5 if tag == "in_window" else inwin[teL] <= 0.5
                i_tr, i_te = trL[s_tr], teL[s_te]
                out["n_%s" % tag].append(int(i_te.size))
                if i_tr.size < MIN_STRATUM or i_te.size < MIN_STRATUM:
                    out["margin_%s" % tag].append(float("nan"))
                    continue
                a, _ = ridge_spearman(Xa[i_tr], y[i_tr], Xa[i_te], y[i_te])
                g, _ = ridge_spearman(gc[i_tr], y[i_tr], gc[i_te], y[i_te])
                out["margin_%s" % tag].append(round(float(a) - float(g), 4))
        m = {k: float(np.nanmedian(np.asarray(v, dtype=float))) for k, v in out.items()}
        res = {
            "label": label,
            "atlas_median": round(m["atlas"], 4),
            "gc_median": round(m["gc"], 4),
            "frame_alone_median": round(m["frame_alone"], 4),
            "gc_plus_frame_median": round(m["gc_frame"], 4),
            "margin_over_gc": round(m["atlas"] - m["gc"], 4),
            "margin_over_gc_plus_frame": round(m["atlas"] - m["gc_frame"], 4),
            "frame_absorption": round((m["atlas"] - m["gc"]) - (m["atlas"] - m["gc_frame"]), 4),
            "baseline_gain_from_frame": round(m["gc_frame"] - m["gc"], 4),
            "margin_in_window_median": round(m["margin_in_window"], 4),
            "margin_out_window_median": round(m["margin_out_window"], 4),
            "min_stratum_margin": round(min(m["margin_in_window"], m["margin_out_window"]), 4),
            "median_n_in_window_test": int(m["n_in_window"]),
            "median_n_out_window_test": int(m["n_out_window"]),
            "n_located": int(hf.sum()),
            "per_seed": out,
        }
        return res

    # ---- the TRUE annotation, recomputed here so the comparison is within-run ----------------
    F_t, hf_t, inwin_t, covered_t, nwin_t = build(cds_true)
    true_arm = measure(F_t, hf_t, inwin_t, "true_annotation")
    true_arm["n_windows"] = nwin_t
    true_arm["window_covered_bases"] = covered_t
    if verbose:
        print("TRUE   absorption %+.4f | in-window %.4f out-of-window %.4f | frame_alone %.4f  (%.0fs)"
              % (true_arm["frame_absorption"], true_arm["margin_in_window_median"],
                 true_arm["margin_out_window_median"], true_arm["frame_alone_median"],
                 time.time() - t0), flush=True)

    # ---- the rotated annotations -------------------------------------------------------------
    shifted = []
    for k in range(N_OFFSETS):
        d = sample_shift_offset(cds_true, glen, seed=k)
        cds_s = circular_shift_cds(cds_true, glen, d)
        rep = shift_match_report(cds_true, cds_s, glen)
        if not (rep["gene_count_match"] and rep["gene_length_multiset_match"]
                and rep["strand_composition_match"] and rep["cds_covered_bases_match"]):
            raise SystemExit(
                "offset %d did not preserve the matched quantities: %r -- the control is not "
                "matched by construction and must not be reported as if it were" % (d, rep)
            )
        F_s, hf_s, inwin_s, covered_s, nwin_s = build(cds_s)
        arm_s = measure(F_s, hf_s, inwin_s, "rotated_offset_%d" % d)
        arm_s["offset"] = d
        arm_s["n_windows"] = nwin_s
        arm_s["window_covered_bases"] = covered_s
        arm_s["match_report"] = rep
        shifted.append(arm_s)
        if verbose:
            print("ROT %d  delta %7d absorption %+.4f | in %.4f out %.4f | frame_alone %.4f "
                  "| windows %d covered %d  (%.0fs)"
                  % (k, d, arm_s["frame_absorption"], arm_s["margin_in_window_median"],
                     arm_s["margin_out_window_median"], arm_s["frame_alone_median"],
                     nwin_s, covered_s, time.time() - t0), flush=True)

    abs_s = [a["frame_absorption"] for a in shifted]
    gain_s = [a["baseline_gain_from_frame"] for a in shifted]
    fa_s = [a["frame_alone_median"] for a in shifted]
    in_s = [a["margin_in_window_median"] for a in shifted]
    out_s = [a["margin_out_window_median"] for a in shifted]
    med_abs_s = float(np.median(np.asarray(abs_s, dtype=float)))
    attributable = true_arm["frame_absorption"] - med_abs_s

    rec = {
        "schema": SCHEMA,
        "date": "2026-11-14",
        "status": "REPORTED, not graded",
        "question": (
            "Would an annotation of the SAME STRUCTURE, placed anywhere on the chromosome, absorb "
            "as much of the atlas margin as the real annotation does?"
        ),
        "not_graded_because": (
            "The frozen acceptance bar does not name this control. Adding a control after seeing "
            "the result it bears on, and letting it select a verdict, is the move the bar exists "
            "to prevent. No threshold is read here and no verdict is computed; the published "
            "FRAME_PARTIAL and its resolved=false status are untouched. These numbers are for the "
            "user to decide on."
        ),
        "why_this_control_replaces_the_proposed_one": (
            "The proposed shuffled-window control was measured invalid first "
            "(wiki/promoteratlas_frame_column_split_2026-11-05.json): frame_features takes window "
            "information only as the union mask, so relocating windows leaves the two CDS-derived "
            "columns byte-identical, and those two columns ALONE match the full five-column table "
            "(0.3183 against 0.3170; at least the full table on 6 of 10 seeds). A circular "
            "rotation of the whole annotation regenerates BOTH halves, so nothing is invariant."
        ),
        "anchor": {
            "committed_native150_median": round(committed_atlas, 4),
            "reproduced_full_cohort_median": round(anchor_med, 4),
            "abs_delta": round(anchor_delta, 4),
            "tolerance": ANCHOR_TOL,
            "ok": True,
            "per_seed": anchor,
        },
        "locate": dataclasses.asdict(lrep) if dataclasses.is_dataclass(lrep) else lrep,
        "design": {
            "n_split_seeds": N_SEEDS,
            "n_rotation_offsets": N_OFFSETS,
            "held_out_fraction": FRAC,
            "splits_shared": (
                "The same fixed folds are used for the true and every rotated arm, so the "
                "comparison is paired and no part of it is explained by different splits."
            ),
            "matched_by_construction": [
                "gene count", "gene-length multiset", "strand composition",
                "total CDS-covered bases", "coding density",
                "intergenic spacing distribution (rotated)",
            ],
            "destroyed_on_purpose": "the phase between the annotation and the measured tiles",
        },
        "true_annotation": true_arm,
        "rotated_annotations": shifted,
        "comparison": {
            "true_frame_absorption": true_arm["frame_absorption"],
            "rotated_frame_absorption": _spread(abs_s),
            "attributable_absorption_vs_rotated_median": round(attributable, 4),
            "true_baseline_gain_from_frame": true_arm["baseline_gain_from_frame"],
            "rotated_baseline_gain_from_frame": _spread(gain_s),
            "true_frame_alone": true_arm["frame_alone_median"],
            "rotated_frame_alone": _spread(fa_s),
            "true_margin_in_window": true_arm["margin_in_window_median"],
            "rotated_margin_in_window": _spread(in_s),
            "true_margin_out_window": true_arm["margin_out_window_median"],
            "rotated_margin_out_window": _spread(out_s),
            "true_stratum_collapse": round(
                true_arm["margin_in_window_median"] - true_arm["margin_out_window_median"], 4),
            "rotated_stratum_collapse": _spread(
                [a["margin_in_window_median"] - a["margin_out_window_median"] for a in shifted]),
            "reading_guide": (
                "attributable_absorption_vs_rotated_median near zero means the absorption is a "
                "property of annotation-SHAPED columns rather than of true positions. A rotated "
                "stratum collapse near the true one means the in/out collapse is an artifact of "
                "stratifying by this spatial structure rather than a fact about corpus windows."
            ),
        },
        "honest_limits": [
            "NOT GRADED. No threshold is read and no verdict is computed; this does not resolve "
            "FRAME_PARTIAL versus FRAME_EXPLAINS and is not an amendment to the frozen bar.",
            "%d rotation offsets is a small sample of the rotation group; the spread across them "
            "is reported rather than summarised to a point, and a median of %d values carries "
            "little precision." % (N_OFFSETS, N_OFFSETS),
            "A rotation preserves structure but not biology: rotated windows sit upstream of "
            "nothing real, so this controls for annotation SHAPE and placement, not for the "
            "possibility that some other real genomic feature co-varies with corpus windows.",
            "The out-of-window stratum remains range-restricted and fatter-tailed in the true "
            "arm, which is a property of the substrate that no control removes.",
            "Window count and covered bases can differ by a few from the true arm, because "
            "corpus_windows skips genes within one regulatory-region length of the linear "
            "coordinate boundary and a rotation moves which genes those are. The per-arm counts "
            "are reported so the size of that effect is visible rather than assumed negligible.",
        ],
        "seconds": round(time.time() - t0, 1),
    }
    OUT.write_text(json.dumps(rec, indent=2) + "\n", encoding="utf-8")
    return rec


def main() -> int:
    rec = run()
    c = rec["comparison"]
    print("")
    print("STATUS   %s" % rec["status"])
    print("ABSORB   true %+.4f | rotated median %+.4f (min %+.4f max %+.4f n=%d)"
          % (c["true_frame_absorption"], c["rotated_frame_absorption"]["median"],
             c["rotated_frame_absorption"]["min"], c["rotated_frame_absorption"]["max"],
             c["rotated_frame_absorption"]["n"]))
    print("         attributable vs rotated median %+.4f" % c["attributable_absorption_vs_rotated_median"])
    print("ALONE    true %.4f | rotated median %.4f"
          % (c["true_frame_alone"], c["rotated_frame_alone"]["median"]))
    print("STRATA   true in %.4f out %.4f (collapse %.4f)"
          % (c["true_margin_in_window"], c["true_margin_out_window"], c["true_stratum_collapse"]))
    print("         rot  in %.4f out %.4f (collapse median %.4f, min %.4f max %.4f)"
          % (c["rotated_margin_in_window"]["median"], c["rotated_margin_out_window"]["median"],
             c["rotated_stratum_collapse"]["median"], c["rotated_stratum_collapse"]["min"],
             c["rotated_stratum_collapse"]["max"]))
    print("wrote %s (%.0fs)" % (OUT, rec["seconds"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
