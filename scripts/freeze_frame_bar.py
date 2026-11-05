"""Freeze the frame-conditioning acceptance bar. Read-only on the experiment; emits the bar.

**Why this is a script and not a hand-written JSON.** Specifying a threshold before asking what the
substrate permits produced four amended thresholds in a single session in this repo, and the
opposite failure -- choosing a threshold after seeing the answer -- is the one a frozen bar exists
to prevent. So this runs AFTER feasibility is established and BEFORE any conditioning number
exists: it measures the window count, the unique-locate rate and the variance of every frame
feature, writes them into the bar as the evidence that the experiment is answerable, and computes
NO ridge, NO atlas number and NO margin. Grep it: nothing here imports `atlas_embed`, `expression`
or sklearn.

**It refuses to overwrite an existing bar.** Re-freezing after a result is visible is exactly the
move the artifact guards against, so a second run needs `--force` and says so loudly.

**The margin is INHERITED, not restated.** It is read out of
`wiki/promoteratlas_arm_acceptance_bar.json` and the script FAILS if that file no longer carries
0.02 -- a bar that restated the number in prose would silently decouple from its parent.

**The check order is IMPORTED from the grader**, not duplicated here, so the two cannot drift. The
bar's test asserts the equality holds in the committed artifact as well.

Run:  uv run python -m scripts.freeze_frame_bar [--force]
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
PARENT_ARM_BAR = REPO / "wiki/promoteratlas_arm_acceptance_bar.json"
PARENT_EXPOSURE_BAR = REPO / "wiki/promoteratlas_exposure_acceptance_bar.json"
ARM = REPO / "wiki/promoteratlas_tile_arm_2026-11-01.json"
FNA = REPO / "data/cache/refseq/GCF_000005845.2/genome.fna"
GFF = REPO / "data/cache/refseq/GCF_000005845.2/annotations.gff3"
SUBSTRATE_CACHE = REPO / "data/mpra"
LOCATE_CACHE = REPO / "data/promoteratlas/frame_locate_cache.json"
OUT = REPO / "wiki/promoteratlas_frame_acceptance_bar.json"

SCHEMA = "promoteratlas-frame-acceptance-bar-v1"
FROZEN_DATE = "2026-11-05"
TILE_LEN = 150

#: The granularity graded as PRIMARY, fixed here so the 66%/19% choice cannot be made after seeing
#: which one gives a better answer.
PRIMARY_GRANULARITY = "in_window_any"
REPORTED_ALSO = "in_window_full"


def _read_inherited_margin() -> float:
    d = json.loads(PARENT_ARM_BAR.read_text(encoding="utf-8"))
    try:
        value = float(d["margin"]["value"])
    except (KeyError, TypeError, ValueError) as e:
        raise SystemExit(
            "REFUSED: %s no longer carries margin.value; the frame bar inherits that number and "
            "may not restate it (%s)" % (PARENT_ARM_BAR.name, e)
        ) from e
    if abs(value - 0.02) > 1e-12:
        raise SystemExit(
            "REFUSED: the parent arm's margin is %r, not 0.02. The frame bar must be re-derived "
            "against whatever the parent now says rather than silently inheriting a stale number."
            % value
        )
    return value


def measure_feasibility() -> dict:
    """Window count, locate buckets and per-feature variance. NO conditioning number."""
    from dna_decode.data.annotations import parse_gff3
    from dna_decode.glm import genomewide as gw
    from dna_decode.glm.corpus_window import (
        cds_rows_from_annotation_table,
        corpus_windows,
        window_coverage,
    )
    from dna_decode.glm.frame_features import FEATURE_NAMES, frame_features
    from dna_decode.glm.tile_locate import (
        LOCATE_FLOOR,
        locate_report,
        locate_sequences,
        read_fasta_single,
        save_locate_cache,
    )

    contig, genome = read_fasta_single(FNA)
    cds = cds_rows_from_annotation_table(parse_gff3(GFF))
    w_sorted = corpus_windows(cds, len(genome), order="sorted")
    w_file = corpus_windows(cds, len(genome), order="file")
    mask = window_coverage(w_sorted, len(genome))

    rows, _ = gw.load_peak_tiles(SUBSTRATE_CACHE, allow_download=False)
    tiles = [t for t in rows if t.category == "tile" and len(t.seq) == TILE_LEN]
    located = locate_sequences([t.seq for t in tiles], genome)
    rep = locate_report(located, tile_length=TILE_LEN, genome_length=len(genome))
    save_locate_cache(LOCATE_CACHE, located, report=rep)

    table = frame_features(located, mask, cds, tile_length=TILE_LEN)

    variance = {}
    for name in FEATURE_NAMES:
        col = table.column(name)
        mean = sum(col) / len(col)
        variance[name] = {
            "mean": round(mean, 6),
            "n_distinct_values": len(set(col)),
            "is_constant": table.is_constant(name),
        }

    return {
        "contig": contig,
        "genome_length": len(genome),
        "cds_rows_with_locus_tag": len(cds),
        "windows_sorted": len(w_sorted),
        "windows_file_order": len(w_file),
        "bases_covered_sorted": sum(mask),
        "genome_fraction_covered": round(sum(mask) / len(genome), 6),
        "n_tiles_scored_by_the_arm": len(tiles),
        "locate": rep.as_dict(),
        "locate_floor": LOCATE_FLOOR,
        "frame_feature_variance": variance,
        "tiles_in_window_any": int(sum(table.column("in_window_any"))),
        "tiles_in_window_full": int(sum(table.column("in_window_full"))),
        "feature_rows_emitted": len(table.seqs),
        "excluded_accounting": table.accounting(),
    }


def build_bar(margin: float, feas: dict) -> dict:
    from dna_decode.glm.frame_verdict import CHECK_ORDER, LOCATE_FLOOR

    arm = json.loads(ARM.read_text(encoding="utf-8"))
    anchors = {
        "committed_atlas_native150_median": float(arm["arms"]["native150"]["median"]),
        "committed_gc_median": float(arm["gc_anchor"]["median"]),
        "headroom_at_risk_atlas_minus_gc": round(
            float(arm["arms"]["native150"]["median"]) - float(arm["gc_anchor"]["median"]), 4
        ),
    }

    outcomes = {
        "INDETERMINATE_HARNESS_DRIFT": {
            "meaning": (
                "The cached embedding pass no longer reproduces the committed arms.native150, so "
                "no frame claim may rest on it. Diagnose the harness; do NOT widen the tolerance."
            ),
            "predicate": "the local anchor fails to reproduce arms.native150 median within 0.02",
        },
        "INDETERMINATE_LOCATE_RATE": {
            "meaning": (
                "The frame is unknown for too much of the cohort, so the scored subset is no "
                "longer the cohort the arm measured."
            ),
            "predicate": "unique-locate rate over distinct tile sequences < %.2f" % LOCATE_FLOOR,
        },
        "INDETERMINATE_BASELINE_NOT_STRENGTHENED": {
            "meaning": (
                "The baseline never received the frame. Adding a column cannot lower a "
                "least-squares TRAINING fit, so this is a plumbing fault and 'no absorption' "
                "would be a statement about the harness. Checked BEFORE informativeness, because "
                "a feed this broken can itself make the frame look uninformative."
            ),
            "predicate": "median training fit of (gc + frame) < median training fit of gc",
        },
        "INDETERMINATE_FRAME_UNINFORMATIVE": {
            "meaning": (
                "The frame carries less expression signal than the gap under test, so 'adding it "
                "did not absorb the gain' is VACUOUSLY true. Pre-committed response: a better "
                "positional feature, never a weaker test."
            ),
            "predicate": "abs(median held-out Spearman of frame features alone) < %.2f" % margin,
        },
        "INDETERMINATE_ABSORPTION_NOT_ATTRIBUTABLE": {
            "meaning": (
                "A cell fired, but the absorption is confounded. Either the frame took no more of "
                "the margin than an equally strong NON-positional feature did, or the "
                "size-matched random stratification null already collapses a stratum margin by "
                "itself. Shrinkage is not absorption. Pre-committed response: a sharper "
                "positional feature or a better-matched control, never dropping the control."
            ),
            "predicate": (
                "(C2 or C3 fired) AND NOT ((C2 fired AND frame absorption - matched absorption > "
                "0) OR (C3 fired AND random-stratification minimum stratum margin >= %.2f))"
                % margin
            ),
        },
        "FRAME_EXPLAINS": {
            "meaning": (
                "Corpus-window position was carrying the gain. arms.native150 keeps its measured "
                "value; the honest limit becomes a quantified positional discount."
            ),
            "predicate": (
                "(C2: atlas - (gc + frame) < %.2f AND frame absorption exceeds the matched "
                "non-positional control) OR (C3: the weaker stratum's within-stratum margin < "
                "%.2f AND the random-stratification null leaves its minimum margin >= %.2f)"
                % (margin, margin, margin)
            ),
        },
        "FRAME_PARTIAL": {
            "meaning": "A material, quantified positional discount rather than a refutation.",
            "predicate": (
                "frame absorption (margin over gc minus margin over gc+frame) > %.2f AND atlas - "
                "(gc + frame) >= %.2f AND the weaker stratum margin >= %.2f"
                % (margin, margin, margin)
            ),
        },
        "FRAME_EXCLUDED_CONDITIONAL": {
            "meaning": (
                "The gain survives the frame-strengthened baseline and BOTH strata. BOUNDED BY "
                "CONSTRUCTION: this excludes THE CORPUS-WINDOW RULE AS REPLICATED -- 200 bp "
                "upstream of an annotated gene start with a 5 bp intergenic gate, under sorted "
                "ordering, on ONE organism's annotation. A sharper or differently-granular "
                "positional feature is NOT excluded, and it may never be written up as "
                "'position excluded'."
            ),
            "predicate": (
                "frame absorption <= %.2f AND atlas - (gc + frame) >= %.2f AND both strata "
                "margins >= %.2f" % (margin, margin, margin)
            ),
        },
    }

    missing = [k for k in CHECK_ORDER if k not in outcomes]
    extra = [k for k in outcomes if k not in CHECK_ORDER]
    if missing or extra:
        raise SystemExit(
            "REFUSED: bar outcomes do not match the grader's CHECK_ORDER (missing %r, extra %r)"
            % (missing, extra)
        )

    return {
        "schema": SCHEMA,
        "frozen": FROZEN_DATE,
        "arm": "promoteratlas tile rung, native 150 bp",
        "criterion": (
            "Does CORPUS-WINDOW POSITION explain the frozen PromoterAtlas encoder's 0.5076 on the "
            "peak-tile rung?"
        ),
        "question": (
            "The subtractive half removed the control-vs-tile detection direction and the atlas "
            "channel barely moved. This is the ADDITIVE half: hand a cheap baseline the exact "
            "positional rule that generated the pretraining corpus and see whether the margin is "
            "absorbed."
        ),
        "not_a_new_arm": (
            "This measures a LIMIT on the existing measured arm. A positional or fine-tuned "
            "PromoterAtlas variant is a NEW ARM under R6 and is candidate 8, not this. This run "
            "must not move arms.native150 or the ceiling verdict whatever it finds."
        ),
        "inherited_margin": margin,
        "inherited_margin_source": (
            "wiki/promoteratlas_arm_acceptance_bar.json margin.value, read at freeze time and "
            "asserted equal to 0.02 rather than restated. Reached here via "
            "wiki/promoteratlas_exposure_acceptance_bar.json, which inherited the same number."
        ),
        "anchors": anchors,
        "primary_frame_granularity": PRIMARY_GRANULARITY,
        "reported_also": REPORTED_ALSO,
        "granularity_rationale": (
            "in_window_any fires on %d of %d uniquely-located tiles and in_window_full on %d, so "
            "the two differ sharply in variance and neither is obviously the right resolution. "
            "The primary is NAMED here, before any conditioning number exists, so the choice "
            "cannot be made after seeing which one answers better. in_window_full means FULLY "
            "COVERED BY THE UNION of windows, which is the definition the committed figure was "
            "measured with."
            % (feas["tiles_in_window_any"], feas["feature_rows_emitted"],
               feas["tiles_in_window_full"])
        ),
        "in_window_any_threshold_bases": 1,
        "feasibility_established_first": feas,
        "controls_required_before_grading": {
            "matched_strength_series": {
                "what": (
                    "The atlas margin is reported over baselines of increasing standalone "
                    "strength: gc, gc + dinucleotide, gc + kmer4, and gc + frame."
                ),
                "why": (
                    "Adding any informative feature to GC strengthens it, which MECHANICALLY "
                    "shrinks the atlas margin whether or not the advantage is positional. "
                    "Shrinkage is not absorption unless it exceeds what a matched-informativeness "
                    "NON-positional feature produces. The subtractive experiment had a "
                    "random-direction control for exactly this; without this series the additive "
                    "one has no analogue."
                ),
                "feeds_outcome": "INDETERMINATE_ABSORPTION_NOT_ATTRIBUTABLE",
            },
            "size_matched_random_stratification": {
                "what": (
                    "Tiles are split into strata of the SAME SIZES at random and the "
                    "within-stratum margins measured, giving the null for how much stratification "
                    "alone moves the margin."
                ),
                "why": (
                    "The real strata are measurably NOT exchangeable -- in-window interquartile "
                    "range 0.4735 with standard deviation 3.380 against out-of-window 0.2710 with "
                    "5.535, so the shapes differ and not merely the scales. A rank correlation "
                    "inside a compressed stratum is attenuated for BOTH channels, so a smaller "
                    "within-stratum margin is expected for reasons unrelated to absorption."
                ),
                "feeds_outcome": "INDETERMINATE_ABSORPTION_NOT_ATTRIBUTABLE",
            },
            "per_stratum_expression_spread": {
                "what": "Each stratum's expression spread is reported beside its margin.",
                "why": (
                    "Without the spread beside the margin a reader cannot tell range restriction "
                    "from absorption."
                ),
            },
        },
        "outcomes": outcomes,
        "verdict_check_order": list(CHECK_ORDER),
        "check_order_rationale": (
            "VALIDITY FIRST, and the ordering is not cosmetic. Run 1 of the parent arm emitted its "
            "STRONGEST verdict while a validity gate was failing, because that gate's predicate "
            "could never fire. Every indeterminate is evaluated before every substantive outcome, "
            "and the grader's test asserts the reachable set EQUALS this order rather than being "
            "contained in it, so no branch can be dead."
        ),
        "deviation_from_the_technical_plan": {
            "what": (
                "The technical plan names SEVEN outcomes; this bar and its grader carry EIGHT. "
                "The addition is INDETERMINATE_ABSORPTION_NOT_ATTRIBUTABLE."
            ),
            "why": (
                "The plan's own pre-execution review found that shrinkage is not absorption "
                "unless it exceeds a matched non-positional control, and that the bar cannot be "
                "frozen without that control and the random-stratification null. A control that "
                "cannot change the verdict is not a control: with only seven outcomes the grader "
                "would be forced to report FRAME_EXPLAINS over a measurement it knows is "
                "confounded."
            ),
            "status": (
                "Recorded for the user to ratify or reject. Not treated as settled, and the plan "
                "file is left unedited."
            ),
        },
        "forecast": {
            "note": (
                "Scored after the run. My forecast record in this project is 2 of 7 and the "
                "errors have been PESSIMISTIC -- the measured bias is that I talk the project out "
                "of runs that then succeed. Recorded to be scored, not to be deferred to."
            ),
            "p_FRAME_EXCLUDED_CONDITIONAL": 0.44,
            "p_FRAME_PARTIAL": 0.25,
            "p_FRAME_EXPLAINS": 0.08,
            "p_INDETERMINATE_FRAME_UNINFORMATIVE": 0.10,
            "p_INDETERMINATE_ABSORPTION_NOT_ATTRIBUTABLE": 0.08,
            "p_INDETERMINATE_BASELINE_NOT_STRENGTHENED": 0.03,
            "p_INDETERMINATE_HARNESS_DRIFT": 0.015,
            "p_INDETERMINATE_LOCATE_RATE": 0.005,
            "reasoning": (
                "Five crude positional features -- two binaries, a coverage fraction, a signed "
                "distance and a CDS flag -- are being asked to absorb a 0.1964 margin that a "
                "128-dimensional learned representation produces. The subtractive half already "
                "found the detection direction worth only 0.0080 of the atlas channel, so the "
                "gain is not obviously the detection axis in any simple form. The likeliest "
                "outcomes are therefore that the margin survives, or that it takes a measurable "
                "but partial discount. FRAME_EXPLAINS is given real weight anyway because the "
                "encoder was pretrained on exactly these windows and 66% of tiles touch one."
            ),
        },
        "pre_commitments": [
            "The frame must be shown INFORMATIVE before any absorption claim is read, and the "
            "response to an uninformative frame is a better positional feature, never a weaker "
            "test.",
            "FRAME_EXCLUDED_CONDITIONAL is BOUNDED: it excludes the corpus-window rule as "
            "replicated, at 200 bp with a 5 bp intergenic gate, under sorted ordering, on one "
            "organism's annotation. It may NEVER be written as 'position excluded'.",
            "The stratified cell's baseline must be re-measured ON THE SAME STRATUM, never "
            "compared to the full-cohort 0.5076, because stratifying changes both n and the peak "
            "composition.",
            "arms.native150 keeps its measured value whatever the outcome, and the ceiling verdict "
            "is not touched by this run.",
            "Absorption is judged against the matched non-positional series and the "
            "size-matched random stratification null. Neither control may be dropped to reach a "
            "cleaner verdict.",
        ],
        "honest_limits_fixed_in_advance": [
            "The replica is scoped to the corpus rule AS SPECIFIED (200 bp, 5 bp gate, sorted "
            "ordering), never 'as the authors ran it': the upstream function threads its gate in "
            "FILE order over a GenBank file this repo does not have, and the two orderings differ "
            "by 26 windows.",
            "One organism, one annotation build, one encoder, one rung.",
            "Ambiguous (325) and not-found (31) tile sequences are excluded and counted, so the "
            "scored set is the uniquely-located subset and not the full 44,106 rows.",
            "689 of the 44,106 scored rows share a sequence with another row, so a sequence "
            "duplicated across two peaks can land on both sides of a leave-peak-out split. This "
            "is a property of the substrate and of every arm on this rung, reported rather than "
            "fixed here.",
            "The frame features are crude by construction; failing to absorb the margin bounds "
            "THESE features, not every possible positional encoding.",
            "Strata are not exchangeable in outcome distribution, which is why the random "
            "stratification null and the per-stratum spreads are mandatory rather than optional.",
        ],
        "frozen_before": [
            "any ridge was fitted on frame features",
            "any atlas-versus-frame margin was computed",
            "any within-stratum number existed",
        ],
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--force", action="store_true",
                    help="overwrite an existing frozen bar (re-freezing after a result is "
                         "visible is the move this artifact exists to prevent)")
    a = ap.parse_args(argv)

    if OUT.exists() and not a.force:
        print("REFUSED: %s already exists." % OUT.name)
        print("A frozen bar is not rewritten. Re-freezing after a conditioning number exists is")
        print("exactly the move it guards against. Pass --force only if nothing has been measured.")
        return 2

    for p in (PARENT_ARM_BAR, ARM, FNA, GFF):
        if not p.is_file():
            print("REFUSED: %s absent" % p)
            return 3

    t0 = time.time()
    margin = _read_inherited_margin()
    print("inherited margin %.2f from %s" % (margin, PARENT_ARM_BAR.name))
    print("measuring feasibility (no conditioning number is computed) ...", flush=True)
    feas = measure_feasibility()
    bar = build_bar(margin, feas)

    probs = {k: v for k, v in bar["forecast"].items() if k.startswith("p_")}
    total = sum(probs.values())
    if abs(total - 1.0) > 1e-6:
        print("REFUSED: forecast probabilities sum to %.6f, not 1.0" % total)
        return 4

    OUT.write_text(json.dumps(bar, indent=2, sort_keys=True), encoding="utf-8")

    print("windows(sorted) %d over %d bases (%.2f%% of the genome)" % (
        feas["windows_sorted"], feas["bases_covered_sorted"],
        100.0 * feas["genome_fraction_covered"]))
    print("locate: %d unique / %d distinct (%.4f), %d ambiguous, %d not found" % (
        feas["locate"]["n_unique"], feas["locate"]["n_distinct"],
        feas["locate"]["unique_rate"], feas["locate"]["n_ambiguous"],
        feas["locate"]["n_not_found"]))
    print("in_window_any %d / in_window_full %d of %d emitted rows" % (
        feas["tiles_in_window_any"], feas["tiles_in_window_full"],
        feas["feature_rows_emitted"]))
    print("no frame feature is constant: %s" % all(
        not v["is_constant"] for v in feas["frame_feature_variance"].values()))
    print("primary granularity FROZEN as %r (also reporting %r)"
          % (PRIMARY_GRANULARITY, REPORTED_ALSO))
    print("outcomes %d, check order matches the grader" % len(bar["outcomes"]))
    print("-> %s  (%.0fs)" % (OUT, time.time() - t0))
    return 0


if __name__ == "__main__":
    sys.exit(main())
