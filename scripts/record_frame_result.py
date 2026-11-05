"""Record the frame-conditioning verdict on the MACHINE-READABLE surface, additively.

The prose memo is written after this, not before: a reader who trusts the JSON must never find it
behind the narrative. Every figure here is READ OUT of
`wiki/promoteratlas_frame_conditioning_2026-11-01.json` rather than retyped, so the two cannot
disagree.

**`arms{}` is asserted byte-identical across the write.** The pre-commitment is that a measured
metric is not retro-edited to match a later finding, and an assertion is the only form of that
promise which cannot rot. The script refuses to save if the block moved by a single byte.

**It is idempotent.** Re-running replaces the `frame_conditioning` sub-block with an identical one
and touches nothing else, so a repeat run is a no-op rather than an append.

Run:  uv run python -m scripts.record_frame_result [--check]
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
ARM = REPO / "wiki/promoteratlas_tile_arm_2026-11-01.json"
RESULT = REPO / "wiki/promoteratlas_frame_conditioning_2026-11-01.json"
BAR = REPO / "wiki/promoteratlas_frame_acceptance_bar.json"
MEMO = "wiki/promoteratlas_frame_result_2026-11-01.md"


def build_block(res: dict, bar: dict) -> dict:
    c1 = res["c1_does_the_encoder_know_corpus_membership"]
    c2 = res["c2_matched_strength_series"]
    c3 = res["c3_stratified"]
    mc = c2["matched_control_premise"]
    st_out = c3["stability_margin_out_window"]
    st_in = c3["stability_margin_in_window"]

    # score the frozen forecast: a hit is the outcome that was given the HIGHEST probability
    probs = {k[2:]: v for k, v in bar["forecast"].items() if k.startswith("p_")}
    modal = max(probs, key=lambda k: probs[k])
    hit = modal == res["verdict"]

    return {
        "date": res["date"],
        "bar": "wiki/promoteratlas_frame_acceptance_bar.json",
        "artifact": "wiki/promoteratlas_frame_conditioning_2026-11-01.json",
        "memo": MEMO,
        "verdict": res["verdict"],
        "status": (
            "MEASURED_BUT_THE_PARTIAL_VS_EXPLAINS_DISTINCTION_IS_UNRESOLVED"
            if not res["verdict_stability"]["resolved"]
            else "MEASURED"
        ),
        "resolved_by_this_run": res["verdict_stability"]["resolved"],
        "anchor_delta_against_the_committed_metric": res["anchor"]["abs_delta"],
        "measured": {
            "c1_auroc_embedding_predicting_corpus_membership": c1["auroc_embedding_median"],
            "c1_auroc_gc_only": c1["auroc_gc_only_median"],
            "subset_atlas": c2["atlas_median"],
            "subset_gc": c2["gc_median"],
            "subset_gc_plus_frame": c2["gc_plus_frame_median"],
            "frame_alone": c2["frame_alone_median"],
            "margin_over_gc": c2["margin_over_gc"],
            "margin_over_gc_plus_frame": c2["margin_over_gc_plus_frame"],
            "frame_absorption": c2["frame_absorption"],
            "attributable_absorption": c2["attributable_absorption"],
            "margin_in_window": c3["margin_in_window"],
            "margin_out_window": c3["margin_out_window"],
            "margin_out_window_seed_range": [st_out["min"], st_out["max"]],
            "margin_in_window_seed_range": [st_in["min"], st_in["max"]],
            "n_seeds_out_window_below_bar": st_out["n_seeds_below_bar"],
            "n_seeds_out_window_negative": st_out["n_seeds_negative"],
            "random_stratification_null_min_margin": c3["random_stratum_min_margin_median"],
            "n_tiles_with_a_frame": res["accounting"]["n_tiles_with_a_frame"],
            "n_tiles_excluded": res["accounting"]["n_tiles_excluded"],
        },
        "reading": (
            "The encoder plainly CARRIES corpus membership (AUROC %.4f against GC's %.4f). What the "
            "gain is made of is the harder question, and the two cells answer it with different "
            "strength. C2, a LINEAR model on five crude positional features, absorbs %.4f of the "
            "%.4f margin and leaves %.4f standing. C3, which controls for membership BY "
            "CONSTRUCTION and which the plan pre-committed as the STRONGER design, finds the "
            "advantage over GC is %.4f inside corpus windows and %.4f outside -- and the "
            "size-matched random stratification null leaves %.4f, so the collapse outside windows "
            "is not an artifact of stratifying."
            % (c1["auroc_embedding_median"], c1["auroc_gc_only_median"],
               c2["frame_absorption"], c2["margin_over_gc"], c2["margin_over_gc_plus_frame"],
               c3["margin_in_window"], c3["margin_out_window"],
               c3["random_stratum_min_margin_median"])
        ),
        "the_strongest_cell_points_HARDER_at_position_than_the_verdict_does": (
            "Outside corpus windows the atlas advantage over GC is NOT DISTINGUISHABLE FROM ZERO: "
            "%d of %d seeds fall below the %.2f bar and %d are NEGATIVE, over a seed range of "
            "[%.4f, %.4f]. The point median landed %.4f above the bar, which is what selected "
            "FRAME_PARTIAL instead of FRAME_EXPLAINS. That is a coin flip, so this run does not "
            "resolve which of the two is right, and FRAME_PARTIAL may NOT be quoted as though the "
            "stratified cell came out clean."
            % (st_out["n_seeds_below_bar"], st_out["n_seeds"], bar["inherited_margin"],
               st_out["n_seeds_negative"], st_out["min"], st_out["max"],
               c3["min_stratum_margin"] - bar["inherited_margin"])
        ),
        "the_attribution_control_did_not_hold_its_own_premise": (
            "`attributable_absorption` %+.4f is only WEAKLY controlled. The control's premise is "
            "matched informativeness, and the non-positional comparators are WEAKER standalone "
            "than the frame (dinucleotide %.4f, k-mer4 %.4f, against the frame's %.4f), so part of "
            "that figure may measure the frame being a BETTER feature rather than a POSITIONAL "
            "one. What the series does establish is directional and still useful: adding "
            "non-positional composition columns did not shrink the atlas margin at all -- it "
            "WIDENED it, because both k-mer baselines lost held-out ground against plain GC -- "
            "while adding the frame shrank it."
            % (c2["attributable_absorption"], mc["dinucleotide_alone_median"],
               mc["kmer4_alone_median"], mc["frame_alone_median"])
        ),
        "may_NEVER_be_written_as": (
            "'position excluded', and equally not 'position explains it'. The bounded form is: the "
            "corpus-window rule AS REPLICATED -- 200 bp upstream of an annotated gene start with a "
            "5 bp intergenic gate, under sorted ordering, on ONE organism's annotation -- absorbs "
            "a material part of the margin, and the advantage outside those windows is "
            "indistinguishable from zero on a statistic this run could not resolve."
        ),
        "metric_is_NOT_retro_edited": (
            "arms.native150 = 0.5076 keeps its measured value, per the frame bar's pre-commitment "
            "(d). The local anchor reproduced it to %.4f on the committed splits before any new "
            "number was read. What changes is the honest limit."
            % res["anchor"]["abs_delta"]
        ),
        "not_a_new_arm": (
            "This is a LIMIT on an existing measured arm. The ceiling verdict is untouched and a "
            "positional or fine-tuned PromoterAtlas variant remains a NEW ARM under R6."
        ),
        "forecast_scored": {
            "modal_prediction": modal,
            "p_modal": probs[modal],
            "outcome": res["verdict"],
            "p_assigned_to_outcome": probs.get(res["verdict"]),
            "hit": hit,
            "standing_record_before": "2 of 7",
            "standing_record_after": "%d of 8" % (2 + (1 if hit else 0)),
            "error_direction": (
                "The modal prediction was that the margin would survive cleanly "
                "(FRAME_EXCLUDED_CONDITIONAL, p=%.2f). It did not: the frame took a material "
                "discount and the stratified cell nearly collapsed. The standing bias on record is "
                "PESSIMISM about whether a run will succeed; this error runs the other way -- "
                "over-confidence that the arm's gain would hold up -- so it is a new data point "
                "rather than another instance of the recorded bias." % probs[modal]
            ),
        },
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--check", action="store_true",
                    help="verify the block is present and consistent; write nothing")
    a = ap.parse_args(argv)

    for p in (ARM, RESULT, BAR):
        if not p.is_file():
            print("REFUSED: %s absent" % p)
            return 3

    raw_before = ARM.read_text(encoding="utf-8")
    arm = json.loads(raw_before)
    res = json.loads(RESULT.read_text(encoding="utf-8"))
    bar = json.loads(BAR.read_text(encoding="utf-8"))

    arms_before = json.dumps(arm["arms"], sort_keys=True)
    block = build_block(res, bar)

    if a.check:
        have = arm.get("corpus_selection_exposure_resolved", {}).get("frame_conditioning")
        if have is None:
            print("ABSENT: no frame_conditioning sub-block")
            return 1
        same = json.dumps(have, sort_keys=True) == json.dumps(block, sort_keys=True)
        print("PRESENT and %s" % ("CONSISTENT with the result artifact" if same
                                  else "DIVERGENT from the result artifact"))
        return 0 if same else 1

    arm.setdefault("corpus_selection_exposure_resolved", {})["frame_conditioning"] = block

    arms_after = json.dumps(arm["arms"], sort_keys=True)
    if arms_after != arms_before:
        print("REFUSED: arms{} changed across the write. A measured metric is not retro-edited.")
        return 4

    text = json.dumps(arm, indent=2, sort_keys=True)
    ARM.write_text(text, encoding="utf-8")

    # re-read and re-assert, because an exit code is not evidence of a write on this host
    back = json.loads(ARM.read_text(encoding="utf-8"))
    if json.dumps(back["arms"], sort_keys=True) != arms_before:
        print("REFUSED: arms{} differs after the round trip")
        return 5
    if "frame_conditioning" not in back["corpus_selection_exposure_resolved"]:
        print("REFUSED: the sub-block is not in the file after writing")
        return 6

    print("arms{} byte-identical across the write: True")
    print("verdict recorded: %s (resolved by the run: %s)"
          % (block["verdict"], block["resolved_by_this_run"]))
    fs = block["forecast_scored"]
    print("forecast: modal %s p=%.2f -> outcome %s, hit=%s, record %s -> %s"
          % (fs["modal_prediction"], fs["p_modal"], fs["outcome"], fs["hit"],
             fs["standing_record_before"], fs["standing_record_after"]))
    print("-> %s" % ARM)
    return 0


if __name__ == "__main__":
    sys.exit(main())
