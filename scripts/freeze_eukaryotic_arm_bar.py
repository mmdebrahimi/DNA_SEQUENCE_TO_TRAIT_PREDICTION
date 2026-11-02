"""Freeze the acceptance bar for the LAST planned arm on the tile-chain ceiling criterion.

The arm is `substrate-change-constructed-eukaryotic-or-human-panel`, declared 2026-10-09 and the
only `planned` arm left (`ceiling_verdict` = `blocked:untried-arm`, 10 measured / 1 planned). It is
the R8 frame check made into an experiment: every other arm changed the METHOD and the rung topped
out at 0.5076 against a replicate-derived ceiling of 0.941. This arm changes the SUBSTRATE instead
and asks whether the ceiling was ever the method's to raise.

**The comparison quantity is DERIVED, not asserted, and it is NOT the headline number.** Raw
Spearman cannot cross a substrate boundary: two panels have different noise floors, different
assays and different dynamic range, so "0.5076 here vs 0.83 there" compares two things that were
never the same measurement. What DOES cross is a within-substrate contrast -- how much a learned
encoder buys over the dumb baseline ON ITS OWN PANEL, with both sides measured the same way on the
same splits. Both halves of that contrast already exist on the bacterial rung as committed arms:

    conv-3width-13k-params-150bp   0.3234
    gc-single-feature              0.3112
    ------------------------------------- gain  +0.0122

and that +0.0122 is the SAME number the project already reported independently as its encoder gate
(`wiki/glm_encoder_gate_result_2026-10-08.md`: median +0.0122, paired CI includes zero). So the
anchor is not invented for this bar; it is two committed arm metrics and a third artifact agreeing.

**Why the conv encoder is the primary and the frozen atlas is NOT.** PromoterAtlas is pretrained on
gammaproteobacterial promoters. Running it on human sequence is out-of-distribution BY
CONSTRUCTION, so a low human atlas number would be unattributable -- it could mean the substrate is
not the limit, or it could mean this encoder does not transfer, and this arm cannot separate those.
A conv encoder is RETRAINED on whichever panel it is pointed at, which removes that confound
entirely. The atlas is still measured, and reported, as a secondary with the OOD caveat attached.

**The published incumbents are REFUSED as a verdict component, and that is pre-committed here.**
Table S6 publishes per-fold values for Enformer, MPRALegNet, Sei, MPRAnn and a biochemical feature
set -- and never states the metric. A head-to-head against an unnamed metric would be an assumption
wearing a number's clothes. They are recorded verbatim with `metric_unstated: true` and they inform
nothing in the predicate. What is legitimate is the INTERNAL contrast, where both sides are ours and
the metric is named.
"""

from __future__ import annotations

import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "wiki/eukaryotic_arm_acceptance_bar.json"

#: Committed arm metrics on the bacterial rung. Read from soraya_runs/arms-ledger.json, not retyped
#: from memory -- `verify()` below re-reads them and refuses on drift.
BACT_CONV = 0.3234
BACT_GC = 0.3112
BACT_BEST = 0.5076
BACT_CEILING = 0.941

CRITERION = ("ceiling: best achievable prediction of bacterial promoter activity from sequence, "
             "reported against the replicate-derived 0.941 ceiling")
ARM = "substrate-change-constructed-eukaryotic-or-human-panel-UNRUN"


def bacterial_anchor() -> dict:
    return {
        "conv_minus_gc": round(BACT_CONV - BACT_GC, 4),
        "conv": BACT_CONV,
        "gc": BACT_GC,
        "best_arm": BACT_BEST,
        "replicate_ceiling": BACT_CEILING,
        "fraction_of_ceiling_captured_by_best": round(BACT_BEST / BACT_CEILING, 4),
        "headroom_best_minus_gc": round(BACT_BEST - BACT_GC, 4),
        "corroborating_artifact": "wiki/glm_encoder_gate_result_2026-10-08.md",
        "why_this_is_the_anchor": (
            "A within-substrate contrast is the only quantity that crosses a substrate boundary. "
            "Raw Spearman does not: the panels have different noise floors and different assays."
        ),
    }


def verify(ledger_path: Path) -> list[str]:
    """Re-read the committed arm metrics and refuse on drift. A bar anchored to retyped numbers is
    anchored to nothing -- this is the same discipline as reading a figure from its artifact."""
    problems = []
    data = json.loads(ledger_path.read_text(encoding="utf-8"))
    arms = data["criteria"].get(CRITERION)
    if arms is None:
        return ["criterion key absent from the arms ledger: %r" % CRITERION]
    for name, expected in (("conv-3width-13k-params-150bp", BACT_CONV),
                           ("gc-single-feature", BACT_GC),
                           ("promoteratlas-base-embeddings-frozen-150bp-padded-to-200", BACT_BEST)):
        got = (arms.get(name) or {}).get("metric")
        if got is None:
            problems.append("arm %s has no metric on disk" % name)
        elif abs(got - expected) > 1e-9:
            problems.append("arm %s: disk %.4f != bar %.4f" % (name, got, expected))
    if ARM not in arms:
        problems.append("the arm this bar is for is not declared: %s" % ARM)
    elif arms[ARM].get("status") != "planned":
        problems.append("arm status is %r, expected 'planned' -- a bar frozen after a run is not a "
                        "bar" % arms[ARM].get("status"))
    return problems


def build() -> dict:
    anchor = bacterial_anchor()
    gain = anchor["conv_minus_gc"]
    return {
        "schema": "eukaryotic-arm-acceptance-bar-v1",
        "frozen": "2026-11-01",
        "criterion": CRITERION,
        "arm": ARM,
        "question": (
            "Was the bacterial rung's ceiling the METHOD's or the SUBSTRATE's? Every measured arm "
            "changed the method and the rung topped out at %.4f of a %.3f replicate ceiling. Run "
            "the same learned-encoder-vs-dumb-baseline contrast on a CONSTRUCTED human panel and "
            "compare the CONTRAST, not the raw number."
            % (anchor["fraction_of_ceiling_captured_by_best"], BACT_CEILING)
        ),
        "primary_measurement": {
            "what": "conv_human_minus_gc_human",
            "definition": (
                "Spearman(predicted, measured activity) for a conv encoder RETRAINED on the human "
                "panel, minus the same for a GC-only ridge, on IDENTICAL published folds."
            ),
            "why_conv_and_not_the_atlas": (
                "PromoterAtlas is pretrained on gammaproteobacterial promoters, so on human "
                "sequence it is out-of-distribution BY CONSTRUCTION and a low number would be "
                "unattributable between 'substrate is not the limit' and 'this encoder does not "
                "transfer'. A retrained conv has no such confound."
            ),
            "splits": "the panel's own published fold assignment (Table S9 `fold` column)",
            "uncertainty": "paired bootstrap over folds; the CI is the verdict input, not the point",
        },
        "secondary_measurements": [
            {"what": "atlas_human_minus_gc_human",
             "status": "reported, NEVER a verdict input",
             "caveat": "frozen bacterial encoder on human sequence = out-of-distribution"},
            {"what": "fraction_of_ceiling_captured",
             "status": "reported only if the human panel supplies a replicate-derived reliability "
                       "estimate; refused otherwise rather than substituted"},
        ],
        "anchors": anchor,
        "outcomes": {
            "SUBSTRATE_LIMITED": {
                "predicate": "paired-bootstrap CI lower bound of conv_human_minus_gc_human > %.4f"
                             % gain,
                "meaning": (
                    "A learned encoder buys MORE on a constructed eukaryotic panel than it bought "
                    "on bacterial tiles. The bacterial ceiling was the substrate's, and the rung's "
                    "ten method arms were competing for headroom that was never there."
                ),
            },
            "METHOD_LIMITED": {
                "predicate": "CI upper bound < %.4f" % gain,
                "meaning": (
                    "A learned encoder buys no more -- or less -- on the new substrate. Changing "
                    "substrate does not rescue the approach, and the limit travels with the method."
                ),
            },
            "SUBSTRATE_NEUTRAL": {
                "predicate": "CI contains %.4f" % gain,
                "meaning": (
                    "No detectable difference. The arm is MEASURED and discriminates nothing, which "
                    "is a real outcome and not a failure -- it retires the arm without moving the "
                    "ceiling story either way."
                ),
            },
            "INDETERMINATE_JOIN_RATE": {
                "predicate": "sequence-to-activity join rate < 0.95 (the loader's MIN_JOIN_RATE)",
                "meaning": "the substrate did not assemble; nothing is scored",
            },
            "INDETERMINATE_NO_FOLDS": {
                "predicate": "the published fold assignment does not resolve onto the loaded rows",
                "meaning": (
                    "scoring on splits of my own invention would silently change what the "
                    "incumbents' folds mean; refuse instead"
                ),
            },
            "INDETERMINATE_DEGENERATE_TARGET": {
                "predicate": "activity mode-share > 0.25 or fewer than 20 distinct levels",
                "meaning": (
                    "the shipped `assay_degeneracy` bar, imported not restated. CcdB posted this "
                    "project's best-looking number purely because 79.3% of its values sat at a "
                    "ceiling, so this screen runs BEFORE any correlation is read."
                ),
            },
        },
        "reachable_outcomes_equal_check_order": True,
        "refused_comparisons": {
            "published_incumbents_head_to_head": {
                "why": (
                    "Table S6 publishes per-fold values for Enformer, MPRALegNet, Sei, MPRAnn and a "
                    "biochemical feature set and NEVER STATES THE METRIC. Comparing our named "
                    "Spearman to an unnamed quantity would be an assumption wearing a number's "
                    "clothes."
                ),
                "what_is_done_instead": (
                    "they are recorded verbatim with metric_unstated: true, as context for what "
                    "this panel supports, and they inform no predicate"
                ),
            },
        },
        "honest_limits_fixed_in_advance": [
            "THE BIGGEST ONE: this is ONE substrate change, not a substrate sweep. A human "
            "lentiMPRA panel differs from bacterial tiles in assay, in element length (230 bp vs "
            "150 bp), in kingdom and in how the library was constructed. The arm cannot say WHICH "
            "of those differences carries any effect it finds.",
            "The bacterial anchor %.4f is treated as a FIXED POINT; its own uncertainty is not "
            "propagated into the CI, so a verdict that lands near the boundary is weaker than the "
            "predicate makes it look." % gain,
            "The conv encoder's hyperparameters were chosen on the BACTERIAL rung. Reusing them on "
            "human sequence is the honest apples-to-apples choice for a method-class comparison, "
            "but it is NOT a tuned human model and must never be reported as one.",
            "A constructed panel is the regime this project already knows learned methods work in "
            "(constructed -> molecular, Spearman 0.761 on TEM-1). So a SUBSTRATE_LIMITED result "
            "CORROBORATES an existing finding rather than discovering something; its value is that "
            "it bounds the tile-chain ceiling, not that it is news about constructed panels.",
            "Nothing here licenses a claim about eukaryotic g->p generally. The project's "
            "organism-level eukaryotic embedding arm is a CLOSED NEGATIVE on natural populations "
            "(Arabidopsis, de-confounded) and this arm does not reopen it: constructed variation "
            "is a different regime, which is the whole reason this arm is scoped to a constructed "
            "panel.",
        ],
        "forecast": {
            "p_SUBSTRATE_LIMITED": 0.55,
            "p_SUBSTRATE_NEUTRAL": 0.2,
            "p_METHOD_LIMITED": 0.1,
            "p_any_INDETERMINATE": 0.15,
            "reasoning": (
                "Published methods reach ~0.73-0.83 on this panel against a bacterial rung whose "
                "BEST is 0.5076, and constructed-molecular is the one regime where this project's "
                "learned arms have measured wins. So headroom probably exists. Against that: the "
                "conv is tiny (13k params), untuned for human sequence, and GC may itself be "
                "stronger on human regulatory elements."
            ),
            "my_record": (
                "2 of 6 in this project. The last bar put p=0.08 on the outcome that happened. "
                "Do NOT defer to these priors."
            ),
        },
        "frozen_before": [
            "any GC-on-human number existed",
            "any conv-on-human number existed",
            "any atlas-on-human number existed",
            "the fold assignment had been resolved onto loaded rows",
        ],
        "known_before_freezing_and_deliberately_not_a_predicate": [
            "the Table S6 incumbent values (published prior art, not a result of this arm)",
            "that the substrate is on disk and MD5-verified (wiki/lentimpra_acquisition_2026-10-27.md)",
        ],
    }


def main() -> int:
    problems = verify(REPO / "soraya_runs/arms-ledger.json")
    if problems:
        print("REFUSING to freeze -- the anchors do not match the committed ledger:")
        for p in problems:
            print("  -", p)
        return 2
    rec = build()
    OUT.write_text(json.dumps(rec, indent=2, sort_keys=True), encoding="utf-8")
    print("anchor conv_minus_gc = %.4f (verified against the arms ledger)"
          % rec["anchors"]["conv_minus_gc"])
    print("outcomes: %s" % ", ".join(sorted(rec["outcomes"])))
    print("-> %s" % OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
