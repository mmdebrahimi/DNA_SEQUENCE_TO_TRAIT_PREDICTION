"""The ten rejection gates, runnable. A DECISION-SUPPORT screen that REFUSES rather than guesses.

WHY THIS IS CODE NOW, AND NOT BEFORE. The gates were written as a prose memo
(`wiki/negative_results_map_2026-06-13.md`) and applied by hand. Coding them earlier would have been
trust-layer theatre — a polished artifact encoding subjective judgements as booleans, built with no real
candidate to teach it what the schema needs. Two candidates have now been screened by hand (PEAR
2026-08-31, HBV 2026-09-01), which is the precondition: enough worked examples to derive the shape,
few enough that it is still grounded.

THE SPLIT THAT MAKES IT HONEST. Eight gates carry a countable rule already ("<20/class = trip",
"majority-censored = trip", "<~3 effective lineages"). TWO DO NOT:

  G1  is the label wet-lab/clinical, or produced by a genomic tool the decoder would compete against?
  G3  is the label an assay reading, or a description of where/why the isolate was collected?

Those are readings of a dataset's methods section, not computations. So they take a HUMAN EVIDENCE
string, and a screen missing it returns `NEEDS_HUMAN_EVIDENCE` and the overall verdict REFUSES. It never
defaults to pass. That is the whole difference between a screen and theatre: the gates that cannot be
computed must block the verdict rather than quietly clear it.

INTENDED LAYER COMES FIRST, and that is the PEAR lesson. The same dataset screens differently depending
on what it is FOR: as an L4 constructed-molecular substrate PEAR's continuous growth readout makes G6
(phenotype censoring at a clinical breakpoint) inapplicable, while as an L1 R/S label source the same
readout must be binarised and G6 is live. Screening without declaring the layer scores one dataset under
incompatible criteria.

A PASSING SCREEN IS NOT A BUILD RECOMMENDATION. PEAR clears every gate that applies to it and is still
not buildable here — its processed data ships as serialized plot objects. The gates bound whether a
usable LABEL exists; they say nothing about whether the artifact is reachable, the work is worth doing,
or the regime has a measured positive (that is `eval/regime.py`).

Pure; no I/O, no network. Every threshold is sourced to the memo or to a measurement, never invented.
"""
from __future__ import annotations

from dataclasses import dataclass, field

# --- verdicts -------------------------------------------------------------------------------------
PASS = "pass"                              # measured, and the gate does not trip
TRIP = "trip"                              # measured, and the gate trips -> the candidate is rejected
NOT_APPLICABLE = "not_applicable"          # this gate cannot apply given the intended layer / design
NEEDS_HUMAN_EVIDENCE = "needs_human_evidence"   # G1/G3 only: a judgement no computation can make
INSUFFICIENT_DATA = "insufficient_data"    # a mechanical gate whose inputs were not supplied

L1_AMR_RS = "L1_AMR_RS"                    # a deterministic R/S decoder label source
L4_FORWARD_CONTINUOUS = "L4_forward_continuous"   # a constructed-variation continuous-fitness substrate
LAYERS = (L1_AMR_RS, L4_FORWARD_CONTINUOUS)

# Thresholds. Each is SOURCED -- to the memo's own wording, or to a measurement that set it.
MIN_PER_CLASS = 20            # G4: memo, "<20/class = trip"
MIN_EFFECTIVE_LINEAGES = 3    # G8: memo, "<~3 effective lineages = trip"
MAJORITY = 0.50               # G6/G9/G10: the memo's "majority = trip" convention
DOMINANT_SOURCE_SHARE = 0.60  # G2: NOT in the memo as a number. Taken from the source-concentration
                              # arm, which refuses a cell when one source exceeds 60%
                              # (scripts/source_diverse_validate.py). Stated so it is arguable.
# G6's L4 form. Copied from the SHIPPED gate `forward_inverse_roundtrip.assay_degeneracy`, which is
# where these were set and measured (CcdB: 79.3% at its ceiling, 8 levels -> excluded).
MAX_MODE_SHARE = 0.25         # >25% of variants at one value -> the target grid collapses onto it
MIN_DISTINCT_VALUES = 20      # fewer distinct levels than targets -> "percentile" is a step function


@dataclass
class GateResult:
    gate: str
    verdict: str
    reason: str
    evidence: dict = field(default_factory=dict)

    def as_dict(self) -> dict:
        return {"gate": self.gate, "verdict": self.verdict, "reason": self.reason,
                "evidence": dict(self.evidence)}


def _num(ev: dict, key: str):
    v = ev.get(key)
    return v if isinstance(v, (int, float)) else None


# --- the two JUDGEMENT gates: they require a human reading, and refuse without one -----------------

def g1_circular_label(ev: dict) -> GateResult:
    """Is the phenotype produced by a genomic tool the decoder would compete against?"""
    prov = (ev.get("label_provenance_evidence") or "").strip()
    if not prov:
        return GateResult("G1", NEEDS_HUMAN_EVIDENCE,
                          "no label-provenance reading supplied. Whether a label is a wet-lab measurement "
                          "or a genome->label tool's output is read from the dataset's methods, not "
                          "computed. Supply `label_provenance_evidence` + `label_is_measured`.")
    measured = ev.get("label_is_measured")
    if measured is None:
        return GateResult("G1", NEEDS_HUMAN_EVIDENCE,
                          "provenance text supplied but `label_is_measured` not asserted",
                          {"label_provenance_evidence": prov})
    if measured:
        return GateResult("G1", PASS, f"label recorded as a measurement: {prov}",
                          {"label_provenance_evidence": prov})
    return GateResult("G1", TRIP,
                      f"label is tool-derived, so scoring a rule against it scores a rule against a "
                      f"rule: {prov}", {"label_provenance_evidence": prov})


def g3_sampling_defined_label(ev: dict) -> GateResult:
    """Is the label an assay reading, or a description of where/why the isolate was collected?"""
    eviden = (ev.get("label_semantics_evidence") or "").strip()
    if not eviden:
        return GateResult("G3", NEEDS_HUMAN_EVIDENCE,
                          "no label-semantics reading supplied. 'assay reading vs sampling context' is a "
                          "judgement about what the label MEANS. Supply `label_semantics_evidence` + "
                          "`label_is_assay_reading`.")
    is_assay = ev.get("label_is_assay_reading")
    if is_assay is None:
        return GateResult("G3", NEEDS_HUMAN_EVIDENCE,
                          "semantics text supplied but `label_is_assay_reading` not asserted",
                          {"label_semantics_evidence": eviden})
    if is_assay:
        return GateResult("G3", PASS, f"label is an assay reading: {eviden}",
                          {"label_semantics_evidence": eviden})
    return GateResult("G3", TRIP, f"label IS the sampling context, not a measurement: {eviden}",
                      {"label_semantics_evidence": eviden})


# --- the eight MECHANICAL gates -------------------------------------------------------------------

def g2_study_equals_class(ev: dict) -> GateResult:
    # APPLICABILITY BEFORE MEASUREMENT. Asking for the input first makes an inapplicable gate report
    # `insufficient_data`, which reads as "go measure this" when the honest answer is "this cannot apply".
    if ev.get("variation_is_constructed"):
        return GateResult("G2", NOT_APPLICABLE,
                          "variation is constructed, so there is no source-vs-class confound to have")
    share = _num(ev, "largest_source_share")
    if share is None:
        return GateResult("G2", INSUFFICIENT_DATA, "no `largest_source_share` supplied")
    v = TRIP if share > DOMINANT_SOURCE_SHARE else PASS
    return GateResult("G2", v, f"largest source holds {share:.0%} of the cohort "
                               f"(bar {DOMINANT_SOURCE_SHARE:.0%})", {"largest_source_share": share})


def g4_surveillance_domination(ev: dict) -> GateResult:
    if ev.get("variation_is_constructed"):
        return GateResult("G4", NOT_APPLICABLE,
                          "not a surveillance corpus; there is no ecosystem to exclude")
    n = _num(ev, "non_ecosystem_min_class_n")
    if n is None:
        return GateResult("G4", INSUFFICIENT_DATA, "no `non_ecosystem_min_class_n` supplied")
    v = TRIP if n < MIN_PER_CLASS else PASS
    return GateResult("G4", v, f"smallest non-surveillance class n={n} (bar {MIN_PER_CLASS})",
                      {"non_ecosystem_min_class_n": n})


def g5_assembly_attrition(ev: dict) -> GateResult:
    if ev.get("genotype_defined_by_construction"):
        return GateResult("G5", NOT_APPLICABLE,
                          "genotypes are known by construction; no assembly needs fetching")
    n = _num(ev, "n_fetchable_assemblies")
    if n is None:
        return GateResult("G5", INSUFFICIENT_DATA, "no `n_fetchable_assemblies` supplied")
    v = TRIP if n < MIN_PER_CLASS else PASS
    return GateResult("G5", v, f"{n} records have a fetchable assembly (bar {MIN_PER_CLASS})",
                      {"n_fetchable_assemblies": n})


def g6_phenotype_censoring(ev: dict, layer: str) -> GateResult:
    """LAYER-DISPATCHED, and this is the subtlest gate.

    G6's letter is MIC interval-censoring at a clinical breakpoint. Its SPIRIT is wider: *the quantitative
    readout cannot separate where it matters*. Under L1 that is breakpoint censoring. Under L4 there is no
    breakpoint at all, but the same failure arrives as an ASSAY FLOOR -- every dead variant scores the same
    value, so the target grid collapses onto the pile-up. That is why the PEAR screen scored G6 OPEN rather
    than n/a on a continuous assay, and it was right to: CcdB posted the forward/inverse sweep's BEST
    number precisely because 79.3% of it was tied at the ceiling.

    So L4 evaluates G6 with DEGENERACY inputs, on the shipped `assay_degeneracy` thresholds. A screen that
    returned n/a here would clear a candidate on the exact failure that has already bitten this repo once.
    """
    if layer == L4_FORWARD_CONTINUOUS:
        share, distinct = _num(ev, "mode_share"), _num(ev, "n_distinct_values")
        if share is None or distinct is None:
            return GateResult("G6", INSUFFICIENT_DATA,
                              "continuous readout: G6 takes its L4 (assay-degeneracy) form and needs "
                              "`mode_share` + `n_distinct_values`. Neither breakpoint censoring nor "
                              "degeneracy has been screened, so this gate is OPEN, not passed.")
        v = TRIP if (share > MAX_MODE_SHARE or distinct < MIN_DISTINCT_VALUES) else PASS
        return GateResult("G6", v, f"assay degeneracy: {share:.1%} of values at the mode over {distinct:.0f} "
                                   f"distinct levels (bars {MAX_MODE_SHARE:.0%} / {MIN_DISTINCT_VALUES})",
                          {"mode_share": share, "n_distinct_values": distinct})
    frac = _num(ev, "censored_fraction")
    if frac is None:
        return GateResult("G6", INSUFFICIENT_DATA, "no `censored_fraction` supplied")
    v = TRIP if frac > MAJORITY else PASS
    return GateResult("G6", v, f"{frac:.0%} of in-class values are interval-censored (bar {MAJORITY:.0%})",
                      {"censored_fraction": frac})


def g7_provenance_not_separable(ev: dict) -> GateResult:
    if ev.get("variation_is_constructed"):
        return GateResult("G7", NOT_APPLICABLE,
                          "no provenance-disjoint split is needed; the split is by position/variant")
    frac = _num(ev, "provenance_field_populated_fraction")
    if frac is None:
        return GateResult("G7", INSUFFICIENT_DATA, "no `provenance_field_populated_fraction` supplied")
    v = PASS if frac > MAJORITY else TRIP
    return GateResult("G7", v, f"submitter/center/collection populated on {frac:.0%} of records",
                      {"provenance_field_populated_fraction": frac})


def g8_dedup_collapses_balance(ev: dict) -> GateResult:
    if ev.get("variation_is_constructed"):
        return GateResult("G8", NOT_APPLICABLE,
                          "ancestry is randomised by construction, so there is no clonality to collapse")
    n = _num(ev, "min_effective_lineages")
    if n is None:
        return GateResult("G8", INSUFFICIENT_DATA, "no `min_effective_lineages` supplied")
    v = TRIP if n < MIN_EFFECTIVE_LINEAGES else PASS
    return GateResult("G8", v, f"smallest class has {n} effective lineages "
                               f"(bar {MIN_EFFECTIVE_LINEAGES})", {"min_effective_lineages": n})


def g9_causal_variant_unrecorded(ev: dict) -> GateResult:
    frac = _num(ev, "loci_without_recorded_variant_fraction")
    if frac is None:
        return GateResult("G9", INSUFFICIENT_DATA, "no `loci_without_recorded_variant_fraction` supplied")
    v = TRIP if frac > MAJORITY else PASS
    return GateResult("G9", v, f"{frac:.0%} of loci record no concrete causal variant "
                               f"(bar {MAJORITY:.0%})", {"loci_without_recorded_variant_fraction": frac})


def g10_variant_class_off_panel(ev: dict) -> GateResult:
    frac = _num(ev, "off_panel_variant_fraction")
    if frac is None:
        return GateResult("G10", INSUFFICIENT_DATA, "no `off_panel_variant_fraction` supplied")
    v = TRIP if frac > MAJORITY else PASS
    return GateResult("G10", v, f"{frac:.0%} of recorded causal variants are indel/structural "
                                f"(bar {MAJORITY:.0%})", {"off_panel_variant_fraction": frac})


LABEL_GATES = ("G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8")   # is a usable LABEL available?
DECODER_GATES = ("G9", "G10")                                    # is the RULE scoreable at all?


@dataclass
class ScreenResult:
    candidate: str
    intended_layer: str
    gates: list
    verdict: str
    reason: str

    def as_dict(self) -> dict:
        return {"candidate": self.candidate, "intended_layer": self.intended_layer,
                "verdict": self.verdict, "reason": self.reason,
                "gates": [g.as_dict() for g in self.gates],
                "contract": ("A PASS bounds only whether a usable LABEL exists and whether the rule is "
                             "scoreable. It is NOT a build recommendation: reachability of the artifact, "
                             "regime fit (eval/regime.py) and worth-doing are separate questions.")}


def screen_candidate(candidate: str, intended_layer: str, evidence: dict) -> ScreenResult:
    """Run all ten gates. REFUSES a verdict while any judgement gate lacks its human reading."""
    if intended_layer not in LAYERS:
        return ScreenResult(candidate, intended_layer, [], "REFUSED",
                            f"intended_layer must be one of {LAYERS}; got {intended_layer!r}. The layer "
                            "decides which gates even apply, so it cannot be inferred.")
    gates = [
        g1_circular_label(evidence), g2_study_equals_class(evidence),
        g3_sampling_defined_label(evidence), g4_surveillance_domination(evidence),
        g5_assembly_attrition(evidence), g6_phenotype_censoring(evidence, intended_layer),
        g7_provenance_not_separable(evidence), g8_dedup_collapses_balance(evidence),
        g9_causal_variant_unrecorded(evidence), g10_variant_class_off_panel(evidence),
    ]
    tripped = [g.gate for g in gates if g.verdict == TRIP]
    if tripped:
        return ScreenResult(candidate, intended_layer, gates, "REJECTED",
                            f"gate(s) {', '.join(tripped)} tripped")
    needs = [g.gate for g in gates if g.verdict == NEEDS_HUMAN_EVIDENCE]
    if needs:
        return ScreenResult(candidate, intended_layer, gates, "REFUSED",
                            f"{', '.join(needs)} require a human reading of the dataset's methods and "
                            "none was supplied. A screen that defaulted these to pass would be theatre.")
    missing = [g.gate for g in gates if g.verdict == INSUFFICIENT_DATA]
    if missing:
        return ScreenResult(candidate, intended_layer, gates, "INCOMPLETE",
                            f"no measurement supplied for {', '.join(missing)}")
    return ScreenResult(candidate, intended_layer, gates, "CLEARS",
                        "every applicable gate passes. This bounds the LABEL question only — not "
                        "artifact reachability, regime fit, or whether the work is worth doing.")
