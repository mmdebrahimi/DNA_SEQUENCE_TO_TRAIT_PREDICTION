"""F1 of the broaden-horizon dataset hunt: a DETERMINISTIC candidate scorecard.

Scores a genotype->phenotype dataset candidate against the 8 quality gates distilled from this project's
negative-results map + embedding-niche 3-part test + horse-coat joinability lesson (see
plans/Broaden_Horizon_Dataset_Hunt_Plan_2026-07-02.md). Pure logic, no I/O — the research (F2) fills in the
per-gate verdicts; this module turns them into an overall verdict + a decoder-type routing.

Gate verdict values: "pass" / "fail" / "unknown". A gate that is "unknown" is NOT a pass (fail-closed for the
overall PASS, but tracked separately so F2 knows what to go verify).
"""
from __future__ import annotations

from dataclasses import dataclass, field

GATES = {
    "G1_accessible": "free + fetchable, no DUA/paywall",
    "G2_non_circular": "phenotype measured independently, not derived from a tool the decoder competes with",
    "G3_sampling_independent": "phenotype not confounded with sampling context (site/source/study/date)",
    "G4_unit_joinable": "per-individual genotype joined to per-individual phenotype (not aggregate)",
    "G5_provenance_separable": "a leakage-free split exists (temporal/accession/cohort)",
    "G6_depth_or_catalog": ">=~100 same-organism units (learned) OR a curated determinant catalog exists",
    "G7_genotype_fetchable": "actual sequence/variants per unit are downloadable",
    "G8_label_not_censored": "quantitative labels are tierable (not all operator-censored at one bound)",
}
GATE_KEYS = list(GATES)
_VALID = {"pass", "fail", "unknown"}


@dataclass
class Candidate:
    name: str
    creature: str
    phenotype: str
    gates: dict            # gate_key -> "pass"/"fail"/"unknown"
    has_curated_catalog: bool = False
    depth_estimate: int | None = None     # # same-organism units, if known
    notes: str = ""
    sources: list = field(default_factory=list)

    def __post_init__(self):
        bad = {k: v for k, v in self.gates.items() if v not in _VALID}
        if bad:
            raise ValueError(f"{self.name}: invalid gate verdicts {bad} (must be pass/fail/unknown)")
        missing = set(GATE_KEYS) - set(self.gates)
        if missing:
            raise ValueError(f"{self.name}: missing gate verdicts {sorted(missing)}")


def decoder_type(c: Candidate) -> str:
    """Which product paradigm applies. Curated catalog -> the VALIDATED deterministic path. Else, if it has
    depth + non-circular + sampling-independent, the LEARNED niche (high bar: embeddings 0-for-4)."""
    if c.has_curated_catalog:
        return "deterministic"          # the project's validated product; G6 satisfied by the catalog
    deep = (c.depth_estimate or 0) >= 100
    if deep and c.gates["G2_non_circular"] == "pass" and c.gates["G3_sampling_independent"] == "pass":
        return "learned-niche"          # eligible, but embeddings are 0-for-4 de-confounded -> high bar
    return "neither"                    # no catalog + not learned-eligible -> not a substrate yet


def score(c: Candidate) -> dict:
    """Overall verdict. PASS = all 8 gates pass (catalog relaxes G6). fail_closed: any fail. unknowns tracked."""
    gates = dict(c.gates)
    if c.has_curated_catalog:
        gates["G6_depth_or_catalog"] = "pass"            # catalog satisfies G6 by definition
    fails = [k for k, v in gates.items() if v == "fail"]
    unknowns = [k for k, v in gates.items() if v == "unknown"]
    dt = decoder_type(c)
    if fails:
        verdict = "REJECT"
    elif unknowns:
        verdict = "VERIFY"                                # nothing fails, but gaps remain to confirm
    elif dt == "neither":
        verdict = "REJECT"                               # all gates pass but no viable decoder paradigm
    else:
        verdict = "PASS"
    n_pass = sum(1 for v in gates.values() if v == "pass")
    return {"name": c.name, "creature": c.creature, "verdict": verdict, "decoder_type": dt,
            "n_pass": n_pass, "fails": fails, "unknowns": unknowns,
            "depth_estimate": c.depth_estimate, "has_curated_catalog": c.has_curated_catalog}


def rank(cands: list[Candidate]) -> list[dict]:
    """PASS > VERIFY > REJECT; then more gates passed; then greater depth; then name (deterministic)."""
    order = {"PASS": 0, "VERIFY": 1, "REJECT": 2}
    scored = [score(c) for c in cands]
    scored.sort(key=lambda s: (order[s["verdict"]], -s["n_pass"], -(s["depth_estimate"] or 0), s["name"]))
    return scored
