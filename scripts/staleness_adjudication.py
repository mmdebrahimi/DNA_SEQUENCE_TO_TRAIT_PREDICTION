"""The hand-adjudication record for the full-corpus staleness flagging pass, as DATA not prose.

WHY THIS IS A FILE. The flagging pass is only half the instrument -- every flag has to be checked against
its artifact by hand, and that verdict is the actual product. Keeping the verdicts as a committed table
(rather than a paragraph in a memo) means the precision number can be RE-DERIVED instead of trusted, and a
later run can be scored against the same rubric.

THE HEADLINE THIS FILE EXISTS TO MAKE HONEST. The benchmark measured 4/5 precision on a 50%-positive test
set. The real corpus ran at 2/9. That is not a contradiction and not a regression -- it is the base-rate
effect, and `precision_from_base_rate` computes it, so the claim is arithmetic rather than hand-waving.
A screen with excellent specificity still returns mostly false positives when the thing it screens for is
rare, which is exactly the regime a documentation corpus is in.

Run: uv run python scripts/staleness_adjudication.py
"""
from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Adjudicated:
    artifact: str
    model_evidence: str
    verdict: str          # "true_positive" | "false_positive"
    why: str


# Every flag the 80-item pass raised, checked against its artifact by hand on 2026-08-27.
ADJUDICATIONS: tuple[Adjudicated, ...] = (
    Adjudicated(
        "wiki/negative_results_map_2026-06-13.md",
        "lists 10 rejection gates (G1-G10) including G9/G10 added 2026-08-26, contradicting 'exactly 8'",
        "true_positive",
        "CLAUDE.md said '8 reusable rejection GATES ... G1-G8'. I added G9/G10 the previous day and never "
        "updated the summarising bullet. The cited file exists and every path resolves, so no mechanical "
        "check could reach it -- only the claim ABOUT the file went stale. Fixed.",
    ),
    Adjudicated(
        "tests/test_models_cache.py",
        "the tests do not reference verify_complete / CompletenessReport / the 4 status buckets",
        "true_positive",
        "RIGHT FILE, WRONG REASON -- `verify_complete` appears 19 times, so the model's stated evidence is "
        "false. But the claim said '8 unit tests' and the file now has 42: a genuine count drift, the same "
        "class already fixed in README. Counted as a true positive because the claim IS stale, with the "
        "caveat that the evidence would not have survived adjudication on its own.",
    ),
    Adjudicated(
        "scripts/bvbrc_strict_mic_4drug_census.py",
        "implements gentamicin breakpoints, contradicting 'gent substrate infeasible'",
        "false_positive",
        "The capability-vs-finding category error, now confirmed on an unseen item (it also produced the "
        "benchmark's single FP). Code existing refutes a CAPABILITY claim ('X is not built'); it SUPPORTS "
        "a FINDING claim, because the script is the instrument that produced the finding.",
    ),
    Adjudicated(
        "plans/Trait_Decoding_Roadmap.md",
        "labelled 'DRAFT 2026-05-26', contradicting 'shipped 2026-05-26'",
        "false_positive",
        "The header is a title convention; 'shipped' refers to the artifact landing, not to a status field.",
    ),
    Adjudicated(
        "project_state/eukaryotic-trait-decoding-cycle-2026-06-07.md",
        "'Verdict: PASS' contradicts the claim that the deterministic decoder is the validated artifact",
        "false_positive",
        "The ledger says CYCLE COMPLETE / H2 FALSIFIED / embedding arm closed -- verbatim what the claim "
        "says. The model latched onto a PASS string belonging to a DIFFERENT gate (G1) in the same file.",
    ),
    Adjudicated(
        "plans/EP8_PathB_PreStage_Manifest.md",
        "the harness is actively implemented with detailed steps, contradicting the claim",
        "false_positive",
        "The claim says the harness 'remains as now-historical scaffolding' -- i.e. it EXISTS and is "
        "historical. Existing is what the claim asserts, not a contradiction of it.",
    ),
    Adjudicated(
        "dna_decode/forward/prosst_scorer.py",
        "the artifact states the forward pass ran on a Kaggle T4, contradicting 'not on Kaggle'",
        "false_positive",
        "wiki/prosst_lift_2026-07-18.md says 'Ran entirely LOCALLY on CPU -- no Kaggle' and has a section "
        "titled 'Local, not Kaggle'. The model asserted the opposite of the text.",
    ),
    Adjudicated(
        "wiki/provdisjoint_census_results.json",
        "contains 'powered': false entries, contradicting 'the NOT_CENSUSED bucket is empty'",
        "false_positive",
        "Two different fields. The report card has 0 NOT_CENSUSED cells (verified: 27 cells, none in that "
        "state); `powered: false` maps to UNDERPOWERED, which the claim itself reports as 3.",
    ),
    Adjudicated(
        "wiki/decoder_v0_ux_and_success_criterion.md",
        "lists 'strict-MIC training labels' as an explicit non-criterion, contradicting the claim",
        "false_positive",
        "The heading is 'Explicit non-criteria for v0' and the claim says v0 trains on CATEGORICAL labels. "
        "The artifact AGREES; the model read the heading correctly and inverted the conclusion.",
    ),
)

N_SCORED = 80          # verdicts recovered before the OOM
N_SUPPORTED = 65
N_UNCLEAR = 6
N_UNPARSEABLE = 4


def tally() -> dict[str, int]:
    tp = sum(a.verdict == "true_positive" for a in ADJUDICATIONS)
    return {"flags": len(ADJUDICATIONS), "true_positive": tp,
            "false_positive": len(ADJUDICATIONS) - tp}


def precision() -> float:
    t = tally()
    return t["true_positive"] / t["flags"]


def specificity() -> float:
    """Share of genuinely-not-stale items the auditor did NOT flag.

    This is the number that should be compared against the benchmark, NOT precision -- it is the one the
    base rate does not move.
    """
    t = tally()
    negatives = N_SCORED - t["true_positive"]
    return (negatives - t["false_positive"]) / negatives


def precision_from_base_rate(base_rate: float, sens: float, spec: float) -> float:
    """Expected precision of a screen. PURE -- the arithmetic behind the headline.

    Makes the benchmark-vs-corpus gap a calculation instead of an excuse: the SAME detector posts high
    precision on a 50%-positive benchmark and low precision on a 2%-positive corpus.
    """
    tp = base_rate * sens
    fp = (1 - base_rate) * (1 - spec)
    return tp / (tp + fp) if (tp + fp) else 0.0


def main() -> int:
    t = tally()
    print(f"flags adjudicated: {t['flags']}  "
          f"true positives: {t['true_positive']}  false positives: {t['false_positive']}")
    print(f"precision on the real corpus: {precision():.2f}")
    print(f"specificity: {specificity():.3f}   (base rate of stale claims: "
          f"{t['true_positive']/N_SCORED:.3f})\n")

    print("Benchmark precision was 4/5 = 0.80 on a 50%-positive set. Same detector, both regimes:")
    for label, br in (("benchmark (50% positive)", 0.50), ("real corpus", t["true_positive"] / N_SCORED)):
        exp = precision_from_base_rate(br, sens=1.0, spec=specificity())
        print(f"  base rate {br:5.3f} -> expected precision {exp:.2f}   [{label}]")
    print("\nThe gap is the base rate, not a regression in the detector.\n")

    for a in ADJUDICATIONS:
        mark = "TP" if a.verdict == "true_positive" else "FP"
        print(f"[{mark}] {a.artifact}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
