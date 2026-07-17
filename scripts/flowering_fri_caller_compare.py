"""Does calling FRI ourselves beat the source's own column? -- and is that question even answerable here?

Compares three FRI sources through the SAME honesty gates as the main scorer (null baseline, ancestry
weighting, attrition), plus a circularity control that decides whether the headline is quotable at all:

  s3_column  -- Table S3's `deleterious_allele`, i.e. the source's own answer (what the v0 cell was handed)
  putative   -- OUR caller over Table S1's variants, annotation-only. Reproducing s3_column is its PASS.
  curated    -- putative + L276R/L294F, the two substitutions Zhang 2020 proved non-functional by transgenic
                assay. The paper's own summary column calls their carriers FUNCTIONAL.

THE CIRCULARITY CONTROL IS THE POINT. The paper says, in its own words, that a016/a060/a089 were
"selected because their association with early flowering in Figure 1" -- Figure 1 being built on the SAME
FT16 phenotype we score against. So L294F/L276R were chosen BECAUSE their carriers are early HERE.
Adding them and then scoring on FT16 is circular: the gain is guaranteed by construction, not earned.

What makes this more than a caveat is that the selection was NOT self-confirming. Three alleles were picked
the same way; the transgenic assay (an independent experiment, common null background) CONFIRMED two and
REFUTED one -- FRI-Van-0/a060 delays flowering, so it is functional despite its early carriers. So:

    `circular` -- the naive data-driven rule "any early-associated allele must be LoF" = curated + a060.
                  It is MECHANISTICALLY WRONG (the paper disproved a060) and it should still SCORE BETTER,
                  because a060's carriers are early here by selection. If it does, that is a live
                  demonstration that "improves the score" is not a valid criterion on this dataset, and
                  that the curated rule's support is the EXPERIMENT, never this number.

So the deliverable is a REFUSAL to quote curated-vs-putative as validation, backed by a control that shows
what quoting it would have cost.

Run:  uv run --with openpyxl python scripts/flowering_fri_caller_compare.py
Exit: 0 = the caller reproduces the source column (its real PASS); 1 = it does not.
"""
from __future__ import annotations

import argparse
import json
from datetime import date
from pathlib import Path

from dna_decode.organism_rules.arabidopsis_fri_caller import (
    VERIFIED_LOF_SUBSTITUTIONS,
    call_all_accessions,
    reference_integrity_ok,
)
from scripts.flowering_tables3_score import (
    TABLE_S3,
    _has_phenotype,
    _is_deleterious,
    load_table_s3,
    score,
)

REPO = Path(__file__).resolve().parents[1]
TABLE_S1 = REPO / "data" / "arabidopsis" / "zhang2020" / "tpj14716-sup-0010-TableS1.xlsx"

# The allele the paper selected the same way as a016/a089 and then REFUTED by transgenic assay.
# "the FRI-Van-0 (a060) allele delayed flowering time in the transgenic lines, indicating that this
#  allele is functional" -- so calling it LoF is WRONG, however early its carriers are here.
REFUTED_EARLY_ASSOCIATED_ALLELE = "a060"


def _fri_calls(rule: str, s3_rows: list[dict]) -> dict[str, str]:
    calls = {a: c.status for a, c in call_all_accessions(TABLE_S1, rule if rule != "circular" else "curated").items()}
    if rule == "circular":
        # The naive move: trust the data association instead of the experiment.
        for r in s3_rows:
            if r["allele_group"] == REFUTED_EARLY_ASSOCIATED_ALLELE and calls.get(r["accession_id"]) == "functional":
                calls[r["accession_id"]] = "lof"
    return calls


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=REPO / "wiki")
    args = ap.parse_args()

    if not reference_integrity_ok(TABLE_S1):
        print("[fri-compare] REFUSING: caller integrity guard FAILED (variant/LoF arithmetic does not "
              "reconcile to the paper's stated 171 variants / 26 LoF).")
        return 1

    rows = load_table_s3()
    out: dict[str, dict] = {"s3_column": score(rows)}
    for rule in ("putative", "curated", "circular"):
        out[rule] = score(rows, fri_calls=_fri_calls(rule, rows))

    # --- the caller's REAL pass: does `putative` reproduce the source's own column? -----------------------
    calls = call_all_accessions(TABLE_S1, "putative")
    agree = dis = 0
    for r in rows:
        c = calls.get(r["accession_id"])
        if c is None or c.status == "unknown":
            continue
        src = "lof" if _is_deleterious(r) else "functional"
        agree += c.status == src
        dis += c.status != src
    reproduces = dis == 0

    # --- what the curated rule actually changes ---------------------------------------------------------
    cur = call_all_accessions(TABLE_S1, "curated")
    flipped = [r for r in rows
               if cur.get(r["accession_id"]) and cur[r["accession_id"]].status == "lof" and not _is_deleterious(r)]
    flipped_ph = [r for r in flipped if _has_phenotype(r)]

    def head(k: str) -> dict:
        p, g = out[k]["pooled"], out[k]["group_weighted"]
        return {"n": p["n"], "abstained": out[k]["n_abstained_unknown_fri"],
                "accuracy": p["accuracy"], "null": p["null_accuracy"],
                "sensitivity": p["sensitivity"], "specificity": p["specificity"],
                "tp": p["tp"], "fp": p["fp"], "tn": p["tn"], "fn": p["fn"],
                "group_weighted_accuracy": g["mean_accuracy"], "group_weighted_null": g["mean_null_accuracy"]}

    rep = {
        "schema": "flowering-fri-caller-compare-v1",
        "date": date.today().isoformat(),
        "caller": "dna_decode/organism_rules/arabidopsis_fri_caller.py",
        "caller_reproduces_source_column": {
            "verdict": "PASS" if reproduces else "FAIL",
            "agree": agree, "disagree": dis,
            "claim": ("the `putative` rule, run over Table S1's variant matrix, reproduces Table S3's own "
                      "`deleterious_allele` column exactly -- so the caller is verified against the source's "
                      "answer, and the cell can now call FRI from variants instead of being handed it"),
        },
        "sources": {k: head(k) for k in out},
        "curated_delta": {
            "verified_substitutions": sorted(VERIFIED_LOF_SUBSTITUTIONS),
            "n_accessions_flipped": len(flipped),
            "n_flipped_phenotyped": len(flipped_ph),
            "flipped": [{"name": r["name"], "allele": r["allele_group"], "ft16": r["FT16_mean"]} for r in flipped],
            "accuracy_delta": out["curated"]["pooled"]["accuracy"] - out["putative"]["pooled"]["accuracy"],
            "specificity_delta": out["curated"]["pooled"]["specificity"] - out["putative"]["pooled"]["specificity"],
            "sensitivity_delta": out["curated"]["pooled"]["sensitivity"] - out["putative"]["pooled"]["sensitivity"],
        },
        "circularity": {
            "verdict": "CURATED_GAIN_IS_NOT_QUOTABLE_AS_VALIDATION",
            "why": ("the paper states a016/a060/a089 were 'selected because their association with early "
                    "flowering in Figure 1', and Figure 1 is built on the SAME FT16 phenotype scored here. "
                    "L294F/L276R were therefore chosen BECAUSE their carriers are early in this data: the "
                    "curated rule's gain on FT16 is guaranteed by construction, not earned."),
            "independent_support_for_curated": ("the transgenic assay in a common null background (Figs 3-6: "
                                                "flowering time + FRI/FLC expression + GFP localisation + "
                                                "western blot) -- an experiment independent of FT16. THAT is "
                                                "the curated rule's justification; this score is not."),
            "control": {
                "rule": "circular = curated + a060 (Van-0), the allele selected the SAME way and then REFUTED",
                "mechanistically_correct": False,
                "paper_verdict_on_a060": ("'the FRI-Van-0 (a060) allele delayed flowering time in the "
                                          "transgenic lines, indicating that this allele is functional'"),
                "accuracy": out["circular"]["pooled"]["accuracy"],
                "accuracy_vs_curated": out["circular"]["pooled"]["accuracy"] - out["curated"]["pooled"]["accuracy"],
                "scores_better_than_curated": (out["circular"]["pooled"]["accuracy"]
                                               > out["curated"]["pooled"]["accuracy"]),
                "lesson": ("if the mechanistically-WRONG rule scores BETTER, then score-on-this-data cannot "
                           "adjudicate the curation -- which is exactly why the experiment, not the number, "
                           "is the warrant. 1 of the paper's 3 data-selected candidates was refuted by "
                           "experiment: the selection was not self-confirming, and neither is ours."),
            },
        },
    }

    stem = f"flowering_fri_caller_compare_{rep['date']}"
    (args.out_dir / f"{stem}.json").write_text(json.dumps(rep, indent=2), encoding="utf-8")

    print(f"[fri-compare] caller reproduces S3's own column: {rep['caller_reproduces_source_column']['verdict']} "
          f"(agree={agree}, disagree={dis})")
    print(f"{'source':11s} {'N':>4s} {'abst':>4s} {'acc':>6s} {'spec':>6s} {'sens':>6s} {'grp-wtd':>7s}")
    for k in ("s3_column", "putative", "curated", "circular"):
        h = head(k)
        print(f"{k:11s} {h['n']:4d} {h['abstained']:4d} {h['accuracy']:6.3f} {h['specificity']:6.3f} "
              f"{h['sensitivity']:6.3f} {h['group_weighted_accuracy']:7.3f}")
    c = rep["circularity"]["control"]
    print(f"\n  curated vs putative: acc {rep['curated_delta']['accuracy_delta']:+.4f} / "
          f"spec {rep['curated_delta']['specificity_delta']:+.4f} / "
          f"sens {rep['curated_delta']['sensitivity_delta']:+.4f} "
          f"({rep['curated_delta']['n_flipped_phenotyped']} accessions flipped)")
    print(f"  CIRCULARITY CONTROL: the mechanistically-WRONG rule scores {c['accuracy']:.4f} "
          f"({c['accuracy_vs_curated']:+.4f} vs curated) -> better={c['scores_better_than_curated']}")
    print(f"  -> {rep['circularity']['verdict']}")
    print(f"  -> {args.out_dir / (stem + '.json')}")
    return 0 if reproduces else 1


if __name__ == "__main__":
    raise SystemExit(main())
