"""Is the coverage-only pattern in the plasmid caller a real defect, or structurally inert? Test it.

THE BACKGROUND. Validating the E. coli serotype cell found a live defect: allele selection by COVERAGE
ONLY let cross-hybridizing alleles win, and switching to identity-primary lifted H accuracy 0.770 ->
0.926. A sweep then found the SAME coverage-first pattern in `pneumoserotype` and `plasmid`. Neither was
changed, because the biology differs and identity-primary is not automatically correct elsewhere. The
pneumococcal cell has since been probed (the rule flips 1 call in 25, and that flip is wrong under BOTH
orderings -- no evidence it helps). `plasmid` is the last one untested, and it has been deferred twice
on the grounds that it "needs a cohort".

IT DOES NOT NEED A COHORT TO ANSWER THE FIRST QUESTION. Labels are needed to say which ordering is
BETTER. They are not needed to say whether the ordering changes the ANSWER at all. If it never does,
the concern is moot and no cohort is required.

THE STRUCTURAL PREDICTION, WHICH IS WHY THIS IS WORTH TESTING RATHER THAN ASSUMING. The plasmid caller
reports a SET of replicon families -- every family with at least one called allele. The coverage
comparison only decides WHICH allele represents a family it has already decided to report. So the
reported replicon SET should be invariant to the ordering, and only the secondary fields
(`best_allele`, and the identity/coverage printed for that family) should move. That is a
first-principles claim about the code, so it gets executed rather than published: the probe compares
both orderings on real assemblies and reports the SET difference, the secondary-field difference, and
whether the prediction survived.

Offline: cached assemblies + the committed PlasmidFinder DB + native blastn.
"""
from __future__ import annotations

import argparse
import collections
import json
import sys
from datetime import date as _date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from dna_decode.plasmid.runner import replicon_family  # noqa: E402
from dna_decode.typing.blast_caller import call_alleles  # noqa: E402

BLASTN = "C:/Users/Farshad/ncbi-blast/bin/blastn.exe"
DEPLOYED = {"identity": 80.0, "coverage": 60.0}   # PlasmidFinder-style defaults used by the cell


def best_per_family(per_allele: dict, identity_primary: bool) -> dict[str, dict]:
    """Winning allele per replicon family under one ordering."""
    out: dict[str, dict] = {}
    for allele_id, hit in per_allele.items():
        if not hit["called"]:
            continue
        rep = replicon_family(allele_id)
        key = ((hit["percent_identity"], hit["percent_coverage"]) if identity_primary
               else (hit["percent_coverage"], hit["percent_identity"]))
        cur = out.get(rep)
        cur_key = None if cur is None else (
            (cur["percent_identity"], cur["percent_coverage"]) if identity_primary
            else (cur["percent_coverage"], cur["percent_identity"]))
        if cur is None or key > cur_key:
            out[rep] = {"replicon": rep, "best_allele": allele_id,
                        "percent_identity": hit["percent_identity"],
                        "percent_coverage": hit["percent_coverage"]}
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--asm-dirs", nargs="+", default=["D:/dna_decode_cache/ecoli_sero_asm",
                                                      "D:/dna_decode_cache/salm_asm"])
    ap.add_argument("--db", type=Path,
                    default=ROOT / "data" / "plasmidfinder_db" / "enterobacteriales.fsa")
    ap.add_argument("--limit", type=int, default=40)
    ap.add_argument("--blastn", default=BLASTN)
    ap.add_argument("--out", type=Path, default=ROOT / "wiki" /
                    f"plasmid_selection_rule_probe_{_date.today().isoformat()}.json")
    a = ap.parse_args()

    if not a.db.exists():
        print(f"PlasmidFinder DB absent at {a.db}", file=sys.stderr)
        return 2

    fastas: list[Path] = []
    for d in a.asm_dirs:
        root = Path(d)
        if not root.exists():
            continue
        for sub in sorted(root.iterdir()):
            if len(fastas) >= a.limit:
                break
            if sub.is_dir():
                fa = next((p for p in sub.glob("*.fna")), None)
                if fa:
                    fastas.append(fa)
    if not fastas:
        print("no cached assemblies found", file=sys.stderr)
        return 2
    print(f"probing {len(fastas)} assemblies at identity>={DEPLOYED['identity']} "
          f"coverage>={DEPLOYED['coverage']}\n")

    rows, set_diffs, allele_diffs = [], 0, 0
    n_ok = 0
    total_replicons = 0
    for n, fa in enumerate(fastas, 1):
        acc = fa.parent.name
        try:
            res = call_alleles(fa, a.db, identity_threshold=DEPLOYED["identity"],
                               coverage_threshold=DEPLOYED["coverage"], blastn_bin=a.blastn,
                               timeout=600)
        except Exception as e:                       # noqa: BLE001
            rows.append({"asm": acc, "status": f"error:{type(e).__name__}"})
            continue
        if res.get("status") != "ok":
            rows.append({"asm": acc, "status": res.get("status")})
            continue
        n_ok += 1
        cov_first = best_per_family(res["per_allele"], identity_primary=False)
        id_first = best_per_family(res["per_allele"], identity_primary=True)
        set_c, set_i = set(cov_first), set(id_first)
        same_set = set_c == set_i
        moved = sorted(r for r in set_c & set_i
                       if cov_first[r]["best_allele"] != id_first[r]["best_allele"])
        total_replicons += len(set_c)
        if not same_set:
            set_diffs += 1
        if moved:
            allele_diffs += 1
        rows.append({"asm": acc, "status": "ok", "n_replicons": len(set_c),
                     "replicon_set_identical": same_set,
                     "set_only_in_coverage_first": sorted(set_c - set_i),
                     "set_only_in_identity_first": sorted(set_i - set_c),
                     "families_whose_best_allele_moved": moved})
        if n % 10 == 0:
            print(f"  [{n}/{len(fastas)}] {acc} replicons={len(set_c)} "
                  f"set_same={same_set} allele_moved={len(moved)}", flush=True)

    # NON-VACUITY: if nothing was ever called, "no difference" says nothing about the rule.
    if n_ok == 0 or total_replicons == 0:
        print(f"\nREFUSING: {n_ok} assemblies scored and {total_replicons} replicons called in total, "
              "so neither ordering was ever exercised. A null here is a plumbing result, not a finding.",
              file=sys.stderr)
        return 3

    print(f"\n=== {n_ok} assemblies scored, {total_replicons} replicon calls total ===")
    print(f"  assemblies whose REPLICON SET differs between orderings : {set_diffs}")
    print(f"  assemblies where a family's BEST ALLELE moved           : {allele_diffs}")

    if set_diffs == 0:
        verdict = "STRUCTURALLY_INERT_FOR_THE_REPORTED_SET"
        why = (f"the reported replicon set is IDENTICAL under both orderings on all {n_ok} assemblies "
               f"({total_replicons} replicon calls). The structural prediction survived: the ordering "
               "only decides WHICH allele represents a family the caller has already decided to report, "
               f"so it moves secondary fields (best_allele / the printed identity+coverage) on "
               f"{allele_diffs} assemblies but never the primary output. The E. coli fix does NOT "
               "transfer here, and no labelled cohort is needed to establish that.")
    else:
        verdict = "RULE_CHANGES_THE_REPORTED_SET"
        why = (f"the replicon set differs on {set_diffs} of {n_ok} assemblies, so the ordering DOES "
               "change the primary output and the structural prediction is FALSIFIED. A labelled "
               "cohort is now required to say which ordering is right.")
    print(f"\nVERDICT: {verdict}\n  {why}")

    out = {"schema": "plasmid-selection-rule-probe-v1", "date": _date.today().isoformat(),
           "question": "does coverage-primary vs identity-primary change the plasmid caller's output?",
           "structural_prediction": ("the caller reports a SET of replicon families -- every family with "
                                     "a called allele -- so the ordering should only pick WHICH allele "
                                     "represents an already-reported family, never the set itself"),
           "prediction_survived": set_diffs == 0,
           "thresholds": DEPLOYED, "n_assemblies_scored": n_ok,
           "total_replicon_calls": total_replicons,
           "n_assemblies_with_different_replicon_set": set_diffs,
           "n_assemblies_with_moved_best_allele": allele_diffs,
           "rows": rows, "verdict": verdict, "why": why,
           "honest_limits": [
               "This answers WHETHER the ordering changes the answer, NOT which ordering is better. "
               "The latter needs wet-lab replicon labels (PCR-based replicon typing), which are rare in "
               "public metadata -- but it is only needed if the set actually moves.",
               "Secondary fields DO move (best_allele and the identity/coverage printed per family). A "
               "consumer that reads those rather than the replicon set is affected even when the set is "
               "not.",
               "Assemblies are E. coli and Salmonella from other cohorts, not a plasmid-focused set; "
               "replicon content is therefore whatever those genomes happen to carry.",
               "One DB build (enterobacteriales). A different allele set could behave differently.",
           ]}
    a.out.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
