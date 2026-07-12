"""HIV blind-spot position-novelty flag — Family D (risk-flagged) of the genome world-model plan (2026-07-11).

QUESTION: the deployed NNRTI catalog misses some phenotypically-RESISTANT isolates that carry no catalogued
DRM (the documented blind spot; a LEARNED/likelihood rescue is a CLOSED negative — ESM below chance, mutation
burden fails, and the blind spot is pocket-mediated per the ddG probe). So D asks a DETERMINISTIC question
instead: do those catalog-missed resistant isolates carry an UN-catalogued mutation at a KNOWN DRM POSITION
(which a cheap position-novelty flag could catch), or is the resistance at other/pocket-distal positions
(confirming the closed negative)? BOTH answers are valid, honest world-model self-awareness findings.

FLAG (deterministic, no model): for an isolate, `position_novel` = it carries a non-wild-type amino acid at a
catalogued DRM position whose SPECIFIC substitution is NOT itself catalogued (i.e. removing the catalogued
DRMs, a mutation still sits at a resistance-associated position). This is a "the catalog call may be
incomplete here" self-awareness flag — NOT a resistance prediction.

MEASURED (on the free Stanford HIVDB PhenoSense fold-change; per NNRTI drug):
  * catalog-negative set = isolates the catalog calls S (`call_from_observed_substitutions` != R).
  * blind-spot = catalog-negative AND truly R (fold >= ILLUSTRATIVE_FOLD_CUTOFF=3).
  * flag_sens_on_blindspot = fraction of blind-spot isolates the position-novelty flag fires on.
  * flag_fp_on_catneg_S = fraction of catalog-negative-AND-truly-S isolates it (wrongly) fires on.
  * lift = (true-R rate among FLAGGED catalog-negative) / (true-R base rate among catalog-negative) — does the
    flag ENRICH for out-of-catalog resistance?

PRE-REGISTERED: FLAG_RECOVERS_BLINDSPOT iff median(flag_sens_on_blindspot) >= SENS_MIN (=0.30) AND
median(lift) >= LIFT_MIN (=1.2) across the powered NNRTI drugs; else BLINDSPOT_NOT_POSITION_LOCAL (the
out-of-catalog resistance is NOT explained by uncatalogued mutations at known positions — consistent with the
pocket-mediated closed negative). Honest, bounded: this is a diagnostic flag, not a decoder. Frozen surface
READ-only; catalog used as-is (no tuning).
"""
from __future__ import annotations

import argparse
import json
import statistics
import sys
from datetime import date as _date
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from dna_decode.data.hiv_amr import call_from_observed_substitutions  # noqa: E402
from scripts.hiv_nnrti_validate import (  # noqa: E402
    DEFAULT_DATA, DRUGS, ILLUSTRATIVE_FOLD_CUTOFF, _DRUG_COL, _observed_rt_mutations, _parse_fold, load_rows,
)

SENS_MIN = 0.30
LIFT_MIN = 1.2
CUTOFF = ILLUSTRATIVE_FOLD_CUTOFF


def _is_catalogued(drug: str, mut: str) -> bool:
    return call_from_observed_substitutions(drug, {"RT": {mut}}).prediction == "R"


def analyse_drug(rows, drug: str):
    col = _DRUG_COL[drug]
    catneg_R = catneg_S = 0                  # catalog-negative, split by true label
    flag_in_bs = flag_in_catneg_S = 0        # position-novel flag fires
    flagged_catneg = flagged_catneg_R = 0    # for lift
    for row in rows:
        fold = _parse_fold(row.get(col, ""))
        if fold is None or fold <= 0:
            continue
        obs = _observed_rt_mutations(row)
        catalog_R = call_from_observed_substitutions(drug, {"RT": obs}).prediction == "R"
        if catalog_R:
            continue                         # only the catalog-negative set matters for the blind spot
        true_R = fold >= CUTOFF
        # position-novel = a non-WT AA at a catalogued position whose specific substitution isn't catalogued
        position_novel = any(not _is_catalogued(drug, m) for m in obs)
        if true_R:
            catneg_R += 1
            if position_novel:
                flag_in_bs += 1
        else:
            catneg_S += 1
            if position_novel:
                flag_in_catneg_S += 1
        if position_novel:
            flagged_catneg += 1
            if true_R:
                flagged_catneg_R += 1
    n_catneg = catneg_R + catneg_S
    base_R = catneg_R / n_catneg if n_catneg else None
    flag_R_rate = flagged_catneg_R / flagged_catneg if flagged_catneg else None
    lift = (flag_R_rate / base_R) if (flag_R_rate is not None and base_R) else None
    return {
        "drug": drug, "n_catalog_negative": n_catneg,
        "n_blindspot_true_R": catneg_R, "n_catneg_true_S": catneg_S,
        "flag_sens_on_blindspot": round(flag_in_bs / catneg_R, 3) if catneg_R else None,
        "flag_fp_on_catneg_S": round(flag_in_catneg_S / catneg_S, 3) if catneg_S else None,
        "base_R_rate_catneg": round(base_R, 3) if base_R is not None else None,
        "flag_R_rate": round(flag_R_rate, 3) if flag_R_rate is not None else None,
        "lift": round(lift, 3) if lift is not None else None,
        "powered": bool(catneg_R >= 5),
    }


def run(path: Path = DEFAULT_DATA):
    rows = load_rows(path)
    per_drug = {d: analyse_drug(rows, d) for d in DRUGS}
    powered = [m for m in per_drug.values() if m["powered"]]
    sens_list = [m["flag_sens_on_blindspot"] for m in powered if m["flag_sens_on_blindspot"] is not None]
    lift_list = [m["lift"] for m in powered if m["lift"] is not None]
    med_sens = statistics.median(sens_list) if sens_list else None
    med_lift = statistics.median(lift_list) if lift_list else None
    if not powered:
        verdict = "NO_POWERED_DRUGS"
    elif med_sens is not None and med_lift is not None and med_sens >= SENS_MIN and med_lift >= LIFT_MIN:
        verdict = "FLAG_RECOVERS_BLINDSPOT"
    else:
        verdict = "BLINDSPOT_NOT_POSITION_LOCAL"
    return {
        "artifact": "hiv_blindspot_position_novelty",
        "schema": "hiv-blindspot-position-novelty-v1",
        "question": "Do catalog-MISSED resistant HIV isolates carry an uncatalogued mutation at a KNOWN DRM "
                    "position (a deterministic position-novelty flag can catch), or is the resistance elsewhere?",
        "flag": "position_novel = a non-WT AA at a catalogued DRM position whose specific substitution is NOT "
                "catalogued (a 'catalog call may be incomplete' self-awareness flag, not a resistance predictor)",
        "label_source": "Stanford HIVDB PhenoSense fold-change (free, independent wet-lab)",
        "prereg": {"SENS_MIN": SENS_MIN, "LIFT_MIN": LIFT_MIN, "fold_cutoff": CUTOFF,
                   "verdict_rule": "FLAG_RECOVERS_BLINDSPOT iff median(sens)>=SENS_MIN AND median(lift)>=LIFT_MIN"},
        "verdict": verdict, "median_flag_sens_on_blindspot": med_sens, "median_lift": med_lift,
        "honest_caveats": [
            "RISK-FLAGGED family: the LEARNED/likelihood blind-spot rescue is a CLOSED negative (ESM below "
            "chance; blind spot is pocket-mediated per the ddG probe). This tests only the DETERMINISTIC "
            "position-novelty angle.",
            "A BLINDSPOT_NOT_POSITION_LOCAL verdict CONFIRMS the closed negative (out-of-catalog resistance is "
            "not at known positions) — a valid finding, not a failure to build.",
            "The flag is a self-awareness DIAGNOSTIC ('catalog may be incomplete here'), NOT a resistance call.",
            "fold>=3 is the illustrative NNRTI cutoff (not a per-drug clinical breakpoint).",
        ],
        "citation": "Rhee 2003 Nucleic Acids Res 31:298-303",
        "per_drug": per_drug,
    }


def render_md(res, generated):
    L = [f"# HIV blind-spot position-novelty flag — is out-of-catalog resistance at known positions? ({generated})", "",
         f"**Verdict: {res['verdict']}** — median flag sensitivity on the blind spot = "
         f"{res['median_flag_sens_on_blindspot']}, median lift = {res['median_lift']} "
         f"(bars: sens>={res['prereg']['SENS_MIN']}, lift>={res['prereg']['LIFT_MIN']}).", "",
         f"{res['question']}", "",
         f"Flag: {res['flag']}. Label = {res['label_source']}.", "",
         "| drug | catalog-neg | blind-spot (true R) | **flag sens** | flag FP (on S) | base R | flag R | lift |",
         "|---|---|---|---|---|---|---|---|"]
    for m in res["per_drug"].values():
        L.append(f"| {m['drug']} | {m['n_catalog_negative']} | {m['n_blindspot_true_R']} | "
                 f"**{m['flag_sens_on_blindspot']}** | {m['flag_fp_on_catneg_S']} | {m['base_R_rate_catneg']} | "
                 f"{m['flag_R_rate']} | {m['lift']} |")
    L += ["", "## Honest caveats"] + [f"- {c}" for c in res["honest_caveats"]]
    L += ["", f"Citation: {res['citation']}."]
    return "\n".join(L)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data", type=Path, default=DEFAULT_DATA)
    ap.add_argument("--out", type=Path, default=None)
    a = ap.parse_args(argv)
    if not a.data.exists():
        print(f"ERROR: HIV NNRTI dataset absent at {a.data} (gitignored)", file=sys.stderr)
        return 2
    today = _date.today().isoformat()
    res = run(a.data)
    out = a.out or (REPO / "wiki" / f"hiv_blindspot_position_novelty_{today}.json")
    out.write_text(json.dumps(res, indent=2), encoding="utf-8")
    out.with_suffix(".md").write_text(render_md(res, today), encoding="utf-8")
    print(render_md(res, today))
    print(f"\n[wrote {out} + .md]  verdict={res['verdict']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
