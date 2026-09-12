"""Does the genotype agree better with the ECOFF anchor or the clinical breakpoint? Measured.

STEP 5 of plans/ECOFF_Anchored_Tiering_Evaluation_Only_Technical_Plan, against the bar frozen in
wiki/ecoff_tiering_acceptance_bar.json BEFORE any carriage number existed.

THE TEST. A determinant-based decoder predicts MECHANISM PRESENCE, not treatment outcome. So among
isolates the CLINICAL breakpoint calls susceptible, the ones the ECOFF calls NON-WILD-TYPE should carry
the decoder's own determinants more often than the ones it calls wild-type. That is the only stratum
where the two anchors disagree, and therefore the only one that can discriminate between them.

THE DETERMINANT DEFINITION IS THE DEPLOYED RULE'S, AND THAT CHOICE IS LOAD-BEARING. The frozen
gentamicin rule keys on a SUBSTRING of the AMRFinder Subclass ('GENTAMICIN'), so compound subclasses
like GENTAMICIN/KANAMYCIN/TOBRAMYCIN count. Keying instead on Class=AMINOGLYCOSIDE would pull in 3,298
STREPTOMYCIN rows (aph(3'')-Ib, aph(6)-Id, aadA) and produce a result about streptomycin co-carriage
that still looked like a gentamicin finding.

WHY A PERMUTATION NULL RATHER THAN AN EFFECT-SIZE THRESHOLD. The discriminating stratum is 25 isolates
against ~2,656. Any effect size picked in advance would have been arbitrary -- the mis-specification
that sank two earlier bars in this project. Shuffling the WT/NWT label with both stratum sizes held
fixed gives a control that self-calibrates to the sizes actually present.

Read-only. No network. MIC = 2**upper (scripts/oxford_score.py:14). Frozen surface READ, never written.
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import date as _date
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from dna_decode.data.ecoff_catalog import ecoff_for, entry_for  # noqa: E402
from dna_decode.data.mic_tiers import breakpoints_for  # noqa: E402  (FROZEN -- read only)
from dna_decode.eval.amr_rules import DRUG_RULE  # noqa: E402  (FROZEN -- read only)
from dna_decode.eval.error_rates import vme_me  # noqa: E402

MIN_STRATUM = 20
SUPPORTED, WEAK, FALSIFIED, INDETERMINATE = (
    "SUPPORTED", "WEAK_DIRECTIONAL", "FALSIFIED", "INDETERMINATE")


def carriers_for_drug(amr: pd.DataFrame, drug: str) -> set[str]:
    """Genomes the DEPLOYED rule counts as carrying >=1 determinant for this drug.

    Replicates `curated_determinants_from_main`'s matching exactly: a row is drug-relevant if any of
    the drug's AMRFinder class tokens appears in Class OR Subclass, and -- when the rule carries a
    `subclass_any` refinement -- it is kept only if one of those tokens appears (as a SUBSTRING) in the
    Subclass, unless `symbol_rescue` matches the element symbol.
    """
    import re
    from dna_decode.data.mic_tiers import amrfinder_classes_for
    rule = DRUG_RULE[drug]
    classes = {c.upper() for c in amrfinder_classes_for(drug)}
    refine = {t.upper() for t in (rule.get("subclass_any") or ())} or None
    rescue = re.compile(rule["symbol_rescue"], re.I) if rule.get("symbol_rescue") else None

    cls = amr["Class"].astype(str).str.upper()
    sub = amr["Subclass"].astype(str).str.upper()
    relevant = pd.Series(False, index=amr.index)
    for c in classes:
        relevant |= cls.str.contains(c, regex=False) | sub.str.contains(c, regex=False)
    keep = relevant
    if refine is not None:
        refined = pd.Series(False, index=amr.index)
        for t in refine:
            refined |= sub.str.contains(t, regex=False)
        if rescue is not None:
            rescued = amr["Gene symbol"].astype(str).str.strip().str.match(rescue)
            refined |= rescued
        keep = relevant & refined
    return set(amr.loc[keep, "Name"].astype(str))


def permutation_gap(labels: np.ndarray, carries: np.ndarray, rng, n_perm: int = 1000) -> dict:
    """Null for the ANCHOR LABEL: shuffle which isolates are NWT, holding both stratum sizes fixed."""
    obs = float(carries[labels == 1].mean() - carries[labels == 0].mean())
    gaps = []
    for _ in range(n_perm):
        perm = rng.permutation(labels)
        gaps.append(float(carries[perm == 1].mean() - carries[perm == 0].mean()))
    g = np.array(gaps)
    return {"observed_gap": obs, "n_perm": n_perm, "perm_mean": float(g.mean()),
            "perm_max": float(g.max()), "perm_p95": float(np.percentile(g, 95)),
            "exceeds_perm_max": bool(obs > g.max())}


def verdict_from_bar(nwt_rate: float | None, wt_rate: float | None, perm: dict | None,
                     stratum_n: int, min_stratum: int = MIN_STRATUM) -> str:
    """The FROZEN rule, applied mechanically (wiki/ecoff_tiering_acceptance_bar.json)."""
    if stratum_n < min_stratum or nwt_rate is None or wt_rate is None or perm is None:
        return INDETERMINATE
    if nwt_rate <= wt_rate:
        return FALSIFIED
    return SUPPORTED if perm["exceeds_perm_max"] else WEAK


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--drug", default="gentamicin")
    ap.add_argument("--column", default="Gentamicin")
    ap.add_argument("--oxford-dir", type=Path, default=ROOT / "data" / "raw" / "oxford")
    ap.add_argument("--n-perm", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", type=Path,
                    default=ROOT / "wiki" / f"ecoff_tiering_result_{_date.today()}.json")
    a = ap.parse_args(argv)

    mic_df = pd.read_csv(a.oxford_dir / "main_data.csv", low_memory=False)
    amr = pd.read_csv(a.oxford_dir / "amrfinder.tsv", sep="\t", low_memory=False)
    scanned = set(amr["Name"].astype(str))
    carriers = carriers_for_drug(amr, a.drug)

    mic = 2.0 ** pd.to_numeric(mic_df[f"{a.column}_upper"], errors="coerce")
    df = pd.DataFrame({"guuid": mic_df["guuid"].astype(str), "mic": mic}).dropna()
    df = df[df["guuid"].isin(scanned)]          # only genomes the decoder could actually be run on
    ecoff = ecoff_for(a.drug)
    s_bp = breakpoints_for(a.drug)["clsi_s"]
    r_bp = breakpoints_for(a.drug)["clsi_r"]
    df["carries"] = df["guuid"].isin(carriers).astype(int)
    df["clinical"] = np.where(df["mic"] >= r_bp, "R", np.where(df["mic"] <= s_bp, "S", "I"))
    df["wt"] = np.where(df["mic"] > ecoff, "NWT", "WT")

    # THE discriminating stratum: the two anchors disagree only among clinically-S isolates.
    s_only = df[df["clinical"] == "S"]
    nwt, wt = s_only[s_only["wt"] == "NWT"], s_only[s_only["wt"] == "WT"]
    stratum_n = int(len(nwt))
    nwt_rate = float(nwt["carries"].mean()) if len(nwt) else None
    wt_rate = float(wt["carries"].mean()) if len(wt) else None

    perm = None
    if stratum_n >= MIN_STRATUM and len(wt):
        labels = np.concatenate([np.ones(len(nwt), int), np.zeros(len(wt), int)])
        carries = np.concatenate([nwt["carries"].to_numpy(), wt["carries"].to_numpy()])
        perm = permutation_gap(labels, carries, np.random.default_rng(a.seed), a.n_perm)

    # VME/ME under each anchor, treating the determinant call as the prediction.
    def errs(truth_pos: pd.Series) -> dict | None:
        pred = df["carries"] == 1
        tp = int((pred & truth_pos).sum()); fp = int((pred & ~truth_pos).sum())
        tn = int((~pred & ~truth_pos).sum()); fn = int((~pred & truth_pos).sum())
        e = vme_me(tp, fp, tn, fn)
        return None if e is None else {"tp": tp, "fp": fp, "tn": tn, "fn": fn,
                                       "vme": e.vme, "vme_n": e.vme_n, "me": e.me, "me_n": e.me_n}

    doc = {
        "schema": "ecoff-tiering-result-v1",
        "date": str(_date.today()),
        "drug": a.drug,
        "bar": "wiki/ecoff_tiering_acceptance_bar.json",
        "cohort": "Oxford E. coli bacteraemia deposit",
        "ecoff_mg_L": ecoff, "ecoff_source_url": entry_for(a.drug).source_url,
        "ecoff_is_tentative": "TENTATIVE" in (entry_for(a.drug).note or ""),
        "clsi_susceptible": s_bp, "clsi_resistant": r_bp,
        "n_joined": int(len(df)), "n_determinant_carriers": int(df["carries"].sum()),
        "determinant_rule": {k: (sorted(v) if isinstance(v, frozenset) else v)
                             for k, v in DRUG_RULE[a.drug].items() if k != "validated"},
        "strata": {
            "clinically_S_and_NWT": {"n": stratum_n, "carriage_rate": nwt_rate},
            "clinically_S_and_WT": {"n": int(len(wt)), "carriage_rate": wt_rate},
            "clinically_R": {"n": int((df["clinical"] == "R").sum()),
                             "carriage_rate": float(df[df["clinical"] == "R"]["carries"].mean())
                             if (df["clinical"] == "R").any() else None},
        },
        "permutation": perm,
        "error_rates_vs_clinical_anchor": errs(df["clinical"] == "R"),
        "error_rates_vs_ecoff_anchor": errs(df["wt"] == "NWT"),
        "verdict": verdict_from_bar(nwt_rate, wt_rate, perm, stratum_n),
        "verdict_is_mechanical": "verdict_from_bar applied to the frozen rule; not authored",
        "honest_limits": [
            "One organism, one UK region, 2008-2018, one collapsed MIC panel. Oxford is single-source "
            "(trips G7, n/a on G2) -- wiki/oxford_gate_screen_2026-09-11.md.",
            "The discriminating stratum is small. A FALSIFIED or WEAK_DIRECTIONAL verdict at this size "
            "is weak evidence against the claim, not a refutation.",
            "This compares CARRIAGE, not accuracy: which anchor the genotype agrees with better. It "
            "does not establish that either anchor predicts clinical outcome.",
            "A SUPPORTED verdict is NOT grounds to edit mic_tiers.py -- that is sha256-pinned in the "
            "v2 prospective lock and changing it restarts the prospective clock. User-authority call.",
            "The v2 rmt/npmA rescue is INERT on this cohort (zero rmt rows), so it cannot have "
            "influenced the result in either direction.",
        ],
    }
    a.out.write_text(json.dumps(doc, indent=2) + "\n", encoding="utf-8")
    s = doc["strata"]
    print(f"{a.drug}  ECOFF {ecoff}  CLSI-S <= {s_bp}")
    print(f"  clinically-S & NWT : n={s['clinically_S_and_NWT']['n']:>5}  "
          f"carriage {s['clinically_S_and_NWT']['carriage_rate']}")
    print(f"  clinically-S & WT  : n={s['clinically_S_and_WT']['n']:>5}  "
          f"carriage {s['clinically_S_and_WT']['carriage_rate']}")
    print(f"  clinically-R       : n={s['clinically_R']['n']:>5}  "
          f"carriage {s['clinically_R']['carriage_rate']}")
    if perm:
        print(f"  gap {perm['observed_gap']:.4f} vs perm max {perm['perm_max']:.4f}")
    print(f"  -> {doc['verdict']}\n-> {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
