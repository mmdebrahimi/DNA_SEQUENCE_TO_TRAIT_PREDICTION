"""ProteinGym-NATIVE three-way: Site-Independent vs BLOSUM vs AlphaMissense on IDENTICAL mutant rows.

The row-312 conservation result cited ProteinGym's aggregated per-assay Spearman for Site-Independent while this
project's BLOSUM/AM were computed on MAVEDB rows -- a coverage-mixing caveat (a pre-build review flagged it).
This closes it: on ProteinGym's OWN assays (its DMS_score + its Site-Independent per-mutant zero-shot scores),
BLOSUM + AlphaMissense are re-scored on the EXACT SAME ProteinGym mutant rows, so all three predictors are
compared on identical variants. ProteinGym-native sidesteps the MAVEDB<->ProteinGym coordinate landmine
entirely (the review's preferred route).

Alignment: ProteinGym DMS_score is already fitness-oriented (higher = fitter). Site-Independent + BLOSUM are
fitness-oriented (higher = more tolerated) -> positive Spearman vs DMS. AlphaMissense is a PATHOGENICITY score
(higher = damaging) -> its Spearman vs DMS is negated so POSITIVE always means "predicts fitness". Single
substitutions only. INTERSECTION = mutants covered by all predictors present for the assay; FULL = each
predictor on its own coverage. Data (ProteinGym DMS + zero-shot Site-Independent + AlphaMissense) on D:.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

_SINGLE = re.compile(r"^([A-Z])(\d+)([A-Z])$")
_LEARNED_SELECTIONS = ("Activity", "Binding", "Stability", "Expression", "OrganismalFitness")


def _blosum():
    from Bio.Align import substitution_matrices
    return substitution_matrices.load("BLOSUM62")


def _bscore(mat, a, b):
    try:
        return float(mat[a, b])
    except (KeyError, IndexError):
        return float(mat[b, a])


def load_am(am_tsv: Path) -> dict:
    lut = {}
    if not am_tsv or not Path(am_tsv).exists():
        return lut
    with open(am_tsv, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) >= 3:
                try:
                    lut[(p[0], p[1])] = float(p[2])
                except ValueError:
                    pass
    return lut


def _si_column(df: pd.DataFrame) -> str | None:
    for c in df.columns:                                    # prefer the named Site-Independent column
        if "site" in c.lower() and "indep" in c.lower():
            return c
    for c in df.columns:                                    # fallback: first column with real numeric content
        if c.lower() in ("mutant", "mutated_sequence", "dms_score", "dms_score_bin"):
            continue
        if pd.to_numeric(df[c], errors="coerce").notna().sum() > 0:
            return c
    return None


def assay_threeway(dms_csv: Path, si_csv: Path | None, accession: str | None, mat, am: dict) -> dict | None:
    dms = pd.read_csv(dms_csv)
    if "mutant" not in dms.columns or "DMS_score" not in dms.columns:
        return None
    dms = dms[~dms["mutant"].astype(str).str.contains(":")]                      # single substitutions only
    dms["DMS_score"] = pd.to_numeric(dms["DMS_score"], errors="coerce")
    dms = dms.dropna(subset=["DMS_score"])
    parsed = dms["mutant"].astype(str).str.extract(_SINGLE)
    dms = dms.assign(_wt=parsed[0], _to=parsed[2]).dropna(subset=["_wt", "_to"])
    if len(dms) < 30:
        return None
    dms["blosum"] = [_bscore(mat, a, b) for a, b in zip(dms["_wt"], dms["_to"])]
    # Site-Independent
    si_map = {}
    if si_csv and Path(si_csv).exists():
        s = pd.read_csv(si_csv)
        col = _si_column(s)
        if col:
            si_map = dict(zip(s["mutant"].astype(str), pd.to_numeric(s[col], errors="coerce")))
    dms["si"] = dms["mutant"].astype(str).map(si_map)
    # AlphaMissense (negated later): key (accession, mutant)
    dms["am"] = [am.get((accession, m)) for m in dms["mutant"].astype(str)] if accession else np.nan

    def sp(col, negate=False):
        d = dms.dropna(subset=[col, "DMS_score"])
        if len(d) < 20:
            return None, 0
        r = float(spearmanr(d[col], d["DMS_score"])[0])
        return round(-r if negate else r, 4), len(d)

    # intersection rows: covered by DMS + BLOSUM(always) + SI + AM (when present)
    cov = ["blosum"]
    if dms["si"].notna().any():
        cov.append("si")
    if dms["am"].notna().any():
        cov.append("am")
    inter = dms.dropna(subset=cov + ["DMS_score"])
    out = {"n_single": len(dms), "n_intersection": len(inter)}
    for name, col, neg in (("blosum", "blosum", False), ("site_independent", "si", False), ("alphamissense", "am", True)):
        full, nf = sp(col, neg)
        out[f"{name}_full"] = full
        out[f"{name}_n"] = nf
        if len(inter) >= 20 and col in inter and inter[col].notna().all():
            ri = float(spearmanr(inter[col], inter["DMS_score"])[0])
            out[f"{name}_inter"] = round(-ri if neg else ri, 4)
        else:
            out[f"{name}_inter"] = None
    return out


def run(reference_csv: Path, dms_dir: Path, si_dir: Path | None, acc_map: Path, am_tsv: Path) -> dict:
    ref = pd.read_csv(reference_csv)
    ref = ref[ref["taxon"].astype(str).str.contains("Human", case=False, na=False)]
    acc = json.loads(Path(acc_map).read_text(encoding="utf-8")) if Path(acc_map).exists() else {}
    mat = _blosum()
    am = load_am(am_tsv)
    selc = "coarse_selection_type" if "coarse_selection_type" in ref.columns else "selection_type"
    rows = []
    for _, r in ref.iterrows():
        did = r["DMS_id"]
        dcsv = Path(dms_dir) / f"{did}.csv"
        if not dcsv.exists():
            continue
        scsv = (Path(si_dir) / f"{did}.csv") if si_dir else None
        res = assay_threeway(dcsv, scsv, acc.get(str(r["UniProt_ID"])), mat, am)
        if res:
            res.update({"DMS_id": did, "uniprot": str(r["UniProt_ID"]), "selection": r[selc]})
            rows.append(res)
    df = pd.DataFrame(rows)
    by_sel = {}
    for st in _LEARNED_SELECTIONS:
        g = df[df["selection"] == st]
        gi = g.dropna(subset=["site_independent_inter", "blosum_inter", "alphamissense_inter"])
        if len(gi) >= 3:
            by_sel[st] = {"n": len(gi),
                          "blosum": round(float(gi["blosum_inter"].median()), 4),
                          "site_independent": round(float(gi["site_independent_inter"].median()), 4),
                          "alphamissense": round(float(gi["alphamissense_inter"].median()), 4)}
    allg = df.dropna(subset=["site_independent_inter", "blosum_inter", "alphamissense_inter"])
    overall = {"n_assays": len(allg),
               "blosum": round(float(allg["blosum_inter"].median()), 4) if len(allg) else None,
               "site_independent": round(float(allg["site_independent_inter"].median()), 4) if len(allg) else None,
               "alphamissense": round(float(allg["alphamissense_inter"].median()), 4) if len(allg) else None}
    # UniProt-median (average per-protein first, then median) -- controls protein over-representation
    up = allg.groupby("uniprot")[["blosum_inter", "site_independent_inter", "alphamissense_inter"]].mean()
    overall_uniprot = {"n_proteins": len(up),
                       "blosum": round(float(up["blosum_inter"].median()), 4) if len(up) else None,
                       "site_independent": round(float(up["site_independent_inter"].median()), 4) if len(up) else None,
                       "alphamissense": round(float(up["alphamissense_inter"].median()), 4) if len(up) else None}
    act = by_sel.get("Activity")
    beats_am = [st for st, v in by_sel.items() if v["site_independent"] >= v["alphamissense"]]
    verdict = "NO_FUNCTION_ROW"
    if act:
        si, amv = act["site_independent"], act["alphamissense"]
        if si is not None and amv is not None:
            if si >= amv - 0.03:
                verdict = "DETERMINISTIC_MATCHES_LEARNED_ON_FUNCTION"
            elif beats_am or si >= amv - 0.08:
                verdict = "DETERMINISTIC_LARGELY_COMPETES"     # beats AM on >=1 modality, or ~function-parity
            else:
                verdict = "DETERMINISTIC_TRAILS_LEARNED_ON_FUNCTION"
    return {"n_assays_threeway": len(allg), "identical_rows": True,
            "by_selection_intersection": by_sel, "overall_assay_median": overall,
            "overall_uniprot_median": overall_uniprot,
            "site_independent_beats_alphamissense_on": beats_am, "verdict": verdict,
            "assays": rows}


def main(argv=None) -> int:
    import argparse
    DG = Path("D:/dna_decode_cache/proteingym")
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--reference", type=Path, default=DG / "pg_reference.csv")
    ap.add_argument("--dms-dir", type=Path, default=DG / "pg_dms" / "DMS_ProteinGym_substitutions")
    ap.add_argument("--si-dir", type=Path, default=DG / "pg_zeroshot" / "Site_Independent")
    ap.add_argument("--acc-map", type=Path, default=DG / "pg_entryname_to_accession.json")
    ap.add_argument("--am", type=Path, default=DG / "am_pg.tsv")
    ap.add_argument("--out", type=Path, default=REPO / "wiki" / "pg_native_threeway_scores.json")
    a = ap.parse_args(argv)
    res = run(a.reference, a.dms_dir, a.si_dir, a.acc_map, a.am)
    a.out.write_text(json.dumps(res, indent=2), encoding="utf-8")
    print(f"three-way assays={res['n_assays_threeway']} | verdict={res['verdict']}")
    print(f"{'selection':16} {'n':>3} {'BLOSUM':>7} {'Site-Ind':>8} {'AlphaMiss':>9}")
    for st, v in res["by_selection_intersection"].items():
        print(f"{st:16} {v['n']:>3} {v['blosum']:>7} {v['site_independent']:>8} {v['alphamissense']:>9}")
    print(f"overall (assay-median): {res['overall_assay_median']}")
    print(f"overall (uniprot-median): {res['overall_uniprot_median']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
