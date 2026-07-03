"""HIV within-subtype de-confounding check — extends the NRTI within-subtype discipline to NNRTI / PI / INSTI.

The project's core de-confounding rail (within-lineage, NOT overall AUROC — overall conflates
lineage/subtype structure with mechanism): does a class catalog's signal hold WITHIN subtype B (mechanism),
or does the class-level number ride subtype structure? The NNRTI/PI/INSTI validators report a
subtype-MIXED cutoff-free AUC (`hiv_nnrti_validate` / `hiv_targetsite_validate`); NRTI already has this check
(`hiv_nrti_within_subtype`, 2026-06-21 — the NRTI catalog HOLDS within B). This closes the gap for the other
three classes.

Method (three groups, apples-to-apples on the SAME rows):
  * Restrict to Method=PhenoSense AND Type=Clinical (mirrors the NRTI within-subtype check).
  * Split each drug's isolates into all / B / non-B by the Subtype column (the .Full datasets carry it; the
    regular filtered sets do NOT — that is why the NRTI check needed the .Full variant).
  * Per group, recompute the SAME cutoff-free AUC the class validators use:
        AUC = P(fold of a called-R isolate > fold of a called-S isolate)   (Mann-Whitney, ties=0.5)
    A within-B AUC materially > 0.5 means the deterministic call orders isolates by the INDEPENDENT lab
    phenotype INSIDE a single subtype -> mechanism, not subtype structure. The (pooled - within-B) delta is
    the direct estimate of the subtype-structure contribution to the class-mixed number.

Cutoff-free BY DESIGN (like the class validators) -> no per-drug clinical breakpoint sourced -> no cutoff
wall. Label = Stanford HIVDB PhenoSense fold (independent wet-lab IC50; NEVER Sierra = circular). Caller =
the FROZEN `dna_decode.data.hiv_amr` dispatch (NNRTI mutant-level; PI/INSTI position-based). An illustrative
fold>=3 sens/spec is reported as a secondary (NOT a per-drug clinical breakpoint).

DATA (gitignored): the .Full HIVDB datasets, which carry the Subtype column:
  data/raw/hiv/{NNRTI,PI,INI}_DataSet.Full.txt
  (download https://hivdb.stanford.edu/download/GenoPhenoDatasets/<CLASS>_DataSet.Full.txt). Cite Rhee 2003.
"""
from __future__ import annotations

import json
import statistics
import sys
from datetime import date as _date
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from dna_decode.data.hiv_amr import (  # noqa: E402
    _RT_WT, INSTI_CLASS, PI_CLASS, call_hiv_observed,
)
from scripts.hiv_nnrti_validate import _auc_rank, _parse_fold, load_rows  # noqa: E402

ILLUSTRATIVE_FOLD_CUTOFF = 3.0            # NOT a per-drug clinical breakpoint — a sensitivity check only
MIN_GROUP_N = 15                          # below this a group is under-powered (mirrors the NRTI check)
HOLD_AUC = 0.60                           # within-B median AUC >= this = materially above chance
SUBTYPE_INFLATION_DELTA = 0.10            # pooled - within-B above this = the mixed number was subtype-inflated

# Per-class config. `mode`: "mutant" (NNRTI, needs real consensus-B WT to match the mutant catalog) or
# "position" (PI/INSTI position-based — the caller parses the position from substitution[1:], so a
# placeholder WT letter is sufficient and correct). Fold columns are the .Full header codes.
_NNRTI_DRUGS = {"efavirenz": "EFV", "nevirapine": "NVP", "etravirine": "ETR",
                "rilpivirine": "RPV", "doravirine": "DOR"}
_PI_DRUGS = {"fosamprenavir": "FPV", "atazanavir": "ATV", "indinavir": "IDV", "lopinavir": "LPV",
             "nelfinavir": "NFV", "saquinavir": "SQV", "tipranavir": "TPV", "darunavir": "DRV"}
_INSTI_DRUGS = {"raltegravir": "RAL", "elvitegravir": "EVG", "dolutegravir": "DTG",
                "bictegravir": "BIC", "cabotegravir": "CAB"}

CLASSES: dict[str, dict] = {
    "NNRTI": {"data": "NNRTI_DataSet.Full.txt", "gene": "RT", "mode": "mutant",
              "positions": tuple(sorted(_RT_WT)), "drugs": _NNRTI_DRUGS},
    "PI": {"data": "PI_DataSet.Full.txt", "gene": "PR", "mode": "position",
           "positions": tuple(PI_CLASS.positions), "drugs": _PI_DRUGS},
    "INSTI": {"data": "INI_DataSet.Full.txt", "gene": "IN", "mode": "position",
              "positions": tuple(INSTI_CLASS.positions), "drugs": _INSTI_DRUGS},
}

_HIVDB_URL = "https://hivdb.stanford.edu/download/GenoPhenoDatasets"


def _group(subtype: str) -> str:
    """B vs non-B vs unknown (blank/'unknown' excluded from the split, like the NRTI check)."""
    s = (subtype or "").strip()
    return "B" if s == "B" else ("non-B" if s and s.lower() != "unknown" else "unknown")


def _observed(row: dict[str, str], gene: str, positions, mode: str) -> dict[str, set[str]]:
    """Build {gene: {substitutions}} from the P-columns at the class's catalogued positions.

    '-'/blank = consensus (no mutation). mutant mode uses the real consensus-B WT (`_RT_WT`) so the string
    matches the NNRTI mutant catalog; position mode uses a placeholder WT ('X') — the position-based caller
    parses the position from substitution[1:] and never reads the WT letter (verified)."""
    out: set[str] = set()
    for pos in positions:
        cell = (row.get(f"P{pos}") or "").strip()
        if cell in ("", "-", ".", "NA"):
            continue
        wt = _RT_WT[pos] if mode == "mutant" else "X"
        for aa in cell:
            if aa.isalpha() and aa != wt and aa != "-":
                out.add(f"{wt}{pos}{aa}")
    return {gene: out}


def _group_metrics(folds_calls: list[tuple[float, bool]]) -> dict:
    """Cutoff-free AUC (call separates fold) + illustrative fold>=3 sens/spec for one (drug, subtype) group."""
    n = len(folds_calls)
    if n < MIN_GROUP_N:
        return {"n": n, "note": "under-powered (<%d)" % MIN_GROUP_N}
    fold_R = [f for f, c in folds_calls if c]
    fold_S = [f for f, c in folds_calls if not c]
    auc = _auc_rank(fold_R, fold_S)
    tp = sum(1 for f, c in folds_calls if c and f >= ILLUSTRATIVE_FOLD_CUTOFF)
    fp = sum(1 for f, c in folds_calls if c and f < ILLUSTRATIVE_FOLD_CUTOFF)
    tn = sum(1 for f, c in folds_calls if not c and f < ILLUSTRATIVE_FOLD_CUTOFF)
    fn = sum(1 for f, c in folds_calls if not c and f >= ILLUSTRATIVE_FOLD_CUTOFF)
    n_R = sum(1 for f, _ in folds_calls if f >= ILLUSTRATIVE_FOLD_CUTOFF)
    return {
        "n": n, "n_label_R_f3": n_R, "n_called_R": len(fold_R), "n_called_S": len(fold_S),
        "auc": round(auc, 4) if auc is not None else None,
        "median_fold_R": round(statistics.median(fold_R), 2) if fold_R else None,
        "median_fold_S": round(statistics.median(fold_S), 2) if fold_S else None,
        "sens_f3": round(tp / (tp + fn), 3) if (tp + fn) else None,
        "spec_f3": round(tn / (tn + fp), 3) if (tn + fp) else None,
    }


def _class_verdict(per_drug: dict) -> dict:
    """FROZEN pure verdict from the per-drug within-B AUCs vs the pooled AUCs.

    HOLDS_WITHIN_SUBTYPE : median within-B AUC >= HOLD_AUC AND not subtype-inflated -> mechanism, not structure
    SUBTYPE_INFLATED     : median(pooled - within-B) > SUBTYPE_INFLATION_DELTA -> the mixed number rode subtype
    WITHIN_B_UNDERPOWERED: too few drugs have a powered within-B group to decide
    MIXED                : within-B present but neither clearly holds nor clearly inflates
    """
    b_aucs, deltas = [], []
    for m in per_drug.values():
        b, allg = m.get("B", {}), m.get("all", {})
        if isinstance(b.get("auc"), (int, float)):
            b_aucs.append(b["auc"])
            if isinstance(allg.get("auc"), (int, float)):
                deltas.append(allg["auc"] - b["auc"])
    if len(b_aucs) < 2:
        return {"verdict": "WITHIN_B_UNDERPOWERED", "n_powered_drugs": len(b_aucs),
                "median_within_b_auc": round(statistics.median(b_aucs), 4) if b_aucs else None}
    med_b = statistics.median(b_aucs)
    med_delta = statistics.median(deltas) if deltas else 0.0
    if med_delta > SUBTYPE_INFLATION_DELTA:
        verdict = "SUBTYPE_INFLATED"
    elif med_b >= HOLD_AUC:
        verdict = "HOLDS_WITHIN_SUBTYPE"
    else:
        verdict = "MIXED"
    return {"verdict": verdict, "n_powered_drugs": len(b_aucs),
            "median_within_b_auc": round(med_b, 4),
            "median_pooled_minus_within_b_auc": round(med_delta, 4)}


def run_class(cls_name: str, path: Path) -> dict:
    cfg = CLASSES[cls_name]
    rows = [r for r in load_rows(path)
            if r.get("Method", "").strip() == "PhenoSense" and r.get("Type", "").strip() == "Clinical"]
    subtype_counts: dict[str, int] = {}
    for r in rows:
        k = r.get("Subtype", "").strip() or "(blank)"
        subtype_counts[k] = subtype_counts.get(k, 0) + 1

    per_drug: dict[str, dict] = {}
    for drug, col in cfg["drugs"].items():
        groups: dict[str, list[tuple[float, bool]]] = {"all": [], "B": [], "non-B": []}
        for r in rows:
            fold = _parse_fold(r.get(col, ""))
            if fold is None or fold <= 0:
                continue
            obs = _observed(r, cfg["gene"], cfg["positions"], cfg["mode"])
            call_R = call_hiv_observed(drug, obs).prediction == "R"
            groups["all"].append((fold, call_R))
            g = _group(r.get("Subtype", ""))
            if g in ("B", "non-B"):
                groups[g].append((fold, call_R))
        per_drug[drug] = {g: _group_metrics(fc) for g, fc in groups.items()}

    verdict = _class_verdict(per_drug)
    nonb = sum(v for k, v in subtype_counts.items() if k not in ("B", "Unknown", "(blank)"))
    return {
        "artifact": "hiv_within_subtype_transfer", "schema": "hiv-within-subtype-v0",
        "drug_class": cls_name, "gene": cfg["gene"], "call_mode": cfg["mode"],
        "catalog": f"frozen dna_decode.data.hiv_amr ({cfg['mode']}-level, consensus-B numbering)",
        "label_source": "Stanford HIVDB PhenoSense fold (independent wet-lab IC50; NOT Sierra)",
        "metric": "cutoff-free AUC = P(fold|called-R > fold|called-S), per subtype group",
        "filter": "Method=PhenoSense AND Type=Clinical",
        "dataset": str(path), "n_clinical_phenosense": len(rows),
        "subtype_mix": {"B": subtype_counts.get("B", 0), "non_B_pooled": nonb},
        "subtype_counts": dict(sorted(subtype_counts.items(), key=lambda kv: -kv[1])),
        "class_verdict": verdict,
        "honest_caveats": [
            "the free HIVDB gp data is ~96% subtype B -> the non-B arm is UNDER-POWERED (a free-data limit, "
            "reported not hidden); within-B is the well-powered arm and is the de-confounding test that matters",
            "cutoff-free AUC needs no clinical breakpoint; the fold>=3 sens/spec is illustrative only",
            "PI/INSTI are POSITION-BASED v0 (deliberate over-call at major positions) -> their AUC has a "
            "built-in ceiling below a mutant-specific catalog; the within-vs-pooled DELTA (not the level) is "
            "the de-confounding readout",
            "a within-B AUC ~ the pooled AUC means the class-mixed number was NOT subtype-inflated; it does "
            "NOT prove non-B generalisation at scale (non-B is small)",
        ],
        "citation": "Rhee et al. 2003 Nucleic Acids Res 31:298-303; dataset CC public per HIVDB Terms of Use",
        "per_drug": per_drug,
    }


def render_md(result: dict, generated: str) -> str:
    v = result["class_verdict"]
    L = [f"# HIV {result['drug_class']} within-subtype de-confounding check ({generated})", ""]
    L.append(f"**Verdict: `{v['verdict']}`** (median within-B AUC = {v.get('median_within_b_auc')}; "
             f"median pooled−within-B = {v.get('median_pooled_minus_within_b_auc')}; "
             f"{v.get('n_powered_drugs')} powered drugs).")
    L.append("")
    L.append(f"Catalog = {result['catalog']}. Label = {result['label_source']}. "
             f"Metric = {result['metric']}. Filter = {result['filter']}; "
             f"N = {result['n_clinical_phenosense']}.")
    L.append(f"**Subtype mix:** B = {result['subtype_mix']['B']}, non-B (pooled) = "
             f"{result['subtype_mix']['non_B_pooled']} -> non-B under-powered; within-B is the test arm.")
    L.append("")
    L.append("| Drug | all AUC (n) | **B AUC (n)** | non-B AUC (n) | pooled−B |")
    L.append("|---|---|---|---|---|")
    for drug, g in result["per_drug"].items():
        def cell(x):
            if not isinstance(x.get("auc"), (int, float)):
                return f"— ({x.get('n', 0)}{'; '+x['note'] if 'note' in x else ''})"
            return f"{x['auc']} ({x['n']})"
        a, b, nb = g.get("all", {}), g.get("B", {}), g.get("non-B", {})
        delta = (round(a["auc"] - b["auc"], 3)
                 if isinstance(a.get("auc"), (int, float)) and isinstance(b.get("auc"), (int, float)) else "—")
        L.append(f"| {drug} | {cell(a)} | **{cell(b)}** | {cell(nb)} | {delta} |")
    L.append("")
    L.append("## What the verdict means")
    L.append("- **`HOLDS_WITHIN_SUBTYPE`** — the deterministic call orders isolates by the independent lab "
             "phenotype INSIDE subtype B (AUC materially > 0.5) and the pooled number is not subtype-inflated. "
             "The catalog is mechanism, not subtype structure — the same rail that NRTI cleared (2026-06-21).")
    L.append("- **`SUBTYPE_INFLATED`** — the pooled AUC exceeds the within-B AUC by more than "
             f"{SUBTYPE_INFLATION_DELTA}: the class-mixed number was riding subtype structure.")
    L.append("")
    L.append("## Honest caveats")
    for c in result["honest_caveats"]:
        L.append(f"- {c}")
    L.append("")
    L.append(f"Citation: {result['citation']}.")
    return "\n".join(L)


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--class", dest="cls", choices=[*CLASSES, "all"], default="all")
    ap.add_argument("--data-dir", type=Path, default=REPO / "data" / "raw" / "hiv")
    ap.add_argument("--out-dir", type=Path, default=REPO / "wiki")
    args = ap.parse_args(argv)
    today = _date.today().isoformat()
    targets = list(CLASSES) if args.cls == "all" else [args.cls]
    rc = 0
    for cls_name in targets:
        path = args.data_dir / CLASSES[cls_name]["data"]
        if not path.exists():
            print(f"ERROR: {cls_name} .Full dataset not found at {path}\n"
                  f"  download: curl -L -o {path} {_HIVDB_URL}/{CLASSES[cls_name]['data']}", file=sys.stderr)
            rc = 2
            continue
        result = run_class(cls_name, path)
        stem = f"hiv_{cls_name.lower()}_within_subtype_{today}"
        (args.out_dir / f"{stem}.md").write_text(render_md(result, today), encoding="utf-8")
        (args.out_dir / f"{stem}.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
        print(render_md(result, today))
        print(f"\n[wrote {args.out_dir / stem}.md + .json]\n" + "=" * 80)
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
