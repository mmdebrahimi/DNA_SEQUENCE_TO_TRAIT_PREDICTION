"""HIV v0.2 absolute-cutoff calibration — PI + NNRTI at the DRMcv.R per-drug clinical lower cutoffs.

The PI/INSTI/NNRTI validators reported an ILLUSTRATIVE uniform fold>=3 sens/spec (the CLAUDE.md flagged
"per-drug PI clinical cutoffs would upgrade to absolute -> a v0.2 item, not fabricated"). This SOURCES those
cutoffs — NOT fabricates them — from Stanford HIVDB's own `DRMcv.R` (`cutoffmat`, the SAME authoritative
script the NRTI cell already used for `NRTI_LOWER_CUTOFF`), and re-scores the frozen catalog at each drug's
real lower cutoff, all + within-subtype-B.

Feasibility gate (verified 2026-07-03 against the fetched DRMcv.R lines 165-182):
  * PI     : all 8 lower cutoffs present (FPV/ATV/IDV/SQV=3, LPV=9, NFV=3, TPV=2, DRV=10) -> CALIBRATED.
  * NNRTI  : EFV/NVP/ETR/RPV=3 present; DOR (doravirine) POSTDATES DRMcv.R -> CUTOFF_UNAVAILABLE for DOR.
  * INSTI  : NOT in DRMcv.R at all (integrase inhibitors postdate the script) -> CUTOFF_UNAVAILABLE (a free
             Monogram INSTI cutoff exists in the literature but is NOT in this canonical free source ->
             external sourcing, deferred; reported as a wall, never guessed).

Honest note: for NNRTI every DRMcv.R lower cutoff IS 3 -> the prior illustrative fold>=3 was already the
clinical cutoff (a nice confirmation, not a change). The genuine v0.2 content is PI, where LPV(9)/TPV(2)/
DRV(10) differ from the illustrative 3, giving clinically-meaningful absolute sens/spec.

Label = PhenoSense fold (independent wet-lab; NOT Sierra). Caller = frozen dna_decode.data.hiv_amr dispatch
(PI position-based v0 -> deliberate over-call -> low spec is EXPECTED; the mutant-specific v0.1 catalogs
lift it). DATA: the .Full HIVDB datasets (carry Subtype). Cite Rhee 2003 + Stanford DRMcv.R.
"""
from __future__ import annotations

import json
import sys
from datetime import date as _date
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from dna_decode.data.hiv_amr import call_hiv_observed  # noqa: E402
from scripts.hiv_within_subtype import (  # noqa: E402
    CLASSES, _group, _observed, _parse_fold, load_rows,
)

# Stanford HIVDB DRMcv.R per-drug LOWER clinical cutoffs (fold), cutoffmat lines 165-182.
# The lower cutoff = susceptible / partially-susceptible boundary used to binarize the fold label.
# DOR + all INSTI are ABSENT from DRMcv.R (postdate it) -> deliberately NOT in this dict (a wall, not a guess).
DRMCV_LOWER_CUTOFF: dict[str, float] = {
    # PI (protease) — DRMcv.R lines 165-172
    "fosamprenavir": 3.0, "atazanavir": 3.0, "indinavir": 3.0, "lopinavir": 9.0,
    "nelfinavir": 3.0, "saquinavir": 3.0, "tipranavir": 2.0, "darunavir": 10.0,
    # NNRTI (RT) — DRMcv.R lines 179-182
    "efavirenz": 3.0, "nevirapine": 3.0, "etravirine": 3.0, "rilpivirine": 3.0,
}
MIN_GROUP_N = 15


def _confusion_at_cutoff(pairs: list[tuple[float, bool]], cutoff: float) -> dict:
    """Absolute sens/spec/balacc at `cutoff` (label_R = fold >= cutoff) for (fold, call_R) pairs."""
    n = len(pairs)
    if n < MIN_GROUP_N:
        return {"n": n, "note": "under-powered (<%d)" % MIN_GROUP_N}
    tp = sum(1 for f, c in pairs if c and f >= cutoff)
    fp = sum(1 for f, c in pairs if c and f < cutoff)
    tn = sum(1 for f, c in pairs if not c and f < cutoff)
    fn = sum(1 for f, c in pairs if not c and f >= cutoff)
    sens = tp / (tp + fn) if (tp + fn) else None
    spec = tn / (tn + fp) if (tn + fp) else None
    balacc = round((sens + spec) / 2, 3) if (sens is not None and spec is not None) else None
    return {"n": n, "n_label_R": tp + fn, "cutoff": cutoff,
            "sens": round(sens, 3) if sens is not None else None,
            "spec": round(spec, 3) if spec is not None else None,
            "balacc": balacc, "confusion": {"tp": tp, "fp": fp, "tn": tn, "fn": fn}}


def run_class(cls_name: str, path: Path) -> dict:
    cfg = CLASSES[cls_name]
    rows = [r for r in load_rows(path)
            if r.get("Method", "").strip() == "PhenoSense" and r.get("Type", "").strip() == "Clinical"]
    per_drug: dict[str, dict] = {}
    calibrated = 0
    for drug, col in cfg["drugs"].items():
        cutoff = DRMCV_LOWER_CUTOFF.get(drug)
        if cutoff is None:
            per_drug[drug] = {"status": "CUTOFF_UNAVAILABLE",
                              "reason": "no DRMcv.R clinical cutoff (drug postdates the script); "
                                        "external sourcing required — not guessed"}
            continue
        calibrated += 1
        allg: list[tuple[float, bool]] = []
        bg: list[tuple[float, bool]] = []
        for r in rows:
            fold = _parse_fold(r.get(col, ""))
            if fold is None or fold <= 0:
                continue
            call_R = call_hiv_observed(drug, _observed(r, cfg["gene"], cfg["positions"], cfg["mode"])).prediction == "R"
            allg.append((fold, call_R))
            if _group(r.get("Subtype", "")) == "B":
                bg.append((fold, call_R))
        per_drug[drug] = {"status": "CALIBRATED", "cutoff": cutoff,
                          "all": _confusion_at_cutoff(allg, cutoff),
                          "within_B": _confusion_at_cutoff(bg, cutoff)}
    return {
        "artifact": "hiv_absolute_cutoff_calibration", "schema": "hiv-absolute-cutoff-v0.2",
        "drug_class": cls_name, "gene": cfg["gene"], "call_mode": cfg["mode"],
        "cutoff_source": "Stanford HIVDB DRMcv.R cutoffmat (per-drug clinical lower cutoff) — SOURCED, not fabricated",
        "label_source": "Stanford HIVDB PhenoSense fold (independent wet-lab IC50; NOT Sierra)",
        "filter": "Method=PhenoSense AND Type=Clinical", "dataset": str(path),
        "n_clinical_phenosense": len(rows),
        "n_drugs_calibrated": calibrated, "n_drugs_total": len(cfg["drugs"]),
        "honest_caveats": [
            "cutoffs SOURCED from DRMcv.R (the script the NRTI cell used); DOR + all INSTI are ABSENT from it "
            "(postdate integrase inhibitors / doravirine) -> CUTOFF_UNAVAILABLE, reported as a wall not guessed",
            "PI/INSTI are POSITION-BASED v0 (deliberate over-call at major positions) -> LOW spec at the "
            "cutoff is EXPECTED; the mutant-specific v0.1 catalogs (hiv_pi_mutant_catalog) lift it",
            "for NNRTI every DRMcv.R lower cutoff is 3 == the prior illustrative fold>=3 -> the illustrative "
            "choice already matched the clinical cutoff (a confirmation, not a change)",
            "within-B is the powered de-confound arm (~96% B data)",
        ],
        "citation": "Rhee 2003 Nucleic Acids Res 31:298-303; cutoffs from Stanford HIVDB DRMcv.R",
        "per_drug": per_drug,
    }


def render_md(result: dict, generated: str) -> str:
    L = [f"# HIV {result['drug_class']} v0.2 absolute-cutoff calibration ({generated})", ""]
    L.append(f"**Cutoffs SOURCED from {result['cutoff_source']}.** "
             f"{result['n_drugs_calibrated']}/{result['n_drugs_total']} drugs calibrated; "
             f"the rest CUTOFF_UNAVAILABLE (external wall).")
    L.append(f"Catalog = frozen dna_decode.data.hiv_amr ({result['call_mode']}-level). "
             f"Label = {result['label_source']}. Filter = {result['filter']}; N = {result['n_clinical_phenosense']}.")
    L.append("")
    L.append("| Drug | cutoff | all sens/spec/**balacc** (n) | within-B sens/spec/**balacc** (n) |")
    L.append("|---|---|---|---|")
    for drug, m in result["per_drug"].items():
        if m.get("status") == "CUTOFF_UNAVAILABLE":
            L.append(f"| {drug} | — | _CUTOFF_UNAVAILABLE (external)_ | — |")
            continue

        def cell(x):
            if "balacc" not in x or x.get("balacc") is None:
                return f"— ({x.get('n', 0)}{'; '+x['note'] if 'note' in x else ''})"
            return f"{x['sens']}/{x['spec']}/**{x['balacc']}** ({x['n']})"
        L.append(f"| {drug} | {m['cutoff']} | {cell(m['all'])} | {cell(m['within_B'])} |")
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
    ap.add_argument("--class", dest="cls", choices=["NNRTI", "PI", "INSTI", "all"], default="all")
    ap.add_argument("--data-dir", type=Path, default=REPO / "data" / "raw" / "hiv")
    ap.add_argument("--out-dir", type=Path, default=REPO / "wiki")
    args = ap.parse_args(argv)
    today = _date.today().isoformat()
    targets = ["NNRTI", "PI", "INSTI"] if args.cls == "all" else [args.cls]
    rc = 0
    for cls_name in targets:
        path = args.data_dir / CLASSES[cls_name]["data"]
        if not path.exists():
            print(f"ERROR: {cls_name} .Full dataset not found at {path}", file=sys.stderr)
            rc = 2
            continue
        result = run_class(cls_name, path)
        stem = f"hiv_{cls_name.lower()}_absolute_cutoff_{today}"
        (args.out_dir / f"{stem}.md").write_text(render_md(result, today), encoding="utf-8")
        (args.out_dir / f"{stem}.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
        print(render_md(result, today))
        print(f"\n[wrote {args.out_dir / stem}.md + .json]\n" + "=" * 80)
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
