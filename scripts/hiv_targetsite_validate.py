"""Validate the v0 HIV PI / INSTI / CAI position-based cells vs the Stanford HIVDB genotype-phenotype data.

Generalises `scripts/hiv_nnrti_validate.py` to the three position-based target-site classes added alongside
the genome-mode caller:
  - PI    -> protease   (PI_DataSet.txt;  drugs FPV/ATV/IDV/LPV/NFV/SQV/TPV/DRV)
  - INSTI -> integrase  (INI_DataSet.txt; drugs RAL/EVG/DTG/BIC/CAB)
  - CAI   -> capsid     (CAI_DataSet.txt; drug LEN / lenacapavir)

Contract (circularity-safe, IDENTICAL to the NNRTI/NRTI cells): the LABEL is the **PhenoSense fold-decreased
susceptibility** (independent in-vitro IC50) — NOT HIVDB's own Sierra/GRT-IS interpretation (rule-vs-rule =
circular). The PREDICTION is `dna_decode.data.hiv_amr.call_hiv_observed` (the v0 position-based catalog).
Primary metric is cutoff-free: per drug, AUC = P(fold of a called-R isolate > fold of a called-S isolate).
A fold>=3 sens/spec is an illustrative secondary (NOT a per-drug clinical breakpoint).

DATA (gitignored at data/raw/hiv/; download from https://hivdb.stanford.edu/download/GenoPhenoDatasets/;
cite Rhee 2003 Nucleic Acids Res 31:298-303). Columns: SeqID, <drug fold-change cols>, P1..Pn = amino acid
at protein position i ('-' = consensus). The HQ-filtered set already excludes redundant same-patient viruses
+ mixtures at major positions (the built-in de-confound).

v0 SCOPE (honest): POSITION-BASED (like NRTI) -> DELIBERATELY over-calls non-resistant polymorphisms/
revertants at a major position; the per-drug AUC + spec quantify it. Mutant-specific deconfounded v0.1
mirrors the NRTI arc. Class-level (each drug in a class shares the major-position set) -> per-drug
differential resistance over-calls drugs a mutation spares (e.g. 2nd-gen INSTIs vs the Q148 pathway).
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
    CAI_CLASS, INSTI_CLASS, PI_CLASS, HIVTargetClass, call_hiv_observed,
)
from scripts.hiv_nnrti_validate import _auc_rank, _parse_fold, load_rows  # noqa: E402 (DRY)

ILLUSTRATIVE_FOLD_CUTOFF = 3.0

# class -> (spec, dataset filename, {full drug name: dataset fold-change column})
_CLASS_SPECS: dict[str, tuple[HIVTargetClass, str, dict[str, str]]] = {
    "PI": (PI_CLASS, "PI_DataSet.txt", {
        "fosamprenavir": "FPV", "atazanavir": "ATV", "indinavir": "IDV", "lopinavir": "LPV",
        "nelfinavir": "NFV", "saquinavir": "SQV", "tipranavir": "TPV", "darunavir": "DRV"}),
    "INSTI": (INSTI_CLASS, "INI_DataSet.txt", {
        "raltegravir": "RAL", "elvitegravir": "EVG", "dolutegravir": "DTG",
        "bictegravir": "BIC", "cabotegravir": "CAB"}),
    "CAI": (CAI_CLASS, "CAI_DataSet.txt", {"lenacapavir": "LEN"}),
}


def _observed_mutations(row: dict[str, str], cls: HIVTargetClass) -> set[str]:
    """Build <WT><pos><MUT> substitutions at the catalogued major positions from the P-columns."""
    out: set[str] = set()
    for pos in cls.positions:
        cell = (row.get(f"P{pos}") or "").strip()
        if cell in ("", "-", ".", "NA"):
            continue
        wt = cls.wt[pos]
        for aa in cell:
            if aa.isalpha() and aa != wt:
                out.add(f"{wt}{pos}{aa}")
    return out


def validate(label: str, path: Path) -> dict:
    cls, _, drug_cols = _CLASS_SPECS[label]
    rows = load_rows(path)
    per_drug = {}
    for drug, col in drug_cols.items():
        fold_R, fold_S = [], []
        tp = fp = tn = fn = n_used = 0
        for row in rows:
            fold = _parse_fold(row.get(col, ""))
            if fold is None:
                continue
            n_used += 1
            observed = _observed_mutations(row, cls)
            call = call_hiv_observed(drug, {cls.gene: observed}).prediction
            (fold_R if call == "R" else fold_S).append(fold)
            label_R = fold >= ILLUSTRATIVE_FOLD_CUTOFF
            if call == "R" and label_R:
                tp += 1
            elif call == "R" and not label_R:
                fp += 1
            elif call == "S" and not label_R:
                tn += 1
            else:
                fn += 1
        auc = _auc_rank(fold_R, fold_S)
        sens = tp / (tp + fn) if (tp + fn) else None
        spec = tn / (tn + fp) if (tn + fp) else None
        per_drug[drug] = {
            "n_isolates_with_fold": n_used,
            "n_called_R": len(fold_R), "n_called_S": len(fold_S),
            "auc_call_separates_fold": round(auc, 4) if auc is not None else None,
            "median_fold_called_R": round(statistics.median(fold_R), 2) if fold_R else None,
            "median_fold_called_S": round(statistics.median(fold_S), 2) if fold_S else None,
            "illustrative_fold3_sens": round(sens, 3) if sens is not None else None,
            "illustrative_fold3_spec": round(spec, 3) if spec is not None else None,
            "illustrative_fold3_confusion": {"tp": tp, "fp": fp, "tn": tn, "fn": fn},
        }
    return {
        "artifact": f"hiv_{label.lower()}_v0_validation",
        "schema": "hiv-targetsite-validation-v0",
        "drug_class": label, "gene": cls.gene,
        "catalog_positions": list(cls.positions),
        "label_source": "Stanford HIVDB PhenoSense fold-change (independent wet-lab; NOT Sierra interpretation)",
        "caller": f"dna_decode.data.hiv_amr v0 ({label} position-based major-DRM catalog)",
        "catalog_source": cls.source,
        "dataset": str(path), "n_isolates": len(rows),
        "primary_metric": "auc_call_separates_fold = P(fold|called-R > fold|called-S), cutoff-free",
        "call_mode": "mutant-level" if cls.major_drms is not None else "position-based",
        "honest_caveats": ([
            f"v0 is POSITION-BASED -> over-calls non-resistant polymorphisms/revertants at a major position "
            f"(the {label} spec quantifies it); mutant-specific deconfounded v0.1 mirrors the NRTI arc"
        ] if cls.major_drms is None else [
            f"v0 is MUTANT-LEVEL ({len(cls.major_drms)} catalogued substitutions); the {label} dataset is "
            f"resistance-enriched (treatment-selected isolates) so the susceptible contrast arm is small"
        ]) + [
            "class-level catalog -> per-drug differential resistance over-calls drugs a mutation spares",
            "fold>=3 sens/spec is illustrative, NOT a per-drug clinical breakpoint",
            "in-distribution vs HIVDB-PhenoSense (NOT provenance-disjoint external validation)",
        ],
        "citation": "Rhee et al. 2003 Nucleic Acids Res 31:298-303; dataset public per HIVDB Terms of Use",
        "per_drug": per_drug,
    }


def render_md(result: dict, generated: str) -> str:
    L = [f"# HIV {result['drug_class']} v0 cell — validation vs Stanford HIVDB PhenoSense ({generated})", ""]
    L.append(f"Gene **{result['gene']}**; position-based catalog at {result['catalog_positions']}.")
    L.append(f"Label = {result['label_source']}. Caller = {result['caller']}.")
    L.append(f"Dataset = {result['n_isolates']} isolates (`{Path(result['dataset']).name}`).")
    L.append(f"Catalog source: {result['catalog_source']}.")
    L.append("")
    L.append("Primary metric (cutoff-free): **AUC = P(fold of a called-R isolate > fold of a called-S isolate)**.")
    L.append("")
    L.append("| Drug | n (fold) | called R / S | **AUC** | median fold R | median fold S | sens@f3 | spec@f3 |")
    L.append("|---|---|---|---|---|---|---|---|")
    for drug, m in result["per_drug"].items():
        L.append(f"| {drug} | {m['n_isolates_with_fold']} | {m['n_called_R']}/{m['n_called_S']} | "
                 f"**{m['auc_call_separates_fold']}** | {m['median_fold_called_R']} | "
                 f"{m['median_fold_called_S']} | {m['illustrative_fold3_sens']} | "
                 f"{m['illustrative_fold3_spec']} |")
    L.append("")
    L.append("## Honest caveats")
    L += [f"- {c}" for c in result["honest_caveats"]]
    L.append("")
    L.append(f"Citation: {result['citation']}.")
    return "\n".join(L)


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--class", dest="cls", required=True, choices=sorted(_CLASS_SPECS))
    ap.add_argument("--data", type=Path, default=None)
    ap.add_argument("--out-md", type=Path, default=None)
    args = ap.parse_args(argv)
    _, fname, _ = _CLASS_SPECS[args.cls]
    data = args.data or (REPO / "data" / "raw" / "hiv" / fname)
    if not data.exists():
        print(f"ERROR: dataset not found at {data}\n"
              f"  download: curl -L -o {data} "
              f"https://hivdb.stanford.edu/download/GenoPhenoDatasets/{fname}", file=sys.stderr)
        return 2
    today = _date.today().isoformat()
    result = validate(args.cls, data)
    out_md = args.out_md or (REPO / "wiki" / f"hiv_{args.cls.lower()}_v0_validation_{today}.md")
    out_md.write_text(render_md(result, today), encoding="utf-8")
    out_md.with_suffix(".json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(render_md(result, today))
    print(f"\n[wrote {out_md} + .json]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
