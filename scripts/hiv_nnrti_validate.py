"""Validate the v0 HIV NNRTI cell against the Stanford HIVDB genotype-phenotype dataset (Wave B).

The FIRST validation of a viral cell against a FREE, independent, wet-lab genotype-phenotype label.

Contract (circularity-safe): the LABEL is the **PhenoSense fold-decreased susceptibility** (Monogram, an
independent in-vitro IC50 measurement) — NOT HIVDB's own Sierra/GRT-IS interpretation (rule-vs-rule =
circular). The PREDICTION is `dna_decode.data.hiv_amr.call_from_observed_substitutions` (the v0 class-level
NNRTI major-DRM catalog). Primary metric is **cutoff-free**: per drug, AUC = P(fold of a called-R isolate >
fold of a called-S isolate) — how well the binary genotypic call orders isolates by the lab phenotype (no
clinical-breakpoint sourcing needed). A fold>=3 sens/spec is reported as an illustrative secondary (NOT a
per-drug clinical breakpoint).

DATA: the Stanford HIVDB high-quality-filtered NNRTI dataset (gitignored at data/raw/hiv/NNRTI_DataSet.txt;
download https://hivdb.stanford.edu/download/GenoPhenoDatasets/NNRTI_DataSet.txt; cite Rhee 2003 Nucleic
Acids Res 31:298-303). Columns: SeqID, {EFV,NVP,ETR,RPV,DOR} fold-change, P1..Pn = amino acid at RT
position i ('-' = consensus). The filtered set already EXCLUDES redundant same-patient viruses + mixtures
at major positions (a built-in de-confound vs over-representation — the clonality analog).

v0 SCOPE (honest): class-level catalog (per-drug differential resistance is the documented v0.1 gap — e.g.
K103N is EFV/NVP but spares ETR/RPV, so v0 will OVER-CALL ETR/RPV; the per-drug AUC quantifies it). The
filtered set carries no Subtype column -> the per-subtype transfer check (is the catalog B-only?) is a v0.1
step needing the UNFILTERED dataset. The Stanford R-script least-squares regression is the v0.1
"validate-vs-underlying-tool" baseline (not run here).
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

from dna_decode.data.hiv_amr import (
    _RT_WT, call_from_observed_substitutions, supported_hiv_drugs,
)

DEFAULT_DATA = REPO / "data" / "raw" / "hiv" / "NNRTI_DataSet.txt"
NNRTI_POSITIONS = sorted(_RT_WT)           # the catalogued RT positions (100,101,103,106,181,188,190,230)
DRUGS = ["efavirenz", "nevirapine", "etravirine", "rilpivirine", "doravirine"]
_DRUG_COL = {"efavirenz": "EFV", "nevirapine": "NVP", "etravirine": "ETR",
             "rilpivirine": "RPV", "doravirine": "DOR"}
ILLUSTRATIVE_FOLD_CUTOFF = 3.0             # NOT a per-drug clinical breakpoint — a sensitivity check only


def _parse_fold(v: str) -> float | None:
    """Parse a fold-change cell; None on NA/blank. Censored '>'/'<' kept at the numeric bound (v0)."""
    s = (v or "").strip()
    if s in ("", "NA", "na", "."):
        return None
    s = s.lstrip("><=~")
    try:
        return float(s)
    except ValueError:
        return None


def _observed_rt_mutations(row: dict[str, str]) -> set[str]:
    """Build <WT><pos><MUT> substitutions at the catalogued NNRTI positions from the P-columns.

    A cell is the amino acid(s) at that RT position; '-' = consensus (no mutation). The filtered set has
    no mixtures at major positions, but handle multi-letter cells defensively (each letter that differs
    from the consensus-B wild-type becomes a substitution)."""
    out: set[str] = set()
    for pos in NNRTI_POSITIONS:
        cell = (row.get(f"P{pos}") or "").strip()
        if cell in ("", "-", ".", "NA"):
            continue
        wt = _RT_WT[pos]
        for aa in cell:
            if aa.isalpha() and aa != wt:
                out.add(f"{wt}{pos}{aa}")
    return out


def _auc_rank(pos_scores: list[float], neg_scores: list[float]) -> float | None:
    """AUC = P(a random pos score > a random neg score), ties=0.5 (Mann-Whitney U, average ranks)."""
    n_pos, n_neg = len(pos_scores), len(neg_scores)
    if n_pos == 0 or n_neg == 0:
        return None
    combined = sorted([(v, 1) for v in pos_scores] + [(v, 0) for v in neg_scores], key=lambda x: x[0])
    ranks = [0.0] * len(combined)
    i = 0
    while i < len(combined):
        j = i
        while j < len(combined) and combined[j][0] == combined[i][0]:
            j += 1
        avg = (i + j - 1) / 2.0 + 1.0  # 1-based average rank for the tie block
        for k in range(i, j):
            ranks[k] = avg
        i = j
    rank_sum_pos = sum(ranks[k] for k in range(len(combined)) if combined[k][1] == 1)
    u = rank_sum_pos - n_pos * (n_pos + 1) / 2.0
    return u / (n_pos * n_neg)


def load_rows(path: Path) -> list[dict[str, str]]:
    with open(path, encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")
        rows = []
        for line in f:
            vals = line.rstrip("\n").split("\t")
            if len(vals) < len(header):
                vals += [""] * (len(header) - len(vals))
            rows.append(dict(zip(header, vals)))
    return rows


def validate(path: Path = DEFAULT_DATA) -> dict:
    rows = load_rows(path)
    per_drug = {}
    for drug in DRUGS:
        col = _DRUG_COL[drug]
        fold_R, fold_S = [], []          # folds of called-R vs called-S isolates
        tp = fp = tn = fn = 0            # at the illustrative fold>=3 cutoff
        n_used = 0
        for row in rows:
            fold = _parse_fold(row.get(col, ""))
            if fold is None:
                continue
            n_used += 1
            observed = _observed_rt_mutations(row)
            call = call_from_observed_substitutions(drug, {"RT": observed}).prediction
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
        "artifact": "hiv_nnrti_v0_validation",
        "schema": "hiv-nnrti-validation-v0",
        "label_source": "Stanford HIVDB PhenoSense fold-change (independent wet-lab; NOT Sierra interpretation)",
        "caller": "dna_decode.data.hiv_amr v0 (class-level NNRTI major-DRM catalog)",
        "dataset": str(path), "n_isolates": len(rows),
        "primary_metric": "auc_call_separates_fold = P(fold|called-R > fold|called-S), cutoff-free",
        "honest_caveats": [
            "v0 is CLASS-LEVEL -> per-drug differential resistance over-calls ETR/RPV (the AUC gap quantifies it)",
            "filtered dataset has no Subtype column -> per-subtype transfer check is v0.1 (needs the unfiltered set)",
            "Stanford R-script least-squares regression is the v0.1 'validate-vs-underlying-tool' baseline (not run here)",
            "fold>=3 sens/spec is illustrative, NOT a per-drug clinical breakpoint",
        ],
        "citation": "Rhee et al. 2003 Nucleic Acids Res 31:298-303; dataset CC public per HIVDB Terms of Use",
        "per_drug": per_drug,
    }


def render_md(result: dict, generated: str) -> str:
    lines = [f"# HIV NNRTI v0 cell — validation vs Stanford HIVDB PhenoSense ({generated})", ""]
    lines.append(f"**The first validated viral cell.** Label = {result['label_source']}.")
    lines.append(f"Caller = {result['caller']}. Dataset = {result['n_isolates']} isolates "
                 f"(`{Path(result['dataset']).name}`).")
    lines.append("")
    lines.append("Primary metric (cutoff-free): **AUC = P(fold of a called-R isolate > fold of a called-S "
                 "isolate)** — how well the genotypic call orders isolates by the independent lab phenotype.")
    lines.append("")
    lines.append("| Drug | n (fold) | called R / S | **AUC** | median fold R | median fold S | sens@f3 | spec@f3 |")
    lines.append("|---|---|---|---|---|---|---|---|")
    for drug, m in result["per_drug"].items():
        lines.append(f"| {drug} | {m['n_isolates_with_fold']} | {m['n_called_R']}/{m['n_called_S']} | "
                     f"**{m['auc_call_separates_fold']}** | {m['median_fold_called_R']} | "
                     f"{m['median_fold_called_S']} | {m['illustrative_fold3_sens']} | "
                     f"{m['illustrative_fold3_spec']} |")
    lines.append("")
    lines.append("## Honest caveats")
    for c in result["honest_caveats"]:
        lines.append(f"- {c}")
    lines.append("")
    lines.append(f"Citation: {result['citation']}.")
    return "\n".join(lines)


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data", type=Path, default=DEFAULT_DATA)
    ap.add_argument("--out-md", type=Path, default=None)
    args = ap.parse_args(argv)
    if not args.data.exists():
        print(f"ERROR: dataset not found at {args.data}\n"
              f"  download: curl -L -o {args.data} "
              f"https://hivdb.stanford.edu/download/GenoPhenoDatasets/NNRTI_DataSet.txt", file=sys.stderr)
        return 2
    today = _date.today().isoformat()
    result = validate(args.data)
    out_md = args.out_md or (REPO / "wiki" / f"hiv_nnrti_v0_validation_{today}.md")
    out_md.write_text(render_md(result, today), encoding="utf-8")
    out_md.with_suffix(".json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(render_md(result, today))
    print(f"\n[wrote {out_md} + .json]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
