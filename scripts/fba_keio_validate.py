"""Validate the FBA cell's gene-KO essentiality against the free Keio gold-standard label.

Genome-wide single-gene-deletion essentiality on iML1515 (glucose M9 aerobic) vs experimental
Keio-collection mutant fitness (Bernstein 2023 method: fitness < -2 = essential-on-glucose).

Emits `wiki/fba_keio_validation_<date>.{md,json}` with the full metric panel. Because gene
essentiality is HIGHLY class-imbalanced (~5% essential), we report accuracy + MCC + PR-AUC
alongside ROC-AUC (the literature warns ROC-AUC over-flatters on this imbalance; Bernstein 2023).

Exit 0 iff accuracy >= --min-accuracy (default 0.85; the iML1515-vs-Keio literature value is ~0.93).
"""
from __future__ import annotations

import argparse
import json
import sys
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.fba.keio import (  # noqa: E402
    KEIO_FITNESS_URL,
    confusion,
    metrics_from_confusion,
    parse_keio_fitness,
)
from dna_decode.fba.model import gene_essentiality, load_model, wildtype_growth  # noqa: E402


def _auc(scores: list[float], labels: list[bool]) -> tuple[float, float]:
    """ROC-AUC + PR-AUC (average precision). Uses sklearn if present, else a pure fallback."""
    try:
        from sklearn.metrics import average_precision_score, roc_auc_score

        return float(roc_auc_score(labels, scores)), float(average_precision_score(labels, scores))
    except Exception:
        # pure ROC-AUC (Mann-Whitney); skip PR-AUC in the fallback
        pos = [s for s, y in zip(scores, labels) if y]
        neg = [s for s, y in zip(scores, labels) if not y]
        if not pos or not neg:
            return float("nan"), float("nan")
        wins = sum((p > n) + 0.5 * (p == n) for p in pos for n in neg)
        return wins / (len(pos) * len(neg)), float("nan")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--min-accuracy", type=float, default=0.85)
    ap.add_argument("--carbon", default="D-Glucose")
    ap.add_argument("--frac", type=float, default=0.01)
    ap.add_argument("--date", default="2026-08-03")
    ap.add_argument("--fitness-tsv", default=None, help="local Keio fitness TSV (skip network)")
    a = ap.parse_args(argv)

    # 1. experimental labels
    if a.fitness_tsv:
        text = Path(a.fitness_tsv).read_text(encoding="utf-8", errors="replace")
    else:
        text = urllib.request.urlopen(
            urllib.request.Request(KEIO_FITNESS_URL, headers={"User-Agent": "dna-decode"}), timeout=60
        ).read().decode("utf-8", "replace")
    labels_full = parse_keio_fitness(text, carbon=a.carbon)
    exp = {b: d["essential"] for b, d in labels_full.items()}

    # 2. FBA genome-wide essentiality
    model = load_model()
    wt = wildtype_growth(model)
    fba = gene_essentiality(model, frac=a.frac)  # {gid: (growth, is_essential)}
    pred = {g: v[1] for g, v in fba.items()}
    growth = {g: v[0] for g, v in fba.items()}

    # 3. join + metrics (intersection of labelled + modelled genes)
    cm = confusion(exp, pred)
    met = metrics_from_confusion(cm)
    keys = sorted(set(exp) & set(pred))
    scores = [-(growth[g]) for g in keys]  # lower growth -> more essential -> higher score
    ylab = [exp[g] for g in keys]
    roc, prauc = _auc(scores, ylab)

    # 4. corroboration: FBA-essential genes ABSENT from the mutant table (no viable mutant)
    fba_ess = {g for g, e in pred.items() if e}
    absent = sorted(fba_ess - set(exp))

    result = {
        "record": "fba-keio-validation-v1",
        "date": a.date,
        "model": "iML1515",
        "label_source": "Keio BW25113 RB-TnSeq mutant fitness (Bernstein 2023); fitness<-2 = essential-on-glucose",
        "carbon_source": a.carbon,
        "wildtype_growth_per_h": round(wt, 4),
        "n_model_genes": len(pred),
        "n_labelled_genes": len(exp),
        "n_scored_intersection": cm["n"],
        "essential_prevalence": round(met["essential_prevalence"], 4),
        "confusion": cm,
        "metrics": {k: round(v, 4) for k, v in met.items()},
        "roc_auc": None if roc != roc else round(roc, 4),
        "pr_auc": None if prauc != prauc else round(prauc, 4),
        "fba_essential_absent_from_mutant_table": len(absent),
        "note_absent": (
            f"{len(absent)} FBA-essential genes have NO viable mutant in the Keio pool "
            "(an independent essentiality signal, excluded from the assayable-set metrics above)"
        ),
        "caveats": [
            "METABOLIC traits only; glucose M9 aerobic medium.",
            "Metrics computed on the ASSAYABLE gene set (genes with a measurable mutant); "
            "absolutely-essential genes without mutants are corroborated separately.",
            "Gene essentiality is highly class-imbalanced -> MCC + PR-AUC are more meaningful than ROC-AUC.",
            "In-distribution vs a published knowledge baseline; not an independent-lab claim.",
        ],
        "min_accuracy_gate": a.min_accuracy,
        "passed": met["accuracy"] >= a.min_accuracy,
    }

    outdir = Path(__file__).resolve().parent.parent / "wiki"
    outdir.mkdir(exist_ok=True)
    (outdir / f"fba_keio_validation_{a.date}.json").write_text(
        json.dumps(result, indent=2), encoding="utf-8"
    )
    md = _render_md(result)
    (outdir / f"fba_keio_validation_{a.date}.md").write_text(md, encoding="utf-8")
    # Windows consoles are cp1252 -- never let a stdout encoding error mask a completed run.
    sys.stdout.buffer.write(md.encode("utf-8", "replace"))
    sys.stdout.buffer.write(b"\n")
    return 0 if result["passed"] else 1


def _render_md(r: dict) -> str:
    m = r["metrics"]
    cm = r["confusion"]
    return f"""# FBA cell -> Keio essentiality validation ({r['date']})

**Claim tested:** the FBA metabolic cell predicts gene-KO essentiality (a cell-level trait) for ANY
iML1515 gene, validated against the free Keio mutant-fitness gold standard.

- Model: **{r['model']}** ({r['organism'] if 'organism' in r else 'E. coli K-12'}); WT growth **{r['wildtype_growth_per_h']} /h** ({r['carbon_source']} minimal aerobic)
- Label: {r['label_source']}
- Genes scored (model AND labelled): **{r['n_scored_intersection']}**  (essential prevalence {r['essential_prevalence']:.1%} -- highly imbalanced)

## Metrics (assayable gene set)

| accuracy | MCC | precision | recall | F1 | ROC-AUC | PR-AUC |
|---|---|---|---|---|---|---|
| **{m['accuracy']:.3f}** | {m['mcc']:.3f} | {m['precision']:.3f} | {m['recall']:.3f} | {m['f1']:.3f} | {r['roc_auc']} | {r['pr_auc']} |

Confusion (positive = essential): TP {cm['tp']} · FP {cm['fp']} · TN {cm['tn']} · FN {cm['fn']} (n={cm['n']})

**Corroboration:** {r['note_absent']}.

## Caveats (honest scope)
""" + "\n".join(f"- {c}" for c in r["caveats"]) + f"""

**Gate:** accuracy {m['accuracy']:.3f} {'>=' if r['passed'] else '<'} {r['min_accuracy_gate']} -> **{'PASS' if r['passed'] else 'FAIL'}**
"""


if __name__ == "__main__":
    raise SystemExit(main())
