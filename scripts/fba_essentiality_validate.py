"""Per-organism FBA gene-essentiality validation vs a free experimental gold standard.

Generalizes the E. coli Keio validation (`fba_keio_validate.py`) across organisms: load the organism's
genome-scale model, run genome-wide single-gene-deletion essentiality, join a per-organism experimental
essential-gene set, and emit the full metric panel (accuracy/MCC/precision/recall/ROC-AUC/PR-AUC).

    uv run python scripts/fba_essentiality_validate.py --organism yeast

HONEST reporting (R2): gene essentiality is highly class-imbalanced, so ACCURACY is flattered by the
majority (viable) class -- **MCC is the discrimination signal**, reported as strong (>=0.5) / moderate
(0.3-0.5) / weak (<0.3). The deliverable is the honest NUMBER, not a pass/fail on an asserted threshold.
E. coli lives in `fba_keio_validate.py` (its Keio-fitness label has its own parser); this script covers the
other organisms. Walled organisms (no fetchable+keyed gold standard) emit an honest wall artifact, no number.

Exit 0 = a real metric panel was produced; 1 = label-walled / no gold standard.
"""
from __future__ import annotations

import argparse
import json
import sys
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from dna_decode.fba.essentiality_labels import (  # noqa: E402
    ESSENTIALITY_LABEL_CONDITION,
    ESSENTIALITY_LABEL_SOURCES,
    LABEL_WALLED,
    MODEL_WALLED,
    parse_essential,
)
from dna_decode.fba.keio import confusion, metrics_from_confusion  # noqa: E402
from dna_decode.fba.medium import apply_rich_medium  # noqa: E402
from dna_decode.fba.model import gene_essentiality, load_model, resolve_model_id, wildtype_growth  # noqa: E402


def _discrimination(mcc: float) -> str:
    return "strong" if mcc >= 0.5 else ("moderate" if mcc >= 0.3 else "weak")


def _auc(scores, labels):
    try:
        from sklearn.metrics import average_precision_score, roc_auc_score
        return float(roc_auc_score(labels, scores)), float(average_precision_score(labels, scores))
    except Exception:
        return float("nan"), float("nan")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organism", required=True,
                    help="ecoli(->use fba_keio_validate) | yeast | saureus | salmonella | pputida "
                         "(paeruginosa -> MODEL_WALLED: no BiGG reconstruction exists)")
    ap.add_argument("--frac", type=float, default=0.01)
    ap.add_argument("--date", default="2026-08-03")
    ap.add_argument("--label-file", default=None, help="local gold-standard file (skip network)")
    ap.add_argument("--medium", choices=("label_matched", "default", "rich"), default="label_matched",
                    help="growth medium: label_matched (default) honours how the gold standard was "
                         "measured; `default` uses the reconstruction's own medium")
    a = ap.parse_args(argv)

    key = a.organism.strip().lower().replace(" ", "_").replace(".", "").replace("-", "_")
    outdir = Path(__file__).resolve().parent.parent / "wiki"
    outdir.mkdir(exist_ok=True)

    if key in ("escherichia_coli", "ecoli", "e_coli"):
        print("E. coli essentiality is validated by scripts/fba_keio_validate.py (Keio-fitness label).")
        return 1

    if key in MODEL_WALLED:
        wall = {
            "record": "fba-essentiality-validation-v1",
            "organism": a.organism,
            "model": None,
            "status": "MODEL_WALLED",
            "wall_kind": "external",
            "reason": MODEL_WALLED[key],
            "note": "there is NO genome-scale model for this organism to validate; the blocker is the "
                    "MODEL, not the label. Refusing to substitute another organism's model.",
        }
        (outdir / f"fba_essentiality_{key}_{a.date}.json").write_text(json.dumps(wall, indent=2), encoding="utf-8")
        print(f"{a.organism}: MODEL_WALLED -- {wall['reason']}")
        return 1

    if key in LABEL_WALLED:
        wall = {
            "record": "fba-essentiality-validation-v1",
            "organism": a.organism,
            "model": resolve_model_id(a.organism),
            "status": "LABEL_WALLED",
            "wall_kind": "external",
            "reason": LABEL_WALLED[key],
            "note": "the FBA engine RUNS on this organism; only a fetchable+model-gene-keyed essential-gene "
                    "gold standard is missing. Not code-closable without a crosswalk / a new label source.",
        }
        (outdir / f"fba_essentiality_{key}_{a.date}.json").write_text(json.dumps(wall, indent=2), encoding="utf-8")
        print(f"{a.organism}: LABEL_WALLED -- {wall['reason']}")
        return 1

    if key not in ESSENTIALITY_LABEL_SOURCES:
        print(f"no essentiality gold standard registered for '{a.organism}'", file=sys.stderr)
        return 1

    kind, url = ESSENTIALITY_LABEL_SOURCES[key]
    text = (Path(a.label_file).read_text(encoding="utf-8", errors="replace") if a.label_file
            else urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "dna-decode"}),
                                        timeout=90).read().decode("utf-8", "replace"))
    ess = parse_essential(kind, text)

    model = load_model(organism=a.organism)
    mid = resolve_model_id(a.organism)
    # Match the medium to how the LABELS were measured. Essentiality is medium-dependent, so scoring a
    # minimal-medium model against rich-medium labels charges the model for biology (yeast/iMM904 vs SGD:
    # MCC 0.2524 -> 0.3773, FP 67 -> 13). `--medium default` restores the old model-default behaviour.
    condition = (ESSENTIALITY_LABEL_CONDITION.get(key, "default") if a.medium == "label_matched"
                 else a.medium)
    supplements_opened: list[str] = []
    if condition == "rich":
        supplements_opened = apply_rich_medium(model)
    wt = wildtype_growth(model)
    fba = gene_essentiality(model, frac=a.frac)          # {gid: (growth, is_essential)}
    pred = {g: v[1] for g, v in fba.items()}
    growth = {g: v[0] for g, v in fba.items()}
    exp = {g: (g in ess) for g in pred}                  # experimental label over the model gene set

    cm = confusion(exp, pred)
    met = metrics_from_confusion(cm)
    keys = sorted(pred)
    roc, pr = _auc([-growth[g] for g in keys], [exp[g] for g in keys])

    result = {
        "record": "fba-essentiality-validation-v1",
        "date": a.date,
        "organism": a.organism,
        "model": mid,
        "status": "SCORED",
        "label_source": f"{kind}: {url}",
        "medium_mode": a.medium,
        "medium_condition": condition,
        "n_rich_supplements_opened": len(supplements_opened),
        "wildtype_growth_per_h": round(wt, 4),
        "n_gold_standard_essential": len(ess),
        "n_model_genes_scored": cm["n"],
        "essential_prevalence": round(met["essential_prevalence"], 4),
        "confusion": cm,
        "metrics": {k: round(v, 4) for k, v in met.items()},
        "roc_auc": None if roc != roc else round(roc, 4),
        "pr_auc": None if pr != pr else round(pr, 4),
        "discrimination": _discrimination(met["mcc"]),
        "caveats": [
            f"METABOLIC-gene essentiality only; medium = {condition}.",
            "Highly class-imbalanced -> MCC (not accuracy) is the discrimination signal.",
            "In-distribution vs a published knowledge baseline; not an independent-lab claim.",
            f"Medium: {condition} (mode={a.medium}). Essentiality is medium-dependent -- the SGD labels "
            "come from YPD (rich), so a minimal-medium score charges the model for biology; measured "
            "effect on yeast/iMM904 was MCC 0.2524 -> 0.3773 with FP 67 -> 13.",
        ],
    }
    (outdir / f"fba_essentiality_{key}_{a.date}.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    md = f"""# FBA essentiality validation: {a.organism} ({mid}) ({a.date})

- WT growth **{wt:.4f} /h** (default medium); genes scored **{cm['n']}**; essential prevalence {met['essential_prevalence']:.1%}
- Gold standard: {kind} ({len(ess)} experimental essential genes)

| accuracy | MCC | precision | recall | ROC-AUC | PR-AUC |
|---|---|---|---|---|---|
| {met['accuracy']:.3f} | **{met['mcc']:.3f}** | {met['precision']:.3f} | {met['recall']:.3f} | {result['roc_auc']} | {result['pr_auc']} |

Confusion (positive=essential): TP {cm['tp']} FP {cm['fp']} TN {cm['tn']} FN {cm['fn']}

**Discrimination: {result['discrimination'].upper()}** (MCC {met['mcc']:.3f}). Accuracy is flattered by the
imbalanced majority class -- MCC is the honest signal.

## Caveats
""" + "\n".join(f"- {c}" for c in result["caveats"]) + "\n"
    (outdir / f"fba_essentiality_{key}_{a.date}.md").write_text(md, encoding="utf-8")
    sys.stdout.buffer.write(md.encode("utf-8", "replace"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
